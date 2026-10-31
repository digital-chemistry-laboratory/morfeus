"""Multiwfn (http://sobereva.com/multiwfn/) interface code.

Uses pexpect for interactive control.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import functools
import math
from pathlib import Path
import re
import shutil
import tempfile
import time
from types import TracebackType
import typing
from typing import Any, Iterable

if typing.TYPE_CHECKING:
    import pexpect

from morfeus.utils import (
    build_execution_env,
    Import,
    requires_dependency,
    requires_executable,
)

RealFunctionSetting = int | tuple[int, str]
VectorFunctionSetting = tuple[int, int, int] | tuple[tuple[int, int, int], str]

NUMBER_PATTERN = r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?"
ORBITAL_LINE_PATTERN = re.compile(
    rf"""
    ^\s*
    (?:Orb:\s*)?
    (?P<index>\d+)
    (?:\s+\(\s*\d+\s*\))?
    \s+
    (?:Ene|E)\(au/eV\):
    \s*{NUMBER_PATTERN}\s+(?P<energy_ev>{NUMBER_PATTERN})
    \s+Occ:\s*(?P<occupation>{NUMBER_PATTERN})
    \s+Typ(?:e)?\s*:\s*(?P<orbital_type>\S+)
    """,
    flags=re.MULTILINE | re.VERBOSE,
)
SINGLE_DETERMINANT_WFN_ERROR_TEXT = (
    "Only closed-shell single-determinant wavefunction is supported by this function"
)
# Multiwfn variants differ slightly in spacing/case/hyphenation; match both forms.
SINGLE_DETERMINANT_WFN_ERROR_PATTERN = (
    rf"{SINGLE_DETERMINANT_WFN_ERROR_TEXT}"
    r"|Only\s+closed-shell\s+single\s+determinant\s+wavefunction\s+"
    r"is\s+supported\s+by\s+this\s+function"
    r"|single[-\s]?determinant\s+wavefunction\s+is\s+supported\s+by\s+this\s+function"
)
SINGLE_DETERMINANT_WFN_ERROR_RE = re.compile(
    SINGLE_DETERMINANT_WFN_ERROR_PATTERN, flags=re.IGNORECASE
)


REAL_SPACE_FUNCTIONS: dict[str, str] = {
    "rho": "1 Electron density",
    "rho_grad_norm": "2 Gradient norm of rho",
    "rho_lapl": "3 Laplacian of rho",
    "spin_density": "5 Electron spin density",
    "kr": "6 Hamiltonian kinetic energy density",
    "gr": "7 Lagrangian kinetic energy density",
    "esp_nuclear": "8 Electrostatic potential from nuclear charges",
    "elf": "9 Electron localization function",
    "lol": "10 Localized orbital locator",
    "lie": "11 Local information entropy",
    "esp_total": "12 Total electrostatic potential",
    "rdg": "13 Reduced density gradient",  # only for grid
    "rdg_pro": "14 RDG with promolecular approximation",
    "sign": "15 Sign",
    "sign_pro": "16 Sign",
    "correlation_alpha": "17 Correlation hole for alpha",
    "alie": "18 Average local ionization energy",
    "deltag_pro": "22 Delta-g",
    "deltag_hirsh": "23 Delta-g",
    "iri": "24 Interaction region indicator",
    "vdw": "25 van der Waals potential",
}


USER_REAL_FUNCTIONS: dict[str, RealFunctionSetting] = {
    "alpha_density": 1,  # open
    "beta_density": 2,  # open
    "spatial_extent_integrand_r2rho": 3,
    "weizsaecker_potential": 4,
    "weizsaecker_functional_integrand": 5,
    "radial_distribution_rdf_4pi_r2_rho": 6,
    "local_temperature_kb_units": (7, "PNAS, 81, 8028"),
    "avg_local_esp": (8, "J. Chem. Phys., 72, 3027"),
    "shape_function_rho_over_N": 9,
    "virial_field_potential_energy": 10,
    "electron_energy": 11,
    "local_nuclear_attraction_energy": 12,
    "kinetic_energy_per_electron_G_over_rho": 13,
    "electron_esp_Vele": 14,
    "bond_metallicity": (15, "J. Phys.: Condens. Matter, 14, 10251"),
    "dimensionless_bond_metallicity": (16, "Chem. Phys. Lett., 471, 174"),
    "energy_density_per_electron": (17, "J. Chem. Phys., 117, 5529 (2002)"),
    "region_of_slow_electrons_RoSE": (18, "Chem. Phys. Lett., 582, 144 (2013)"),
    "SEDD": (19, "J. Chem. Theory Comput., 10, 3745 (2014)"),
    "DORI": (20, "J. Chem. Theory Comput., 10, 3745 (2014)"),
    "dft_linear_response_kernel": (24, "Phys. Chem. Chem. Phys., 14, 3960"),  # closed
    "electronic_momentum_fluctuation_magnitude": (25, "Theor. Chim. Acc., 127, 393"),
    "thomas_fermi_functional_integrand": (26, "J. Mol. Model., 9, 342"),
    "local_electron_affinity_EAL": (27, "J. Mol. Model., 9, 342"),
    "local_mulliken_electronegativity": (28, "J. Mol. Model., 9, 342"),
    "local_hardness": (29, "J. Mol. Model., 9, 342"),
    "ellipticity_of_density_hessian": 30,
    "eta_index": (
        31,
        "Angew. Chem. Int. Ed., 53, 2766 (2014); J. Phys. Chem. A, 114, 552 (2010)",
    ),
    "modified_eta_index_tian_lu": 32,
    "paem_pair_density_muller_xc": (33, "J. Comput. Chem., 35, 965 (2014)"),
    "paem_dft_xc_closedshell": 34,  # closed
    "G_over_V_ratio_at_bcp": (35, "J. Chem. Phys., 117, 5529 (2002)"),
    "on_top_pair_density_pi_rr": (36, "Int. J. Quantum Chem., 61, 197 (1995)"),
    # "angle_hessvec2_vs_plane_normal": ("38", "J. Phys. Chem. A, 115, 12512 (2011)"),
    "electrostatic_potential_excluding_nucleus_K": (
        39,
        "J. Phys. Chem. A, 118, 1697 (2014)",
    ),
    "steric_energy_density_weizsaecker_integrand": 40,
    "steric_potential": (41, "J. Chem. Phys., 126, 244103"),
    "steric_charge": (42, "J. Chem. Phys., 126, 244103"),
    "steric_force_magnitude": (43, "J. Chem. Phys., 126, 244103"),
    "damped_steric_potential": 44,
    "steric_force_from_damped_steric_potential": 45,
    "directly_damped_steric_force": 46,
    "shannon_entropy": 50,
    "fisher_information": 51,
    "second_fisher_information": (52, "J. Chem. Phys., 126, 191107"),
    "ghosh_entropy_GBP_local_entropy": (53, "PNAS, 81, 8028"),
    "ghosh_entropy_GBP_alt_kinetic_t": (54, "PNAS, 81, 8028"),
    "renyi_entropy_quadratic_integrand_rho2": 55,
    "renyi_entropy_cubic_integrand_rho3": 56,
    "pauli_potential": (60, "Comput. Theor. Chem., 1006, 92–99"),  # closed
    "pauli_force_magnitude": 61,  # closed
    "pauli_charge": 62,  # closed
    "quantum_potential": 63,
    "quantum_force_magnitude": 64,
    "quantum_charge": 65,
    "electrostatic_force_magnitude": 66,
    "electrostatic_charge": 67,
    "energy_density_electrostatic_sbl": 68,
    "energy_density_quantum_sbl_hamiltonian": 69,
    "energy_density_quantum_sbl_lagrangian": -69,
    "phase_space_fisher_information_PS_FID": (
        70,
        "Chem. Phys., 435, 49 (2014)",
    ),
    "electron_linear_momentum_magnitude": 74,
    "magnetic_dipole_moment_magnitude": (78, "Theoret. Chim. Acta, 6, 341"),
    "grad_norm_energy_density": 79,
    "lapl_energy_density": 80,
    "local_electron_correlation": 87,
    "local_dyn_correlation": 88,
    "local_nondyn_correlation": 89,
    "vdw": 92,
    "repulsion": 93,
    "dispersion": 94,
    "fplus": 95,
    "fminus": 96,
    "fnull": 97,
    "fdual": 98,
    "iri": 99,
    "disequilibrium_rho2": (100, "Int. J. Quantum Chem., 113, 2589 (2013)"),
    "positive_part_of_ESP": 101,
    "negative_part_of_ESP": 102,
    "magnitude_electric_field": 103,
    "ultrastrong_interaction": 819,
    "bni": 820,
    "shubin_rho_over_grad": 821,
    "local_hf_exchange_energy": 999,
    "dft_xc": 1000,
    "dft_xc_closed": 1100,  # closed
    "dft_xc_alpha": 1101,  # open
    "dft_xc_beta": 1102,  # open
}

VECTORS: dict[str, VectorFunctionSetting] = {
    "coordinates": (900, 901, 902),
    "dipole_integrand_minus": (21, 22, 23),
    "electron_linear_momentum": (71, 72, 73),
    "magnetic_dipole_moment": ((75, 76, 77), "Theoret. Chim. Acta, 6, 341"),
    "hamiltonian_kinetic_energy": (81, 82, 83),
    "lagrangian_kinetic_energy": (84, 85, 86),
}

REAL_FUNCTIONS_EXPENSIVE: set[str] = {
    "electron_esp_Vele",
    "avg_local_esp",
    "electrostatic_potential_excluding_nucleus_K",
    "paem_pair_density_muller_xc",
    "paem_dft_xc_closedshell",
    "pauli_force_magnitude",
    "pauli_charge",
    "pauli_potential",
    "energy_density_electrostatic_sbl",
    "quantum_potential",
    "quantum_force_magnitude",
    "quantum_charge",
    "electrostatic_force_magnitude",
    "electrostatic_charge",
    "positive_part_of_ESP",
    "negative_part_of_ESP",
    "magnitude_electric_field",
    "steric_charge",
    "esp_total",
}

ALL_FUNCTIONS: set[str] = (
    set(REAL_SPACE_FUNCTIONS.keys())
    | set(USER_REAL_FUNCTIONS.keys())
    | set(VECTORS.keys())
)
FAST_FUNCTIONS: set[str] = ALL_FUNCTIONS - REAL_FUNCTIONS_EXPENSIVE

CHARGE_MODELS: dict[str, str] = {
    "hirshfeld": "1 Hirshfeld atomic",
    "vdd": "2 Voronoi deformation density",
    "mulliken": "5 Mulliken atom & basis function",
    "adch": "11 Atomic dipole corrected Hirshfeld",
    "cm5": "16 CM5 atomic charge",
    "12cm5": "-16 Generate 1.2",
}

CHARGE_MODEL_CITATIONS: dict[str, str] = {
    "hirshfeld": "Theor. Chim. Acta. (Berl), 44, 129-138 (1977)",
    "vdd": "J. Comput. Chem., 25, 189.",
    "adch": (
        "Tian Lu, Feiwu Chen, Atomic dipole moment corrected Hirshfeld population "
        "method, J. Theor. Comput. Chem., 11, 163 (2012)"
    ),
}

SURFACE_MODELS: dict[str, str] = {
    "esp": "1 Electrostatic potential",
    "alie": "2 Average local ionization energy",
    "lea": "4 Local electron affinity",
    "leae": "-4 Local electron attachment energy",
    "electron": "11 Electron density",
    "sign": "12 Sign",
}

BOND_ORDER_MODELS: dict[str, str] = {
    "mayer": "1 Mayer bond order analysis",
    "wiberg": "3 Wiberg bond order analysis in Lowdin orthogonalized basis",
    "mulliken": "4 Mulliken bond order",
    "fuzzy": "7 Fuzzy bond order analysis",
    "laplacian": "8 Laplacian bond order",
}

MISSING_BASIS_FUNCTION_INFORMATION_ERROR = (
    "The input file you used does not contain basis function information!"
)
WFN_UNSUPPORTED_CHARGE_MODELS: set[str] = {"mulliken"}
WFN_UNSUPPORTED_BOND_ORDER_MODELS: set[str] = {"mayer", "wiberg", "mulliken"}

GRID_QUALITIES: dict[str, str] = {
    "low": "1 Low[- ]quality grid",
    "medium": "2 Medium[- ]quality grid",
    "high": "3 High[- ]quality grid",
}

FUZZY_PARTITIONING_SCHEMES: dict[str, str] = {
    "becke": "1 Becke",
    "hirshfeld": "3 Hirshfeld",
    "mbis": "5 MBIS",
}


@dataclass
class CommandStep:
    """Command step with optional flag for conditional execution.

    Args:
        cmd: Command to send
        expect: Pattern to expect before sending
        optional: If True, only execute if expect pattern is found in output;
                  if False, always execute regardless
    """

    cmd: str | None = None
    expect: str | None = None
    optional: bool = False


@dataclass
class WaitProgress:
    """Wait for a progress bar to complete."""

    max_wait: float = 30.0


@dataclass
class MultiwfnRunResult:
    """Result from a single Multiwfn run."""

    stdout: str
    workdir: Path


@dataclass
class MultiwfnResults:
    """Cached Multiwfn results."""

    charges: dict[str, dict[int, float]] = field(default_factory=dict)
    surfaces: dict[str, dict[str, Any]] = field(default_factory=dict)
    atomic_descriptors: dict[str, dict[int, float]] = field(default_factory=dict)
    grid_descriptors: dict[tuple[str, int], dict[int, float]] = field(
        default_factory=dict
    )
    electric_moments: dict[str, float] | None = None
    orbitals: list[dict[str, Any]] | None = None
    bond_orders: dict[str, dict[tuple[int, int], float]] = field(default_factory=dict)
    fukui: dict[str, dict[int, float]] | None = None
    superdelocalizabilities: dict[str, dict[int, float]] | None = None
    localization: dict[int, float] | None = None
    atomic_vector_descriptors: dict[str, dict[int, tuple[float, float, float]]] = field(
        default_factory=dict
    )
    citations: set[str] = field(default_factory=set)


@requires_dependency([Import("pexpect")], globals())
class _PexpectSession:
    """Pexpect session manager for Multiwfn."""

    def __init__(
        self,
        base_command: str,
        workdir: Path,
        debug: bool,
        env: dict[str, str] | None = None,
        expect_timeout: float = 3.0,
    ) -> None:
        self._child = pexpect.spawn(
            base_command,
            cwd=str(workdir),
            encoding="utf-8",
            env=env,
        )
        self._timeout = 10
        self._expect_timeout = expect_timeout
        self._debug = debug
        self._transcript: list[str] = []
        self._last_command_pos = 0

    def _append_transcript(self, chunk: Any) -> None:
        """Append text emitted by the child process."""
        if chunk is not None:
            self._transcript.append(str(chunk))

    def read_available(self) -> None:
        """Read available output without blocking."""
        try:
            chunk = self._child.read_nonblocking(size=4096, timeout=0.1)
            if chunk:
                self._append_transcript(chunk)
        except (pexpect.TIMEOUT, pexpect.EOF):
            pass

    def wait_for_progress(self, max_wait: float = 30.0) -> None:
        """Read progress output until completion, EOF, or the idle timeout."""
        last_activity = time.monotonic()
        saw_progress = False
        while True:
            try:
                chunk = self._child.read_nonblocking(size=8192, timeout=0.5)
            except pexpect.EOF:
                return
            except pexpect.TIMEOUT:
                timeout = max_wait if saw_progress else min(1.0, max_wait)
                if time.monotonic() - last_activity >= timeout:
                    return
                continue
            if not chunk:
                continue
            self._append_transcript(chunk)
            last_activity = time.monotonic()
            has_progress, percentage = self._has_progress_bar(chunk)
            if self._debug:
                print(f"[DEBUG] Progress output: {chunk!r}")
            if (percentage is not None and percentage >= 100) or (
                saw_progress and not has_progress
            ):
                return
            saw_progress = saw_progress or has_progress

    def send(self, cmd: str) -> None:
        """Send command to process.

        Args:
            cmd: Command to send
        """
        self._child.sendline(cmd)
        if self._debug:
            print(f"[DEBUG] Sent: {cmd!r}")
        self._last_command_pos = len(self._transcript)

    def _has_progress_bar(self, text: str) -> tuple[bool, float | None]:
        """Check if text contains a progress bar and extract percentage.

        Args:
            text: Text chunk to inspect.

        Returns:
            Tuple of (has_progress_bar, percentage).
        """
        bar_match = re.search(r"\[([#=\-]+)\]", text)
        if bar_match:
            bar = bar_match.group(1)
            if "-" not in bar:
                return True, 100.0
        # Match patterns like: "Progress: [###---] 14.3 %" or "14.3 %"
        progress_pattern = r"Progress:\s*\[[#-]+\]\s+([\d.]+)\s*%"
        match = re.search(progress_pattern, text)
        if match:
            try:
                percentage = float(match.group(1))
                return True, percentage
            except ValueError:
                return True, None
        return False, None

    def expect(self, pattern: str, timeout: float | None = None) -> bool:
        """Wait for a pattern, returning False on timeout or process exit."""
        if self._debug:
            print(f"[DEBUG] Expecting: {pattern!r}")
        try:
            self._child.expect(
                pattern, timeout=self._expect_timeout if timeout is None else timeout
            )
            self._append_transcript(self._child.before)
            self._append_transcript(self._child.after)
            return True
        except pexpect.TIMEOUT:
            self.read_available()
        except pexpect.EOF:
            self._append_transcript(self._child.before)
        if self._debug:
            print(f"[DEBUG] Pattern not found: {pattern!r}")
        return False

    def load(self, file_path: Path, reset: bool = False) -> None:
        """Load an input after optionally resetting an existing process."""
        self._transcript.clear()
        self._last_command_pos = 0
        if reset:
            self.send("r")
        if not self.expect("Now input file path"):
            raise RuntimeError("Multiwfn did not reach the input-file prompt.")
        self.send(str(file_path))
        self.wait_for_main_menu()

    def wait_for_main_menu(self) -> None:
        """Require the complete main menu before reusing the process."""
        pattern = "300 Other functions"
        if not self._has_pattern_in_buffered_output(pattern) and not self.expect(
            pattern
        ):
            raise RuntimeError("Multiwfn did not return to the main menu.")

    def close(self, graceful: bool = False) -> None:
        """Release the child process, optionally quitting from the main menu."""
        try:
            if graceful:
                self.send("q")
                self.wait_for_exit()
        finally:
            self._child.close(force=True)

    def wait_for_exit(self) -> None:
        """Wait for process to exit."""
        try:
            self._child.expect(pexpect.EOF, timeout=self._timeout)
        except pexpect.TIMEOUT:
            if self._debug:
                print("[WARN] Timeout waiting for EOF")
        self._append_transcript(self._child.before)
        if self._child.isalive():
            self._child.close()

    def get_output_since_last_command(self) -> str:
        """Get output since last command."""
        return "".join(self._transcript[self._last_command_pos :])

    @property
    def stdout(self) -> str:
        """Full session transcript."""
        return "".join(self._transcript)

    def _has_pattern_in_buffered_output(self, pattern: str) -> bool:
        """Check output since the last command for an already consumed prompt."""
        return re.search(pattern, self.get_output_since_last_command()) is not None

    def _process_command(
        self, item: str | CommandStep | WaitProgress, use_robust: bool
    ) -> None:
        """Send a command after checking its required or optional prompt."""
        if isinstance(item, WaitProgress):
            self.wait_for_progress(max_wait=item.max_wait)
            return
        if isinstance(item, str):
            item = CommandStep(item)
        if item.expect and (item.optional or (use_robust and item.cmd is not None)):
            matched = self._has_pattern_in_buffered_output(item.expect) or self.expect(
                item.expect, timeout=1.0 if item.optional else None
            )
            if not matched:
                if item.optional:
                    return
                if SINGLE_DETERMINANT_WFN_ERROR_RE.search(self.stdout):
                    raise RuntimeError(SINGLE_DETERMINANT_WFN_ERROR_TEXT)
                raise RuntimeError(
                    f"Expected pattern not found before command {item.cmd!r}: {item.expect!r}"
                )
        if item.cmd is not None:
            self.send(item.cmd)
            if not use_robust:
                self.read_available()


class Multiwfn:
    """Multiwfn interface using pexpect for interactive control.

    Args:
        file_path: Path to molden/wfn file from xtb or other QM program.
        run_path: Optional archive directory. Each instance creates a unique
            subdirectory whose results are retained. Otherwise all intermediate
            files are temporary and removed when the instance closes.
        robust_mode: If True, wait for expect patterns between commands (slower but
            safer).
            If False, send commands without waiting (faster but may desync).
        debug: If True, print debug info including sent commands and expect patterns.
        env_variables: Environment variables to use for the Multiwfn process.
        max_wait: Maximum time to wait for progress updates in seconds.
        has_spin: Spin state. If None, detect from file for wavefunction inputs
            (True=open-shell, False=closed-shell). For .cub/.grd inputs, spin
            state is left undefined.
        expect_timeout: Time to wait for expect patterns in seconds.
        persistent: Reuse the process and loaded wavefunction between analyses.
            Use a context manager or close() to release the process.
            Do not share an instance across concurrent tasks.
    """

    def __init__(
        self,
        file_path: str | Path,
        run_path: str | Path | None = None,
        robust_mode: bool = True,
        debug: bool = False,
        env_variables: dict[str, str] | None = None,
        max_wait: int = 30,
        has_spin: bool | None = None,
        expect_timeout: float = 3.0,
        persistent: bool = False,
    ) -> None:
        self._persistent = persistent
        self._session: _PexpectSession | None = None
        self._loaded_file: Path | None = None
        self._session_dir: tempfile.TemporaryDirectory[str] | None = None
        self._file_path = Path(file_path).resolve()
        if not self._file_path.exists():
            raise FileNotFoundError(f"File not found: {self._file_path}")

        if run_path is not None:
            archive = Path(run_path).resolve()
            archive.mkdir(parents=True, exist_ok=True)
            self._run_path = Path(tempfile.mkdtemp(prefix="morfeus_", dir=archive))
            self._temp_dir: tempfile.TemporaryDirectory[str] | None = None
        else:
            self._temp_dir = tempfile.TemporaryDirectory(prefix="morfeus_multiwfn_")
            self._run_path = Path(self._temp_dir.name).resolve()
        self._robust_mode = robust_mode
        self._debug = debug
        self._env_variables = env_variables
        self._results = MultiwfnResults()
        self._settings_ini_path: Path | None = None
        self._max_wait = max_wait
        self._expect_timeout = expect_timeout

        self._results.citations = {
            # Both following papers ***MUST BE CITED IN MAIN TEXT*** if Multiwfn is used:
            # See "How to cite Multiwfn.pdf" in Multiwfn binary package for more
            # information.
            (
                "Tian Lu, Feiwu Chen, J. Comput. Chem., 33, 580 (2012) DOI: "
                "10.1002/jcc.22885"
            ),
            "Tian Lu, J. Chem. Phys., 161, 082503 (2024) DOI: 10.1063/5.0216272",
        }
        self._has_spin = has_spin
        if self._has_spin is None and self._file_path.suffix.lower() not in {
            ".cub",
            ".grd",
        }:
            self._has_spin = self._detect_spin_state()

    def _detect_spin_state(self) -> bool:
        """Load the molden file and parse the multiplicity."""
        result = self._run_commands([])
        mult_pattern = r"Expected multiplicity:\s*(\d+)"

        match = re.search(mult_pattern, result.stdout)
        if not match:
            self._close_session()
            raise RuntimeError("Could not detect multiplicity from the file.")

        mult = int(match.group(1))
        return mult > 1

    def load_settingini(self, file_path: str | Path) -> None:
        """Load custom settings.ini file for all future runs.

        Args:
            file_path: Path to custom settings.ini file.

        Raises:
            FileNotFoundError: If settings.ini does not exist.
        """
        settings_path = Path(file_path).resolve()
        if not settings_path.exists() or settings_path.name != "settings.ini":
            raise FileNotFoundError(f"Settings file not found: {settings_path}")
        self._settings_ini_path = settings_path
        self._close_session()
        self._results = MultiwfnResults(citations=set(self._results.citations))

    def load_settingsini(self, file_path: str | Path) -> None:
        """Backward-compatible alias for load_settingini()."""
        self.load_settingini(file_path)

    def _close_session(self, graceful: bool = True) -> None:
        """Discard a reusable process and its private output directory."""
        session = self._session
        self._session = None
        self._loaded_file = None
        try:
            if session is not None:
                session.close(graceful=graceful)
        finally:
            if self._session_dir is not None:
                self._session_dir.cleanup()
                self._session_dir = None

    def close(self) -> None:
        """Release the process and clean up the temporary working directory."""
        try:
            self._close_session()
        finally:
            if self._temp_dir is not None:
                self._temp_dir.cleanup()
                self._temp_dir = None

    def __enter__(self) -> "Multiwfn":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        del exc_type, exc_value, traceback
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _copy_settings_ini(self, workdir: str | Path) -> None:
        """Copy configured settings.ini into a working directory."""
        if self._settings_ini_path is None:
            return

        target_dir = Path(workdir)
        target_dir.mkdir(parents=True, exist_ok=True)
        target = target_dir / "settings.ini"
        if target.resolve() != self._settings_ini_path:
            shutil.copy2(self._settings_ini_path, target)

    def _commands_change_iuserfunc(self, function: int = -1) -> list[CommandStep]:
        """Adjust the iuserfunc setting.

        Refer to a complete list (ctrl F: !Below functions can be selected by
        "iuserfunc" parameter):
        https://github.com/foxtran/MultiWFN/blob/master/src/function.f90

        Args:
            function: User-defined function index.

        Returns:
            Command steps to apply iuserfunc.
        """
        return [
            CommandStep("1000", expect="300 Other functions"),
            CommandStep("2", expect="2 Set iuserfunc"),
            CommandStep(
                str(function), expect="Input index of the user-defined function"
            ),
        ]

    def _parse_custom_commands(self, commands: list) -> list[CommandStep]:
        custom_commands = []
        for step in commands:
            if not isinstance(step, CommandStep):
                cmd = None
                expect = None
                optional = False
                if isinstance(step, str):
                    cmd = step
                elif isinstance(step, tuple) and len(step) in (2, 3):
                    cmd = step[0]
                    expect = step[1]
                    if len(step) == 3:
                        optional = step[2]
                else:
                    raise ValueError(f"Invalid custom command: {step}")
                step = CommandStep(cmd, expect, optional)
            custom_commands.append(step)

        return custom_commands

    def _commands_fuzzy_partitioning(
        self, partitioning: str = "becke"
    ) -> list[CommandStep]:
        """Enter fuzzy analysis and select its atomic partitioning scheme."""
        partition_cmd = FUZZY_PARTITIONING_SCHEMES[partitioning].split(" ")[0]
        return [
            CommandStep("15", expect="15 Fuzzy atomic space analysis"),
            CommandStep("-1", expect="-1 Select method for partitioning atomic spaces"),
            CommandStep(partition_cmd, expect="Select atomic space partition method"),
        ]

    def _build_fuzzy_integration_commands(
        self,
        menu_cmd: str,
        menu_pattern: str,
        prefix_commands: list[CommandStep] | None = None,
        partitioning: str = "becke",
    ) -> list[CommandStep | WaitProgress]:
        """Build command sequence for fuzzy atomic space integration."""
        commands: list[CommandStep | WaitProgress] = list(prefix_commands or [])
        commands.extend(
            [
                *self._commands_fuzzy_partitioning(partitioning),
                CommandStep("1", expect="1 Perform integration in fuzzy atomic spaces"),
                CommandStep(menu_cmd, expect=menu_pattern),
                WaitProgress(self._max_wait),
                CommandStep("0", expect="0 Return"),
            ]
        )
        return commands

    def _setup_workdir(self, subdir: str | None) -> Path:
        """Create a fresh result directory for an uncached calculation."""
        label = Path(subdir) if subdir else Path("analysis")
        if label.is_absolute() or ".." in label.parts:
            raise ValueError("subdir must be relative to the run directory.")
        parent = self._run_path / label
        parent.mkdir(parents=True, exist_ok=True)
        return Path(tempfile.mkdtemp(prefix="run-", dir=parent))

    @staticmethod
    def _normalize_option(value: str) -> str:
        """Normalize option keys to lower-case."""
        return value.strip().lower()

    def _is_wfn_input(self) -> bool:
        """Return True when the current input file is a .wfn file."""
        return self._file_path.suffix.lower() == ".wfn"

    def list_options(self) -> dict[str, list[str]]:
        """List available option names for common Multiwfn analyses."""
        return {
            "charges": sorted(CHARGE_MODELS),
            "surface": sorted(SURFACE_MODELS),
            "bond_order": sorted(BOND_ORDER_MODELS),
            "grid_quality": sorted(GRID_QUALITIES),
            "descriptors": sorted(ALL_FUNCTIONS),
            "descriptors_fast": sorted(FAST_FUNCTIONS),
        }

    def get_citations(self) -> list[str]:
        """Get all collected citations for the current object."""
        return sorted(self._results.citations)

    def run_commands(
        self,
        commands: Iterable[str | CommandStep | WaitProgress],
        subdir: str | None = None,
    ) -> MultiwfnRunResult:
        """Run custom commands that return to the main menu.

        Custom commands invalidate cached results and the loaded input because
        they may modify the wavefunction or analysis settings.

        Args:
            commands: Commands returning to the main menu, without a final quit.
            subdir: Subdirectory within run_path to store results.

        Returns:
            MultiwfnRunResult with this run's stdout and workdir.
        """
        commands = list(commands)
        if commands:
            last = commands[-1]
            if (last.cmd if isinstance(last, CommandStep) else last) == "q":
                commands.pop()
        self._loaded_file = None
        try:
            result: MultiwfnRunResult = self._run_commands(commands, subdir)
            return result
        finally:
            self._loaded_file = None
            self._results = MultiwfnResults(citations=set(self._results.citations))

    @requires_executable(["Multiwfn"])
    def _run_commands(
        self,
        commands: Iterable[str | CommandStep | WaitProgress],
        subdir: str | None = None,
    ) -> MultiwfnRunResult:
        """Run a descriptor from main menu to main menu, managing process lifetime."""
        workdir = self._setup_workdir(subdir)
        try:
            new_session = self._session is None
            if self._session is None:
                session_path = workdir
                if self._persistent:
                    self._session_dir = tempfile.TemporaryDirectory(
                        prefix=".multiwfn-session-", dir=self._run_path
                    )
                    session_path = Path(self._session_dir.name)
                self._copy_settings_ini(session_path)
                self._session = _PexpectSession(
                    "Multiwfn",
                    session_path,
                    self._debug,
                    env=build_execution_env(env_variables=self._env_variables),
                    expect_timeout=self._expect_timeout,
                )
            session = self._session
            if self._loaded_file != self._file_path:
                session.load(self._file_path, reset=not new_session)
                self._loaded_file = self._file_path
                start = 0
            else:
                # Retain the last menu for the first command's prompt check.
                session._transcript = session._transcript[session._last_command_pos :]
                session._last_command_pos = 0
                start = len(session._transcript)
            for item in commands:
                session._process_command(item, self._robust_mode)
            session.wait_for_main_menu()
            if self._persistent:
                self._collect_session_files(workdir)
            return MultiwfnRunResult(
                stdout="".join(session._transcript[start:]), workdir=workdir
            )
        except BaseException:
            self._close_session(graceful=False)
            raise
        finally:
            if not self._persistent:
                self._close_session()

    def _collect_session_files(self, workdir: Path) -> None:
        """Move this run's output files into its public result directory."""
        if self._session_dir is None:
            raise RuntimeError("Missing persistent session directory.")
        for source in Path(self._session_dir.name).iterdir():
            if source.name == "settings.ini":
                continue
            if not source.is_file():
                raise RuntimeError(f"Unsupported Multiwfn output directory: {source}")
            shutil.copy2(source, workdir / source.name)
            source.unlink()

    def get_charges(self, model: str = "adch") -> dict[int, float]:
        """Calculate atomic charges.

        Args:
            model: Charge model to use.
                - 'hirshfeld': Standard Hirshfeld charges
                - 'adch': Atomic dipole moment corrected Hirshfeld charges
                - 'vdd': Voronoi deformation density
                - 'mulliken': Mulliken charges
                - 'cm5': CM5 charges
                - '12cm5': 1.2*CM5 charges

        Raises:
            ValueError: If given charge model is not available.
            RuntimeError: If the charge output file is missing or cannot be parsed.

        Returns:
            Atomic charges indexed by atom number (1-indexed)
        """
        normalized_model = self._normalize_option(model)
        if normalized_model in self._results.charges:
            return self._results.charges[normalized_model]

        if normalized_model not in CHARGE_MODELS:
            choices = ", ".join(CHARGE_MODELS.keys())
            raise ValueError(
                f"Charge model {model!r} not supported. Choose between {choices}."
            )
        if self._is_wfn_input() and normalized_model in WFN_UNSUPPORTED_CHARGE_MODELS:
            raise ValueError(MISSING_BASIS_FUNCTION_INFORMATION_ERROR)

        menu_pattern = CHARGE_MODELS[normalized_model]
        menu_cmd = menu_pattern.split(" ")[0]

        commands = [
            CommandStep("7", expect="7 Population analysis"),
            CommandStep(menu_cmd, expect=menu_pattern),
            CommandStep(
                "1",
                expect=(
                    "1 Use buil[dt]-in sphericalized atomic"
                    if normalized_model != "mulliken"
                    else " 1 Output Mulliken population and atomic"
                ),
            ),
            WaitProgress(self._max_wait),
            CommandStep("y", expect="y/n"),
        ]
        if normalized_model == "mulliken":
            commands += [CommandStep("0", expect="Mulliken population")]
        commands.append(CommandStep("0", expect="0 Return"))

        citation = CHARGE_MODEL_CITATIONS.get(normalized_model, None)
        if citation is not None:
            self._results.citations.add(citation)

        result = self._run_commands(commands, subdir=normalized_model)

        chg_files = sorted(result.workdir.glob("*.chg"))
        if not chg_files:
            raise RuntimeError(f"Multiwfn produced no charge file in {result.workdir}.")
        charges = self._parse_chg_file(chg_files[0])
        if not charges:
            raise RuntimeError(f"Could not parse atomic charges from {chg_files[0]}.")

        self._results.charges[normalized_model] = charges

        return charges

    def get_surface(self, model: str = "esp") -> dict[str, Any]:
        """Calculate molecular surface properties.

        Args:
            model: Surface analysis model to use.
                - 'esp': Electrostatic Potential
                - 'alie': Average Local Ionization Energy
                - 'lea': Local Electron Affinity
                - 'leae': Local Electron Attachment Energy
                - 'electron': Electron Density
                - 'sign': Sign(lambda2)*rho

        Raises:
            ValueError: If given model is not available.

        Returns:
            Dictionary with 'atomic' and 'global' keys:
                - 'atomic': Dict mapping atom index (1-based) to atomic properties
                    - area_total, area_positive, area_negative: Surface areas (A^2)
                    - min_value, max_value: Min/max function values
                    - avg_all, avg_positive, avg_negative: Average values
                    - var_all, var_positive, var_negative: Variance values
                - 'global': Dict with global surface statistics
                    - surface_area, volume, average, variance, minimum, maximum
        """
        normalized_model = self._normalize_option(model)
        if normalized_model in self._results.surfaces:
            return self._results.surfaces[normalized_model]

        if normalized_model not in SURFACE_MODELS:
            choices = ", ".join(SURFACE_MODELS.keys())
            raise ValueError(
                f"Surface model {model!r} not supported. Choose between {choices}."
            )
        menu_pattern = SURFACE_MODELS[normalized_model]
        menu_cmd = menu_pattern.split(" ")[0]

        self._results.citations.add(
            "Tian Lu, Feiwu Chen, Quantitative analysis of molecular surface based "
            "on improved Marching Tetrahedra algorithm, J. Mol. Graph. Model., 38, "
            "314-323 (2012)"
        )
        if normalized_model == "esp":
            self._results.citations.add("Phys. Chem. Chem. Phys., 23, 20323 (2021)")

        commands = [
            CommandStep("12", expect="12 Quantitative analysis of molecular surface"),
            CommandStep("2", expect="2 Select mapped function"),
            CommandStep(menu_cmd, expect=menu_pattern),
            CommandStep("0", expect="0 Start analysis"),
            WaitProgress(self._max_wait),
            CommandStep("11", expect="11 Output surface properties"),
            CommandStep("n", expect="surface facets to locsurf"),
            CommandStep("-1", expect="-1 Return to upper level"),
            CommandStep("-1", expect="-1 Return to main menu"),
        ]

        result = self._run_commands(commands, subdir=normalized_model)

        atomic_props = self._parse_atomic_surface_properties(result.stdout)
        global_props = self._parse_global_surface_properties(result.stdout)
        surface_result = {"atomic": atomic_props, "global": global_props}

        self._results.surfaces[normalized_model] = surface_result

        return surface_result

    def get_bond_order(self, model: str) -> dict[tuple[int, int], float]:
        """Calculate bond orders.

        Args:
            model: Bond order model to use.

        Returns:
            Bond order matrix as {(i, j): bond_order}.

        Raises:
            ValueError: If given bond order model is not available.
        """
        normalized_model = self._normalize_option(model)
        if normalized_model in self._results.bond_orders:
            return self._results.bond_orders[normalized_model]

        if normalized_model not in BOND_ORDER_MODELS:
            choices = ", ".join(sorted(BOND_ORDER_MODELS.keys()))
            raise ValueError(
                "Bond order type " f"{model!r} not supported. Choose between {choices}."
            )
        if (
            self._is_wfn_input()
            and normalized_model in WFN_UNSUPPORTED_BOND_ORDER_MODELS
        ):
            raise ValueError(MISSING_BASIS_FUNCTION_INFORMATION_ERROR)

        menu_pattern = BOND_ORDER_MODELS[normalized_model]
        menu_cmd = menu_pattern.split(" ")[0]

        if normalized_model == "laplacian":
            self._results.citations.add(
                "Tian Lu and Feiwu Chen, Bond Order Analysis Based on the Laplacian of "
                "Electron Density in Fuzzy Overlap Space, J. Phys. Chem. A, 117, "
                "3100-3108 (2013)"
            )
        commands = [
            CommandStep("9", expect="9 Bond order analysis"),
            CommandStep(menu_cmd, expect=menu_pattern),
            WaitProgress(self._max_wait),
            CommandStep("n", expect="current folder"),
            CommandStep("0", expect="0 Return"),
        ]

        if normalized_model in {"fuzzy", "laplacian"}:
            commands = [
                *self._commands_fuzzy_partitioning(),
                CommandStep("0", expect="0 Return"),
                *commands,
            ]
        result = self._run_commands(commands, subdir=normalized_model)
        bond_order_matrix = self._parse_bond_orders(result.stdout)

        self._results.bond_orders[normalized_model] = bond_order_matrix

        return bond_order_matrix

    def grid_to_descriptors(
        self, grid_path: str | Path | None, userfunction: int = -1
    ) -> dict[int, float]:
        """Calculate atomic densities from integration in fuzzy atomic spaces.

        Args:
            grid_path: Path to a grid file in .cub or .grd. If None, use the
                object's file_path.
            userfunction: User-defined function index for integration.

        Returns:
            Dictionary mapping atom index (1-based) to integrated density value.

        Raises:
            ValueError: If input file is not a cube/grid file.
            FileNotFoundError: If a provided grid_path does not exist.
        """
        grids = {".cub", ".grd"}
        grid_file = Path(grid_path).resolve() if grid_path is not None else None
        if grid_file is not None:
            if grid_file.suffix.lower() not in grids:
                raise ValueError(
                    "Density analysis requires a .cub or .grd file as input."
                )
            if not grid_file.exists():
                raise FileNotFoundError(f"Grid file not found: {grid_file}")
            active_grid = grid_file
        else:
            if self._file_path.suffix.lower() not in grids:
                raise ValueError(
                    "Density analysis requires a .cub or .grd file as input."
                )
            active_grid = self._file_path

        grid_key = (str(active_grid), userfunction)
        if grid_key in self._results.grid_descriptors:
            return self._results.grid_descriptors[grid_key]

        file_path = self._file_path
        try:
            if active_grid != self._file_path:
                self._file_path = active_grid

            commands = self._build_fuzzy_integration_commands(
                "100",
                "100 User-defined function",
                prefix_commands=self._commands_change_iuserfunc(function=userfunction),
            )
            result = self._run_commands(commands, subdir="grid")
        finally:
            self._file_path = file_path

        grid_descriptors = self._parse_atomic_values(result.stdout)

        self._results.grid_descriptors[grid_key] = grid_descriptors
        return grid_descriptors

    def get_grid(
        self,
        descriptor: str,
        grid_quality: list | str,
        output_path: str | Path,
    ) -> Path:
        """Generate a cube and export it to an explicit, unused destination.

        Args:
            descriptor: Descriptor to calculate on grid.
            grid_quality: Grid quality, one of 'low', 'medium', 'high' or custom commands.
            output_path: Destination for the cube, retained after close().
                Parent directories are created if needed; existing paths are rejected.

        Returns:
            Path to the generated cube file.

        Raises:
            ValueError: If the descriptor is not supported.
            FileNotFoundError: If no cube file was generated.
            FileExistsError: If the output destination already exists.
            RuntimeError: If the generated cube is ambiguous.
            BaseException: If export fails; the incomplete output is removed.
        """
        output_path = Path(output_path).absolute()
        if output_path.exists() or output_path.is_symlink():
            raise FileExistsError(f"Output already exists: {output_path}")
        menu_cmd, menu_pattern, commands = self._get_descriptor_function(descriptor)

        if (
            isinstance(grid_quality, str)
            and self._normalize_option(grid_quality) in GRID_QUALITIES
        ):
            normalized_grid_quality = self._normalize_option(grid_quality)
            grid_pattern = GRID_QUALITIES[normalized_grid_quality]
            grid_cmd = grid_pattern.split(" ")[0]
            grid_quality_commands = [CommandStep(grid_cmd, expect=grid_pattern)]
        elif isinstance(grid_quality, list):
            grid_quality_commands = self._parse_custom_commands(grid_quality)
        else:
            choices = ", ".join(GRID_QUALITIES.keys())
            raise ValueError(
                f"Grid quality {grid_quality!r} not supported. Choose between {choices}."
            )

        grid_commands: list[CommandStep | WaitProgress] = [
            *commands,
            CommandStep(
                "5",
                expect="5 Output and plot specific property within a spatial region",
            ),
            CommandStep(menu_cmd, expect=menu_pattern),
            *grid_quality_commands,
            WaitProgress(self._max_wait),
            CommandStep("2", expect="2 Export data to a Gaussian-type cube file"),
            CommandStep("0", expect="0 Return to main menu"),
        ]
        result = self._run_commands(grid_commands, subdir=descriptor)

        cubes = list(result.workdir.glob("*.cub"))
        if not cubes:
            raise FileNotFoundError("No .cub file was generated by Multiwfn.")
        if len(cubes) != 1:
            raise RuntimeError("Multiwfn produced more than one cube file.")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        target = output_path.open("xb")
        try:
            with target, cubes[0].open("rb") as source:
                shutil.copyfileobj(source, target)
        except BaseException:
            output_path.unlink()
            raise
        return output_path

    def get_fukui(self) -> dict[str, dict[int, float]]:
        """Calculate Fukui functions.

        Returns:
            Dictionary with all Fukui descriptors

        Raises:
            ValueError: If the system is open-shell or spin is undefined.
        """
        if self._has_spin is None:
            raise ValueError(
                "Fukui function analysis requires known spin state. "
                "Initialize Multiwfn with has_spin=True/False."
            )
        if self._has_spin:
            raise ValueError("Fukui function analysis requires closed-shell.")

        if self._results.fukui is not None:
            return self._results.fukui

        commands = [
            CommandStep("22", expect="22 Conceptual DFT"),
            CommandStep("6", expect="6 Calculate condensed OW Fukui"),
            WaitProgress(self._max_wait),
            CommandStep(
                "",
                expect=SINGLE_DETERMINANT_WFN_ERROR_PATTERN,
                optional=True,
            ),
            CommandStep(
                "0",
                expect=r"(Fukui potential|0 Return)",
            ),
        ]
        result = self._run_commands(commands, subdir="fukui")
        self._raise_if_single_determinant_wavefunction_error(
            result.stdout, analysis="Fukui function analysis"
        )
        fukui_results = self._parse_fukui(result.stdout)

        self._results.fukui = fukui_results

        return fukui_results

    def get_superdelocalizabilities(self) -> dict[str, dict[int, float]]:
        """Calculate superdelocalizabilities.

        Returns:
            Dictionary with all superdelocalizability descriptors.

        Raises:
            ValueError: If the system is open-shell or spin is undefined.
            RuntimeError: If Multiwfn reports unsupported wavefunction type or no
                superdelocalizability table can be parsed from output.
        """
        if self._has_spin is None:
            raise ValueError(
                "Superdelocalizability analysis requires known spin state. "
                "Initialize Multiwfn with has_spin=True/False."
            )
        if self._has_spin:
            raise ValueError("Superdelocalizability analysis requires closed-shell.")

        if self._results.superdelocalizabilities is not None:
            return self._results.superdelocalizabilities

        commands = [
            CommandStep("22", expect="22 Conceptual DFT"),
            CommandStep("8", expect="8 Calculate nucleophilic and electrophilic"),
            WaitProgress(self._max_wait),
            CommandStep(
                "",
                expect=SINGLE_DETERMINANT_WFN_ERROR_PATTERN,
                optional=True,
            ),
            # Prompt text after option 8 varies across Multiwfn builds; send "0"
            # unconditionally to return to the previous menu.
            CommandStep("0"),
        ]
        result = self._run_commands(commands, subdir="superdeloc")
        self._raise_if_single_determinant_wavefunction_error(
            result.stdout, analysis="Superdelocalizability analysis"
        )
        superdeloc_matrix = self._parse_superdelocalizabilities(result.stdout)
        if not superdeloc_matrix["d_n"]:
            nonempty_lines = [
                line.strip() for line in result.stdout.splitlines() if line.strip()
            ]
            excerpt = "\n".join(nonempty_lines[-20:])
            if excerpt:
                excerpt = f"\nOutput excerpt:\n{excerpt}"
            raise RuntimeError(
                f"{SINGLE_DETERMINANT_WFN_ERROR_TEXT}. "
                "Superdelocalizability analysis returned no atomic values; "
                "Multiwfn output format may be unsupported."
                f"{excerpt}"
            )

        self._results.superdelocalizabilities = superdeloc_matrix

        return superdeloc_matrix

    def get_vector(self, descriptor: str) -> dict[int, tuple[float, float, float]]:
        """Calculate atomic 3D descriptors.

        Args:
            descriptor: Descriptor name

        Returns:
            Dictionary of atomic descriptors as 3D vectors.

        Raises:
            ValueError: for invalid descriptor input.
        """
        if descriptor not in VECTORS:
            raise ValueError(f"Vector descriptor {descriptor!r} not supported")

        if descriptor in self._results.atomic_vector_descriptors:
            return self._results.atomic_vector_descriptors[descriptor]

        menu_cmd = "100"
        menu_pattern = "100 User-defined function"
        functions = VECTORS[descriptor]
        if len(functions) == 2 and isinstance(functions[1], str):
            functions, citation = functions
            self._results.citations.add(citation)

        vectors: dict[int, list[float]] = {}
        for component_idx, (func, axis) in enumerate(zip(functions, ("x", "y", "z"))):
            self._check_function_selection(func, descriptor)
            label = f"{descriptor}_{axis}"
            result = self._run_fuzzy_integration(
                menu_cmd,
                menu_pattern,
                self._commands_change_iuserfunc(function=func),
                label,
            )
            for atom_idx, value in result.items():
                components = vectors.setdefault(atom_idx, [0.0, 0.0, 0.0])
                components[component_idx] = value

        vectors_3d = {
            atom_idx: (components[0], components[1], components[2])
            for atom_idx, components in vectors.items()
        }

        self._results.atomic_vector_descriptors[descriptor] = vectors_3d
        return vectors_3d

    def _check_function_selection(self, func: int, descriptor: str) -> None:
        closed_shell_only = {24, 34, 60, 61, 62, 1100}
        open_shell_only = {1, 2, 1101, 1102}

        if self._has_spin is None:
            if func in closed_shell_only or func in open_shell_only:
                raise ValueError(
                    f"Descriptor {descriptor!r} requires known spin state. "
                    "Initialize Multiwfn with has_spin=True/False."
                )
            return

        if self._has_spin:  # open-shell
            if func in closed_shell_only:
                raise ValueError(
                    f"Descriptor {descriptor!r} not available for open-shell systems."
                )
        else:  # closed-shell
            if func in open_shell_only:
                raise ValueError(
                    f"Descriptor {descriptor!r} not available for closed-shell systems."
                )

    def _get_descriptor_function(
        self, descriptor: str
    ) -> tuple[str, str, list[CommandStep]]:
        if descriptor not in ALL_FUNCTIONS:
            raise ValueError(f"Descriptor {descriptor!r} not supported.")
        if descriptor in VECTORS:
            raise ValueError(
                f"Descriptor {descriptor!r} is vector-valued. Use get_vector() instead."
            )

        commands: list[CommandStep] = []

        if descriptor in REAL_SPACE_FUNCTIONS:
            menu_pattern = REAL_SPACE_FUNCTIONS[descriptor]
            menu_cmd = menu_pattern.split(" ")[0]

        else:
            function_setting = USER_REAL_FUNCTIONS[descriptor]
            func, citation = (
                function_setting
                if isinstance(function_setting, tuple)
                else (function_setting, None)
            )
            if citation is not None:
                self._results.citations.add(citation)
            menu_cmd = "100"
            menu_pattern = "100 User-defined function"

            self._check_function_selection(func, descriptor)
            commands.extend(self._commands_change_iuserfunc(function=func))

        return menu_cmd, menu_pattern, commands

    def get_descriptor(
        self, descriptor: str, partitioning: str = "becke"
    ) -> dict[int, float]:
        """Calculate atomic densities from integration in fuzzy atomic spaces.

        Args:
            descriptor: Descriptor name.
            partitioning: Atomic space partitioning scheme to use. Supported values
                are 'becke', 'hirshfeld', and 'mbis'.

        Returns:
            Dictionary mapping atom index (1-based) to integrated density value.

        Raises:
            ValueError: If the descriptor is unsupported, vector-valued, requires grid
                data, is incompatible with the current spin state, or the
                partitioning scheme is unsupported.
        """
        if descriptor == "rdg":
            raise ValueError(f"Descriptor {descriptor!r} requires grid data!")
        normalized_partitioning = self._normalize_option(partitioning)
        if normalized_partitioning not in FUZZY_PARTITIONING_SCHEMES:
            choices = ", ".join(FUZZY_PARTITIONING_SCHEMES.keys())
            raise ValueError(
                f"Partitioning scheme {partitioning!r} not supported. "
                f"Choose between {choices}."
            )
        cache_key = (
            descriptor
            if normalized_partitioning == "becke"
            else f"{descriptor}:{normalized_partitioning}"
        )
        if cache_key in self._results.atomic_descriptors:
            return self._results.atomic_descriptors[cache_key]

        menu_cmd, menu_pattern, commands = self._get_descriptor_function(descriptor)
        descriptor_result = self._run_fuzzy_integration(
            menu_cmd,
            menu_pattern,
            commands,
            descriptor,
            normalized_partitioning,
        )
        self._results.atomic_descriptors[cache_key] = descriptor_result
        return descriptor_result

    def _run_fuzzy_integration(
        self,
        menu_cmd: str,
        menu_pattern: str,
        commands: list[CommandStep],
        descriptor: str,
        partitioning: str = "becke",
    ) -> dict[int, float]:
        command_sequence = self._build_fuzzy_integration_commands(
            menu_cmd,
            menu_pattern,
            prefix_commands=commands,
            partitioning=partitioning,
        )
        subdir = (
            descriptor if partitioning == "becke" else f"{descriptor}_{partitioning}"
        )
        result = self._run_commands(command_sequence, subdir=subdir)
        descriptor_result = self._parse_atomic_values(result.stdout)
        return descriptor_result

    def get_electric_moments(self) -> dict[str, float]:
        """Calculate electric dipole/multipole moments."""
        if self._results.electric_moments is not None:
            return self._results.electric_moments

        commands = [
            CommandStep("300", expect="300 Other functions"),
            CommandStep("5", expect="5 Calculate electric dipole/multipole moments"),
            CommandStep("0", expect="0 Return"),
        ]
        result = self._run_commands(commands, subdir="electric_moments")
        moments = self._parse_electric_moments(result.stdout)

        self._results.electric_moments = moments

        return moments

    def get_orbitals(self) -> list[dict[str, Any]]:
        """Return orbital records in their reported order.

        Energies are read directly from Multiwfn's eV column.
        Spin labels A+B, A, and B map to restricted, alpha, and beta, respectively.

        Returns:
            Records with original indices, spin labels, occupations and energies in eV.
        """
        if self._results.orbitals is None:
            commands = [
                CommandStep("6", expect="6 Check & modify wavefunction"),
                CommandStep("3", expect="3 List all orbitals"),
                CommandStep("-1", expect="-1 Return"),
            ]
            result = self._run_commands(commands, subdir="orbitals")
            self._results.orbitals = self._parse_orbital_list(result.stdout)
        return self._results.orbitals

    def get_frontier_gaps(
        self,
        occupation_tolerance: float = 1e-6,
        energy_tolerance_ev: float = 1e-5,
    ) -> dict[str, dict[str, Any]]:
        """Compute same-spin frontier orbital-energy differences.

        Restricted occupations must be near 0 or 2; alpha/beta occupations near
        0 or 1. Fractional and restricted open-shell occupations are unsupported.
        Physical interpretation requires appropriate source orbital energies;
        these differences are not excitation energies.

        Args:
            occupation_tolerance: Absolute occupation tolerance, in [0, 0.5).
            energy_tolerance_ev: Nonnegative absolute tolerance for reporting
                degenerate frontier indices. Does not change the frontier energies.

        Returns:
            Per-channel frontier energies, degenerate index lists and signed gap
            in eV. Keys are 'restricted', or 'alpha' and 'beta'.

        Raises:
            ValueError: If tolerances or orbital data are invalid, or occupations
                and spin channels do not support automatic classification.
        """
        if (
            not math.isfinite(occupation_tolerance)
            or not 0 <= occupation_tolerance < 0.5
        ):
            raise ValueError("occupation_tolerance must be finite and in [0, 0.5).")
        if not math.isfinite(energy_tolerance_ev) or energy_tolerance_ev < 0:
            raise ValueError("energy_tolerance_ev must be finite and nonnegative.")
        orbitals = self.get_orbitals()
        spins = {orbital["spin"] for orbital in orbitals}
        if spins not in ({"restricted"}, {"alpha", "beta"}):
            raise ValueError("Orbital spin channels are empty, mixed, or incomplete.")
        if any(
            not math.isfinite(orbital[key])
            for orbital in orbitals
            for key in ("energy_ev", "occupation")
        ):
            raise ValueError("Orbital energies and occupations must be finite.")
        identities = {(orbital["spin"], orbital["index"]) for orbital in orbitals}
        if len(identities) != len(orbitals) or any(
            orbital["index"] < 1 for orbital in orbitals
        ):
            raise ValueError("Invalid or duplicate orbital identity.")
        return {
            spin: self._frontier_gap(
                [orbital for orbital in orbitals if orbital["spin"] == spin],
                occupation_tolerance,
                energy_tolerance_ev,
            )
            for spin in sorted(spins)
        }

    @staticmethod
    def _frontier_gap(
        orbitals: list[dict[str, Any]],
        occupation_tolerance: float,
        energy_tolerance_ev: float,
    ) -> dict[str, Any]:
        """Classify occupations and compute frontier extrema within one channel."""
        spin = orbitals[0]["spin"]
        full_occupation = 2 if spin == "restricted" else 1
        occupied = []
        unoccupied = []
        for orbital in orbitals:
            occupation = orbital["occupation"]
            if (
                occupation < -occupation_tolerance
                or occupation > full_occupation + occupation_tolerance
            ):
                raise ValueError(
                    f"Occupation outside the valid range for {spin} orbitals."
                )
            if abs(occupation) <= occupation_tolerance:
                unoccupied.append(orbital)
            elif abs(occupation - full_occupation) <= occupation_tolerance:
                occupied.append(orbital)
            elif spin == "restricted" and abs(occupation - 1) <= occupation_tolerance:
                raise ValueError(
                    "Restricted open-shell frontier classification is unsupported."
                )
            else:
                raise ValueError(
                    "Fractional occupations detected; use get_orbitals() for explicit analysis."
                )
        if not occupied or not unoccupied:
            raise ValueError(
                f"Missing occupied or unoccupied orbitals in the {spin} channel."
            )
        homo = max(orbital["energy_ev"] for orbital in occupied)
        lumo = min(orbital["energy_ev"] for orbital in unoccupied)
        return {
            "homo_energy_ev": homo,
            "homo_indices": sorted(
                orbital["index"]
                for orbital in occupied
                if abs(orbital["energy_ev"] - homo) <= energy_tolerance_ev
            ),
            "lumo_energy_ev": lumo,
            "lumo_indices": sorted(
                orbital["index"]
                for orbital in unoccupied
                if abs(orbital["energy_ev"] - lumo) <= energy_tolerance_ev
            ),
            "gap_ev": lumo - homo,
        }

    def get_localization(self) -> dict[int, float]:
        """Calculate localization indices using Becke atomic spaces."""
        if self._results.localization is not None:
            return self._results.localization

        commands = [
            *self._commands_fuzzy_partitioning(),
            CommandStep("4", expect="4 Calculate localization"),
            WaitProgress(self._max_wait),
            CommandStep("n", expect="y/n"),
            CommandStep("0", expect="0 Return"),
        ]
        result = self._run_commands(commands, subdir="localization")
        localization = self._parse_localization_values(result.stdout)

        self._results.localization = localization
        return localization

    def _table_lines(
        self, stdout: str, headers: tuple[str, ...], end: str | None = None
    ) -> Iterable[str]:
        """Yield lines after a matching header and before an optional end marker."""
        lines = iter(stdout.splitlines())
        for line in lines:
            if all(re.search(pattern, line) for pattern in headers):
                break
        for line in lines:
            if end is not None and re.search(end, line):
                break
            yield line

    def _raise_if_single_determinant_wavefunction_error(
        self, stdout: str, analysis: str
    ) -> None:
        """Raise a clear error for unsupported multi-determinant wavefunctions."""
        if SINGLE_DETERMINANT_WFN_ERROR_RE.search(stdout):
            raise RuntimeError(
                f"{analysis} failed: {SINGLE_DETERMINANT_WFN_ERROR_TEXT}."
            )

    def _parse_fukui(self, stdout: str) -> dict[str, dict[int, float]]:
        """Parse Fukui function descriptors from stdout."""
        keys = ("f_plus", "f_minus", "f_zero", "dd")
        result: dict[str, dict[int, float]] = {key: {} for key in keys}
        pattern = (
            rf"\s*(\d+)\([A-Za-z][a-z]?\s*\)?\s+({NUMBER_PATTERN})\s+"
            rf"({NUMBER_PATTERN})\s+({NUMBER_PATTERN})\s+({NUMBER_PATTERN})"
        )
        for line in self._table_lines(stdout, ("Atom index", r"OW f\+")):
            match = re.search(pattern, line)
            if match:
                index = int(match.group(1))
                for key, value in zip(keys, match.groups()[1:]):
                    result[key][index] = float(value)
        return result

    def _parse_superdelocalizabilities(
        self, stdout: str
    ) -> dict[str, dict[int, float]]:
        """Parse two- or four-column tables with D_N/D_E or DN/DE headers."""
        keys = ("d_n", "d_e", "d_n_0", "d_e_0")
        result: dict[str, dict[int, float]] = {key: {} for key in keys}
        headers = (
            r"(?i)\bAtom(?:\s+index)?\b",
            r"(?i)\bD[_\s]?N\b",
            r"(?i)\bD[_\s]?E\b",
        )
        for line in self._table_lines(stdout, headers, r"(?i)sum of"):
            match = re.search(r"\s*(\d+)\([A-Za-z][a-z]?\s*\)?\s+(.*)$", line)
            if match:
                index = int(match.group(1))
                values = re.findall(NUMBER_PATTERN, match.group(2))
                count = 4 if len(values) >= 4 else 2
                if len(values) >= count:
                    for key, value in zip(keys[:count], values):
                        result[key][index] = float(value)
        return result

    def _parse_electric_moments(self, stdout: str) -> dict[str, float]:
        """Parse electric moments from stdout."""
        patterns = {
            "dipole_magnitude_au": (
                rf"Magnitude of dipole moment:\s*({NUMBER_PATTERN})\s+a\.u\."
            ),
            "quadrupole_traceless_magnitude": (
                rf"Magnitude of the traceless quadrupole moment tensor:\s*({NUMBER_PATTERN})"
            ),
            "quadrupole_spherical_magnitude": (
                rf"Magnitude: \|Q_2\|=\s*({NUMBER_PATTERN})"
            ),
            "octopole_spherical_magnitude": (
                rf"Magnitude: \|Q_3\|=\s*({NUMBER_PATTERN})"
            ),
            "electronic_spatial_extent": (
                rf"Electronic spatial extent <r\^2>:\s*({NUMBER_PATTERN})"
            ),
        }

        moments: dict[str, float] = {}
        for key, pattern in patterns.items():
            match = re.search(pattern, stdout)
            if match:
                moments[key] = self._parse_float(match.group(1))

        return moments

    def _parse_orbital_list(self, stdout: str) -> list[dict[str, Any]]:
        """Parse orbital rows without silently discarding malformed records."""
        spins = {"A+B": "restricted", "A": "alpha", "B": "beta"}
        orbitals = []
        for line in stdout.splitlines():
            if not re.search(r"(?:Ene|E)\(au/eV\):", line):
                continue
            match = ORBITAL_LINE_PATTERN.match(line)
            if match is None or match["orbital_type"] not in spins:
                raise RuntimeError(f"Could not parse orbital row: {line!r}.")
            orbitals.append(
                {
                    "index": int(match["index"]),
                    "spin": spins[match["orbital_type"]],
                    "occupation": self._parse_float(match["occupation"]),
                    "energy_ev": self._parse_float(match["energy_ev"]),
                }
            )
        if not orbitals:
            raise RuntimeError("Could not parse orbital listing from stdout.")
        return orbitals

    def _parse_chg_file(self, chg_path: Path) -> dict[int, float]:
        """Parse .chg file to {atom_idx: charge}."""
        charges = {}
        with chg_path.open("r", encoding="utf-8", errors="replace") as f:
            for idx, line in enumerate(f, start=1):
                parts = line.split()
                if len(parts) != 5:
                    continue
                try:
                    charges[idx] = float(parts[4])
                except ValueError:
                    continue
        return charges

    def _parse_atomic_values(self, stdout: str) -> dict[int, float]:
        """Parse fuzzy atomic space integration values from stdout."""
        atomic_values: dict[int, float] = {}
        pattern = rf"\s*(\d+)\([A-Z][a-z]?\s*\)\s+({NUMBER_PATTERN})"
        for line in self._table_lines(stdout, ("Atomic space", "Value"), "Summing up"):
            match = re.search(pattern, line)
            if match:
                atomic_values[int(match.group(1))] = float(match.group(2))
        if not atomic_values:
            raise RuntimeError("Could not parse the Multiwfn atomic integration table.")
        return atomic_values

    def _parse_localization_values(self, stdout: str) -> dict[int, float]:
        """Parse wrapped localization index output from stdout."""
        localization: dict[int, float] = {}
        in_section = False
        atom_value_pattern = re.compile(
            r"(\d+)\([A-Z][a-z]?\s*\)\s*:\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)"
        )

        for line in stdout.splitlines():
            if not in_section:
                if "Localization index" in line:
                    in_section = True
                continue

            if not line.strip():
                if localization:
                    break
                continue

            matches = list(atom_value_pattern.finditer(line))
            if not matches:
                if localization:
                    break
                continue

            for match in matches:
                atom_idx = int(match.group(1))
                localization[atom_idx] = self._parse_float(match.group(2))

        return localization

    def _parse_bond_orders(self, stdout: str) -> dict[tuple[int, int], float]:
        """Parse bond order matrix from Multiwfn output."""
        bond_orders: dict[tuple[int, int], float] = {}
        atom_pattern = re.compile(r"(\d+)\([A-Z][a-z]?\s*\)")
        float_pattern = re.compile(NUMBER_PATTERN)

        for line in stdout.splitlines():
            atom_matches = list(atom_pattern.finditer(line))
            if len(atom_matches) < 2:
                continue
            atom_i = int(atom_matches[0].group(1))
            atom_j = int(atom_matches[1].group(1))
            float_matches = float_pattern.findall(line)
            if not float_matches:
                continue
            bond_orders[(atom_i, atom_j)] = float(float_matches[-1])

        if bond_orders:
            return bond_orders

        for line in stdout.splitlines():
            if "(" not in line or ")" not in line:
                continue
            tokens = line.replace(":", " ").split()
            atom_tokens = [t for t in tokens if "(" in t and ")" in t]
            if len(atom_tokens) < 2:
                continue
            try:
                atom_i = int(atom_tokens[0].split("(")[0])
                atom_j = int(atom_tokens[1].split("(")[0])
                bond_orders[(atom_i, atom_j)] = float(tokens[-1])
            except (ValueError, IndexError):
                continue

        return bond_orders

    def _parse_atomic_surface_properties(
        self, stdout: str
    ) -> dict[int, dict[str, float]]:
        """Parse the separate ESP tables or the single atomic surface table."""
        if "All/Positive/Negative area" in stdout:
            sections = [
                (
                    "All/Positive/Negative area",
                    (
                        "area_total",
                        "area_positive",
                        "area_negative",
                        "min_value",
                        "max_value",
                    ),
                    "All/Positive/Negative average",
                    True,
                ),
                (
                    "All/Positive/Negative average",
                    (
                        "avg_all",
                        "avg_positive",
                        "avg_negative",
                        "var_all",
                        "var_positive",
                        "var_negative",
                    ),
                    r"^\s*$",
                    True,
                ),
                (r"Atom#\s+Pi", ("pi", "nu", "nu_sigma2"), r"^\s*$", False),
            ]
        else:
            sections = [
                (
                    r"Atom#\s+Area\(Ang\^2\)",
                    ("area_total", "min_value", "max_value", "avg_all", "var_all"),
                    r"^\s*$",
                    False,
                )
            ]
        atomic_props: dict[int, dict[str, float]] = {}
        for header, keys, end, skip_invalid in sections:
            for line in self._table_lines(stdout, (header,), end):
                parts = line.split()
                if len(parts) < len(keys) + 1 or not parts[0].isdigit():
                    continue
                index = int(parts[0])
                if keys[0] != "area_total" and index not in atomic_props:
                    continue
                try:
                    values = {key: float(value) for key, value in zip(keys, parts[1:])}
                except ValueError:
                    if skip_invalid:
                        continue
                    raise
                atomic_props.setdefault(index, {}).update(values)
        return atomic_props

    def _parse_global_surface_properties(self, stdout: str) -> dict[str, float]:
        """Parse global surface properties from stdout."""
        props: dict[str, float] = {}
        for line in self._table_lines(
            stdout,
            ("================= Summary of surface analysis",),
            "Surface analysis finished!",
        ):
            if ":" not in line:
                continue
            if not self._parse_summary_minmax(line, props):
                self._parse_summary_key_value(line, props)
        return props

    def _parse_summary_minmax(self, line: str, props: dict[str, float]) -> bool:
        """Parse a summary line containing both minimal and maximal values."""
        if "Minimal value:" not in line or "Maximal value:" not in line:
            return False
        for label, value in re.findall(
            r"(Minimal value|Maximal value)\s*:\s*([-\d.+Ee]+)", line
        ):
            key = re.sub(r"\s+", "_", label.strip().lower())
            props[key] = self._parse_float(value)
        return True

    def _parse_summary_key_value(self, line: str, props: dict[str, float]) -> None:
        """Parse a single summary key/value line."""
        key_raw, rest = line.split(":", 1)
        key = re.sub(r"[()\[\]=]", "", key_raw.strip().lower())
        key = re.sub(r"\s+", "_", key)
        value = self._parse_first_number(rest)
        if value is not None:
            props[key] = value
        elif re.search(r"\bnan\b", rest, flags=re.IGNORECASE):
            props[key] = float("nan")

    @staticmethod
    def _parse_first_number(text: str) -> float | None:
        match = re.search(NUMBER_PATTERN, text)
        if match:
            return float(match.group(0))
        return None

    @staticmethod
    def _parse_float(token: str) -> float:
        if token.lower() == "nan":
            return float("nan")
        return float(token)

    def parse_surfanalysis(self, filepath: Path) -> dict[str, Any]:
        """Parse surfanalysis.txt for global surface statistics."""
        data: dict[str, Any] = {"raw": "", "statistics": {}}

        if not filepath.exists():
            return data

        with filepath.open("r", encoding="utf-8", errors="replace") as f:
            content = f.read()
            data["raw"] = content

            minima_match = re.search(r"Number of surface minima:\s*(\d+)", content)
            maxima_match = re.search(r"Number of surface maxima:\s*(\d+)", content)

            if minima_match:
                data["statistics"]["num_minima"] = int(minima_match.group(1))
            if maxima_match:
                data["statistics"]["num_maxima"] = int(maxima_match.group(1))

            extrema_pattern = re.compile(
                (
                    rf"^\s*\*?\s*\d+\s+({NUMBER_PATTERN})\s+{NUMBER_PATTERN}\s+"
                    rf"{NUMBER_PATTERN}\b"
                ),
                flags=re.MULTILINE,
            )
            matches = extrema_pattern.findall(content)

            if matches:
                extrema_values = [float(m) for m in matches]
                data["statistics"]["num_extrema"] = len(extrema_values)
                data["statistics"]["extrema_min"] = min(extrema_values)
                data["statistics"]["extrema_max"] = max(extrema_values)
                data["statistics"]["extrema_average"] = sum(extrema_values) / len(
                    extrema_values
                )

        return data


def cli(file: str) -> Any:
    """CLI entry point for Multiwfn."""
    partial_func = functools.partial(Multiwfn, file)
    return partial_func


def molden2aim(
    path: str,
    molden_name: str = "xtb.molden",
    molden2aim_executable_path: str = "molden2aim.exe",
) -> str:
    """Normalize molden file using molden2aim.

    Args:
        path: Path to working directory containing the molden2aim executable.
        molden_name: Name of the molden file to convert.
        molden2aim_executable_path: Name of the molden2aim executable.

    Returns:
        Name of the converted molden file.
    """
    session = _PexpectSession(
        base_command=f"./{molden2aim_executable_path} -i {molden_name}",
        workdir=Path(path),
        debug=False,
        env=build_execution_env(),
    )
    commands = ["Yes", "No", "No", "No"]
    for item in commands:
        session.send(item)
    session.wait_for_exit()

    return molden_name.split(".molden")[0] + "_new.molden"
