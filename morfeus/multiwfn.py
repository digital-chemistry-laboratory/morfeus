"""Multiwfn (http://sobereva.com/multiwfn/) interface code.

Uses pexpect for interactive control.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import functools
from pathlib import Path
import re
import shutil
import tempfile
import time
from types import TracebackType
import typing
from typing import Any, cast, Iterable

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
HOMO_LINE_PATTERN = re.compile(
    rf"Orbital\s+(\d+)\s+is\s+HOMO,.*?({NUMBER_PATTERN})\s+eV\b"
)
LUMO_LINE_PATTERN = re.compile(
    rf"Orbital\s+(\d+)\s+is\s+LUMO,.*?({NUMBER_PATTERN})\s+eV\b"
)
HOMO_LUMO_GAP_PATTERN = re.compile(rf"HOMO-LUMO gap:.*?({NUMBER_PATTERN})\s+eV\b")
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
    \s+Typ(?:e)?\s*:\s*(?P<orbital_type>A\+B|A|B)\b
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
    "low": "1 Low quality grid",
    "medium": "2 Medium quality grid",
    "high": "3 High quality grid",
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
class ProgressState:
    """Track progress bar waiting state."""

    last_activity: float
    saw_progress: bool
    saw_please_wait: bool
    idle_timeout: float


@dataclass
class MultiwfnRunResult:
    """Result from a single Multiwfn run."""

    stdout: str
    workdir: Path


@dataclass
class MultiwfnResults:
    """Cached Multiwfn results."""

    charges: dict[str, dict[int, float]] | None = None
    surfaces: dict[str, dict[str, Any]] | None = None
    atomic_descriptors: dict[str, dict[int, float]] | None = None
    grid_descriptors: dict[str, dict[int, float]] | None = None
    electric_moments: dict[str, float] | None = None
    gap: dict[str, dict[str, Any]] | None = None
    bond_orders: dict[str, dict[tuple[int, int], float]] | None = None
    fukui: dict[str, dict[int, float]] | None = None
    superdelocalizabilities: dict[str, dict[int, float]] | None = None
    localization: dict[int, float] | None = None
    atomic_vector_descriptors: (
        dict[str, dict[int, tuple[float, float, float]]] | None
    ) = None
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
    ) -> None:
        self._child = pexpect.spawn(
            base_command,
            cwd=str(workdir),
            encoding="utf-8",
            env=env,
        )
        self._timeout = 10
        self._expect_timeout = 3
        self._debug = debug
        self._transcript: list[str] = []
        self._last_command_pos = 0

    @staticmethod
    def _coerce_output_chunk(chunk: Any) -> str:
        """Normalize pexpect output chunks to text for transcript storage."""
        if chunk is None:
            return ""
        if isinstance(chunk, str):
            return chunk
        return str(chunk)

    def _append_transcript(self, chunk: Any) -> None:
        """Append normalized output chunk to transcript."""
        self._transcript.append(self._coerce_output_chunk(chunk))

    def read_available(self) -> None:
        """Read available output without blocking."""
        try:
            chunk = self._child.read_nonblocking(size=4096, timeout=0.1)
            if chunk:
                self._append_transcript(chunk)
        except (pexpect.TIMEOUT, pexpect.EOF):
            pass

    def wait_for_progress(self, max_wait: float = 30.0) -> None:
        """Wait for a progress bar to appear and complete."""
        state = self._init_progress_state(max_wait)
        if self._debug:
            print("[DEBUG] Waiting for progress bar " f"(max_idle_wait={max_wait}s)")
        while True:
            chunk = self._read_progress_chunk()
            if chunk is None:
                if self._should_stop_on_timeout(state):
                    return
                continue
            if chunk == "EOF":
                return
            if chunk == "":
                continue
            if self._handle_progress_chunk(state, chunk):
                return

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

    def expect(self, pattern: str) -> bool:
        """Wait for pattern. Returns True on match, False on timeout."""
        if self._debug:
            print(f"[DEBUG] Expecting: {pattern!r}")

        try:
            self._child.expect(pattern, timeout=self._expect_timeout)
            self._append_transcript(self._child.before)
            self._append_transcript(self._child.after)
            if self._debug:
                print("[DEBUG] Got match")
            return True
        except pexpect.TIMEOUT:
            self.read_available()
            recent_output = self.get_output_since_last_command()
            if self._debug:
                print(
                    f"[WARN] Timeout ({self._expect_timeout}s) waiting for: {pattern!r}"
                )
                print(f"Recent output: {recent_output[-500:]}")
            return False
        except pexpect.EOF:
            self._append_transcript(self._child.before)
            self._append_transcript(self._child.after)
            if self._debug:
                print(f"[WARN] EOF while waiting for: {pattern!r}")
            return False

    def try_expect(self, pattern: str, timeout: float = 1.0) -> bool:
        """Try to match a pattern with short timeout. Silent on timeout.

        Does not extend timeout for progress bars (use expect() for that).

        Args:
            pattern: Regular expression to match.
            timeout: Timeout in seconds.

        Returns:
            True if the pattern is matched.
        """
        if self._debug:
            print(f"[DEBUG] Try expecting: {pattern!r} (timeout={timeout}s)")
        try:
            self._child.expect(pattern, timeout=timeout)
            self._append_transcript(self._child.before)
            self._append_transcript(self._child.after)
            if self._debug:
                print("[DEBUG] Got match")
            return True
        except pexpect.TIMEOUT:
            self.read_available()
            if self._debug:
                print("[DEBUG] Pattern not found (optional)")
            return False
        except pexpect.EOF:
            before = self._coerce_output_chunk(self._child.before)
            after = self._coerce_output_chunk(self._child.after)
            self._append_transcript(before)
            self._append_transcript(after)
            if self._debug:
                print("[DEBUG] EOF while checking optional pattern")
            return re.search(pattern, f"{before}{after}") is not None

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

    def get_recent_output(self, n_chunks: int = 3, max_chars: int = 1000) -> str:
        """Get recent output for debugging."""
        recent = (
            self._transcript[-n_chunks:]
            if len(self._transcript) >= n_chunks
            else self._transcript
        )
        output = "".join(recent)
        return output[-max_chars:] if len(output) > max_chars else output

    def get_output_since_last_command(self) -> str:
        """Get output since last command."""
        if self._last_command_pos < len(self._transcript):
            return "".join(self._transcript[self._last_command_pos :])
        return ""

    @property
    def stdout(self) -> str:
        """Full session transcript."""
        return "".join(self._transcript)

    def _execute_command(self, cmd: str, expect: str | None, use_robust: bool) -> None:
        """Execute command with optional expect."""
        if use_robust and expect:
            pattern = expect
            matched = self.expect(pattern)
            if not matched:
                # Some Multiwfn builds emit the single-determinant warning and then
                # change prompts abruptly. Surface this domain error instead of a
                # generic prompt-mismatch RuntimeError.
                transcript = self.stdout
                buffered = self.get_output_since_last_command()
                if SINGLE_DETERMINANT_WFN_ERROR_RE.search(buffered) or (
                    transcript and SINGLE_DETERMINANT_WFN_ERROR_RE.search(transcript)
                ):
                    raise RuntimeError(SINGLE_DETERMINANT_WFN_ERROR_TEXT)
                raise RuntimeError(
                    f"Expected pattern not found before command {cmd!r}: {pattern!r}"
                )
        self.send(cmd)
        if not use_robust:
            self.read_available()

    def _has_pattern_in_buffered_output(
        self, pattern: str, debug_message: str | None = None
    ) -> bool:
        """Return True if pattern is present in buffered output."""
        buffered_output = self.get_output_since_last_command()
        pattern_found = re.search(pattern, buffered_output) is not None
        if self._debug and pattern_found and debug_message:
            print(debug_message)
        return pattern_found

    def _resolve_step_execution(
        self, item: CommandStep, use_robust: bool
    ) -> tuple[bool, str | None]:
        """Resolve whether a command step should run and which expect to use."""
        if not item.expect:
            return True, None

        if item.optional:
            pattern_found = self.try_expect(item.expect, timeout=1.0)
            if not pattern_found:
                pattern_found = self._has_pattern_in_buffered_output(
                    item.expect,
                    "[DEBUG] Optional pattern found in buffered output",
                )
            if self._debug:
                status = "found" if pattern_found else "not found"
                print(f"[DEBUG] Optional pattern {status}: {item.expect!r}")
            if not pattern_found:
                return False, None
            # Avoid matching the same already-consumed pattern twice.
            return True, None

        if use_robust and self._has_pattern_in_buffered_output(
            item.expect, "[DEBUG] Required pattern found in buffered output"
        ):
            # Avoid expecting the same prompt again if it was already consumed.
            return True, None
        return True, item.expect

    def _process_command(
        self, item: str | CommandStep | WaitProgress, use_robust: bool
    ) -> None:
        """Process command step."""
        if isinstance(item, str):
            self._execute_command(item, None, use_robust)
            return
        if isinstance(item, WaitProgress):
            if self._debug:
                print(f"[DEBUG] Executing WaitProgress (max_wait={item.max_wait}s)")
            self.wait_for_progress(max_wait=item.max_wait)
            return

        should_execute, command_expect = self._resolve_step_execution(item, use_robust)
        if not should_execute:
            return

        if item.cmd is not None:
            self._execute_command(item.cmd, command_expect, use_robust)

    def _init_progress_state(self, max_wait: float) -> ProgressState:
        """Initialize progress waiting state.

        Args:
            max_wait: Maximum idle time in seconds before timing out.

        Returns:
            Initialized progress waiting state.
        """
        now = time.monotonic()
        return ProgressState(
            last_activity=now,
            saw_progress=False,
            saw_please_wait=False,
            idle_timeout=max_wait,
        )

    def _read_progress_chunk(self) -> str | None:
        """Read output chunk while waiting for progress.

        Returns:
            Chunk string, "EOF" on end-of-file, or None on timeout.
        """
        try:
            chunk = self._child.read_nonblocking(size=8192, timeout=0.5)
            return cast(str, chunk)
        except pexpect.TIMEOUT:
            if self._debug:
                print("[DEBUG] No new output while waiting for progress")
            return None
        except pexpect.EOF:
            if self._debug:
                print("[DEBUG] EOF reached while waiting for progress")
            return "EOF"

    def _should_stop_on_timeout(self, state: ProgressState) -> bool:
        """Check if progress waiting should stop due to inactivity.

        Args:
            state: Current progress waiting state.

        Returns:
            True if waiting should stop.
        """
        timeout = (
            state.idle_timeout if state.saw_progress else min(1.0, state.idle_timeout)
        )
        if (time.monotonic() - state.last_activity) >= timeout:
            if self._debug:
                if state.saw_progress:
                    print("[DEBUG] No progress updates within idle timeout; done")
                else:
                    print("[DEBUG] No progress observed within timeout; done")
            return True
        return False

    def _handle_progress_chunk(self, state: ProgressState, chunk: str) -> bool:
        """Handle a progress chunk and return True when waiting is done.

        Args:
            state: Current progress waiting state.
            chunk: Output chunk from the process.

        Returns:
            True when waiting should stop.
        """
        self._append_transcript(chunk)
        if self._debug:
            print(f"[DEBUG] Read chunk while waiting: {chunk!r}")
        state.last_activity = time.monotonic()
        if "please wait" in chunk.lower():
            state.saw_please_wait = True
            if self._debug:
                print("[DEBUG] 'Please wait' seen while waiting")
        has_progress, current_progress = self._has_progress_bar(chunk)
        if not has_progress:
            if state.saw_progress:
                if self._debug:
                    print("[DEBUG] Progress bar no longer detected; done")
                return True
            return False
        state.saw_progress = True
        if self._debug:
            progress_info = (
                f" ({current_progress}%)" if current_progress is not None else ""
            )
            print(f"[DEBUG] Progress bar detected{progress_info}")
        if current_progress is not None and current_progress >= 100.0:
            if self._debug:
                print("[DEBUG] Progress bar completed (>=100%)")
            return True
        return False


class Multiwfn:
    """Multiwfn interface using pexpect for interactive control.

    Args:
        file_path: Path to molden/wfn file from xtb or other QM program.
        run_path: Directory for output files. Defaults to temporary directory
            that is cleaned up when the instance is closed or garbage-collected.
        robust_mode: If True, wait for expect patterns between commands (slower but
            safer).
            If False, send commands without waiting (faster but may desync).
        debug: If True, print debug info including sent commands and expect patterns.
        env_variables: Environment variables to use for the Multiwfn process.
        max_wait: Maximum time to wait for progress updates in seconds.
        has_spin: Spin state. If None, detect from file for wavefunction inputs
            (True=open-shell, False=closed-shell). For .cub/.grd inputs, spin
            state is left undefined.
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
    ) -> None:
        self._file_path = Path(file_path).resolve()
        if not self._file_path.exists():
            raise FileNotFoundError(f"File not found: {self._file_path}")

        if run_path:
            self._run_path = Path(run_path).resolve()
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
        result = self.run_commands([CommandStep("q", expect="gracefully")])
        mult_pattern = r"Expected multiplicity:\s*(\d+)"

        match = re.search(mult_pattern, result.stdout)
        if not match:
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
        self._copy_settings_ini(self._run_path)

    def load_settingsini(self, file_path: str | Path) -> None:
        """Backward-compatible alias for load_settingini()."""
        self.load_settingini(file_path)

    def close(self) -> None:
        """Clean up temporary working directory, if created."""
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

    def _build_fuzzy_integration_commands(
        self,
        menu_cmd: str,
        menu_pattern: str,
        prefix_commands: list[CommandStep] | None = None,
    ) -> list[CommandStep | WaitProgress]:
        """Build command sequence for fuzzy atomic space integration."""
        commands: list[CommandStep | WaitProgress] = list(prefix_commands or [])
        commands.extend(
            [
                CommandStep("15", expect="15 Fuzzy atomic space analysis"),
                CommandStep("1", expect="1 Perform integration in fuzzy atomic spaces"),
                CommandStep(menu_cmd, expect=menu_pattern),
                WaitProgress(self._max_wait),
                CommandStep("0", expect="0 Return"),
                CommandStep("q", expect="gracefully"),
            ]
        )
        return commands

    @staticmethod
    def _parse_user_function(
        function_setting: RealFunctionSetting,
    ) -> tuple[int, str | None]:
        """Parse a user-defined scalar function setting."""
        if isinstance(function_setting, tuple):
            return function_setting
        return function_setting, None

    def _parse_vector_functions(self, descriptor: str) -> tuple[int, int, int]:
        """Parse vector function settings and cache citation if present."""
        vector_setting = VECTORS[descriptor]
        if len(vector_setting) == 2 and isinstance(vector_setting[1], str):
            functions = vector_setting[0]
            self._results.citations.add(vector_setting[1])
            return functions
        return vector_setting

    def _setup_workdir(self, subdir: str | None) -> Path:
        """Get or create working directory."""
        workdir = self._run_path / subdir if subdir else self._run_path
        workdir.mkdir(parents=True, exist_ok=True)

        # Copy a custom settings.ini file from run_path
        self._copy_settings_ini(workdir)
        return workdir

    def _set_env(self) -> dict[str, str]:
        """Set environment variables for Multiwfn execution."""
        return build_execution_env(env_variables=self._env_variables)

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

    @requires_executable(["Multiwfn"])
    def run_commands(
        self,
        commands: Iterable[str | CommandStep | WaitProgress],
        subdir: str | None = None,
    ) -> MultiwfnRunResult:
        """Run Multiwfn with a custom command sequence.

        Args:
            commands: Sequence of commands (strings or CommandStep objects).
            subdir: Subdirectory within run_path to store results.

        Returns:
            MultiwfnRunResult with stdout, workdir
        """
        workdir = self._setup_workdir(subdir)

        env = self._set_env()
        session = _PexpectSession("Multiwfn", workdir, self._debug, env=env)
        session.send(str(self._file_path))
        session.read_available()

        for item in commands:
            session._process_command(item, self._robust_mode)

        session.wait_for_exit()

        return MultiwfnRunResult(
            stdout=session.stdout,
            workdir=workdir,
        )

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

        Returns:
            Atomic charges indexed by atom number (1-indexed)
        """
        normalized_model = self._normalize_option(model)
        if (
            self._results.charges is not None
            and normalized_model in self._results.charges
        ):
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
                    "1 Use build-in sphericalized atomic"
                    if normalized_model != "mulliken"
                    else " 1 Output Mulliken population and atomic"
                ),
            ),
            CommandStep("y", expect="y/n"),
        ]
        if normalized_model == "mulliken":
            commands += [CommandStep("0", expect="Mulliken population")]
        commands += [
            CommandStep("0", expect="0 Return"),
            CommandStep("q", expect="gracefully"),
        ]

        citation = CHARGE_MODEL_CITATIONS.get(normalized_model, None)
        if citation is not None:
            self._results.citations.add(citation)

        result = self.run_commands(commands, subdir=normalized_model)

        chg_files = sorted(result.workdir.glob("*.chg"))
        if not chg_files:
            charges = {}
        else:
            charges = self._parse_chg_file(chg_files[0])

        if self._results.charges is None:
            self._results.charges = {}
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
        if (
            self._results.surfaces is not None
            and normalized_model in self._results.surfaces
        ):
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
            CommandStep("q", expect="Exit program gracefully"),
        ]

        result = self.run_commands(commands, subdir=normalized_model)

        atomic_props = self._parse_atomic_surface_properties(result.stdout)
        global_props = self._parse_global_surface_properties(result.stdout)
        surface_result = {"atomic": atomic_props, "global": global_props}

        if self._results.surfaces is None:
            self._results.surfaces = {}
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
        if (
            self._results.bond_orders is not None
            and normalized_model in self._results.bond_orders
        ):
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
            CommandStep("q", expect="gracefully"),
        ]

        result = self.run_commands(commands, subdir=normalized_model)
        bond_order_matrix = self._parse_bond_orders(result.stdout)

        if self._results.bond_orders is None:
            self._results.bond_orders = {}
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

        grid_key = str(active_grid)
        if (
            self._results.grid_descriptors is not None
            and grid_key in self._results.grid_descriptors
        ):
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
            result = self.run_commands(commands, subdir="grid")
        finally:
            self._file_path = file_path

        grid_descriptors = self._parse_atomic_values(result.stdout)

        if self._results.grid_descriptors is None:
            self._results.grid_descriptors = {}
        self._results.grid_descriptors[grid_key] = grid_descriptors
        return grid_descriptors

    def get_grid(
        self,
        descriptor: str,
        grid_quality: str,
        grid_file_name: str | None = None,
    ) -> Path:
        """Generate a grid (.cub) file for a real space function.

        Args:
            descriptor: Descriptor to calculate on grid.
            grid_quality: Grid quality, one of 'low', 'medium', 'high'.
            grid_file_name: Optional name for output cube file.

        Returns:
            Path to the generated cube file.

        Raises:
            ValueError: If the descriptor is not supported.
            FileNotFoundError: If no cube file was generated.
        """
        menu_cmd, menu_pattern, commands = self._get_descriptor_function(descriptor)

        normalized_grid_quality = self._normalize_option(grid_quality)
        if normalized_grid_quality not in GRID_QUALITIES:
            choices = ", ".join(GRID_QUALITIES.keys())
            raise ValueError(
                f"Grid quality {grid_quality!r} not supported. Choose between {choices}."
            )
        grid_pattern = GRID_QUALITIES[normalized_grid_quality]
        grid_cmd = grid_pattern.split(" ")[0]
        existing_cub_files: set[Path] = set()
        grid_workdir = self._run_path / descriptor
        if grid_workdir.exists():
            existing_cub_files = set(grid_workdir.glob("*.cub"))

        grid_commands: list[CommandStep | WaitProgress] = [
            *commands,
            CommandStep(
                "5",
                expect="5 Output and plot specific property within a spatial region",
            ),
            CommandStep(menu_cmd, expect=menu_pattern),
            CommandStep(grid_cmd, expect=grid_pattern),
            WaitProgress(self._max_wait),
            CommandStep("2", expect="2 Export data to a Gaussian-type cube file"),
            CommandStep("0", expect="0 Return to main menu"),
            CommandStep("q", expect="gracefully"),
        ]
        result = self.run_commands(grid_commands, subdir=descriptor)

        cub_candidates: list[Path] = list(result.workdir.glob("*.cub"))
        if not cub_candidates:
            raise FileNotFoundError("No .cub file was generated by Multiwfn.")
        new_cub_candidates: list[Path] = [
            path for path in cub_candidates if path not in existing_cub_files
        ]
        selection_pool: list[Path] = (
            new_cub_candidates if new_cub_candidates else cub_candidates
        )
        cub_file_path: Path = max(
            selection_pool, key=lambda path: path.stat().st_mtime_ns
        )
        if grid_file_name is not None:  # Rename the file
            new_path: Path = result.workdir / grid_file_name
            cub_file_path.rename(new_path)
            return new_path

        return cub_file_path

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
            CommandStep("q", expect="gracefully"),
        ]
        result = self.run_commands(commands, subdir="fukui")
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
            # Do not gate quit on a specific prompt here; conceptual-DFT submenu
            # prompts vary across Multiwfn builds.
            CommandStep("q"),
        ]
        result = self.run_commands(commands, subdir="superdeloc")
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

        if self._results.atomic_vector_descriptors is None:
            self._results.atomic_vector_descriptors = {}
        elif descriptor in self._results.atomic_vector_descriptors:
            return self._results.atomic_vector_descriptors[descriptor]

        menu_cmd = "100"
        menu_pattern = "100 User-defined function"
        functions = self._parse_vector_functions(descriptor)

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

            if menu_cmd in {"13"}:
                raise ValueError(f"Descriptor {descriptor!r} requires grid data!")

        else:
            function_setting = USER_REAL_FUNCTIONS[descriptor]
            func, citation = self._parse_user_function(function_setting)
            if citation is not None:
                self._results.citations.add(citation)
            menu_cmd = "100"
            menu_pattern = "100 User-defined function"

            self._check_function_selection(func, descriptor)
            commands.extend(self._commands_change_iuserfunc(function=func))

        return menu_cmd, menu_pattern, commands

    def get_descriptor(self, descriptor: str) -> dict[int, float]:
        """Calculate atomic densities from integration in fuzzy atomic spaces.

        Args:
            descriptor: Descriptor name.

        Returns:
            Dictionary mapping atom index (1-based) to integrated density value.

        Raises:
            ValueError: If the descriptor is unsupported, vector-valued, requires grid
                data, or is incompatible with the current spin state.
        """
        if (
            self._results.atomic_descriptors is not None
            and descriptor in self._results.atomic_descriptors
        ):
            return self._results.atomic_descriptors[descriptor]

        try:
            menu_cmd, menu_pattern, commands = self._get_descriptor_function(descriptor)
            descriptor_result = self._run_fuzzy_integration(
                menu_cmd, menu_pattern, commands, descriptor
            )
            if self._results.atomic_descriptors is None:
                self._results.atomic_descriptors = {}
            self._results.atomic_descriptors[descriptor] = descriptor_result
            return descriptor_result
        except ValueError as e:
            raise e

    def _run_fuzzy_integration(
        self,
        menu_cmd: str,
        menu_pattern: str,
        commands: list[CommandStep],
        descriptor: str,
    ) -> dict[int, float]:
        command_sequence = self._build_fuzzy_integration_commands(
            menu_cmd,
            menu_pattern,
            prefix_commands=commands,
        )
        result = self.run_commands(command_sequence, subdir=descriptor)
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
            CommandStep("q", expect="gracefully"),
        ]
        result = self.run_commands(commands, subdir="electric_moments")
        moments = self._parse_electric_moments(result.stdout)

        if self._results.electric_moments is None:
            self._results.electric_moments = {}
        self._results.electric_moments.update(moments)

        return moments

    def get_gaps(
        self,
        n: int = 3,
        occupation_thresholds: tuple[float, float] = (0.5, 1.5),
    ) -> dict[str, Any]:
        """Return orbital energies and HOMO-centered adjacent orbital gaps.

        Args:
            n: Half-window size around HOMO. A window of `n=3` spans
                `homo_index - 3` through `homo_index + 3`.
            occupation_thresholds: `(lower, upper)` thresholds used to classify
                unoccupied, partially occupied and occupied orbitals.

        Returns:
            A dictionary with:
            - `orbitals`: `{index: (energy_ev, occupation, orbital_type)}`
            - `gaps`: HOMO index/energy and `sequence_ev`, where each entry is
              `{(i, j): gap_ev}` for adjacent orbitals `i -> j`.

        Raises:
            ValueError: If spin state is unknown, `n < 0`, or thresholds are invalid.
        """
        if self._has_spin is None:
            raise ValueError(
                "Gap analysis requires known spin state. "
                "Initialize Multiwfn with has_spin=True/False."
            )
        if n < 0:
            raise ValueError("Gap window n must be >= 0.")

        lower, upper = self._validate_occupation_thresholds(occupation_thresholds)
        cache_key = self._gap_cache_key(n=n, lower=lower, upper=upper)

        cached = self._results.gap
        if cached is not None and cache_key in cached:
            return cached[cache_key]

        orbitals = self._get_orbital_listing(subdir="orbital_gaps")
        gap_result = self._build_homo_window_gaps(
            orbitals,
            n=n,
            occupation_thresholds=(lower, upper),
        )

        if cached is None:
            cached = {}
            self._results.gap = cached
        cached[cache_key] = gap_result
        return gap_result

    def _get_homo_lumo(self) -> dict[str, float | int]:
        """Calculate eV-based HOMO/LUMO energies and HOMO-LUMO gap."""
        commands = [
            CommandStep("0", expect="0 Show molecular structure and view orbitals"),
            CommandStep("q", expect="gracefully"),
        ]
        result = self.run_commands(commands, subdir="homo_lumo")
        return self._parse_homo_lumo_gap(result.stdout)

    def _get_somo_gaps(
        self, occupation_thresholds: tuple[float, float] = (0.5, 1.5)
    ) -> dict[str, Any]:
        """Calculate eV-based SOMO gaps to overall HOMO and LUMO."""
        orbitals = self._get_orbital_listing(subdir="somo_gaps")
        return self._compute_somo_gaps(orbitals, occupation_thresholds)

    @staticmethod
    def _orbital_listing_commands() -> list[CommandStep]:
        """Build command sequence for the `List all orbitals` screen."""
        return [
            CommandStep("6", expect="6 Check & modify wavefunction"),
            CommandStep("3", expect="3 List all orbitals"),
            CommandStep("-1", expect="-1 Return"),
            CommandStep("q", expect="gracefully"),
        ]

    def _get_orbital_listing(self, subdir: str) -> list[dict[str, Any]]:
        """Run orbital listing and parse rows."""
        result = self.run_commands(self._orbital_listing_commands(), subdir=subdir)
        return self._parse_orbital_list(result.stdout)

    def _gap_cache_key(self, n: int, lower: float, upper: float) -> str:
        """Build deterministic cache key for gap requests."""
        if self._has_spin is None:
            raise ValueError(
                "Gap analysis requires known spin state. "
                "Initialize Multiwfn with has_spin=True/False."
            )
        spin_label = "open" if self._has_spin else "closed"
        return f"{spin_label}:n:{n}:lower:{lower:.8g}:upper:{upper:.8g}"

    @staticmethod
    def _validate_occupation_thresholds(
        occupation_thresholds: tuple[float, float],
    ) -> tuple[float, float]:
        """Validate and normalize occupation thresholds."""
        lower, upper = occupation_thresholds
        if lower < 0:
            raise ValueError("Lower occupation threshold must be >= 0.")
        if upper <= lower:
            raise ValueError("Upper occupation threshold must be greater than lower.")
        return lower, upper

    @staticmethod
    def _split_orbitals_by_occupation(
        orbitals: list[dict[str, Any]],
        lower: float,
        upper: float,
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
        """Split orbitals into unoccupied, partially occupied and occupied groups."""
        unoccupied = [orb for orb in orbitals if orb["occupation"] < lower]
        partially_occupied = [
            orb for orb in orbitals if lower <= orb["occupation"] < upper
        ]
        occupied = [orb for orb in orbitals if orb["occupation"] >= upper]
        return unoccupied, partially_occupied, occupied

    def get_localization(self) -> dict[int, float]:
        """Calculate localization index from integration in fuzzy atomic spaces."""
        if self._results.localization is not None:
            return self._results.localization

        commands = [
            CommandStep("15", expect="15 Fuzzy atomic space analysis"),
            CommandStep("4", expect="4 Calculate localization"),
            WaitProgress(self._max_wait),
            CommandStep("n", expect="y/n"),
            CommandStep("0", expect="0 Return"),
            CommandStep("q", expect="gracefully"),
        ]
        result = self.run_commands(commands, subdir="localization")
        localization = self._parse_localization_values(result.stdout)

        self._results.localization = localization
        return localization

    def _parse_atomic_table(
        self,
        stdout: str,
        header_keywords: list[str],
        keys: list[str],
        stop_keyword: str | None = None,
    ) -> dict[str, dict[int, float]]:
        """Parse atomic table with 4 float columns from stdout.

        Args:
            stdout: Output to parse
            header_keywords: List of strings that must all appear in the header line
            keys: List of 4 keys for the output dictionary
            stop_keyword: Optional keyword to stop parsing (e.g., "Sum of")

        Returns:
            Dictionary mapping key names to {atom_idx: value}
        """
        lines = stdout.replace("\r", "\n").splitlines()
        result: dict[str, dict[int, float]] = {key: {} for key in keys}
        row_pattern = re.compile(
            (
                rf"\s*(\d+)\([A-Za-z][a-z]?\s*\)?\s+({NUMBER_PATTERN})\s+"
                rf"({NUMBER_PATTERN})\s+({NUMBER_PATTERN})\s+({NUMBER_PATTERN})"
            )
        )

        for i, line in enumerate(lines):
            if all(keyword in line for keyword in header_keywords):
                # Found the header line
                for j in range(i + 1, len(lines)):
                    data_line = lines[j]
                    if stop_keyword and stop_keyword in data_line:
                        break
                    if not data_line.strip():
                        continue
                    # Parse format: "     1(C        val1        val2        val3        val4"
                    match = row_pattern.search(data_line)
                    if match:
                        atom_idx = int(match.group(1))
                        for k, key in enumerate(keys, start=2):
                            result[key][atom_idx] = self._parse_float(match.group(k))
                break

        return result

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
        return self._parse_atomic_table(
            stdout,
            header_keywords=["Atom index", "OW f+"],
            keys=["f_plus", "f_minus", "f_zero", "dd"],
        )

    def _parse_superdelocalizabilities(
        self, stdout: str
    ) -> dict[str, dict[int, float]]:
        """Parse superdelocalizability descriptors from stdout.

        Multiwfn output differs across versions/builds:
        - Header labels may be `D_N`/`D_E` or `DN`/`DE`
        - Rows may contain either 2 columns (D_N, D_E) or 4 columns
          (D_N, D_E, D_N_0, D_E_0)

        Args:
            stdout: Multiwfn output text containing the superdelocalizability table.

        Returns:
            Mapping with keys `d_n`, `d_e`, `d_n_0`, and `d_e_0`, each containing
            atom-indexed values when available.
        """
        result: dict[str, dict[int, float]] = {
            "d_n": {},
            "d_e": {},
            "d_n_0": {},
            "d_e_0": {},
        }
        lines = stdout.replace("\r", "\n").splitlines()

        atom_header_re = re.compile(r"\bAtom(?:\s+index)?\b", flags=re.IGNORECASE)
        dn_header_re = re.compile(r"\bD[_\s]?N\b", flags=re.IGNORECASE)
        de_header_re = re.compile(r"\bD[_\s]?E\b", flags=re.IGNORECASE)
        row_re = re.compile(r"\s*(\d+)\([A-Za-z][a-z]?\s*\)?\s+(.*)$")

        for i, line in enumerate(lines):
            if (
                atom_header_re.search(line)
                and dn_header_re.search(line)
                and de_header_re.search(line)
            ):
                for data_line in lines[i + 1 :]:
                    if "sum of" in data_line.lower():
                        break
                    if not data_line.strip():
                        continue

                    row_match = row_re.search(data_line)
                    if row_match is None:
                        continue

                    atom_idx = int(row_match.group(1))
                    values = re.findall(NUMBER_PATTERN, row_match.group(2))
                    if len(values) < 2:
                        continue

                    result["d_n"][atom_idx] = self._parse_float(values[0])
                    result["d_e"][atom_idx] = self._parse_float(values[1])
                    if len(values) >= 4:
                        result["d_n_0"][atom_idx] = self._parse_float(values[2])
                        result["d_e_0"][atom_idx] = self._parse_float(values[3])
                break

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

    def _parse_homo_lumo_gap(self, stdout: str) -> dict[str, float | int]:
        """Parse eV-based HOMO/LUMO energies and HOMO-LUMO gap from stdout."""
        homo_match = HOMO_LINE_PATTERN.search(stdout)
        lumo_match = LUMO_LINE_PATTERN.search(stdout)
        gap_match = HOMO_LUMO_GAP_PATTERN.search(stdout)

        if homo_match is None or lumo_match is None or gap_match is None:
            raise RuntimeError("Could not parse HOMO/LUMO energies and HOMO-LUMO gap.")

        return {
            "homo_index": int(homo_match.group(1)),
            "homo_energy_ev": self._parse_float(homo_match.group(2)),
            "lumo_index": int(lumo_match.group(1)),
            "lumo_energy_ev": self._parse_float(lumo_match.group(2)),
            "homo_lumo_gap_ev": self._parse_float(gap_match.group(1)),
        }

    def _parse_somo_gaps(
        self, stdout: str, occupation_thresholds: tuple[float, float] = (0.5, 1.5)
    ) -> dict[str, Any]:
        """Parse orbital list and compute eV-based SOMO gaps."""
        orbitals = self._parse_orbital_list(stdout)
        return self._compute_somo_gaps(orbitals, occupation_thresholds)

    def _compute_somo_gaps(
        self,
        orbitals: list[dict[str, Any]],
        occupation_thresholds: tuple[float, float] = (0.5, 1.5),
    ) -> dict[str, Any]:
        """Compute eV-based SOMO-related energy gaps from parsed orbital rows."""
        lower, upper = self._validate_occupation_thresholds(occupation_thresholds)
        unoccupied, somos, occupied = self._split_orbitals_by_occupation(
            orbitals, lower, upper
        )
        if not occupied or not unoccupied:
            raise RuntimeError(
                "Could not identify doubly occupied and unoccupied orbitals."
            )
        if not somos:
            raise RuntimeError("No partially occupied orbitals (SOMOs) were found.")

        overall_homo = max(occupied, key=lambda orb: orb["energy_ev"])
        overall_lumo = min(unoccupied, key=lambda orb: orb["energy_ev"])

        sorted_somos = sorted(somos, key=lambda orb: orb["index"])
        somo_indices = [somo["index"] for somo in sorted_somos]
        somo_to_lumo_gap_ev: dict[int, float] = {}
        somo_to_homo_gap_ev: dict[int, float] = {}

        for somo in sorted_somos:
            somo_idx = somo["index"]
            somo_to_lumo_gap_ev[somo_idx] = (
                overall_lumo["energy_ev"] - somo["energy_ev"]
            )
            somo_to_homo_gap_ev[somo_idx] = (
                somo["energy_ev"] - overall_homo["energy_ev"]
            )

        return {
            "overall_homo_index": int(overall_homo["index"]),
            "overall_homo_energy_ev": overall_homo["energy_ev"],
            "overall_lumo_index": int(overall_lumo["index"]),
            "overall_lumo_energy_ev": overall_lumo["energy_ev"],
            "somo_indices": somo_indices,
            "somo_to_lumo_gap_ev": somo_to_lumo_gap_ev,
            "somo_to_homo_gap_ev": somo_to_homo_gap_ev,
            "occupation_thresholds": {"lower": lower, "upper": upper},
        }

    def _build_homo_window_gaps(
        self,
        orbitals: list[dict[str, Any]],
        n: int = 3,
        occupation_thresholds: tuple[float, float] = (0.5, 1.5),
    ) -> dict[str, Any]:
        """Build HOMO-centered gap sequence with a symmetric index window."""
        lower, upper = self._validate_occupation_thresholds(occupation_thresholds)
        if n < 0:
            raise ValueError("Gap window n must be >= 0.")

        if not orbitals:
            raise RuntimeError("No orbitals available for gap analysis.")

        _, partially_occupied, occupied = self._split_orbitals_by_occupation(
            orbitals, lower, upper
        )
        # Unrestricted orbital listings may contain only 1/0 occupations.
        # In that case, fall back to the partially occupied bucket.
        occupied_candidates = occupied or partially_occupied
        if not occupied_candidates:
            raise RuntimeError("Could not identify occupied orbitals.")

        homo = max(occupied_candidates, key=lambda orb: orb["energy_ev"])
        homo_index = homo["index"]

        orbital_table: dict[int, tuple[float, float, str]] = {}
        for orb in orbitals:
            idx = orb["index"]
            orbital_table[idx] = (
                orb["energy_ev"],
                orb["occupation"],
                orb["type"],
            )

        sorted_indices = sorted(orbital_table)
        start_idx = max(sorted_indices[0], homo_index - n)
        end_idx = min(sorted_indices[-1], homo_index + n)

        sequence_ev: list[dict[tuple[int, int], float]] = []
        for idx in range(start_idx, end_idx):
            if idx not in orbital_table or (idx + 1) not in orbital_table:
                continue
            gap_ev = orbital_table[idx + 1][0] - orbital_table[idx][0]
            sequence_ev.append({(idx, idx + 1): gap_ev})

        return {
            "orbitals": orbital_table,
            "gaps": {
                "homo_index": homo_index,
                "homo_energy_ev": homo["energy_ev"],
                "sequence_ev": sequence_ev,
            },
        }

    def _parse_orbital_list(self, stdout: str) -> list[dict[str, Any]]:
        """Parse orbital listing from `List all orbitals` output."""
        orbitals = [
            {
                "index": int(match["index"]),
                "energy_ev": self._parse_float(match["energy_ev"]),
                "occupation": self._parse_float(match["occupation"]),
                "type": match["orbital_type"],
            }
            for match in ORBITAL_LINE_PATTERN.finditer(stdout)
        ]

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
        lines = stdout.replace("\r", "\n").splitlines()
        atomic_values: dict[int, float] = {}
        value_pattern = NUMBER_PATTERN

        for i, line in enumerate(lines):
            if "Atomic space" in line and "Value" in line:
                for j in range(i + 1, len(lines)):
                    data_line = lines[j]
                    if "Summing up" in data_line:
                        break
                    if not data_line.strip():
                        continue
                    # Parse format: "    1(C )            0.00607663"
                    match = re.search(
                        rf"\s*(\d+)\([A-Z][a-z]?\s*\)\s+({value_pattern})",
                        data_line,
                    )
                    if match:
                        atom_idx = int(match.group(1))
                        value = self._parse_float(match.group(2))
                        atomic_values[atom_idx] = value
                break

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

    def _parse_atomic_areas(
        self, lines: list[str], start_idx: int
    ) -> tuple[dict[int, dict[str, float]], int]:
        """Parse atomic surface areas section from stdout."""
        atomic_props: dict[int, dict[str, float]] = {}
        for i in range(start_idx + 1, len(lines)):
            line = lines[i]
            if "Atom#   All/Positive/Negative average" in line:
                return atomic_props, i
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 6 and parts[0].isdigit():
                try:
                    atom_idx = int(parts[0])
                    atomic_props[atom_idx] = {
                        "area_total": self._parse_float(parts[1]),
                        "area_positive": self._parse_float(parts[2]),
                        "area_negative": self._parse_float(parts[3]),
                        "min_value": self._parse_float(parts[4]),
                        "max_value": self._parse_float(parts[5]),
                    }
                except ValueError:
                    continue
        return atomic_props, len(lines)

    def _parse_atomic_averages(
        self,
        lines: list[str],
        atomic_props: dict[int, dict[str, float]],
        start_idx: int,
    ) -> dict[int, dict[str, float]]:
        """Parse atomic surface averages section from stdout."""
        for i in range(start_idx + 1, len(lines)):
            line = lines[i]
            if not line.strip():
                break
            parts = line.split()
            if len(parts) >= 7 and parts[0].isdigit():
                try:
                    atom_idx = int(parts[0])
                    if atom_idx in atomic_props:
                        atomic_props[atom_idx].update(
                            {
                                "avg_all": self._parse_float(parts[1]),
                                "avg_positive": self._parse_float(parts[2]),
                                "avg_negative": self._parse_float(parts[3]),
                                "var_all": self._parse_float(parts[4]),
                                "var_positive": self._parse_float(parts[5]),
                                "var_negative": self._parse_float(parts[6]),
                            }
                        )
                except ValueError:
                    continue
        return atomic_props

    def _parse_atomic_charge_separation(
        self,
        lines: list[str],
        atomic_props: dict[int, dict[str, float]],
        start_idx: int,
    ) -> dict[int, dict[str, float]]:
        """Parse atomic internal charge separation section from stdout."""
        for i in range(start_idx + 1, len(lines)):
            line = lines[i]
            if not line.strip():
                break
            parts = line.split()
            if len(parts) >= 4 and parts[0].isdigit():
                atom_idx = int(parts[0])
                if atom_idx in atomic_props:
                    atomic_props[atom_idx].update(
                        {
                            "pi": self._parse_float(parts[1]),
                            "nu": self._parse_float(parts[2]),
                            "nu_sigma2": self._parse_float(parts[3]),
                        }
                    )
        return atomic_props

    def _parse_atomic_surface_properties(
        self, stdout: str
    ) -> dict[int, dict[str, float]]:
        """Parse atomic surface properties from stdout."""
        lines = stdout.replace("\r", "\n").splitlines()
        atomic_props: dict[int, dict[str, float]] = {}

        for i, line in enumerate(lines):
            if "All/Positive/Negative area" in line:
                atomic_props, avg_idx = self._parse_atomic_areas(lines, i)
                self._parse_atomic_averages(lines, atomic_props, avg_idx)
                for j in range(avg_idx + 1, len(lines)):
                    if "Atom#           Pi" in lines[j]:
                        self._parse_atomic_charge_separation(lines, atomic_props, j)
                        break
                break
            if "Atom#      Area(Ang^2)" in line:
                atomic_props = self._parse_atomic_simple_table(lines, i)
                break

        return atomic_props

    def _parse_atomic_simple_table(
        self, lines: list[str], start_idx: int
    ) -> dict[int, dict[str, float]]:
        """Parse atomic surface properties from single-table format."""
        atomic_props: dict[int, dict[str, float]] = {}
        for i in range(start_idx + 1, len(lines)):
            line = lines[i]
            if not line.strip():
                break
            parts = line.split()
            if len(parts) >= 6 and parts[0].isdigit():
                atom_idx = int(parts[0])
                atomic_props[atom_idx] = {
                    "area_total": self._parse_float(parts[1]),
                    "min_value": self._parse_float(parts[2]),
                    "max_value": self._parse_float(parts[3]),
                    "avg_all": self._parse_float(parts[4]),
                    "var_all": self._parse_float(parts[5]),
                }
        return atomic_props

    def _parse_global_surface_properties(self, stdout: str) -> dict[str, float]:
        """Parse global surface properties from stdout."""
        lines = stdout.replace("\r", "\n").splitlines()
        props: dict[str, float] = {}
        in_summary = False
        for line in lines:
            if "================= Summary of surface analysis" in line:
                in_summary = True
                continue
            if not in_summary:
                continue
            if "Surface analysis finished!" in line:
                break
            if not line.strip() or ":" not in line:
                continue

            handled = self._parse_summary_minmax(line, props)
            if handled:
                continue
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
