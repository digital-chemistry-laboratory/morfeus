"""xtb interface code."""

from __future__ import annotations

from collections.abc import Iterable
import functools

from typing import Any, cast
import numpy as np
from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import nullcontext
import subprocess
import shutil
from dataclasses import dataclass
import re
import json

from morfeus.data import DEBYE_TO_AU, AU_TO_DEBYE
from morfeus.io import read_geometry, write_xyz, write_xtb_inp
from morfeus.typing import Array1DFloat, Array2DFloat, ArrayLike2D
from morfeus.utils import convert_elements, requires_executable


@dataclass
class XTBResults:
    """Stores xTB descriptors."""

    charges: list[float] | None = None
    bond_orders: list[tuple[int, int, float]] | None = None
    homo: dict[str, float] | None = None  # Stores both Eh and eV units
    lumo: dict[str, float] | None = None  # Stores both Eh and eV units
    gap: float | None = None  # Unit in eV
    atom_dipole_vect: Array2DFloat | None = None  # Unit in a.u.
    dipole_vect: Array1DFloat | None = None  # Unit in a.u.
    dipole_moment: float | None = None  # Unit in debye
    atom_polarizabilities: list[float] | None = None
    mol_polarizability: float | None = None
    g_solv: float | None = None  # Unit in Eh
    g_solv_hb: float | None = None  # Unit in Eh
    atom_hb_strengths: list[float] | None = None
    fod_pop: list[float] | None = None
    ip: float | None = None
    ea: float | None = None
    fukui_plus: list[float] | None = None
    fukui_minus: list[float] | None = None
    fukui_radical: list[float] | None = None


@requires_executable(["xtb"])
class XTB:
    """Calculates electronic properties with the xtb program.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        method: Method to use in xtb. Currently works with:
            - 2: GFN2-xTB (default)
            - 1: GFN1-xTB
            - ptb: PTB
        charge: Molecular charge
        n_unpaired: Number of unpaired electrons
        solvent: Solvent. Uses the ALPB solvation model
        electronic_temperature: Electronic temperature (K)
        n_processes: Number of parallel processes in xtb
        run_path: Folder path to run xTB calculation. If not provided, runs in a temporary folder
    """

    _charge: int
    _coordinates: Array2DFloat
    _electronic_temperature: int | None
    _elements: Array1DFloat
    _n_unpaired: int | None
    _results: XTBResults
    _solvent: str | None
    _method: int | str
    _n_processes: int | None

    _xyz_input_file: str = "xtb.xyz"
    _xtb_input_file: str = "xtb.inp"

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        method: int | str = 2,
        charge: int = 0,
        n_unpaired: int | None = None,
        solvent: str | None = None,
        electronic_temperature: int | None = None,
        n_processes: int | None = None,
        run_path: Path | str | None = None,
    ) -> None:
        # Converting elements to atomic numbers if the are symbols
        self._elements = np.array(convert_elements(elements, output="numbers"))

        # Store settings
        self._coordinates = np.array(coordinates)
        self._method = str(method).lower()
        self._charge = charge
        self._solvent = solvent
        self._n_unpaired = n_unpaired
        self._electronic_temperature = electronic_temperature
        self._n_processes = n_processes

        self._run_path = Path(run_path) if run_path else None

        self._default_xtb_command = (
            f"xtb {XTB._xyz_input_file} --json --chrg {self._charge}"
        )
        if self._method in ["1", "2"]:
            self._default_xtb_command += f" --gfn {int(self._method)}"
        elif self._method == "ptb":
            self._default_xtb_command += " --ptb"
        else:
            raise ValueError(
                f"Method {self._method!r} not supported. Choose between: 2, 1, or ptb."
            )
        if self._solvent is not None:
            if self._method == "ptb":
                raise ValueError(
                    "Solvation is not available with PTB. Remove solvent or use another xTB method."
                )
            self._default_xtb_command += f" --alpb {self._solvent}"
        if self._n_unpaired is not None:
            self._default_xtb_command += f" --uhf {self._n_unpaired}"
        if self._electronic_temperature is not None:
            self._default_xtb_command += f" --etemp {self._electronic_temperature}"
        if self._n_processes is not None:
            self._default_xtb_command += f" --parallel {self._n_processes}"

        self._results = XTBResults()
        self._corrected: bool | None = None

    def get_bond_order(self, i: int, j: int) -> float:
        """Returns bond order between two atoms.

        Args:
            i: Index of atom 1 (1-indexed)
            j: Index of atom 2 (1-indexed)

        Returns:
            bond_order: Bond order

        Raises:
            ValueError: If no bond exists between the given atoms.
        """
        bonds_orders = self.get_bond_orders()
        if (i, j) in bonds_orders:
            bond_order = bonds_orders[(i, j)]
        elif (j, i) in bonds_orders:
            bond_order = bonds_orders[(j, i)]
        else:
            raise ValueError(f"No bond order calculated between atoms {i} and {j}.")

        return bond_order

    def get_bond_orders(self) -> dict[tuple[int, int], float]:
        """Returns bond orders."""

        if self._results.bond_orders is None:
            self._run_xtb("sp")
            self._results.bond_orders = cast(
                list[tuple[int, int, float]], self._results.bond_orders
            )
        bond_orders = {(x[0], x[1]): round(x[2], 12) for x in self._results.bond_orders}

        return bond_orders

    def get_charges(self) -> dict[int, float]:
        """Returns atomic charges."""

        if self._results.charges is None:
            self._run_xtb("sp")
            self._results.charges = cast(list[float], self._results.charges)
        charges = {i: charge for i, charge in enumerate(self._results.charges, start=1)}

        return charges

    def get_homo(self, unit="Eh") -> float:
        """Returns HOMO energy.

        Args:
            unit: 'Eh' or 'eV'

        Returns:
            HOMO energy (Eh or eV)

        Raises:
            ValueError: If unit does not exist
        """

        if self._results.homo is None:
            self._run_xtb("sp")
            self._results.homo = cast(dict[str, float], self._results.homo)

        if unit == "Eh":
            return self._results.homo["Eh"]
        elif unit == "eV":
            return self._results.homo["eV"]
        else:
            raise ValueError("Unit must be either 'Eh' or 'eV'.")

    def get_lumo(self, unit="Eh") -> float:
        """Returns LUMO energy.

        Args:
            unit: 'Eh' or 'eV'

        Returns:
            LUMO energy (Eh or eV)

        Raises:
            ValueError: If unit does not exist
        """

        if self._results.lumo is None:
            self._run_xtb("sp")
            self._results.lumo = cast(dict[str, float], self._results.lumo)

        if unit == "Eh":
            return self._results.lumo["Eh"]
        elif unit == "eV":
            return self._results.lumo["eV"]
        else:
            raise ValueError("Unit must be either 'Eh' or 'eV'.")

    def get_homo_lumo_gap(self) -> float:
        """Returns HOMO-LUMO gap (eV)."""

        if self._results.gap is None:
            self._run_xtb("sp")
            self._results.gap = cast(float, self._results.gap)

        return self._results.gap

    def get_dipole(self) -> Array1DFloat:
        """Returns molecular dipole vector (a.u.)."""

        if self._results.dipole_vect is None:
            self._run_xtb("sp")

        return self._results.dipole_vect

    def get_dipole_moment(self, unit="au") -> float:
        """Returns molecular dipole moment.

        Args:
            unit: 'au' or 'debye'

        Returns:
            Molecular dipole moment (a.u. or debye)

        Raises:
            ValueError: If unit does not exist
        """
        if self._results.dipole_moment is None:
            self._run_xtb("sp")
            self._results.dipole_moment = cast(float, self._results.dipole_moment)

        if unit == "debye":
            return self._results.dipole_moment
        elif unit == "au":
            return round(self._results.dipole_moment * DEBYE_TO_AU, 3)
        else:
            raise ValueError("Unit must be either 'au' or 'debye'.")

    def get_atom_dipoles(self) -> dict[int, Array1DFloat]:
        """Returns atomic dipole vectors (a.u.).

        Raises:
            ValueError: If the chosen method is GFN1-xTB (does not support atomic dipoles)
        """

        if self._method == "1":
            raise ValueError(
                "Atomic dipoles are not available with GFN1-xTB. Choose another xtb method."
            )

        if self._results.atom_dipole_vect is None:
            self._run_xtb("sp")
            self._results.atom_dipole_vect = cast(
                np.ndarray, self._results.atom_dipole_vect
            )
        atom_dipole_vectors = {
            i: row for i, row in enumerate(self._results.atom_dipole_vect, start=1)
        }

        return atom_dipole_vectors

    def get_atom_dipole_moments(self, unit="au") -> dict[int, float]:
        """Returns atomic dipole moments.

        Args:
            unit: 'au' or 'debye'

        Returns:
            Atomic dipole moments (a.u. or debye)

        Raises:
            ValueError: If unit does not exist
        """
        dipole_vectors = self.get_atom_dipoles()
        dipole_moments = {
            i: round(float(np.linalg.norm(vect)), 8)
            for i, vect in dipole_vectors.items()
        }
        if unit == "debye":
            dipole_moments = {
                i: round(norm * AU_TO_DEBYE, 8) for i, norm in dipole_moments.items()
            }
        elif unit != "au":
            raise ValueError("Unit must be either 'au' or 'debye'.")

        return dipole_moments

    def get_atom_polarizabilities(self) -> dict[int, float]:
        """Returns atomic polarizabilities.

        Raises:
            ValueError: If the chosen method is not GFN2-xTB (necessary for polarizability calculations)
        """
        if self._method != "2":
            raise ValueError("Polarizability is only available with GFN2-xTB.")

        if self._results.atom_polarizabilities is None:
            self._run_xtb("sp")
            self._results.atom_polarizabilities = cast(
                list[float], self._results.atom_polarizabilities
            )
        atom_polarizabilities = {
            i: polar
            for i, polar in enumerate(self._results.atom_polarizabilities, start=1)
        }

        return atom_polarizabilities

    def get_molecular_polarizability(self) -> float:
        """Returns molecular polarizability.

        Raises:
            ValueError: If the chosen method is not GFN2-xTB (necessary for polarizability calculations)
        """
        if self._method != "2":
            raise ValueError("Polarizability is only available with GFN2-xTB.")

        if self._results.mol_polarizability is None:
            self._run_xtb("sp")
            self._results.mol_polarizability = cast(
                float, self._results.mol_polarizability
            )

        return self._results.mol_polarizability

    def get_solvation_energy(self) -> float:
        """Returns solvation free energy (Eh).

        Raises:
            ValueError: If no solvent is specified (necessary for solvation calculations)
        """
        if self._solvent is None:
            raise ValueError("Solvation energy is only available with solvent.")

        if self._results.g_solv is None:
            self._run_xtb("sp")
            self._results.g_solv = cast(float, self._results.g_solv)

        return self._results.g_solv

    def get_solvation_h_bond_correction(self) -> float:
        """
        Returns hydrogen bonding correction to the solvation free energy (Eh).
        The hydrogen bonding correction is 0.0 for non-polar solvents.

        Raises:
            ValueError: If no solvent is specified (necessary for solvation calculations)
        """

        if self._solvent is None:
            raise ValueError(
                "Hydrogen bonding correction to solvation is only available with solvent."
            )

        if self._results.g_solv_hb is None:
            self._run_xtb("sp")
            self._results.g_solv_hb = cast(float, self._results.g_solv_hb)

        return self._results.g_solv_hb

    def get_atom_solvation_h_bond_strengths(self) -> dict[int, float]:
        """Returns atomic hydrogen bonding strengths related to solvent.

        Raises:
            ValueError: If no solvent is specified (necessary for solvation calculations)
            ValueError: If the specified solvent is not polar (necessary for hydrogen bonding contribution)
        """
        if self._solvent is None:
            raise ValueError(
                "Atomic hydrogen bonding strengths are only available with (polar) solvent."
            )

        if self._results.atom_hb_strengths is None:
            self._run_xtb("sp")
            self._results.atom_hb_strengths = cast(
                list[float], self._results.atom_hb_strengths
            )

        if self._results.atom_hb_strengths == []:
            raise ValueError(
                f"No hydrogen bonding contribution calculated with {self._solvent!r}. Provide a polar solvent."
            )

        atom_hb_strengths = {
            i: strength
            for i, strength in enumerate(self._results.atom_hb_strengths, start=1)
        }

        return atom_hb_strengths

    def get_fod_population(self) -> dict[int, float]:
        """Returns atomic fractional occupation number weighted density population.
        The FOD calculation is performed by default with an electronic temperature of 5000 K.
        """
        if self._results.fod_pop is None:
            self._run_xtb("fod")
            self._results.fod_pop = cast(list[float], self._results.fod_pop)

        fod_pop = {i: pop for i, pop in enumerate(self._results.fod_pop, start=1)}

        return fod_pop

    def get_nfod(self) -> float:
        """Returns NFOD descriptor.
        NFOD is the integration over all space of the fractional occupation number weighted density (FOD).
        The FOD calculation is performed by default with an electronic temperature of 5000 K.
        """
        fod_pop = self.get_fod_population()
        nfod = sum(fod_pop.values())

        return nfod

    def get_ip(self, corrected: bool = True) -> float:
        """Returns ionization potential.

        Args:
            corrected: Whether to apply empirical correction term

        Returns:
            Ionization potential (eV)
        """
        if self._results.ip is None or self._corrected != corrected:
            self._corrected = corrected
            self._run_xtb("ipea")
            self._results.ip = cast(float, self._results.ip)

        return self._results.ip

    def get_ea(self, corrected: bool = True) -> float:
        """Returns electron affinity.

        Args:
            corrected: Whether to apply empirical correction term

        Returns:
            Electron affinity (eV)
        """
        if self._results.ea is None or self._corrected != corrected:
            self._corrected = corrected
            self._run_xtb("ipea")
            self._results.ea = cast(float, self._results.ea)

        return self._results.ea

    def get_chemical_potential(self, corrected: bool = True) -> float:
        """Returns chemical potential (eV).

        Args:
            corrected: Whether to apply empirical correction term
        """
        if self._results.ip is None or self._corrected != corrected:
            self._corrected = corrected
            self._run_xtb("ipea")
        chem_pot = (
            -(self.get_ip(corrected=corrected) + self.get_ea(corrected=corrected)) / 2
        )

        return chem_pot

    def get_hardness(self, corrected: bool = True) -> float:
        """Returns hardness (eV).

        Args:
            corrected: Whether to apply empirical correction term
        """
        if self._results.ip is None or self._corrected != corrected:
            self._corrected = corrected
            self._run_xtb("ipea")
        hardness = self.get_ip(corrected=corrected) - self.get_ea(corrected=corrected)

        return hardness

    def get_fukui(self, variety: str, corrected: bool = True) -> dict[int, float]:
        """Calculate Fukui coefficients.

        Args:
            variety: Type of Fukui coefficient: 'nucleophilicity' = 'minus', 'electrophilicity' = 'plus',
                'radical', 'dual', 'local_nucleophilicity' or 'local_electrophilicity'.
                Note: 'nucleophilicity' and 'electrophilicity' are synonym for respectively 'minus' and 'plus'.
            corrected: Whether to apply empirical correction term to the inonization potential and electron affinity
                (only applicable for local electrophilicity calculation)

        Returns:
            fukui: Atomic Fukui coefficients

        Raises:
            ValueError: When variety does not exist
        """
        if self._results.fukui_plus is None:
            self._run_xtb("fukui")

        varieties = [
            "minus",
            "nucleophilicity",
            "plus",
            "electrophilicity",
            "radical",
            "dual",
            "local_electrophilicity",
            "local_nucleophilicity",
        ]
        fukui: Array1DFloat
        if variety in ["local_nucleophilicity", "nucleophilicity", "minus"]:
            fukui = self._results.fukui_minus
        elif variety in ["electrophilicity", "plus"]:
            fukui = self._results.fukui_plus
        elif variety == "radical":
            fukui = self._results.fukui_radical
        elif variety == "dual":
            fukui = np.around(
                np.array(self._results.fukui_plus)
                - np.array(self._results.fukui_minus),
                3,
            )
        elif variety == "local_electrophilicity":
            fukui_radical: Array1DFloat = np.array(self._results.fukui_radical)
            fukui_dual: Array1DFloat = np.array(self._results.fukui_plus) - np.array(
                self._results.fukui_minus
            )
            chem_pot = self.get_chemical_potential(corrected=corrected)
            hardness = self.get_hardness(corrected=corrected)
            fukui = (
                -(chem_pot / hardness) * fukui_radical
                + 1 / 2 * (chem_pot / hardness) ** 2 * fukui_dual
            )
        else:
            raise ValueError(
                f"Variety {variety!r} does not exist. "
                f"Choose one of {', '.join(varieties)}."
                f"Note: 'nucleophilicity' and 'electrophilicity' are synonym for respectively 'minus' and 'plus'."
            )

        fukui = {i: float(fukui) for i, fukui in enumerate(fukui, start=1)}

        return fukui

    def get_global_descriptor(self, variety: str, corrected: bool = True) -> float:
        """Calculate global reactivity descriptors.

        Args:
            corrected: Whether to apply empirical correction term
            variety: Type of descriptor: 'electrophilicity', 'nucleophilicity',
                'electrofugality' or 'nucleofugality'

        Returns:
            descriptor: Global reactivity descriptor (eV)

        Raises:
            ValueError: When variety does not exist
        """
        varieties = [
            "electrofugality",
            "electrophilicity",
            "nucleofugality",
            "nucleophilicity",
        ]
        if variety == "electrophilicity":
            descriptor = (
                self.get_ip(corrected=corrected) + self.get_ea(corrected=corrected)
            ) ** 2 / (
                8
                * (self.get_ip(corrected=corrected) - self.get_ea(corrected=corrected))
            )
        elif variety == "nucleophilicity":
            descriptor = -self.get_ip(corrected=corrected)
        elif variety == "electrofugality":
            descriptor = (
                3 * self.get_ip(corrected=corrected) - self.get_ea(corrected=corrected)
            ) ** 2 / (
                8
                * (self.get_ip(corrected=corrected) - self.get_ea(corrected=corrected))
            )
        elif variety == "nucleofugality":
            descriptor = (
                self.get_ip(corrected=corrected) - 3 * self.get_ea(corrected=corrected)
            ) ** 2 / (
                8
                * (self.get_ip(corrected=corrected) - self.get_ea(corrected=corrected))
            )
        else:
            raise ValueError(
                f"Variety {variety!r} does not exist. "
                f"Choose one of {', '.join(varieties)}."
            )

        return descriptor

    def _run_xtb(self, runtype: str) -> None:  # noqa: C901
        """Run xTB calculation and parse results.

        Args:
            runtype: Type of calculation to perform:
                - 'sp': Single point calculation
                - 'ipea': Ionization potential and electron affinity calculation
                - 'fukui': Fukui coefficient calculation
                - 'fod': Fractional occupation density calculation

        Raises:
            ValueError: If runtype does not exist
            ValueError: If the xtb method chosen is not compatible with the requested calculation
            RuntimeError: If the xtb calculation failed
        """
        # Set xtb command
        runtypes = ["sp", "ipea", "fukui", "fod"]
        if runtype == "sp":
            command = self._default_xtb_command
            if self._solvent is not None:
                command += f" --input {XTB._xtb_input_file}"
        elif self._method == "ptb":
            raise ValueError(
                "PTB can only be used for calculations of bond orders, charges, dipole, and HOMO/LUMO energies."
                "\nFor other descriptors, choose another xtb method."
            )
        elif runtype == "ipea":
            command = self._default_xtb_command + " --vipea"
        elif runtype == "fukui":
            command = self._default_xtb_command + " --vfukui"
        elif runtype == "fod":
            command = self._default_xtb_command + " --fod"
        else:
            raise ValueError(
                f"Runtype {runtype!r} does not exist. Choose one of {', '.join(runtypes)}."
            )

        # Create temporary directory or use provided path
        tmp_dir_context = (
            TemporaryDirectory() if self._run_path is None else nullcontext()
        )
        with tmp_dir_context as tmp_dir:
            run_folder = (
                Path(cast(str, tmp_dir))
                if self._run_path is None
                else self._run_path / runtype
            )
            if self._run_path:
                if run_folder.exists():
                    shutil.rmtree(run_folder)
                run_folder.mkdir(parents=True)

            # Write xyz input file
            xyz_file = run_folder / XTB._xyz_input_file
            write_xyz(xyz_file, self._elements, self._coordinates)

            # To
            if self._solvent is not None:
                xtb_inp = run_folder / XTB._xtb_input_file
                write_xtb_inp(xtb_inp, {"write": ["gbsa=true"]})

            # Run xtb
            with open(run_folder / "xtb.out", "w") as stdout, open(
                run_folder / "xtb.err", "w"
            ) as stderr:
                subprocess.run(
                    command.split(), cwd=run_folder, stdout=stdout, stderr=stderr
                )

            # Return error if xtb fails
            with open(run_folder / "xtb.err", "r") as f:
                err_content = f.read()
            if not re.search(r"(?<!ab)normal termination of xtb", err_content):
                with open(run_folder / "xtb.out", "r") as f:
                    out_content = f.read()
                    start_error_idx = out_content.find("######")
                    error = (
                        out_content[start_error_idx:]
                        if start_error_idx != -1
                        else err_content
                    )
                raise RuntimeError(f"xtb calculation failed. Error:\n{error}")

            # Extract results from xtb files
            if runtype == "sp":
                self._parse_json(run_folder / "xtbout.json")
                self._parse_wbo(run_folder / "wbo")
                self._parse_out_sp(run_folder / "xtb.out")
            elif runtype == "ipea":
                self._parse_out_ipea(run_folder / "xtb.out")
            elif runtype == "fukui":
                self._parse_out_fukui(run_folder / "xtb.out")
            elif runtype == "fod":
                self._parse_fod(run_folder / "fod")

    def _parse_json(self, json_file: Path | str) -> None:
        """Parse 'xtbout.json' file."""
        with open(json_file, "r") as f:

            # The code below fixes an error in the json file outputed when running PTB with xtb 6.7.1
            # TODO: Remove/update when new xtb version is released
            lines_without_error = [line for line in f if line.strip() != ","]
            json_fixed = "".join(lines_without_error)
            data = json.loads(json_fixed)

        self._results.charges = data["partial charges"]
        self._results.gap = data["HOMO-LUMO gap / eV"]
        self._results.dipole_vect = np.array(data["dipole / a.u."])
        if self._method != "1":
            self._results.atom_dipole_vect = np.array(data["atomic dipole moments"])

    def _parse_wbo(self, wbo_file: Path | str) -> None:
        """Parse 'wbo' file."""
        wbos = []
        with open(wbo_file, "r") as f:
            for line in f:
                columns = line.split()
                wbos.append((int(columns[0]), int(columns[1]), float(columns[2])))
        self._results.bond_orders = wbos

    def _parse_out_sp(self, out_file: Path | str) -> None:  # noqa: C901
        """Parse 'xtb.out' file from xtb sp calculation."""

        with open(out_file, "r") as f:
            lines = f.readlines()
        homo, lumo, atom_polarizabilities, atom_hb_strengths = {}, {}, [], []
        in_polarizability_block, in_gbsa_block = False, False
        for i, line in enumerate(lines):
            if "(HOMO)" in line:
                homo["Eh"] = float(line.split()[-3])
                homo["eV"] = float(line.split()[-2])
            elif "(LUMO)" in line:
                lumo["Eh"] = float(line.split()[-3])
                lumo["eV"] = float(line.split()[-2])
            elif "dipole" in line:
                if self._method == "2":
                    dipole_line = lines[i + 3].split()
                    dipole_moment = float(dipole_line[-1])
                elif self._method == "1":
                    dipole_line = lines[i + 2].split()
                    dipole_moment = float(dipole_line[-1])
                elif self._method == "ptb":
                    if "Total dipole" in line:
                        dipole_line = lines[i + 1].split()
                        dipole_moment = float(dipole_line[-1])

            elif self._method == "2" and "α(0)" in line:
                if "Mol." in line:
                    mol_polarizability = float(line.split()[-1])
                else:
                    in_polarizability_block = True
            elif in_polarizability_block:
                if line.strip():
                    atom_polarizabilities.append(float(line.split()[-1]))
                else:
                    in_polarizability_block = False

            elif "Gsolv" in line:
                g_solv = float(line.split()[-3])
            elif "Ghb" in line:
                g_solv_hb = float(line.split()[-3])
            elif "H-bond" in line and "correction" not in line:
                in_gbsa_block = True
            elif in_gbsa_block:
                if line.strip():
                    atom_hb_strengths.append(float(line.split()[-1]))
                else:
                    in_gbsa_block = False

        self._results.homo = homo
        self._results.lumo = lumo
        self._results.dipole_moment = dipole_moment
        if self._method == "2":
            self._results.atom_polarizabilities = atom_polarizabilities
            self._results.mol_polarizability = mol_polarizability
        if self._solvent is not None:
            self._results.g_solv = g_solv
            self._results.g_solv_hb = g_solv_hb
            self._results.atom_hb_strengths = atom_hb_strengths

    def _parse_out_ipea(self, out_file: Path | str) -> None:
        """Parse 'xtb.out' file from xtb ipea calculation."""

        with open(out_file, "r") as f:
            ip = ea = shift = None
            for line in f:
                if "delta SCC IP (eV)" in line:
                    ip = float(line.split()[-1])
                elif "delta SCC EA (eV)" in line:
                    ea = float(line.split()[-1])
                elif "empirical IP shift (eV)" in line:
                    shift = float(line.split()[-1])
                if ip and ea and shift:
                    break

        if ip is None or ea is None or shift is None:
            raise ValueError(f"Failed to parse IP and/or EA from {out_file}.")

        if not self._corrected:
            ip += shift
            ea += shift

        self._results.ip = ip
        self._results.ea = ea

    def _parse_out_fukui(self, out_file: Path | str) -> None:
        """Parse 'xtb.out' file from xtb fukui calculation."""

        with open(out_file, "r") as f:
            fukui_plus, fukui_minus, fukui_radical = [], [], []
            in_fukui_block = False
            for line in f:
                if "Fukui functions:" in line:
                    in_fukui_block = True
                if in_fukui_block:
                    if "-------------" in line or not line.strip():
                        in_fukui_block = False
                        break
                    if "Fukui functions" not in line and "#" not in line:
                        columns = line.split()
                        fukui_plus.append(float(columns[-3]))
                        fukui_minus.append(float(columns[-2]))
                        fukui_radical.append(float(columns[-1]))

        if not fukui_plus or not fukui_minus or not fukui_radical:
            raise ValueError(f"Failed to parse Fukui coefficients from {out_file}.")

        self._results.fukui_plus = fukui_plus
        self._results.fukui_minus = fukui_minus
        self._results.fukui_radical = fukui_radical

    def _parse_fod(self, fod_file: Path | str) -> None:
        """Parse 'fod' file."""
        with open(fod_file, "r") as f:
            lines = f.readlines()
        fod = [float(line.strip()) for line in lines]
        self._results.fod_pop = fod


def cli(file: str) -> Any:
    """CLI for XTB.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(XTB, elements, coordinates)
