"""xtb interface code."""

from __future__ import annotations

from collections.abc import Iterable
import functools

# import typing
from typing import Any
import numpy as np
from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import nullcontext
import subprocess
import shutil
from dataclasses import dataclass
import re

from morfeus.data import DEBYE_TO_AU  # , ANGSTROM_TO_BOHR, HARTREE_TO_EV
from morfeus.io import read_geometry, write_xyz
from morfeus.typing import Array1DFloat, Array2DFloat, ArrayLike2D
from morfeus.utils import convert_elements  # , Import, requires_dependency

# if typing.TYPE_CHECKING:
#     import xtb
#     import xtb.interface
#     import xtb.utils

# IPEA_CORRECTIONS = {"1": 5.700, "2": 4.846}


# @requires_dependency(
#     [Import("xtb"), Import("xtb.interface"), Import("xtb.utils")], globals()
# )
@dataclass
class XTBResults:
    """Stores xTB descriptors."""

    charges: list[float] = None
    bond_orders: list[tuple[int, int, float]] = None
    homo: dict[str, float] = None  # Stores both Eh and eV units
    lumo: dict[str, float] = None  # Stores both Eh and eV units
    dipole_vect: Array1DFloat = None  # Unit in a.u.
    dipole_moment: float = None  # Unit in debye
    ip: float = None
    ea: float = None
    fukui_plus: list[float] = None
    fukui_minus: list[float] = None
    fukui_radical: list[float] = None


class XTB:
    """Calculates electronic properties with the xtb program.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        version: Version of xtb to use. Currently works with '1' or '2'.
        charge: Molecular charge
        n_unpaired: Number of unpaired electrons
        solvent: Solvent. Uses the ALPB solvation model
        electronic_temperature: Electronic temperature (K)
        run_path: Folder path to run xTB calculation. If not provided, runs in a temporary folder
    """

    _charge: int
    _coordinates: Array2DFloat
    _electronic_temperature: int | None
    _elements: Array1DFloat
    _n_unpaired: int | None
    _results: Any
    _solvent: str | None
    _version: int

    _xyz_input: str = "xtb.xyz"

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        version: int | str | None = 2,
        charge: int | None = 0,
        n_unpaired: int | None = None,
        solvent: str | None = None,
        electronic_temperature: int | None = None,
        run_path: Path | str | None = None,
    ) -> None:
        # Converting elements to atomic numbers if the are symbols
        self._elements = np.array(convert_elements(elements, output="numbers"))

        # Store settings
        self._coordinates = np.array(coordinates)
        self._version = version
        self._charge = charge
        self._solvent = solvent
        self._n_unpaired = n_unpaired
        self._electronic_temperature = electronic_temperature

        self._run_path = Path(run_path) if run_path else None

        self._default_xtb_command = (
            f"xtb {XTB._xyz_input} --gfn {self._version} --chrg {self._charge}"
        )
        if self._solvent is not None:
            self._default_xtb_command += f" --alpb {self._solvent}"
        if self._n_unpaired is not None:
            self._default_xtb_command += f" --uhf {self._n_unpaired}"
        if self._electronic_temperature is not None:
            self._default_xtb_command += f" --etemp {self._electronic_temperature}"

        self._results = XTBResults()
        self._corrected: bool | None = None

    def get_bond_order(self, i: int, j: int) -> float:
        """Returns bond order between two atoms.

        Args:
            i: Index of atom 1 (1-indexed)
            j: Index of atom 2 (1-indexed)

        Returns:
            bond_order: Bond order
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
        bond_orders = {(x[0], x[1]): x[2] for x in self._results.bond_orders}

        return bond_orders

    def get_charges(self) -> dict[int, float]:
        """Returns atomic charges."""

        if self._results.charges is None:
            self._run_xtb("sp")
        charges = {i: charge for i, charge in enumerate(self._results.charges, start=1)}

        return charges

    def get_homo(self, unit="Eh") -> float:
        """Returns HOMO energy.

        Args:
            unit: 'Eh' or 'eV'

        Returns:
            HOMO energy (Eh or eV)

        Raises:
            ValueError: When unit does not exist
        """

        if self._results.homo is None:
            self._run_xtb("sp")

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
            ValueError: When unit does not exist
        """

        if self._results.lumo is None:
            self._run_xtb("sp")

        if unit == "Eh":
            return self._results.lumo["Eh"]
        elif unit == "eV":
            return self._results.lumo["eV"]
        else:
            raise ValueError("Unit must be either 'Eh' or 'eV'.")

    def get_dipole(self) -> Array1DFloat:
        """Returns molecular dipole vector (a.u.)."""

        if self._results.dipole_vect is None:
            self._run_xtb("sp")

        return self._results.dipole_vect

    def get_dipole_moment(self, unit="deybe") -> float:
        """Returns molecular dipole moment.

        Args:
            unit: 'deybe' or 'au'

        Returns:
            Molecular dipole moment (deybe or a.u.)

        Raises:
            ValueError: When unit does not exist
        """

        if self._results.dipole_moment is None:
            self._run_xtb("sp")

        if unit == "deybe":
            return self._results.dipole_moment
        elif unit == "au":
            return round(
                self._results.dipole_moment * DEBYE_TO_AU,
                len(str(self._results.dipole_moment).split(".")[-1]),
            )
        else:
            raise ValueError("Unit must be either 'deybe' or 'au'.")

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

        return self._results.ea

    def get_chem_pot(self, corrected: bool = True) -> float:
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
            fukui = list(
                np.array(self._results.fukui_plus) - np.array(self._results.fukui_minus)
            )
        elif variety == "local_electrophilicity":
            fukui_radical: Array1DFloat = np.array(self._results.fukui_radical)
            fukui_dual: Array1DFloat = np.array(self._results.fukui_plus) - np.array(
                self._results.fukui_minus
            )
            chem_pot = self.get_chem_pot(corrected=corrected)
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

    def _run_xtb(self, runtype: str) -> None:
        """Run xTB calculation and parse results.

        Args:
            runtype: Type of calculation to perform: 'sp', 'ipea' or 'fukui'
                'sp': Single point calculation
                'ipea': Ionization potential and electron affinity calculation
                'fukui': Fukui coefficient calculation
        """
        # Set xtb command
        runtypes = ["sp", "ipea", "fukui"]
        if runtype == "sp":
            command = self._default_xtb_command
        elif runtype == "ipea":
            command = self._default_xtb_command + " --vipea"
        elif runtype == "fukui":
            command = self._default_xtb_command + " --vfukui"
        else:
            raise ValueError(
                f"Runtype {runtype!r} does not exist. Choose one of {', '.join(runtypes)}."
            )

        # Create temporary directory or use provided path
        tmp_dir_context = (
            TemporaryDirectory() if self._run_path is None else nullcontext()
        )
        with tmp_dir_context as tmp_dir:
            run_folder = Path(tmp_dir) if self._run_path is None else self._run_path
            if self._run_path:
                if run_folder.exists():
                    shutil.rmtree(run_folder)
                run_folder.mkdir(parents=True)

            # Write xyz input file
            xyz_file = run_folder / XTB._xyz_input
            write_xyz(xyz_file, self._elements, self._coordinates)

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
                raise RuntimeError(f"xTB calculation failed. Error:\n{error}")

            # Extract results from xtb files
            if runtype == "sp":
                self._parse_charges(run_folder / "charges")
                self._parse_wbo(run_folder / "wbo")
                self._parse_out_sp(run_folder / "xtb.out")
            elif runtype == "ipea":
                self._parse_out_ipea(run_folder / "xtb.out")
            elif runtype == "fukui":
                self._parse_out_fukui(run_folder / "xtb.out")

    def _parse_charges(self, charges_file: Path | str) -> None:
        """Parse 'charges' file."""
        with open(charges_file, "r") as f:
            lines = f.readlines()
        charges = [float(line.strip()) for line in lines]
        self._results.charges = charges

    def _parse_wbo(self, wbo_file: Path | str) -> None:
        """Parse 'wbo' file."""
        wbos = []
        with open(wbo_file, "r") as f:
            for line in f:
                columns = line.split()
                wbos.append((int(columns[0]), int(columns[1]), float(columns[2])))
        self._results.bond_orders = wbos

    def _parse_out_sp(self, out_file: Path | str) -> None:
        """Parse 'xtb.out' file from xtb sp calculation."""
        with open(out_file, "r") as f:
            lines = f.readlines()
        homo, lumo, dipole_vect, dipole_moment = {}, {}, None, None
        for i, line in enumerate(lines):
            if "(HOMO)" in line:
                homo["Eh"] = float(line.split()[-3])
                homo["eV"] = float(line.split()[-2])
            elif "(LUMO)" in line:
                lumo["Eh"] = float(line.split()[-3])
                lumo["eV"] = float(line.split()[-2])
            elif "dipole" in line:
                if self._version == 2:
                    dipole_line = lines[i + 3].split()
                    dipole_vect = np.array(
                        [
                            float(dipole_line[-4]),
                            float(dipole_line[-3]),
                            float(dipole_line[-2]),
                        ]
                    )
                    dipole_moment = float(dipole_line[-1])
                elif self._version == 1:
                    dipole_line = lines[i + 2].split()
                    dipole_vect = np.array(
                        [
                            float(dipole_line[0]),
                            float(dipole_line[1]),
                            float(dipole_line[2]),
                        ]
                    )
                    dipole_moment = float(dipole_line[-1])
            if homo and lumo and dipole_moment:
                break
        self._results.homo = homo
        self._results.lumo = lumo
        self._results.dipole_vect = dipole_vect
        self._results.dipole_moment = dipole_moment

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
                    if "-------------" in line:
                        in_fukui_block = False
                        break
                    if "Fukui functions" not in line and "#" not in line:
                        columns = line.split()
                        fukui_plus.append(float(columns[-3]))
                        fukui_minus.append(float(columns[-2]))
                        fukui_radical.append(float(columns[-1]))
        self._results.fukui_plus = fukui_plus
        self._results.fukui_minus = fukui_minus
        self._results.fukui_radical = fukui_radical


def cli(file: str) -> Any:
    """CLI for XTB.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(XTB, elements, coordinates)
