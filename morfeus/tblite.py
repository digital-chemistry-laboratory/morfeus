"""Native tblite interface code."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
import typing
from typing import Any, cast

import numpy as np

from morfeus.data import (
    ANGSTROM_TO_BOHR,
    AU_TO_DEBYE,
    HARTREE,
    HARTREE_TO_EV,
    K_B,
)
from morfeus.typing import Array1DFloat, Array2DFloat, ArrayLike2D
from morfeus.utils import convert_elements, Import, requires_dependency

if typing.TYPE_CHECKING:
    from tblite.interface import Calculator


class TBLite:
    """Calculate xTB properties with the native tblite Python interface.

    This class intentionally exposes only the subset of properties provided by
    tblite result containers. Use :class:`morfeus.xtb.XTB` for descriptors that
    depend on xtb CLI workflows such as Fukui functions, IP/EA corrections, FOD,
    or detailed solvation terms.

    Args:
        elements: Elements as atomic symbols or numbers.
        coordinates: Coordinates (A).
        method: xTB method. Currently supports 1/GFN1-xTB, 2/GFN2-xTB,
            and IPEA1-xTB.
        charge: Molecular charge.
        n_unpaired: Number of unpaired electrons.
        solvent: Optional solvent for implicit solvation.
        solvation_model: Solvation model, either "alpb" or "gbsa".
        electronic_temperature: Electronic temperature (K).
    """

    _method_map: dict[str, str] = {
        "1": "GFN1-xTB",
        "gfn1": "GFN1-xTB",
        "gfn1-xtb": "GFN1-xTB",
        "2": "GFN2-xTB",
        "gfn2": "GFN2-xTB",
        "gfn2-xtb": "GFN2-xTB",
        "ipea1": "IPEA1-xTB",
        "ipea1-xtb": "IPEA1-xTB",
    }

    _solvation_models: dict[str, str] = {
        "alpb": "alpb-solvation",
        "gbsa": "gbsa-solvation",
    }

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        method: int | str = 2,
        charge: int = 0,
        n_unpaired: int | None = None,
        solvent: str | None = None,
        solvation_model: str = "alpb",
        electronic_temperature: float | None = None,
    ) -> None:
        self._elements = np.array(convert_elements(elements, output="numbers"))
        self._coordinates = np.array(coordinates)
        self._charge = charge
        self._n_unpaired = n_unpaired
        self._solvent = solvent
        self._electronic_temperature = electronic_temperature

        method_key = str(method).lower()
        try:
            self._method = self._method_map[method_key]
        except KeyError as error:
            supported_methods = ", ".join(sorted(self._method_map))
            raise ValueError(
                f"Method {method!r} not supported. Choose one of: {supported_methods}."
            ) from error

        solvation_key = solvation_model.lower()
        try:
            self._solvation = self._solvation_models[solvation_key]
        except KeyError as error:
            supported_solvation = ", ".join(sorted(self._solvation_models))
            raise ValueError(
                "Solvation model "
                f"{solvation_model!r} not supported. Choose one of: "
                f"{supported_solvation}."
            ) from error

        self._results = TBLiteResults()

    def get_bond_order(self, i: int, j: int) -> float:
        """Return bond order between two atoms.

        Args:
            i: Index of atom 1 (1-indexed).
            j: Index of atom 2 (1-indexed).

        Returns:
            Bond order between atoms i and j.

        Raises:
            ValueError: If no bond order exists between the given atoms.
        """
        bond_orders = self.get_bond_orders()
        if (i, j) in bond_orders:
            bond_order = bond_orders[(i, j)]
        elif (j, i) in bond_orders:
            bond_order = bond_orders[(j, i)]
        else:
            raise ValueError(f"No bond order calculated between atoms {i} and {j}.")

        return bond_order

    def get_bond_orders(self) -> dict[tuple[int, int], float]:
        """Return bond orders."""
        if self._results.bond_orders is None:
            self._run_tblite()
            self._results.bond_orders = cast(Array2DFloat, self._results.bond_orders)

        bond_order_matrix = self._results.bond_orders
        bond_orders = {}
        rows, cols = np.triu_indices_from(bond_order_matrix, k=1)
        for row, col in zip(rows, cols):
            bond_order = float(bond_order_matrix[row, col])
            if not np.isclose(bond_order, 0):
                bond_orders[(row + 1, col + 1)] = round(bond_order, 12)

        return bond_orders

    def get_charges(self) -> dict[int, float]:
        """Return atomic charges (1-indexed)."""
        if self._results.charges is None:
            self._run_tblite()
            self._results.charges = cast(Array1DFloat, self._results.charges)

        return {
            i: float(charge) for i, charge in enumerate(self._results.charges, start=1)
        }

    def get_dipole(self) -> Array1DFloat:
        """Return molecular dipole vector (a.u.)."""
        if self._results.dipole is None:
            self._run_tblite()
            self._results.dipole = cast(Array1DFloat, self._results.dipole)

        return self._results.dipole

    def get_dipole_moment(self, unit: str = "au") -> float:
        """Return molecular dipole moment.

        Args:
            unit: 'au' or 'debye'.

        Returns:
            Molecular dipole moment (a.u. or debye).

        Raises:
            ValueError: If given unit is not supported.
        """
        dipole_moment = float(np.linalg.norm(self.get_dipole()))
        if unit == "au":
            return round(dipole_moment, 8)
        elif unit == "debye":
            return round(dipole_moment * AU_TO_DEBYE, 8)
        else:
            raise ValueError("Unit must be either 'au' or 'debye'.")

    def get_energy(self) -> float:
        """Return total energy (Eh)."""
        if self._results.total_energy is None:
            self._run_tblite()
            self._results.total_energy = cast(float, self._results.total_energy)

        return self._results.total_energy

    def get_homo(self, unit: str = "Eh") -> float:
        """Return HOMO energy."""
        homo, _ = self._get_frontier_orbitals()

        return self._convert_orbital_energy(homo, unit)

    def get_homo_lumo_gap(self, unit: str = "eV") -> float:
        """Return HOMO-LUMO gap."""
        homo, lumo = self._get_frontier_orbitals()
        gap = lumo - homo

        return self._convert_orbital_energy(gap, unit)

    def get_lumo(self, unit: str = "Eh") -> float:
        """Return LUMO energy."""
        _, lumo = self._get_frontier_orbitals()

        return self._convert_orbital_energy(lumo, unit)

    def get_ea(self, corrected: bool = True) -> float:
        """Raise for unsupported electron affinity calculations."""
        raise NotImplementedError("Electron affinity requires the xtb CLI backend.")

    def get_fukui(self, variety: str, corrected: bool = True) -> dict[int, float]:
        """Raise for unsupported Fukui function calculations."""
        raise NotImplementedError("Fukui functions require the xtb CLI backend.")

    def get_ip(self, corrected: bool = True) -> float:
        """Raise for unsupported ionization potential calculations."""
        raise NotImplementedError("Ionization potential requires the xtb CLI backend.")

    def _convert_orbital_energy(self, energy: float, unit: str) -> float:
        """Convert orbital energy from Hartree."""
        if unit == "Eh":
            return energy
        elif unit == "eV":
            return energy * HARTREE_TO_EV
        else:
            raise ValueError("Unit must be either 'Eh' or 'eV'.")

    def _get_frontier_orbitals(self) -> tuple[float, float]:
        """Return HOMO and LUMO energies in Hartree."""
        if (
            self._results.orbital_energies is None
            or self._results.orbital_occupations is None
        ):
            self._run_tblite()
            self._results.orbital_energies = cast(
                Array1DFloat, self._results.orbital_energies
            )
            self._results.orbital_occupations = cast(
                Array1DFloat, self._results.orbital_occupations
            )

        energies = np.ravel(self._results.orbital_energies)
        occupations = np.ravel(self._results.orbital_occupations)
        order = np.argsort(energies)
        energies = energies[order]
        occupations = occupations[order]

        occupied = np.where(occupations > 1e-8)[0]
        unoccupied = np.where(occupations <= 1e-8)[0]
        if len(occupied) == 0 or len(unoccupied) == 0:
            raise ValueError("Could not determine frontier orbitals from tblite results.")

        homo_index = occupied[-1]
        lumo_candidates = unoccupied[unoccupied > homo_index]
        if len(lumo_candidates) == 0:
            raise ValueError("Could not determine LUMO from tblite results.")

        return float(energies[homo_index]), float(energies[lumo_candidates[0]])

    @requires_dependency(
        [Import(module="tblite.interface", item="Calculator")], globals()
    )
    def _run_tblite(self) -> None:
        """Run tblite single-point calculation and cache supported results."""
        calculator = Calculator(
            self._method,
            self._elements,
            self._coordinates * ANGSTROM_TO_BOHR,
            charge=float(self._charge),
            uhf=self._n_unpaired,
        )
        calculator.set("verbosity", 0)
        if self._electronic_temperature is not None:
            calculator.set(
                "temperature", float(self._electronic_temperature) * K_B / HARTREE
            )
        if self._solvent is not None:
            calculator.add(self._solvation, self._solvent)

        result = calculator.singlepoint()
        self._extract_result(result)

    def _extract_result(self, result: Any) -> None:
        """Extract supported quantities from a tblite result container."""
        self._results.total_energy = float(result.get("energy"))
        self._results.charges = np.array(result.get("charges"))
        self._results.bond_orders = np.array(result.get("bond-orders"))
        self._results.dipole = np.array(result.get("dipole"))
        self._results.orbital_energies = np.array(result.get("orbital-energies"))
        self._results.orbital_occupations = np.array(
            result.get("orbital-occupations")
        )


@dataclass
class TBLiteResults:
    """Stores native tblite results."""

    charges: Array1DFloat | None = None
    bond_orders: Array2DFloat | None = None
    dipole: Array1DFloat | None = None
    total_energy: float | None = None
    orbital_energies: Array1DFloat | None = None
    orbital_occupations: Array1DFloat | None = None
