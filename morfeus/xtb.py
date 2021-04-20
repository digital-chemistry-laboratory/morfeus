"""xtb interface code."""

import functools
import typing
from typing import Any, Dict, Iterable, Optional, Union

import numpy as np

from morfeus.data import ANGSTROM_TO_BOHR, HARTREE_TO_EV
from morfeus.io import read_geometry
from morfeus.typing import Array1D, ArrayLike2D
from morfeus.utils import convert_elements, Import, requires_dependency

if typing.TYPE_CHECKING:
    import xtb
    import xtb.interface

IPEA_CORRECTIONS = {"1": 5.700, "2": 4.846}


@requires_dependency([Import("xtb"), Import("xtb.interface")], globals())
class XTB:
    """Calculates electronic properties with the xtb-python package.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        version: Version of xtb to use. Currently works with '1' or '2'.
        charge: Molecular charge
        n_unpaired: Number of unpaired electrons
        solvent: Solvent. See xtb-python documentation
        electronic_temperature: Electronic temperature (K)
    """

    _charge: int
    _coordinates: np.ndarray
    _electronic_temperature: Optional[int]
    _elements: np.ndarray
    _n_unpaired: Optional[int]
    _results: Any
    _solvent: Optional[str]
    _version: str

    def __init__(
        self,
        elements: Union[Iterable[int], Iterable[str]],
        coordinates: ArrayLike2D,
        version: str = "2",
        charge: int = 0,
        n_unpaired: Optional[int] = None,
        solvent: Optional[str] = None,
        electronic_temperature: Optional[int] = None,
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

        # Set up results dictionary
        self._results = {
            -1: None,
            0: None,
            1: None,
        }

        self._params = {
            "1": xtb.interface.Param.GFN1xTB,
            "2": xtb.interface.Param.GFN2xTB,
        }

    def get_bond_order(self, i: int, j: int) -> float:
        """Returns bond order between two atoms.

        Args:
            i: Index of atom 1 (1-indexed)
            j: Index of atom 2 (1-indexed)

        Returns:
            bond_order: Bond order
        """
        bo_matrix = self.get_bond_orders()
        bond_order: float = bo_matrix[i - 1, j - 1]

        return bond_order

    def get_bond_orders(self, charge_state: int = 0) -> Array1D:
        """Returns bond orders.

        Args:
            charge_state: Charge state relative to original charge

        Returns:
            bond_orders: Bond orders
        """
        self._check_results(charge_state)
        bond_orders: np.ndarray = self._results[charge_state].get_bond_orders()

        return bond_orders

    def _get_charges(self, charge_state: int = 0) -> Array1D:
        """Returns atomic charges."""
        self._check_results(charge_state)
        charges: np.ndarray = self._results[charge_state].get_charges()

        return charges

    def get_charges(self, charge_state: int = 0) -> Dict[int, float]:
        """Returns atomic charges.

        Args:
            charge_state: Charge state relative to original charge

        Returns:
            charges: Atomic charges
        """
        self._check_results(charge_state)
        charges = self._results[charge_state].get_charges()
        charges = {i: charge for i, charge in enumerate(charges, start=1)}

        return charges

    def get_dipole(self, charge_state: int = 0) -> Array1D:
        """Calculate dipole vector (a.u.).

        Args:
            charge_state: Charge state relative to original charge

        Returns:
            dipole: Dipole vector
        """
        self._check_results(charge_state)
        dipole: np.ndarray = self._results[charge_state].get_dipole()

        return dipole

    def get_ea(self, corrected: bool = False) -> float:
        """Calculate electron affinity.

        Args:
            corrected: Whether to apply correction term

        Returns:
            ea: Electron affinity (eV)
        """
        # Calculate energies
        energy_neutral = self._get_energy(0)
        energy_anion = self._get_energy(-1)

        # Calculate electron affinity
        ea = (energy_neutral - energy_anion) * HARTREE_TO_EV
        if corrected:
            ea -= IPEA_CORRECTIONS[self._version]

        return ea

    def get_fukui(self, variety: str) -> Dict[int, float]:
        """Calculate Fukui coefficients.

        Args:
            variety: Type of Fukui coefficient: 'nucleophilicity', 'electrophilicity',
                'radical', 'dual', 'local_nucleophilicity' or 'local_electrophilicity'.

        Returns:
            fukui: Atomic Fukui coefficients

        Raises:
            ValueError: When variety does not exist
        """
        varieties = [
            "dual",
            "electrophilicity",
            "local_electrophilicity",
            "local_nucleophilicity",
            "nucleophilicity",
            "radical",
        ]
        fukui: np.ndarray
        if variety in ["local_nucleophilicity", "nucleophilicity"]:
            fukui = self._get_charges(0) - self._get_charges(1)
        elif variety == "electrophilicity":
            fukui = self._get_charges(-1) - self._get_charges(0)
        elif variety == "radical":
            fukui = (self._get_charges(-1) - self._get_charges(0)) / 2
        elif variety == "dual":
            fukui = (
                2 * self._get_charges(0) - self._get_charges(1) - self._get_charges(-1)
            )
        elif variety == "local_electrophilicity":
            fukui_radical = np.array(list(self.get_fukui("radical").values()))
            fukui_dual = np.array(list(self.get_fukui("dual").values()))
            chem_pot = -(self.get_ip() + self.get_ea()) / 2
            hardness = self.get_ip() - self.get_ea()
            fukui = (
                -(chem_pot / hardness) * fukui_radical
                + 1 / 2 * (chem_pot / hardness) ** 2 * fukui_dual
            )
        else:
            raise ValueError(
                f"Variety '{variety}' does not exist. "
                f"Choose one of {', '.join(varieties)}."
            )

        fukui = {i: fukui for i, fukui in enumerate(fukui, start=1)}

        return fukui

    def get_global_descriptor(self, variety: str, corrected: bool = False) -> float:
        """Calculate global reactivity descriptors.

        Args:
            corrected: Whether to apply correction term
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
                f"Variety '{variety}' does not exist. "
                f"Choose one of {', '.join(varieties)}."
            )

        return descriptor

    def get_homo(self) -> float:
        """Calculate HOMO energy.

        Returns:
            homo_energy: HOMO energy (a.u.)
        """
        eigenvalues = self._get_eigenvalues()
        occupations = self._get_occupations()
        homo_index = int(occupations.sum().round(0) / 2) - 1
        homo_energy: float = eigenvalues[homo_index]

        return homo_energy

    def get_ip(self, corrected: bool = False) -> float:
        """Calculate ionization potential.

        Args:
            corrected: Whether to apply correction term

        Returns:
            ip: Ionization potential (eV)
        """
        # Calculate energies
        energy_neutral = self._get_energy(0)
        energy_cation = self._get_energy(1)

        # Calculate ionization potential
        ip = (energy_cation - energy_neutral) * HARTREE_TO_EV
        if corrected:
            ip -= IPEA_CORRECTIONS[self._version]

        return ip

    def get_lumo(self) -> float:
        """Calculate LUMO energy.

        Returns:
            lumo_energy: LUMO energy (a.u.)
        """
        eigenvalues = self._get_eigenvalues()
        occupations = self._get_occupations()
        homo_index = int(occupations.sum().round(0) / 2) - 1
        lumo_index = homo_index + 1
        lumo_energy: float = eigenvalues[lumo_index]

        return lumo_energy

    def _check_results(self, charge_state: int) -> None:
        """Checks whether results are already calculated and does calculation."""
        if self._results[charge_state] is None:
            self._sp(charge_state)

    def _get_eigenvalues(self) -> Array1D:
        """Get orbital eigenvalues."""
        self._check_results(0)
        eigenvalues: np.ndarray = self._results[0].get_orbital_eigenvalues()
        return eigenvalues

    def _get_energy(self, charge_state: int = 0) -> float:
        """Get total energy."""
        self._check_results(charge_state)
        energy: float = self._results[charge_state].get_energy()
        return energy

    def _get_occupations(self) -> Array1D:
        """Get occupation numbers."""
        self._check_results(0)
        occupations: np.ndarray = self._results[0].get_orbital_occupations()
        return occupations

    def _sp(self, charge_state: int = 0) -> None:
        """Perform single point calculation."""
        # Set up calculator
        calc = xtb.interface.Calculator(
            self._params[self._version],
            self._elements,
            self._coordinates * ANGSTROM_TO_BOHR,
            charge=self._charge + charge_state,
            uhf=self._n_unpaired,
        )
        calc.set_verbosity(xtb.libxtb.VERBOSITY_MUTED)

        # Set solvent
        if self._solvent:
            solvent = xtb.utils.get_solvent(self._solvent)
            if solvent is None:
                raise Exception(f"Solvent '{self._solvent}' not recognized")
            calc.set_solvent(solvent)

        # Set electronic temperature
        if self._electronic_temperature:
            calc.set_electronic_temperature(self._electronic_temperature)

        # Do singlepoint calculation and store the result
        res = calc.singlepoint()
        self._results[charge_state] = res


def cli(file: str) -> Any:
    """CLI for XTB.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(XTB, elements, coordinates)
