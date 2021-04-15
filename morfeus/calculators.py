"""Internal calculators."""

import typing
from typing import Dict, Iterable, List, Union

import numpy as np
from scipy.spatial.distance import cdist

from morfeus.d3_data import c6_reference_data, r2_r4
from morfeus.data import ANGSTROM, BOHR, EV, HARTREE
from morfeus.geometry import Atom
from morfeus.typing import Array1D, ArrayLike2D
from morfeus.utils import convert_elements, get_radii, Import, requires_dependency

if typing.TYPE_CHECKING:
    import ase
    from dftd4.calculators import D3_model, D4_model
    from dftd4.utils import extrapolate_c_n_coeff


@requires_dependency(
    [
        Import("ase"),
        Import(module="dftd4.calculators", item="D3_model"),
        Import(module="dftd4.utils", item="extrapolate_c_n_coeff"),
    ],
    globals(),
)
class D3Grimme:
    """Calculates D3-like Cᴬᴬ coefficients with dftd4.

    Args:
        elements : Elements as atomic symbols or numbers
        coordinates : Coordinates (Å)
        order: Maximum order for the CN coefficients.

    Attributes:
        c_n_coefficients: Cᴬᴬ coefficients (a.u.)
        polarizabilities: Atomic polarizabilities (a.u.)
    """

    c_n_coefficients: Dict[int, Array1D]
    polarizabilities: Array1D
    _atoms: List["ase.Atoms"]

    def __init__(
        self,
        elements: Union[Iterable[int], Iterable[str]],
        coordinates: ArrayLike2D,
        order: int = 6,
    ) -> None:
        # Convert elements to atomic numbers
        elements = convert_elements(elements, output="numbers")

        # Do calculation
        atoms = ase.Atoms(numbers=elements, positions=coordinates)
        calc = D3_model(energy=False, forces=False)
        polarizabilities = (
            calc.get_property("polarizibilities", atoms=atoms) * (ANGSTROM / BOHR) ** 3
        )
        c6_coefficients_all = (
            calc.get_property("c6_coefficients", atoms=atoms)
            * EV
            / HARTREE
            * (ANGSTROM / BOHR) ** 6
        )
        c6_coefficients = np.diag(c6_coefficients_all)
        c_n_coefficients = {}
        c_n_coefficients[6] = c6_coefficients

        # Extrapolate
        for i in range(8, order + 1, 2):
            c_n_coefficients[i] = np.array(
                [
                    extrapolate_c_n_coeff(c6, element, element, i)
                    for c6, element in zip(c6_coefficients, elements)
                ]
            )

        # Store attributes
        self.polarizabilities = polarizabilities
        self.c_n_coefficients = c_n_coefficients
        self._atoms = atoms

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


@requires_dependency(
    [
        Import("ase"),
        Import(module="dftd4.calculators", item="D4_model"),
        Import(module="dftd4.utils", item="extrapolate_c_n_coeff"),
    ],
    globals(),
)
class D4Grimme:
    """Calculates D4 Cᴬᴬ dispersion coefficients with dftd4.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        order: Maximum order for the CN coefficients.
        charge: Molecular charge

    Attributes:
        c_n_coefficients: Cᴬᴬ coefficients (a.u.)
        charges: Atomic charges.
        polarizabilities: Atomic polarizabilities (a.u.)
    """

    c_n_coefficients: Dict[int, Array1D]
    charges: Array1D
    polarizabilities: Array1D
    _atoms: List["ase.Atoms"]

    def __init__(
        self,
        elements: Union[Iterable[int], Iterable[str]],
        coordinates: ArrayLike2D,
        order: int = 8,
        charge: int = 0,
    ) -> None:
        # Convert elements to atomic numbers
        elements = convert_elements(elements, output="numbers")

        # Set up atoms object
        charges = np.zeros(len(elements))
        charges[0] = charge
        atoms = atoms = ase.Atoms(
            numbers=elements, positions=coordinates, charges=charges
        )

        # Do calculation
        calc = D4_model(energy=False, forces=False)
        polarizabilities = (
            calc.get_property("polarizibilities", atoms=atoms) * (ANGSTROM / BOHR) ** 3
        )
        charges = calc.get_property("charges", atoms=atoms)
        c6_coefficients_all = (
            calc.get_property("c6_coefficients", atoms=atoms)
            * EV
            / HARTREE
            * (ANGSTROM / BOHR) ** 6
        )
        c6_coefficients = np.diag(c6_coefficients_all)
        c_n_coefficients = {}
        c_n_coefficients[6] = c6_coefficients

        # Extrapolate
        for i in range(8, order + 1, 2):
            c_n_coefficients[i] = np.array(
                [
                    extrapolate_c_n_coeff(c6, element, element, i)
                    for c6, element in zip(c6_coefficients, elements)
                ]
            )

        # Store attributes
        self.polarizabilities = polarizabilities
        self.charges = charges
        self.c_n_coefficients = c_n_coefficients
        self._atoms = atoms

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


class D3Calculator:
    """Calculates D3 Cᴬᴬ coefficients.

    Procedure as described in J. Chem. Phys. 2010, 132, 154104.

    Args:
        elements : Elements as atomic symbols or numbers
        coordinates : Coordinates (Å)
        order: Maximum order for the CN coefficients.

    Attributes:
        c_n_coefficients: Cᴬᴬ coefficients (a.u.)
        coordination_numbers: Atomic coordination numbers.
    """

    c_n_coefficients: Dict[int, Array1D]
    coordination_numbers: Array1D
    _atoms: List[Atom]

    def __init__(
        self,
        elements: Union[Iterable[int], Iterable[str]],
        coordinates: ArrayLike2D,
        order: int = 8,
    ) -> None:
        # Convert elements to atomic numbers
        elements = convert_elements(elements, output="numbers")
        coordinates = np.array(coordinates)

        # Load the covalent radii
        radii = get_radii(elements, radii_type="pyykko")

        # Set up the atoms objects
        atoms = []
        for i, (element, coordinate, radius) in enumerate(
            zip(elements, coordinates, radii), start=1
        ):
            atom = Atom(element, coordinate, radius, i)
            atoms.append(atom)
        self._atoms = atoms

        # Calculate the coordination numbers according to Grimme's recipe.
        k_1 = 16
        k_2 = 4 / 3
        k_3 = 4
        for cn_atom in atoms:
            # Get coordinates and radii of all other atoms
            other_coordinates = np.array(
                [atom.coordinates for atom in atoms if atom is not cn_atom]
            )
            other_radii = np.array(
                [atom.radius for atom in atoms if atom is not cn_atom]
            )

            # Calculate distances
            dists = cdist(
                np.array(cn_atom.coordinates).reshape(-1, 3),
                other_coordinates.reshape(-1, 3),
            )

            # Calculate coordination number
            coordination_number = np.sum(
                1
                / (
                    1
                    + np.exp(-k_1 * (k_2 * (cn_atom.radius + other_radii) / dists - 1))
                )
            )
            cn_atom.coordination_number = coordination_number

        # Calculate the C_N coefficients
        c_n_coefficients: Dict[int, List[float]] = {
            i: [] for i in range(6, order + 1, 2)
        }
        for atom in atoms:
            # Get the reference data for atom
            reference_data = c6_reference_data[atom.element]
            cn = atom.coordination_number

            # Take out the coordination numbers and c6(aa) values
            c_6_ref = reference_data[:, 0]
            cn_1 = reference_data[:, 1]
            cn_2 = reference_data[:, 2]

            # Calculate c6  and c8 according to the recipe
            r = (cn - cn_1) ** 2 + (cn - cn_2) ** 2
            L = np.exp(-k_3 * r)
            W = np.sum(L)
            Z = np.sum(c_6_ref * L)
            c_6 = Z / W

            temp_coefficients: Dict[int, float] = {
                i: np.nan for i in range(6, order + 1, 2)
            }
            for i in range(6, order + 1, 2):
                if i == 6:
                    temp_coefficients[i] = c_6
                if i == 8:
                    c_8 = 3 * c_6 * r2_r4[atom.element] ** 2
                    temp_coefficients[i] = c_8
                elif i == 10:
                    c_10 = 49.0 / 40.0 * c_8 ** 2 / c_6
                    temp_coefficients[i] = c_10
                elif i > 10:
                    c_n = (
                        temp_coefficients[i - 6]
                        * (temp_coefficients[i - 2] / temp_coefficients[i - 4]) ** 3
                    )
                    temp_coefficients[i] = c_n
            for key, value in temp_coefficients.items():
                c_n_coefficients[key].append(value)

        # Set up attributes
        coordination_numbers = np.array(
            [atom.coordination_number for atom in self._atoms]
        )
        c_n_coefficients = {
            key: np.array(value) for key, value in c_n_coefficients.items()
        }

        self._atoms = atoms
        self.c_n_coefficients = c_n_coefficients
        self.coordination_numbers = coordination_numbers

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"
