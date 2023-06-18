"""Internal calculators."""

from __future__ import annotations

from collections.abc import Iterable
import typing

import numpy as np
from scipy.spatial.distance import cdist

from morfeus.d3_data import c6_reference_data, r2_r4
from morfeus.data import ANGSTROM_TO_BOHR
from morfeus.geometry import Atom
from morfeus.typing import Array1DFloat, Array1DInt, Array2DFloat, ArrayLike2D
from morfeus.utils import convert_elements, get_radii, Import, requires_dependency

if typing.TYPE_CHECKING:
    from dftd4.interface import DispersionModel


@requires_dependency(
    [
        Import(module="dftd4.interface", item="DispersionModel"),
    ],
    globals(),
)
class D4Grimme:
    """Calculates D4 Cᴬᴬ coefficients with dftd4.

    Args:
        elements : Elements as atomic symbols or numbers
        coordinates : Coordinates (Å)
        order: Maximum order for the CN coefficients.
        charge: Molecular charge

    Attributes:
        charges: Partial charges
        c_n_coefficients: Cᴬᴬ coefficients (a.u.)
        coordination_numbers: Coordination numbers
        polarizabilities: Atomic polarizabilities (a.u.)
    """

    charges: dict[int, float]
    c_n_coefficients: dict[int, Array1DFloat]
    coordination_numbers: dict[int, float]
    polarizabilities: dict[int, float]

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        charge: int = 0,
        order: int = 8,
    ) -> None:
        # Convert elements to atomic numbers
        elements: Array1DInt = np.array(convert_elements(elements, output="numbers"))
        coordinates = np.array(coordinates) * ANGSTROM_TO_BOHR

        # Do calculation
        calc = DispersionModel(elements, coordinates, charge=charge)
        properties = calc.get_properties()
        polarizabilities = properties["polarizibilities"]
        charges = properties["partial charges"]
        coordination_numbers = properties["coordination numbers"]
        c6_coefficients_all = properties["c6 coefficients"]
        c6_coefficients: Array1DFloat = np.diag(c6_coefficients_all)
        c_n_coefficients = {}
        c_n_coefficients[6] = c6_coefficients

        # Extrapolate
        for i in range(8, order + 1, 2):
            c_n_coefficients[i] = np.array(
                [
                    extrapolate_c_n(c6, element, element, i)
                    for c6, element in zip(c6_coefficients, elements)
                ]
            )

        # Store attributes
        self.polarizabilities = polarizabilities
        self.charges = charges
        self.c_n_coefficients = c_n_coefficients
        self.coordination_numbers = coordination_numbers

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self.polarizabilities)!r} atoms)"


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

    c_n_coefficients: dict[int, Array1DFloat]
    coordination_numbers: Array1DFloat
    _atoms: list[Atom]

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        order: int = 8,
    ) -> None:
        # Convert elements to atomic numbers
        elements = convert_elements(elements, output="numbers")
        coordinates: Array2DFloat = np.array(coordinates)

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
            other_coordinates: Array2DFloat = np.array(
                [atom.coordinates for atom in atoms if atom is not cn_atom]
            )
            other_radii: Array1DFloat = np.array(
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
        c_n_coefficients: dict[int, list[float]] = {
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

            temp_coefficients: dict[int, float] = {
                i: np.nan for i in range(6, order + 1, 2)
            }
            for i in range(6, order + 1, 2):
                if i == 6:
                    temp_coefficients[i] = c_6
                if i == 8:
                    c_8 = 3 * c_6 * r2_r4[atom.element] ** 2
                    temp_coefficients[i] = c_8
                elif i == 10:
                    c_10 = 49.0 / 40.0 * c_8**2 / c_6
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
        coordination_numbers: Array1DFloat = np.array(
            [atom.coordination_number for atom in self._atoms]
        )
        c_n_coefficients: Array1DFloat = {
            key: np.array(value) for key, value in c_n_coefficients.items()
        }

        self._atoms = atoms
        self.c_n_coefficients = c_n_coefficients
        self.coordination_numbers = coordination_numbers

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


def extrapolate_c_n(c_6: float, i: int, j: int, n: int) -> float:
    """Cᴬᴮ coefficients to order n.

    Following eqs. 6-8 in 10.1063/1.3382344

    Args:
        c_6: C_6 coefficient (a.u.)
        i: Atom number of atom A
        j: Atom number of atom B
        n: Order

    Returns:
        c_n: C_n coefficient (a.u.)

    Raises:
        ValueError: For n < 6 and uneven n.
    """
    if n == 6:
        c_n = c_6
    elif n == 8:
        c_n = 3 * r2_r4[i] * r2_r4[j] * c_6
    elif n == 10:
        c_8 = extrapolate_c_n(c_6, i, j, 8)
        c_n = 49.0 / 40.0 * c_8**2 / c_6
    elif n == 9:
        pass
    elif n > 10 and n % 2 == 0:
        c_n = (
            extrapolate_c_n(c_6, i, j, n - 6)
            * (extrapolate_c_n(c_6, i, j, n - 2) / extrapolate_c_n(c_6, i, j, n - 4))
            ** 3
        )
    else:
        raise ValueError("Only defined for even n >= 6.")
    return c_n
