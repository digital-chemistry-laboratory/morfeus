"""Help classes and functions related to geometry."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
import math

import numpy as np
from scipy.spatial.transform import Rotation

from morfeus.data import ANGSTROM_TO_BOHR
from morfeus.typing import (
    Array1DBool,
    Array1DFloat,
    Array2DFloat,
    ArrayLike1D,
    ArrayLike2D,
)
from morfeus.utils import get_connectivity_matrix


class Atom:
    """Atom common for morfeus calculations.

    Args:
        element: Atomic number (starting from 1)
        coordinates: Coordinates (Å)
        radius: vdW radius (Å)
        index: Atom index (starting from 1)

    Attributes:
        accessible_mask: Boolean mask for accessible points
        accessible_points: Points accessible to solvent (Å)
        area: Solvent-accessible surface area (Å²)
        cone: Cone tangent to atom
        coordinates: Coordinates (Å)
        coordination_number: Coordination number
        element: Atomic number (1-indexed)
        index: Atom index (1-indexed)
        invisible_mask: Boolean mask for invisible points
        occluded_mask: Boolean mask for occluded points
        occluded_points: Points occluded by other atoms (Å)
        p_values: P values
        point_areas: Point areas (Å²)
        point_volumes: Point volumes (Å³)
        point: Points (Å)
        proximal_mask: Boolean mask for proximal points
        radius: vdW radius (Å)
        volume: Volume inside solvent-accessible surface area (Å³)
    """

    accessible_mask: Array1DBool
    accessible_points: Array2DFloat
    area: float
    cone: "Cone"
    coordinates: Array2DFloat
    coordination_number: float
    element: int
    index: int
    invisible_mask: Array1DBool
    occluded_mask: Array1DBool
    occluded_points: Array2DFloat
    p_values: Array1DFloat
    point_areas: Array1DFloat
    point_volumes: Array1DFloat
    points: Array2DFloat
    proximal_mask: Array1DBool
    radius: float
    volume: float

    def __init__(
        self, element: int, coordinates: ArrayLike1D, radius: float, index: int
    ) -> None:
        # Set up initial attributes
        self.coordinates = np.array(coordinates)
        self.element = element
        self.index = index
        self.radius = radius

    def get_cone(self) -> None:
        """Construct cone for atom."""
        vector = self.coordinates
        normal = vector / np.linalg.norm(vector)
        sin_alpha = self.radius / np.linalg.norm(vector)
        alpha = math.asin(sin_alpha)
        cone = Cone(alpha, [self], normal)
        self.cone = cone

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.index!r})"


class Cone:
    """Cone used in cone angle calculations.

    Args:
        angle: Cone angle (rad)
        atoms: Atoms that are tangent to cone (1-indexed)
        normal: Normal vector of cone

    Attributes:
        angle: Cone angle (rad)
        atoms: Atoms that are tangent to cone (1-indexed)
        normal: Normal vector of cone
    """

    angle: float
    atoms: list[Atom]
    normal: Array1DFloat

    def __init__(
        self, angle: float, atoms: Sequence[Atom], normal: ArrayLike1D
    ) -> None:
        self.angle = angle
        self.atoms = list(atoms)
        self.normal = np.array(normal)

    def is_inside(self, atom: Atom) -> bool:
        """Tests if atom lies inside the cone.

        Args:
            atom: Atom to test.

        Returns:
            True if inside, False if outside.
        """
        # Get vertex angle of atom
        beta = atom.cone.angle

        # Calculate angle between cone normal vector and unit vector to atom
        cos = np.dot(atom.cone.normal, self.normal)

        # Take into account numerical problems that sometimes gives a value
        # somewhat above 1
        if 1 - cos > 0 and 1 - cos < 1e-5:
            cos = 1
        angle = math.acos(cos)

        # Check if atom lies inside cone, within numerical reason
        diff = self.angle - (beta + angle)
        if diff > -1e-5:
            return True
        else:
            return False

    def is_inside_points(
        self, points: ArrayLike2D, method: str = "cross"
    ) -> Array1DBool:
        """Test if points are inside cone of atom.

        Args:
            points: Points to check (Å)
            method: Method for testing: 'angle', 'cross' or 'dot'

        Returns:
            is_inside: Boolean array with points marked as inside

        Raises:
            ValueError: When method not supported
        """
        points: Array2DFloat = np.array(points)
        if method in ["cross", "dot"]:
            # Calculate radius of cone at distance of each point
            cone_distances = np.dot(points, self.normal)
            cone_radii = np.tan(self.angle) * cone_distances

            # Calculate orthogonal distance of points to cone normal vector
            if method == "cross":
                orth_distances = np.linalg.norm(np.cross(-points, self.normal), axis=1)
            elif method == "dot":
                orth_distances = np.linalg.norm(
                    -points - np.dot(points, self.normal).reshape(-1, 1) * self.normal,
                    axis=1,
                )

            # Determine if distance is smaller than cone radius.
            is_inside = orth_distances < cone_radii
        elif method == "angle":
            norm_points = points / np.linalg.norm(points, axis=1).reshape(-1, 1)

            # Calculate angle between cone normal vector and unit vector to
            # atom
            cos = np.dot(norm_points, self.normal)

            # Take into account numerical problems that sometimes gives a value
            # somewhat above 1
            cos[np.logical_and(1 - cos > 0, 1 - cos < 1e-5)] = 1
            angle = np.arccos(cos)

            # Check if atom lies inside cone, within numerical reason
            diff = self.angle - angle

            is_inside = diff > -1e-5
        else:
            raise ValueError(f"method={method} not supported.")

        return is_inside

    def __repr__(self) -> str:
        atoms = ", ".join([str(atom.index) for atom in self.atoms])
        return f"{self.__class__.__name__}(Tangent atoms: {atoms})"


class Sphere:
    """Sphere class for creating and holding points on vdW surface.

    Args:
        center: Coordinates for center (Å)
        density: Area per point (Å²) for empty sphere
            and volume per point (Å³) for filled sphere.
        filled: Whether a sphere with internal points should be constructed (works only
            with method='projection')
        method: Method for generating points: 'fibonacci', 'polar' or 'projection'
        radius: Radius (Å)

    Attributes:
        area: Area (Å²)
        center: Coordinates for sphere center (Å)
        circumference: Circumference (Å)
        density: Density of points (Å² or Å³)
        points: Points in/on sphere (Å)
        radius: Radius (Å)
        volume: Volume (Å³)
    """

    area: float
    center: Array1DFloat
    circumference: float
    density: float
    points: Array2DFloat
    radius: float
    volume: float

    def __init__(
        self,
        center: ArrayLike1D,
        radius: float,
        density: float = 0.005,
        method: str = "fibonacci",
        filled: bool = False,
    ) -> None:
        self.center = np.array(center)
        self.radius = radius
        self.circumference = math.pi * radius * 2
        self.area = 4 * radius**2 * math.pi
        self.volume = 4 * radius**3 * math.pi / 3
        self.density = density

        if method == "polar":
            self.points = self._get_points_polar(density=density)
        elif method == "projection":
            self.points = self._get_points_projected(density=density, filled=filled)
        elif method == "fibonacci":
            self.points = self._get_points_fibonacci(density=density)

    @staticmethod
    def _get_cartesian_coordinates(
        r: float, theta: ArrayLike1D, phi: ArrayLike1D
    ) -> Array2DFloat:
        """Converts polar to Cartesian coordinates.

        Args:
            phi: Phi angles (radians)
            r: Radius (Å)
            theta: Theta angles (radians)

        Returns:
            points: Cartesian points (Å)
        """
        # Calculate x, y and z coordinates
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        # Stack coordinates as columns
        points: Array2DFloat = np.column_stack((x, y, z))

        return points

    def _get_points_fibonacci(self, density: float) -> Array2DFloat:
        """Construct points on sphere surface by the Fibonacci golden spiral method.

        Method vectorized from
        https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere

        Args:
            density: Area per point (Å²)

        Returns:
            points: Surface points (Å)
        """
        # Generate points on unit sphere
        rnd = 1
        n = int(round((self.area / density)))
        offset = 2.0 / n
        increment = math.pi * (3.0 - math.sqrt(5.0))

        i = np.arange(n)
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - np.square(y))
        phi = np.mod((i + rnd), n) * increment
        x = np.cos(phi) * r
        z = np.sin(phi) * r

        # Generate points and adjust with radius and center
        points: Array2DFloat = np.column_stack((x, y, z))
        points = points * self.radius
        points = points + self.center

        return points

    def _get_points_polar(self, density: float) -> Array2DFloat:
        """Construct points on sphere by polar coordinate method.

        Args:
            density: Area per point (Å²)

        Returns:
            points: Points on sphere (Å)
        """
        # Calculate number of points
        n = int(round((self.area / density / 2) ** (1 / 2)))

        # Set up points along theta and phi
        theta = np.linspace(0, math.pi, n)
        phi = np.linspace(0, 2 * math.pi, 2 * n)

        # Combine together all the possible combinations of theta and phi
        combined_theta_phi: Array2DFloat = np.dstack(np.meshgrid(theta, phi)).reshape(
            -1, 2
        )
        theta = combined_theta_phi[:, 0]
        phi = combined_theta_phi[:, 1]

        # Get the Cartesian coordinates
        points = self._get_cartesian_coordinates(self.radius, theta, phi)

        # Adjust to sphere center
        points = points + self.center

        return points

    def _get_points_projected(
        self, density: float, filled: bool = False
    ) -> Array2DFloat:
        """Construct points on sphere surface by projection.

        Args:
            density: Area per point (Å²) for empty sphere and volume
                per point (Å³) for filled sphere
            filled: Whether to generate internal points

        Returns:
            points: Array of surface points (Å)
        """
        # Calculate number of points from density of empty or filled sphere.
        if filled:
            n = int(round((self.volume / density * 6 / math.pi) ** (1 / 3)))
        else:
            n = int(round((self.area / density * 6 / math.pi) ** (1 / 3)))

        # Generate points in box
        r = self.radius
        x = np.linspace(-r, r, n)
        y = np.linspace(-r, r, n)
        z = np.linspace(-r, r, n)
        points: Array2DFloat = np.stack(np.meshgrid(x, y, z), -1).reshape(-1, 3)

        # Remove points outside of sphere
        lengths = np.linalg.norm(points, axis=1)
        points = points[lengths <= r]

        # Project points onto sphere surface if sphere is empty
        if not filled:
            points = points / np.linalg.norm(points, axis=1).reshape(-1, 1) * r

        # Adjust with sphere center
        points = points + self.center

        return points

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(center: {self.center}, "
            f"radius: {self.radius})"
        )


class Bond:
    """Bond internal coordinate.

    Args:
        atom_1: Index of atom 1 (1-indexed)
        atom_2: Index of atom 2 (1-indexed)

    Attributes:
        i: Index of atom 1 (1-indexed)
        j: Index of atom 2 (1-indexed)
    """

    i: int
    j: int
    atoms: list[int]

    def __init__(self, atom_1: int, atom_2: int) -> None:
        self.i = atom_1
        self.j = atom_2
        self.atoms = [atom_1, atom_2]

    def get_b_vector(self, coordinates: ArrayLike2D) -> Array1DFloat:
        """Calculate vector of B matrix for internal coordinate.

        Code adapted from pyberny: https://github.com/jhrmnn/pyberny.

        Args:
            coordinates: Coordinates (Å)

        Returns:
            b_vector: Vector of B matrix (a.u.)
        """
        coordinates: Array2DFloat = np.array(coordinates)
        i, j = self.i - 1, self.j - 1
        v = (coordinates[i] - coordinates[j]) * ANGSTROM_TO_BOHR
        r = np.linalg.norm(v)
        grad: Array1DFloat = np.array([v / r, -v / r])

        b_vector = np.zeros(coordinates.size)
        b_vector[i * 3 : i * 3 + 3] = grad[0]
        b_vector[j * 3 : j * 3 + 3] = grad[1]

        return b_vector

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return (self.i == other.i or self.i == other.j) and (
                self.j == other.i or self.j == other.j
            )
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.i, self.j))

    def __repr__(self) -> str:
        atoms = ", ".join(str(i) for i in sorted(self.atoms))
        return f"{self.__class__.__name__}({atoms})"


class Angle:
    """Bond internal coordinate.

    Code adapted from pyberny: https://github.com/jhrmnn/pyberny.

    Args:
        atom_1: Index of atom 1 (1-indexed)
        atom_2: Index of atom 2 (1-indexed)
        atom_3: Index of atom 3 (1-indexed)

    Attributes:
        i: Index of atom 1 (1-indexed)
        j: Index of atom 2 (1-indexed)
        k: Index of atom 3 (1-indexed)
    """

    i: int
    j: int
    k: int
    atoms: list[int]

    def __init__(self, atom_1: int, atom_2: int, atom_3: int) -> None:
        self.i = atom_1
        self.j = atom_2
        self.k = atom_3
        self.atoms = [atom_1, atom_2, atom_3]

    def get_b_vector(self, coordinates: ArrayLike2D) -> Array1DFloat:
        """Calculate vector of B matrix for internal coordinate.

        Code adapted from pyberny: https://github.com/jhrmnn/pyberny.

        Args:
            coordinates: Coordinates (Å)

        Returns:
            b_vector: Vector of B matrix (a.u.)
        """
        coordinates: Array2DFloat = np.array(coordinates)
        i, j, k = self.i - 1, self.j - 1, self.k - 1
        v_1 = (coordinates[i] - coordinates[j]) * ANGSTROM_TO_BOHR
        v_2 = (coordinates[k] - coordinates[j]) * ANGSTROM_TO_BOHR
        dot_product = np.dot(v_1, v_2) / (np.linalg.norm(v_1) * np.linalg.norm(v_2))
        if dot_product < -1:
            dot_product = -1
        elif dot_product > 1:
            dot_product = 1
        phi = np.arccos(dot_product)
        if abs(phi) > np.pi - 1e-6:
            grad = [
                (np.pi - phi) / (2 * np.linalg.norm(v_1) ** 2) * v_1,
                (1 / np.linalg.norm(v_1) - 1 / np.linalg.norm(v_2))
                * (np.pi - phi)
                / (2 * np.linalg.norm(v_1))
                * v_1,
                (np.pi - phi) / (2 * np.linalg.norm(v_2) ** 2) * v_2,
            ]
        else:
            grad = [
                1 / np.tan(phi) * v_1 / np.linalg.norm(v_1) ** 2
                - v_2 / (np.linalg.norm(v_1) * np.linalg.norm(v_2) * np.sin(phi)),
                (v_1 + v_2) / (np.linalg.norm(v_1) * np.linalg.norm(v_2) * np.sin(phi))
                - 1
                / np.tan(phi)
                * (v_1 / np.linalg.norm(v_1) ** 2 + v_2 / np.linalg.norm(v_2) ** 2),
                1 / np.tan(phi) * v_2 / np.linalg.norm(v_2) ** 2
                - v_1 / (np.linalg.norm(v_1) * np.linalg.norm(v_2) * np.sin(phi)),
            ]

        b_vector = np.zeros(coordinates.size)
        b_vector[i * 3 : i * 3 + 3] = grad[0]
        b_vector[j * 3 : j * 3 + 3] = grad[1]
        b_vector[k * 3 : k * 3 + 3] = grad[2]

        return b_vector

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return (
                (self.i == other.i or self.i == other.k)
                and (self.k == other.i or self.k == other.k)
                and self.j == other.j
            )
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.i, self.j, self.k))

    def __repr__(self) -> str:
        atoms = ", ".join(str(i) for i in sorted(self.atoms))
        return f"{self.__class__.__name__}({atoms})"


class Dihedral:
    """Bond internal coordinate.

    Code adapted from pyberny: https://github.com/jhrmnn/pyberny.

    Args:
        atom_1: Index of atom 1 (1-indexed)
        atom_2: Index of atom 2 (1-indexed)
        atom_3: Index of atom 3 (1-indexed)
        atom_4: Index of atom 4 (1-indexed)

    Attributes:
        i: Index of atom 1 (1-indexed)
        j: Index of atom 2 (1-indexed)
        k: Index of atom 3 (1-indexed)
        l: Index of atom 4 (1-indexed)
    """

    i: int
    j: int
    k: int
    l: int
    atoms: list[int]

    def __init__(self, atom_1: int, atom_2: int, atom_3: int, atom_4: int) -> None:
        self.i = atom_1
        self.j = atom_2
        self.k = atom_3
        self.l = atom_4
        self.atoms = [atom_1, atom_2, atom_3, atom_4]

    def get_b_vector(self, coordinates: ArrayLike2D) -> Array1DFloat:
        """Calculate vector of B matrix for internal coordinate.

        Code adapted from pyberny: https://github.com/jhrmnn/pyberny.

        Args:
            coordinates: Coordinates (Å)

        Returns:
            b_vector: Vector of B matrix (a.u.)
        """
        coordinates: Array2DFloat = np.array(coordinates)
        i, j, k, l = self.i - 1, self.j - 1, self.k - 1, self.l - 1
        v_1 = (coordinates[i] - coordinates[j]) * ANGSTROM_TO_BOHR
        v_2 = (coordinates[l] - coordinates[k]) * ANGSTROM_TO_BOHR
        w = (coordinates[k] - coordinates[j]) * ANGSTROM_TO_BOHR
        ew = w / np.linalg.norm(w)
        a_1 = v_1 - np.dot(v_1, ew) * ew
        a_2 = v_2 - np.dot(v_2, ew) * ew
        sgn = np.sign(np.linalg.det(np.array([v_2, v_1, w])))
        sgn = sgn or 1
        dot_product = np.dot(a_1, a_2) / (np.linalg.norm(a_1) * np.linalg.norm(a_2))
        if dot_product < -1:
            dot_product = -1
        elif dot_product > 1:
            dot_product = 1
        phi = np.arccos(dot_product) * sgn
        g: Array1DFloat
        if abs(phi) > np.pi - 1e-6:
            g = np.cross(w, a_1)
            g = g / np.linalg.norm(g)
            A = np.dot(v_1, ew) / np.linalg.norm(w)
            B = np.dot(v_2, ew) / np.linalg.norm(w)
            grad = [
                g / (np.linalg.norm(g) * np.linalg.norm(a_1)),
                -((1 - A) / np.linalg.norm(a_1) - B / np.linalg.norm(a_2)) * g,
                -((1 + B) / np.linalg.norm(a_2) + A / np.linalg.norm(a_1)) * g,
                g / (np.linalg.norm(g) * np.linalg.norm(a_2)),
            ]
        elif abs(phi) < 1e-6:
            g = np.cross(w, a_1)
            g = g / np.linalg.norm(g)
            A = np.dot(v_1, ew) / np.linalg.norm(w)
            B = np.dot(v_2, ew) / np.linalg.norm(w)
            grad = [
                g / (np.linalg.norm(g) * np.linalg.norm(a_1)),
                -((1 - A) / np.linalg.norm(a_1) + B / np.linalg.norm(a_2)) * g,
                ((1 + B) / np.linalg.norm(a_2) - A / np.linalg.norm(a_1)) * g,
                -g / (np.linalg.norm(g) * np.linalg.norm(a_2)),
            ]
        else:
            A = np.dot(v_1, ew) / np.linalg.norm(w)
            B = np.dot(v_2, ew) / np.linalg.norm(w)
            grad = [
                1 / np.tan(phi) * a_1 / np.linalg.norm(a_1) ** 2
                - a_2 / (np.linalg.norm(a_1) * np.linalg.norm(a_2) * np.sin(phi)),
                ((1 - A) * a_2 - B * a_1)
                / (np.linalg.norm(a_1) * np.linalg.norm(a_2) * np.sin(phi))
                - 1
                / np.tan(phi)
                * (
                    (1 - A) * a_1 / np.linalg.norm(a_1) ** 2
                    - B * a_2 / np.linalg.norm(a_2) ** 2
                ),
                ((1 + B) * a_1 + A * a_2)
                / (np.linalg.norm(a_1) * np.linalg.norm(a_2) * np.sin(phi))
                - 1
                / np.tan(phi)
                * (
                    (1 + B) * a_2 / np.linalg.norm(a_2) ** 2
                    + A * a_1 / np.linalg.norm(a_1) ** 2
                ),
                1 / np.tan(phi) * a_2 / np.linalg.norm(a_2) ** 2
                - a_1 / (np.linalg.norm(a_1) * np.linalg.norm(a_2) * np.sin(phi)),
            ]
        b_vector = np.zeros(coordinates.size)
        b_vector[i * 3 : i * 3 + 3] = grad[0]
        b_vector[j * 3 : j * 3 + 3] = grad[1]
        b_vector[k * 3 : k * 3 + 3] = grad[2]
        b_vector[l * 3 : l * 3 + 3] = grad[3]

        return b_vector

    def __eq__(self, other: object) -> bool:
        if isinstance(other, self.__class__):
            return all(
                [
                    self.i == other.i,
                    self.j == other.j,
                    self.k == other.k,
                    self.l == other.l,
                ]
            ) or all(
                [
                    self.i == other.l,
                    self.j == other.k,
                    self.k == other.j,
                    self.l == other.i,
                ]
            )
        else:
            return False

    def __hash__(self) -> int:
        return hash((self.i, self.j, self.k, self.l))

    def __repr__(self) -> str:
        atoms = ", ".join(str(i) for i in sorted(self.atoms))
        return f"{self.__class__.__name__}({atoms})"


class InternalCoordinates:
    """Internal coordinates.

    Attributes:
        internal_coordinates: Internal coordinates
    """

    internal_coordinates: list[Bond | Angle | Dihedral]

    def __init__(self) -> None:
        self.internal_coordinates = []

    def add_bond(self, atom_1: int, atom_2: int) -> None:
        """Add bonds to internal coordinates.

        Args:
            atom_1: Index of atom 1 (1-indexed)
            atom_2: Index of atom 2 (1-indexed)
        """
        bond = Bond(atom_1, atom_2)
        if bond not in self.internal_coordinates:
            self.internal_coordinates.append(bond)

    def add_angle(self, atom_1: int, atom_2: int, atom_3: int) -> None:
        """Add angle to internal coordinates.

        Args:
            atom_1: Index of atom 1 (1-indexed)
            atom_2: Index of atom 2 (1-indexed)
            atom_3: Index of atom 3 (1-indexed)
        """
        angle = Angle(atom_1, atom_2, atom_3)
        if angle not in self.internal_coordinates:
            self.internal_coordinates.append(angle)

    def add_dihedral(self, atom_1: int, atom_2: int, atom_3: int, atom_4: int) -> None:
        """Add dihedral angle to internal coordinates.

        Args:
            atom_1: Index of atom 1 (1-indexed)
            atom_2: Index of atom 2 (1-indexed)
            atom_3: Index of atom 3 (1-indexed)
            atom_4: Index of atom 4 (1-indexed)
        """
        dihedral = Dihedral(atom_1, atom_2, atom_3, atom_4)
        if dihedral not in self.internal_coordinates:
            self.internal_coordinates.append(dihedral)

    def add_internal_coordinate(self, atoms: Sequence[int]) -> None:
        """Add internal coordinate automatically depending on number of atoms.

        Args:
            atoms: Sequence of atom indices (1-index).
        """
        if len(atoms) == 2:
            self.add_bond(*atoms)
        elif len(atoms) == 3:
            self.add_angle(*atoms)
        elif len(atoms) == 4:
            self.add_dihedral(*atoms)

    def detect_bonds(
        self,
        coordinates: ArrayLike2D,
        elements: Iterable[int] | Iterable[str] | None,
        radii: ArrayLike1D | None = None,
        radii_type: str = "pyykko",
        scale_factor: float = 1.2,
    ) -> None:
        """Detect bonds based on covalent radii cutoffs."""
        # Detect bonds based on covalent radii
        connectivity_matrix = get_connectivity_matrix(
            coordinates,
            elements=elements,
            radii=radii,
            radii_type=radii_type,
            scale_factor=scale_factor,
        )
        indices = np.where(connectivity_matrix)

        bonds = set()
        for i, j in zip(*indices):
            bond = frozenset([i + 1, j + 1])
            bonds.add(bond)

        # Add each bond as an internal coordinate
        for bond in bonds:
            i, j = sorted(bond)
            self.add_bond(i, j)

    def get_B_matrix(self, coordinates: ArrayLike2D) -> Array2DFloat:
        """Calculate B matrix for coordinates.

        Args:
            coordinates: Coordinates (Å)

        Returns:
            B_matrix: B matrix
        """
        coordinates: Array2DFloat = np.array(coordinates)
        b_vectors = []
        for internal_coordinate in self.internal_coordinates:
            b_vectors.append(internal_coordinate.get_b_vector(coordinates))
        B_matrix: Array2DFloat = np.vstack(b_vectors)

        return B_matrix

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}({len(self.internal_coordinates)} "
            "coordinates)"
        )


def rotate_coordinates(
    coordinates: ArrayLike2D,
    vector: ArrayLike1D,
    axis: ArrayLike1D,
) -> Array2DFloat:
    """Rotates coordinates by the rotation that aligns vector with axis.

    Args:
        axis: Reference vector to align input vector with
        coordinates: Coordinates to rotate
        vector: Vector to align

    Returns:
        rotated_coordinates: Rotated coordinates.
    """
    coordinates: Array2DFloat = np.array(coordinates)
    vector: Array1DFloat = np.array(vector)
    axis: Array1DFloat = np.array(axis)

    # Get real part of quaternion
    real = np.dot(vector, axis).reshape(1) + 1

    # Get imaginary dimensions and handle case of antiparallel vectors
    if real < 1e-6:
        w: Array1DFloat = np.cross(vector, np.array([0, 0, 1]))
        if np.linalg.norm(w) < 1e-6:
            w = np.cross(vector, np.array([1, 0, 0]))
    else:
        w = np.cross(vector, axis)

    # Form quaternion and normalize
    q = np.concatenate((w, real))
    q = q / np.linalg.norm(q)

    # Rotate atomic coordinates
    rotation = Rotation.from_quat(q)
    rotated_coordinates = rotation.apply(coordinates)

    return rotated_coordinates


def kabsch_rotation_matrix(
    P: ArrayLike2D, Q: ArrayLike2D, center: bool = True
) -> Array2DFloat:
    """Construct Kabsch rotation matrix.

    Constructs the rotation matrix that overlays the points in P with the points in Q
    with minimum RMSD. https://en.wikipedia.org/wiki/Kabsch_algorithm

    Args:
        P: Coordinates to be rotated
        Q: Reference coordinates
        center: Whether to center P and Q at origin.

    Returns:
        R: Rotation matrix
    """
    P: Array2DFloat = np.array(P)
    Q: Array2DFloat = np.array(Q)

    # Calculate centroids and center coordinates
    if center:
        centroid_P = np.mean(P, axis=0)
        centroid_Q = np.mean(Q, axis=0)
        P -= centroid_P
        Q -= centroid_Q

    # Calculate cross-covariance matrix
    H = P.T @ Q

    # Perform SVD
    U, S, V_T = np.linalg.svd(H)

    # Determine correction for right-hand coordinate system
    d = np.sign(np.linalg.det(V_T.T @ U.T))

    # Obtain rotation matrix
    R = V_T.T @ np.array([[1, 0, 0], [0, 1, 0], [0, 0, d]]) @ U.T

    return R


def sphere_line_intersection(
    vector: ArrayLike1D, center: ArrayLike1D, radius: float
) -> Array2DFloat:
    """Get points of intersection between line and sphere.

    Follows the procedure outlined here: http://paulbourke.net/geometry/circlesphere/.

    Args:
        vector: Vector giving direction of line
        center: Center of sphere
        radius: Radius of sphere

    Returns:
        intersection_points: Intersection points
    """
    vector: Array1DFloat = np.array(vector)
    center: Array1DFloat = np.array(center)

    # Set up points
    p_1 = vector
    p_2 = vector * 2

    # Set up terms for quadratic equation
    a = (p_2[0] - p_1[0]) ** 2 + (p_2[1] - p_1[1]) ** 2 + (p_2[2] - p_1[2]) ** 2
    b = 2 * (
        (p_2[0] - p_1[0]) * (p_1[0] - center[0])
        + (p_2[1] - p_1[1]) * (p_1[1] - center[1])
        + (p_2[2] - p_1[2]) * (p_1[2] - center[2])
    )
    c = (
        center[0] ** 2
        + center[1] ** 2
        + center[2] ** 2
        + p_1[0] ** 2
        + p_1[1] ** 2
        + p_1[2] ** 2
        - 2 * (center[0] * p_1[0] + center[1] * p_1[1] + center[2] * p_1[2])
        - radius**2
    )

    # Determine value within the square root and select cases
    within_sqrt = b**2 - 4 * a * c
    if within_sqrt < 0:
        us = []
    elif within_sqrt == 0:
        us = [-b / (2 * a)]
    elif within_sqrt > 0:
        us = [
            (-b + math.sqrt(within_sqrt)) / (2 * a),
            (-b - math.sqrt(within_sqrt)) / (2 * a),
        ]

    # Calculate intersection points.
    intersection_points: Array2DFloat = np.vstack([p_1 + u * (p_2 - p_1) for u in us])

    return intersection_points
