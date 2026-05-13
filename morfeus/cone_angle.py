"""Cone angle code."""

from __future__ import annotations

from collections.abc import Iterable
import functools
import itertools
import math
import typing
from typing import Any
import warnings

import numpy as np

from morfeus.data import atomic_symbols, jmol_colors
from morfeus.geometry import Atom, Cone
from morfeus.io import read_geometry
from morfeus.plotting import get_drawing_cone
from morfeus.typing import Array1DFloat, Array2DFloat, ArrayLike1D, ArrayLike2D
from morfeus.utils import (
    check_distances,
    convert_elements,
    get_radii,
    Import,
    requires_dependency,
)

if typing.TYPE_CHECKING:
    from matplotlib.colors import hex2color
    import pyvista as pv
    from pyvistaqt import BackgroundPlotter


class ConeAngle:
    """Calculates and stores the results of exact cone angle calculation.

    As described in J. Comput. Chem. 2013, 34, 1189.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        atom_1: Index of central atom (1-indexed)
        radii: vdW radii (Å)
        radii_type: Type of vdW radii: 'alvarez', 'bondi', 'crc' or 'truhlar'
        method: Method of calculation: 'internal' or 'libconeangle' (default)

    Attributes:
        cone_angle: Exact cone angle (degrees)
        tangent_atoms: Atoms tangent to cone (1-indexed)

    Raises:
        RunTimeError: If cone angle could not be found by internal algorithm
        ValueError: If atoms within vdW radius of central atom or if exception happened
            with libconeangle or if wrong method chosen
    """

    cone_angle: float
    tangent_atoms: list[int]
    _atoms: list[Atom]
    _max_2_cone: Cone

    def __init__(  # noqa: C901
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        atom_1: int,
        radii: ArrayLike1D | None = None,
        radii_type: str = "crc",
        method: str = "libconeangle",
    ) -> None:
        if atom_1 == 0:
            raise IndexError("Atom indices should not be 0 (1-indexed).")

        # Convert elements to atomic numbers if they are symbols
        elements = convert_elements(elements, output="numbers")
        coordinates: Array2DFloat = np.array(coordinates)

        # Get radii if they are not supplied
        if radii is None:
            radii = get_radii(elements, radii_type=radii_type)
        radii: Array1DFloat = np.array(radii)

        # Check so that no atom is within vdW distance of atom 1
        within = check_distances(elements, coordinates, atom_1, radii=radii)
        if len(within) > 0:
            atom_string = " ".join([str(i) for i in within])
            raise ValueError("Atoms within vdW radius of central atom:", atom_string)

        # Set up coordinate array and translate coordinates
        coordinates -= coordinates[atom_1 - 1]

        # Get list of atoms as Atom objects
        atoms: list[Atom] = []
        for i, (element, coord, radius) in enumerate(
            zip(elements, coordinates, radii), start=1
        ):
            if i != atom_1:
                atom = Atom(element, coord, radius, i)
                atom.get_cone()
                atoms.append(atom)
        self._atoms = atoms

        # Calculate cone angle
        if method == "libconeangle":
            try:
                from libconeangle import cone_angle

                angle, axis, tangent_atoms = cone_angle(coordinates, radii, atom_1 - 1)
                self.cone_angle = angle
                self.tangent_atoms = [i + 1 for i in tangent_atoms]
                atoms = [
                    atom for atom in self._atoms if atom.index in self.tangent_atoms
                ]
                self._cone = Cone(self.cone_angle, atoms, axis)
            except ImportError:
                warnings.warn(
                    "Failed to import libconeangle. Defaulting to method='internal'",
                    stacklevel=2,
                )
                self._cone_angle_internal()
        elif method == "internal":
            self._cone_angle_internal()
        else:
            raise ValueError(
                "Method not implemented. Choose between 'libconeangle' and 'internal'"
            )

    def print_report(self) -> None:
        """Prints report of results."""
        tangent_atoms = [
            atom for atom in self._atoms if atom.index in self.tangent_atoms
        ]
        tangent_labels = [
            f"{atomic_symbols[atom.element]}{atom.index}" for atom in tangent_atoms
        ]
        tangent_string = " ".join(tangent_labels)
        print(f"Cone angle: {self.cone_angle:.1f}")
        print(f"No. tangent atoms: {len(tangent_atoms)}")
        print(f"Tangent to: {tangent_string}")

    def _cone_angle_internal(self) -> None:
        """Calculates cone angle with internal algorithm.

        Raises:
            RuntimeError: If cone cannot be found.
        """
        # Search for cone over single atoms
        cone = self._search_one_cones()

        # Prune out atoms that lie in the shadow of another atom's cone
        if cone is None:
            loop_atoms = list(self._atoms)
            remove_atoms: set[Atom] = set()
            for cone_atom in loop_atoms:
                for test_atom in loop_atoms:
                    if cone_atom is not test_atom:
                        if cone_atom.cone.is_inside(test_atom):
                            remove_atoms.add(test_atom)
            for atom in remove_atoms:
                loop_atoms.remove(atom)
            self._loop_atoms = loop_atoms

        # Search for cone over pairs of atoms
        if cone is None:
            cone = self._search_two_cones()

        # Search for cones over triples of atoms
        if cone is None:
            cone = self._search_three_cones()

        # Check if no cone was found
        if cone is None:
            raise RuntimeError("Cone not found")

        # Set attributes
        self._cone = cone
        self.cone_angle = math.degrees(cone.angle * 2)
        self.tangent_atoms = [atom.index for atom in cone.atoms]

    def _get_upper_bound(self) -> float:
        """Calculates upper bound for apex angle.

        Returns:
            upper_bound: Upper bound to apex angle (radians)
        """
        # Calculate unit vector to centroid
        coordinates: Array2DFloat = np.array([atom.coordinates for atom in self._atoms])
        centroid_vector = np.mean(coordinates, axis=0)
        centroid_unit_vector = centroid_vector / np.linalg.norm(centroid_vector)

        # Getting sums of angle to centroid and vertex angle.
        angle_sums = []
        for atom in self._atoms:
            cone = atom.cone
            cos_angle = np.dot(centroid_unit_vector, cone.normal)
            vertex_angle = math.acos(cos_angle)
            angle_sum = cone.angle + vertex_angle
            angle_sums.append(angle_sum)

        # Select upper bound as the maximum angle
        upper_bound = max(angle_sums)

        return upper_bound

    def _search_one_cones(self) -> Cone | None:
        """Searches over cones tangent to one atom.

        Returns:
            max_1_cone: Largest cone tangent to one atom
        """
        # Get the largest cone
        atoms = self._atoms
        alphas: list[float] = []
        for atom in atoms:
            alphas.append(atom.cone.angle)
        idx = int(np.argmax(alphas))
        max_1_cone = atoms[idx].cone

        # Check if all atoms are contained in cone. If yes, return cone,
        # otherwise, return None.
        in_atoms = []
        test_atoms = [atom for atom in atoms if atom not in max_1_cone.atoms]
        for atom in test_atoms:
            in_atoms.append(max_1_cone.is_inside(atom))
        if all(in_atoms):
            return max_1_cone
        else:
            return None

    def _search_two_cones(self) -> Cone | None:
        """Search over cones tangent to two atoms.

        Returns:
            max_2_cone: Largest cone tangent to two atoms
        """
        # Create two-atom cones
        loop_atoms = self._loop_atoms
        cones = []
        for atom_i, atom_j in itertools.combinations(loop_atoms, r=2):
            cone = _get_two_atom_cone(atom_i, atom_j)
            cones.append(cone)

        # Select largest two-atom cone
        angles = [cone.angle for cone in cones]
        idx = int(np.argmax(angles))
        max_2_cone = cones[idx]
        self._max_2_cone = max_2_cone

        # Check if all atoms are contained in cone. If yes, return cone,
        # otherwise, return None
        in_atoms = []
        for atom in loop_atoms:
            in_atoms.append(max_2_cone.is_inside(atom))

        if all(in_atoms):
            return max_2_cone
        else:
            return None

    def _search_three_cones(self) -> Cone:
        """Search over cones tangent to three atoms.

        Returns:
            min_3_cone: Smallest cone tangent to three atoms
        """
        # Create three-atom cones
        loop_atoms = self._loop_atoms
        cones = []
        for atom_i, atom_j, atom_k in itertools.combinations(loop_atoms, r=3):
            three_cones = _get_three_atom_cones(atom_i, atom_j, atom_k)
            cones.extend(three_cones)

        # Get upper and lower bound to apex angle
        upper_bound = self._get_upper_bound()
        lower_bound = self._max_2_cone.angle

        # Remove cones from consideration which are outside the bounds
        remove_cones = []
        for cone in cones:
            if cone.angle - lower_bound < -1e-5 or upper_bound - cone.angle < -1e-5:
                remove_cones.append(cone)

        for cone in reversed(remove_cones):
            cones.remove(cone)

        # Keep only cones that encompass all atoms
        keep_cones: list[Cone] = []
        for cone in cones:
            in_atoms = []
            for atom in loop_atoms:
                in_atoms.append(cone.is_inside(atom))
            if all(in_atoms):
                keep_cones.append(cone)

        # Take the smallest cone that encompasses all atoms
        angles = [cone.angle for cone in keep_cones]
        idx = int(np.argmin(angles))
        min_3_cone = keep_cones[idx]

        return min_3_cone

    @requires_dependency(
        [
            Import(module="matplotlib.colors", item="hex2color"),
            Import(module="pyvista", alias="pv"),
            Import(module="pyvistaqt", item="BackgroundPlotter"),
        ],
        globals(),
    )
    def draw_3D(
        self,
        atom_scale: float = 1,
        background_color: str = "white",
        cone_color: str = "steelblue",
        cone_opacity: float = 0.75,
    ) -> None:
        """Draw a 3D representation of the molecule with the cone.

        Args:
            atom_scale: Scaling factor for atom size
            background_color: Background color for plot
            cone_color: Cone color
            cone_opacity: Cone opacity
        """
        # Set up plotter
        p = BackgroundPlotter()
        p.set_background(background_color)

        # Draw molecule
        for atom in self._atoms:
            color = hex2color(jmol_colors[atom.element])
            radius = atom.radius * atom_scale
            sphere = pv.Sphere(center=list(atom.coordinates), radius=radius)
            p.add_mesh(sphere, color=color, opacity=1, name=str(atom.index))

        # Determine direction and extension of cone
        angle = math.degrees(self._cone.angle)
        coordinates: Array2DFloat = np.array([atom.coordinates for atom in self._atoms])
        radii: Array1DFloat = np.array([atom.radius for atom in self._atoms])
        if angle > 180:
            normal = -self._cone.normal
        else:
            normal = self._cone.normal
        projected = np.dot(normal, coordinates.T) + np.array(radii)

        max_extension = np.max(projected)
        if angle > 180:
            max_extension += 1

        # Make the cone
        cone = get_drawing_cone(
            center=[0, 0, 0] + (max_extension * normal) / 2,
            direction=-normal,
            angle=angle,
            height=max_extension,
            capping=False,
            resolution=100,
        )
        p.add_mesh(cone, opacity=cone_opacity, color=cone_color)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


def _get_two_atom_cone(atom_i: Atom, atom_j: Atom) -> Cone:
    """Creates a cone tangent to two atoms.

    Args:
        atom_i: First tangent atom
        atom_j: Second tangent atom

    Returns:
        cone: Cone tangent to the two atoms
    """
    # Get the cone angle
    cone_i = atom_i.cone
    cone_j = atom_j.cone
    beta_i = cone_i.angle
    beta_j = cone_j.angle
    beta_ij = math.acos(np.dot(atom_i.cone.normal, atom_j.cone.normal))
    alpha_ij = (beta_ij + beta_i + beta_j) / 2

    # Get the cone normal
    a_ij = (1 / math.sin(beta_ij)) * math.sin(0.5 * (beta_ij + beta_i - beta_j))
    b_ij = (1 / math.sin(beta_ij)) * math.sin(0.5 * (beta_ij - beta_i + beta_j))
    c_ij = 0

    n = a_ij * cone_i.normal + b_ij * cone_j.normal + c_ij
    n = n / np.linalg.norm(n)

    # Create cone
    angle = alpha_ij
    normal = n
    cone = Cone(angle, [atom_i, atom_j], normal)

    return cone


def _get_three_atom_cones(atom_i: Atom, atom_j: Atom, atom_k: Atom) -> list[Cone]:
    """Creates cones tangent to three atoms.

    Args:
        atom_i: First tangent atom
        atom_j: Second tangent atom
        atom_k: Third tangent atom

    Returns:
        cones: Cones tangent to the three atoms
    """
    # Set up vertex angles
    beta_i = atom_i.cone.angle
    beta_j = atom_j.cone.angle
    beta_k = atom_k.cone.angle

    # Set up angles between atom vectors
    beta_ij = math.acos(np.dot(atom_i.cone.normal, atom_j.cone.normal))

    # Set up normal vectors to atoms
    m_i = atom_i.cone.normal
    m_j = atom_j.cone.normal
    m_k = atom_k.cone.normal

    # Setup matrices
    u: Array1DFloat = np.array([math.cos(beta_i), math.cos(beta_j), math.cos(beta_k)])
    v: Array1DFloat = np.array([math.sin(beta_i), math.sin(beta_j), math.sin(beta_k)])
    N: Array2DFloat = np.array(
        [np.cross(m_j, m_k), np.cross(m_k, m_i), np.cross(m_i, m_j)]
    ).T
    P: Array2DFloat = N.T @ N
    gamma = np.dot(m_i, np.cross(m_j, m_k))

    # Set up coefficients of quadratic equation
    A = u @ P @ u
    B = v.T @ P @ v
    C = u.T @ P @ v
    D = gamma**2

    # Solve quadratic equation
    p2 = (A - B) ** 2 + 4 * C**2
    p1 = 2 * (A - B) * (A + B - 2 * D)
    p0 = (A + B - 2 * D) ** 2 - 4 * C**2

    roots = np.roots([p2, p1, p0])
    roots = np.real_if_close(roots, tol=1e10)
    roots[np.isclose(roots, 1, rtol=1e-9, atol=0.0)] = 1
    roots[np.isclose(roots, -1, rtol=1e-9, atol=0.0)] = -1

    cos_roots = [
        math.acos(roots[0]),
        2 * np.pi - math.acos(roots[0]),
        math.acos(roots[1]),
        2 * np.pi - math.acos(roots[1]),
    ]

    # Test roots and keep only those that are physical
    angles = []
    D_tests = []
    for root in cos_roots:
        alpha = root / 2
        test = (
            A * math.cos(alpha) ** 2
            + B * math.sin(alpha) ** 2
            + 2 * C * math.sin(alpha) * math.cos(alpha)
        )
        D_test = abs(test - D)
        angles.append(alpha)
        D_tests.append(D_test)
    angles = np.array(angles)
    D_tests = np.array(D_tests)
    physical_angles = angles[np.argsort(D_tests)][:2]

    # Create cones for physical angles
    cones = []
    for alpha in physical_angles:
        # Calculate normal vector
        a_ij = (
            math.cos(alpha - beta_i) - math.cos(alpha - beta_j) * math.cos(beta_ij)
        ) / math.sin(beta_ij) ** 2
        b_ij = (
            math.cos(alpha - beta_j) - math.cos(alpha - beta_i) * math.cos(beta_ij)
        ) / math.sin(beta_ij) ** 2
        c_ij_squared = 1 - a_ij**2 - b_ij**2 - 2 * a_ij * b_ij * math.cos(beta_ij)
        # Set c_ij_squared to 0 if negative due to numerical precision.
        if c_ij_squared < 0:
            c_ij_squared = 0
        c_ij = math.sqrt(c_ij_squared)
        p = N @ (u * math.cos(alpha) + v * math.sin(alpha)).reshape(-1)
        sign = np.sign(gamma) * np.sign(np.dot(p, np.cross(m_i, m_j)))
        if np.sign(c_ij) != sign:
            c_ij = -c_ij
        n = a_ij * m_i + b_ij * m_j + c_ij * 1 / math.sin(beta_ij) * np.cross(m_i, m_j)

        # Create cone
        cone = Cone(alpha, [atom_i, atom_j, atom_k], n)
        cones.append(cone)

    return cones


def cli(file: str) -> Any:
    """CLI for cone angle.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(ConeAngle, elements, coordinates)
