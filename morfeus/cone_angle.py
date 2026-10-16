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
from morfeus.typing import Array2DFloat, ArrayLike1D, ArrayLike2D
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
        coordinates = np.array(coordinates)

        # Get radii if they are not supplied
        if radii is None:
            radii = get_radii(elements, radii_type=radii_type)
        radii = np.array(radii)

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
        coordinates = np.array([atom.coordinates for atom in self._atoms])
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

        Constructs the same cones as ``_get_three_atom_cones``, for every
        triple of atoms at once with array operations instead of a Python
        loop per triple, and solves the tangency quadratic in closed form
        instead of with ``np.roots``. The mathematics and the selected cone
        are identical to the looped implementation.

        Numerically degenerate candidate cones are not physical and are
        discarded instead of raising from math.acos or a division by zero.

        Returns:
            min_3_cone: Smallest cone tangent to three atoms

        Raises:
            RuntimeError: If no cone encompassing all atoms is found.
        """
        # Set up vertex angles and normal vectors of all atoms and triples
        atoms = self._loop_atoms
        m = np.array([atom.cone.normal for atom in atoms])
        beta = np.array([atom.cone.angle for atom in atoms])
        triples = np.array(
            list(itertools.combinations(range(len(atoms)), 3)), dtype=int
        ).reshape(-1, 3)
        i, j, k = triples.T
        m_i, m_j, m_k = m[i], m[j], m[k]
        beta_i, beta_j, beta_k = beta[i], beta[j], beta[k]

        # Set up angles between atom vectors
        beta_ij = np.arccos(np.clip(np.einsum("ta,ta->t", m_i, m_j), -1, 1))

        # Set up matrices
        u = np.stack([np.cos(beta_i), np.cos(beta_j), np.cos(beta_k)], axis=1)
        v = np.stack([np.sin(beta_i), np.sin(beta_j), np.sin(beta_k)], axis=1)
        N = np.stack(
            [np.cross(m_j, m_k), np.cross(m_k, m_i), np.cross(m_i, m_j)], axis=2
        )
        P = np.einsum("tab,tac->tbc", N, N)
        gamma = np.einsum("ta,ta->t", m_i, np.cross(m_j, m_k))

        # Set up coefficients of quadratic equation
        A = np.einsum("ta,tab,tb->t", u, P, u)
        B = np.einsum("ta,tab,tb->t", v, P, v)
        C = np.einsum("ta,tab,tb->t", u, P, v)
        D = gamma**2

        # Solve quadratic equation in cos(2 * alpha) in closed form
        p2 = (A - B) ** 2 + 4 * C**2
        p1 = 2 * (A - B) * (A + B - 2 * D)
        p0 = (A + B - 2 * D) ** 2 - 4 * C**2
        with np.errstate(divide="ignore", invalid="ignore"):
            sqrt_discriminant = np.sqrt((p1**2 - 4 * p2 * p0).astype(complex))
            roots = np.stack(
                [
                    (-p1 + sqrt_discriminant) / (2 * p2),
                    (-p1 - sqrt_discriminant) / (2 * p2),
                ],
                axis=1,
            )
        cos_2_alpha = roots.real.copy()
        cos_2_alpha[np.isclose(cos_2_alpha, 1, rtol=1e-9, atol=0.0)] = 1
        cos_2_alpha[np.isclose(cos_2_alpha, -1, rtol=1e-9, atol=0.0)] = -1
        # A root marginally outside [-1, 1] is a boundary cone displaced by
        # solver noise (near a double root the error grows with the square
        # root of the coefficient error), e.g. a 180-degree cone. Clamp it to
        # the domain so it stays a candidate; roots further outside are not
        # physical cones and are excluded through `valid` below.
        outside = (np.abs(cos_2_alpha) > 1) & (np.abs(cos_2_alpha) <= 1 + 1e-7)
        cos_2_alpha[outside] = np.sign(cos_2_alpha[outside])
        valid = (
            (np.abs(roots.imag) < 1e10 * np.finfo(float).eps)
            & (np.abs(cos_2_alpha) <= 1)
            & np.isfinite(cos_2_alpha)
        )
        cos_2_alpha = np.clip(cos_2_alpha, -1, 1)

        # Four apex angle candidates per triple: acos(x) / 2 and
        # (2 * pi - acos(x)) / 2 for each root x
        acos_roots = np.arccos(cos_2_alpha)
        angles = (
            np.stack(
                [
                    acos_roots[:, 0],
                    2 * np.pi - acos_roots[:, 0],
                    acos_roots[:, 1],
                    2 * np.pi - acos_roots[:, 1],
                ],
                axis=1,
            )
            / 2
        )

        # Test roots and keep only the two most physical per triple
        cos_angles = np.cos(angles)
        sin_angles = np.sin(angles)
        D_tests = np.abs(
            A[:, None] * cos_angles**2
            + B[:, None] * sin_angles**2
            + 2 * C[:, None] * sin_angles * cos_angles
            - D[:, None]
        )
        D_tests[~np.repeat(valid, 2, axis=1)] = np.inf
        physical = np.argsort(D_tests, axis=1, kind="stable")[:, :2]
        t_idx = np.repeat(np.arange(len(triples)), 2)
        keep_valid = np.repeat(valid, 2, axis=1)[t_idx, physical.ravel()]
        alpha = angles[t_idx, physical.ravel()][keep_valid]
        t_idx = t_idx[keep_valid]

        # Calculate normal vectors of the cones for the physical angles
        m_i, m_j = m_i[t_idx], m_j[t_idx]
        beta_i, beta_j, beta_ij = beta_i[t_idx], beta_j[t_idx], beta_ij[t_idx]
        cross_ij = N[:, :, 2][t_idx]
        with np.errstate(divide="ignore", invalid="ignore"):
            sin_beta_ij_sq = np.sin(beta_ij) ** 2
            a_ij = (
                np.cos(alpha - beta_i) - np.cos(alpha - beta_j) * np.cos(beta_ij)
            ) / sin_beta_ij_sq
            b_ij = (
                np.cos(alpha - beta_j) - np.cos(alpha - beta_i) * np.cos(beta_ij)
            ) / sin_beta_ij_sq
            # Set c_ij squared to 0 if negative due to numerical precision
            c_ij_sq = 1 - a_ij**2 - b_ij**2 - 2 * a_ij * b_ij * np.cos(beta_ij)
            c_ij = np.sqrt(np.clip(c_ij_sq, 0, None))
            p = np.einsum(
                "tab,tb->ta",
                N[t_idx],
                u[t_idx] * np.cos(alpha)[:, None] + v[t_idx] * np.sin(alpha)[:, None],
            )
            sign = np.sign(gamma[t_idx]) * np.sign(np.einsum("ta,ta->t", p, cross_ij))
            c_ij = np.where(np.sign(c_ij) != sign, -c_ij, c_ij)
            normals = (
                a_ij[:, None] * m_i
                + b_ij[:, None] * m_j
                + (c_ij / np.sin(beta_ij))[:, None] * cross_ij
            )

        # Get upper and lower bound to apex angle and remove cones outside
        upper_bound = self._get_upper_bound()
        lower_bound = self._max_2_cone.angle
        in_bounds = (alpha - lower_bound >= -1e-5) & (upper_bound - alpha >= -1e-5)
        alpha = alpha[in_bounds]
        normals = normals[in_bounds]
        t_idx = t_idx[in_bounds]

        # Keep only cones that encompass all atoms (Cone.is_inside, vectorized)
        cos_angle = normals @ m.T
        cos_angle = np.where(
            (1 - cos_angle > 0) & (1 - cos_angle < 1e-5), 1.0, cos_angle
        )
        contains_all = (
            alpha[:, None] - (beta[None, :] + np.arccos(np.clip(cos_angle, -1, 1)))
            > -1e-5
        ).all(axis=1)

        # Take the smallest cone that encompasses all atoms
        keep = np.flatnonzero(contains_all)
        if keep.size == 0:
            raise RuntimeError("Cone not found")
        idx = keep[int(np.argmin(alpha[keep]))]
        min_3_cone = Cone(
            float(alpha[idx]),
            [atoms[t] for t in triples[t_idx[idx]]],
            normals[idx],
        )

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
        coordinates = np.array([atom.coordinates for atom in self._atoms])
        radii = np.array([atom.radius for atom in self._atoms])
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
    u = np.array([math.cos(beta_i), math.cos(beta_j), math.cos(beta_k)])
    v = np.array([math.sin(beta_i), math.sin(beta_j), math.sin(beta_k)])
    N = np.array([np.cross(m_j, m_k), np.cross(m_k, m_i), np.cross(m_i, m_j)]).T
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
    # Clamp real roots marginally outside [-1, 1] to the domain boundary;
    # see the corresponding comment in ConeAngle._search_three_cones.
    outside = np.isreal(roots) & (np.abs(roots) > 1) & (np.abs(roots) <= 1 + 1e-7)
    roots[outside] = np.sign(roots[outside].real)

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
