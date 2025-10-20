"""Pyramidalization code."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
import functools
import itertools
from typing import Any, cast

import numpy as np
import scipy.spatial

from morfeus.io import read_geometry
from morfeus.typing import (
    Array1DBool,
    Array1DFloat,
    Array1DInt,
    Array2DFloat,
    ArrayLike1D,
    ArrayLike2D,
)
from morfeus.utils import get_connectivity_matrix


class Pyramidalization:
    """Calculates and stores results of pyramidalization and alpha angle.

    As described in Struct. Chem. 1991, 2, 107 and alternatively according to bond
    angles as in J. Comput. Chem. 2012, 33 (27), 2173–2179.

    Args:
        coordinates: Coordinates (Å)
        atom_index: Index of pyramidalized atom (1-indexed)
        neighbor_indices: Indices of neighbors to pyramidalized atom
        elements: Elements as atomic symbols or numbers
        radii: Covalent radii used to determine connectivity (Å)
        radii_type: Covalent radii type: 'pyykko'
        excluded_atoms: Indices of atoms to exclude
        method: Method for detecting neighbors: 'connectivity' or 'distance'. Ignored if
            neighbor_indices is given.

    Attributes:
        alpha: Average alpha angle (degrees)
        alphas: Alpha angles for all permutations of neighbors (degrees)
        neighbor_indices: Indices of neighbors to pyramidalized atom
        P_angle: Pyramidalization according to Gavrish
        P: Pyramidalization according to Radhakrishnan
    """

    alpha: float
    alphas: Array1DFloat
    neighbor_indices: list[int]
    P_angle: float
    P: float

    def __init__(  # noqa: C901
        self,
        coordinates: ArrayLike2D,
        atom_index: int,
        neighbor_indices: Sequence[int] | None = None,
        elements: Iterable[int] | Iterable[str] | None = None,
        radii: ArrayLike1D | None = None,
        radii_type: str = "pyykko",
        excluded_atoms: Sequence[int] | None = None,
        method: str = "distance",
        scale_factor: float = 1.2,
    ) -> None:
        if 0 in {
            atom_index,
            *(neighbor_indices or ()),
            *(excluded_atoms or ()),
        }:
            raise IndexError("Atom indices should not be 0 (1-indexed).")

        coordinates: Array2DFloat = np.array(coordinates)
        atom_coordinates = coordinates[atom_index - 1]

        if neighbor_indices is None:
            neighbor_indices = []
        else:
            neighbor_indices = list(neighbor_indices)

        if excluded_atoms is None:
            excluded_atoms = []
        excluded_atoms: Array1DBool = np.array(excluded_atoms, dtype=bool)

        # Get 3 closest neighbors
        if len(neighbor_indices) > 0:
            if len(neighbor_indices) != 3:
                raise Exception(f"Only {len(neighbor_indices)} neighbors.")
            neighbors: Array1DInt = np.array(neighbor_indices) - 1
        elif method == "distance":
            # Generate mask for excluded atoms
            mask: Array1DBool = np.zeros(len(coordinates), dtype=bool)
            mask[excluded_atoms - 1] = True
            mask[atom_index - 1] = True

            # Get three closest atoms not in the excluded atoms
            distances = scipy.spatial.distance.cdist(
                atom_coordinates.reshape(1, -1), coordinates
            ).reshape(-1)
            distances[mask] = np.inf
            neighbors = np.argsort(distances)[:3]
        elif method == "connectivity":
            # Construct connectivity matrix and get closest neighbors.
            if elements is None and radii is None:
                raise Exception("Connectivity requires elements or radii.")
            # if radii is None:
            #    radii = get_radii(elements, radii_type="pyykko")
            connectivity_matrix = get_connectivity_matrix(
                coordinates,
                elements=elements,
                radii=radii,
                radii_type=radii_type,
                scale_factor=scale_factor,
            )
            connected_atoms = np.where(connectivity_matrix[atom_index - 1, :])[0]
            neighbors = connected_atoms[~np.isin(connected_atoms, excluded_atoms - 1)]
            if len(neighbors) != 3:
                raise Exception(f"{len(neighbors)} neighbors. 3 expected.")

        # Get unit vectors between central atom and neighbors
        a = coordinates[neighbors[0]] - atom_coordinates
        a /= np.linalg.norm(a)
        b = coordinates[neighbors[1]] - atom_coordinates
        b /= np.linalg.norm(b)
        c = coordinates[neighbors[2]] - atom_coordinates
        c /= np.linalg.norm(c)

        # Calculate alpha for all permutations
        alphas: list[float] = []
        vectors: list[Array1DFloat] = []
        cos_alphas: list[float] = []
        thetas: list[float] = []
        for v_1, v_2, v_3 in itertools.permutations([a, b, c], 3):
            # Calculate cos_alpha
            normal: Array1DFloat = np.cross(v_1, v_2)
            normal /= np.linalg.norm(normal)
            cos_alpha = np.dot(v_3, normal)

            # Test if normal vector is colinear with v_3
            if cos_alpha < 0:
                continue
            alpha = np.arccos(cos_alpha)

            # Check for "acute" pyramid and correct angle
            v_1_2 = v_1 + v_2
            v_1_2 /= np.linalg.norm(v_1_2)
            cos_angle = np.dot(v_1_2, v_3)
            if cos_angle > 0:
                alpha = -alpha
            alphas.append(alpha)
            cos_alphas.append(cos_alpha)
            vectors.append((v_1, v_2))

            # Calculate theta angle
            cos_theta = np.dot(v_1, v_2)
            theta = np.rad2deg(np.arccos(cos_theta))
            thetas.append(theta)

        # Calculate P
        v_1, v_2 = vectors[0]
        sin_theta = np.linalg.norm(np.cross(v_1, v_2))
        sin_theta = cast(float, sin_theta)
        P = sin_theta * cos_alphas[0]

        # Correct P if pyramid is "acute" on average
        if np.mean(alphas) < 0:
            P = 2 - P

        # Calculate P according to Gavrish method
        P_angle = np.sqrt(360 - sum(thetas))

        # Store attributes
        self.P = P
        self.P_angle = P_angle
        self.alpha = np.rad2deg(np.mean(alphas))
        self.alphas = np.rad2deg(alphas)
        self.neighbor_indices = (neighbors + 1).tolist()

    def print_report(self) -> None:
        """Print report of results."""
        print(f"P: {self.P:.3f}")
        print(f"P_angle: {self.P_angle:.3f}")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({round(self.P, 3)!r})"


def cli(file: str) -> Any:
    """CLI for pyramidalization.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(Pyramidalization, coordinates, elements=elements)
