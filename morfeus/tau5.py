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


class Tau5:
    """Calculates and stores results of Tau5.
    As described by: ....

    Args:
        cooridinates: Coordinates (Å)
        atom_index: Index of the metal center (1-indexed)
        neighbor_indices: Indices of neighbors to metal center
        elements: Elements as atomic symbols or numbers
        radii: Covalent radii used to determine connectivity (Å)
        radii_type: Covalent radii type: 'pyykko'
        excluded_atoms: Indices of atoms to exclude
        method: Method for detecting neighbors: 'connectivity' or 'distance'. Ignored if
            neighbor_indices is given.

    Attributes:
        tau5: Tau5 distortion term
        neighbor_indices: Indices of neighbors to metal center
    """

    og_tau4: float
    imp_tau4: float
    neighbor_indices: list[int]
      
    def __init__(
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
        coordinates: Array2DFloat = np.array(coordinates)
        atom_coordinates = coordinates[atom_index - 1]
        
        if neighbor_indices is None:
            neighbor_indices = []
        else:
            neighbor_indices = list(neighbor_indices)
            
        if excluded_atoms is None:
            excluded_atoms = []
        excluded_atoms: Array1DBool = np.array(excluded_atoms, dtype=bool)
            
        #get 5 closest neighbors
        if len(neighbor_indices) > 0:
            if len(neighbor_indices) != 5:
                raise Exception(f"Only {len(neighbor_indices)} neighbors.")
            neighbors: Array1DInt = np.array(neighbor_indices) - 1
        elif method == "distance":
            # Generate mask for excluded atoms
            mask: Array1DBool = np.zeros(len(coordinates), dtype=bool)
            mask[excluded_atoms - 1] = True
            mask[atom_index - 1] = True
            
            #Get four closest atoms not in the excluded atoms
            distances = scipy.spatial.distance.cdist(
                atom_coordinates.reshape(1, -1), coordinates
                ).reshape(-1)
            
            distances[mask] = np.inf
            neighbors = np.argsort(distances)[:5]
        elif method == "connectivity":
            #construct connectivity matrix and get closest neighbors
            if elements is None and radii is None:
                raise Exception("Connectivity requires elements or radii.")
            connectivity_matrix = get_connectivity_matrix(
                coordinates,
                elements=elements,
                radii=radii,
                radii_type=radii_type,
                scale_factor=scale_factor,
                )
            connected_atoms = np.where(connectivity_matrix[atom_index - 1, :])[0]
            neighbors = connected_atoms[~np.isin(connected_atoms, excluded_atoms - 1)]
            
            if len(neighbors) != 5:
                raise Exception(f"{len(neighbors)} neighbors. 4 expected.")
        
        vectors = coordinates[neighbors] - atom_coordinates
        norm_vectors = vectors / np.linalg.norm(vectors, axis=1).reshape(-1, 1)
        
        angles =    [np.rad2deg(
            np.arctan2(
            np.linalg.norm(np.cross(v_1, v_2)),
                np.dot(v_1, v_2),
            )) for v_1, v_2 in itertools.combinations(norm_vectors, 2)]
        
        beta, alpha = np.sort(angles)[::-1][:2]

        tau5 = (beta - alpha) / 60
        
        #store attributes
        self.tau5 = tau5
        self.neighbor_indices = (neighbors + 1).tolist()
        
    def print_report(self) -> None:
        """Print report of results."""
        print(f"Tau5: {self.tau5:.3f}")
        
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({round(self.tau5, 3)!r})"


def cli(file: str) -> Any:
    """CLI for Tau5.
    
    Args:
        file: Geometry file

    Returns:
        Partially instantiated Tau5 class
    """

    elements, coordinates = read_geometry(file)
    return functools.partial(Tau5, coordinates, elements=elements)