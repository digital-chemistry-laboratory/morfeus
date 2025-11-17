"""Tau4 code."""

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


class Tau4:
    """Calculates and stores results of Tau4.
    Gives both the original Tau4 and the improved Tau4.
    Original Tau4 as described by: ....
    Improved Tau4 as described by: ....

    Args:
        coordinates: Coordinates (Å)
        atom_index: Index of the metal center (1-indexed)
        neighbor_indices: Indices of neighbors to metal center
        elements: Elements as atomic symbols or numbers
        radii: Covalent radii used to determine connectivity (Å)
        radii_type: Covalent radii type: 'pyykko'
        excluded_atoms: Indices of atoms to exclude
        method: Method for detecting neighbors: 'connectivity' or 'distance'. Ignored if
            neighbor_indices is given.

    Attributes:
        tau4: Original Tau4 distortion term
        imp_tau4: Improved Tau4 distortion term
        neighbor_indices: Indices of neighbors to metal center
    """
    
    tau4: float
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
            
        #get 4 closest neighbors
        if len(neighbor_indices) > 0:
            if len(neighbor_indices) != 4:
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
            neighbors = np.argsort(distances)[:4]
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
            
            if len(neighbors) != 4:
                raise Exception(f"{len(neighbors)} neighbors. 4 expected.")
        
        vectors = coordinates[neighbors] - atom_coordinates
        norm_vectors = vectors / np.linalg.norm(vectors, axis=1).reshape(-1, 1)
        
        angles =    [np.rad2deg(
            np.arctan2(
            np.linalg.norm(np.cross(v_1, v_2)),
                np.dot(v_1, v_2),
            )) for v_1, v_2 in itertools.combinations(norm_vectors, 2)]
        
        beta, alpha = np.sort(angles)[::-1][:2]
        
        tetrahedral_angle = np.rad2deg(np.arccos(-1/3))
        
        tau4 = (360 - (alpha + beta)) / (360 - 2 * tetrahedral_angle)
        imp_tau4 = ((beta - alpha) / (360 - tetrahedral_angle) + (180 - beta) / (180 - tetrahedral_angle))
        
        #store attributes
        self.tau4 = tau4
        self.imp_tau4 = imp_tau4
        self.neighbor_indices = (neighbors + 1).tolist()
        
    def print_report(self) -> None:
        """Print report of results."""
        print(f"Tau4: {self.tau4:.3f}")
        print(f"Improved Tau4: {self.imp_tau4:.3f}")
        
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({round(self.tau4, 3)!r})"
        return f"{self.__class__.__name__}({round(self.imp_tau4, 3)!r})"


def cli(file: str) -> Any:
    """CLI for Tau4.
    
    Args:
        file: Geometry file
        
    Returns:
        Partially instantiated Tau4 class
    """

    elements, coordinates = read_geometry(file)
    return functools.partial(Tau4, coordinates, elements=elements)
    
    