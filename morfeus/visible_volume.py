"""Visible volume code."""

from __future__ import annotations

from collections.abc import Iterable
import functools
from typing import Any

import numpy as np

from morfeus.geometry import Cone
from morfeus.io import read_geometry
from morfeus.sasa import SASA
from morfeus.typing import (
    Array1DAny,
    Array1DFloat,
    Array1DInt,
    Array2DFloat,
    ArrayLike1D,
    ArrayLike2D,
)
from morfeus.utils import check_distances, convert_elements, get_radii


class VisibleVolume:
    """Calculates and stores visible volume and area.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        metal_index: Index of metal atom (1-indexed)
        include_hs: Whether to include H atoms in the calculation
        radii: Atomic radii (Å)
        radii_type: Radii type: 'alvarez', 'bondi', 'crc', 'pyykko', 'rahm' or 'truhlar'
        radius: Radius of sphere to divide proximal and distal (Å)
        density: Area per point on atom surface (Å²)

    Attributes:
        distal_area: Distal area (Å²)
        distal_visible_area: Distal visible area (Å²)
        distal_visible_volume: Distal visible volume (Å³)
        distal_volume: Distal volume (Å³)
        invisible_area: Invisible area (Å²)
        invisible_volume: Invisible volume (Å³)
        proximal_area: Proximal area (Å²)
        proximal_visible_area: Proximal visible area (Å²)
        proximal_visible_volume: Proximal visible volume (Å³)
        proximal_volume: Proximal volume (Å³)
        total_area: Total area (Å²)
        total_volume: Total volume (Å³)
        visible_area: Visible area (Å²)
        visible_volume: Visible volume (Å³)
    """

    distal_area: float
    distal_visible_area: float
    distal_visible_volume: float
    distal_volume: float
    invisible_area: float
    invisible_volume: float
    proximal_area: float
    proximal_visible_area: float
    proximal_visible_volume: float
    proximal_volume: float
    total_area: float
    total_volume: float
    visible_area: float
    visible_volume: float

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        metal_index: int,
        include_hs: bool = True,
        radii: ArrayLike1D | None = None,
        radii_type: str = "pyykko",
        radius: float = 3.5,
        density: float = 0.01,
    ) -> None:
        # Set up arrays and get radii
        elements: Array1DInt = np.array(convert_elements(elements, output="numbers"))
        coordinates: Array2DFloat = np.array(coordinates)
        if radii is None:
            radii = get_radii(elements, radii_type=radii_type)
        radii: Array1DFloat = np.array(radii)

        # Check so that no atom is within vdW distance of metal atom
        within = check_distances(elements, coordinates, metal_index, radii=radii)
        if len(within) > 0:
            atom_string = " ".join([str(i) for i in within])
            raise Exception("Atoms within vdW radius of central atom:", atom_string)

        # Center coordinate system around metal and remove it
        coordinates -= coordinates[metal_index - 1]
        elements: Array1DFloat = np.delete(elements, metal_index - 1)
        coordinates: Array2DFloat = np.delete(coordinates, metal_index - 1, axis=0)
        radii: Array1DFloat = np.delete(radii, metal_index - 1)

        # Remove H atoms
        if include_hs is False:
            h_mask = elements == 1
            elements = elements[~h_mask]
            coordinates = coordinates[~h_mask]
            radii = radii[~h_mask]

        # Construct SASA object
        sasa = SASA(
            elements,
            coordinates,
            radii=radii,
            radii_type=radii_type,
            probe_radius=0.0,
            density=density,
        )

        # Set invisible and proximal masks for each atom
        atoms = sasa._atoms
        for atom in atoms:
            atom.invisible_mask = np.zeros(len(atom.accessible_points), dtype=bool)
            atom.proximal_mask = np.linalg.norm(atom.accessible_points, axis=1) < radius

        # Check points on other atoms against cone for each atom
        atom_coordinates: Array2DFloat = np.array([atom.coordinates for atom in atoms])
        atoms: Array1DAny = np.array(atoms)
        for atom in atoms:
            # Calculate distances to other atoms
            atom.get_cone()
            cone = Cone(atom.cone.angle, [atom.index], atom.cone.normal)
            atom_dist = np.linalg.norm(atom.coordinates)
            other_distances = np.dot(atom_coordinates, cone.normal)

            # Check whether points are (1) within cone and (2) beyond atom
            # center
            check_atoms = atoms[other_distances >= (atom_dist - atom.radius)]
            for check_atom in check_atoms:
                if check_atom == atom:
                    continue
                is_inside_mask = cone.is_inside_points(
                    check_atom.accessible_points, method="cross"
                )
                is_beyond_mask = (
                    np.dot(check_atom.accessible_points, cone.normal) >= atom_dist
                )
                invisible_mask = np.logical_and(is_inside_mask, is_beyond_mask)
                check_atom.invisible_mask = np.logical_or(
                    check_atom.invisible_mask, invisible_mask
                )

        # Calculate visible, invisible and proximal_visible volume
        visible_volume = 0
        invisible_volume = 0
        proximal_visible_volume = 0
        proximal_volume = 0
        visible_area = 0
        invisible_area = 0
        proximal_visible_area = 0
        proximal_area = 0
        for atom in atoms:
            point_areas = atom.point_areas[atom.accessible_mask]
            point_volumes = atom.point_volumes[atom.accessible_mask]

            invisible_volume += point_volumes[atom.invisible_mask].sum()
            visible_volume += point_volumes[~atom.invisible_mask].sum()
            proximal_visible_volume += point_volumes[
                np.logical_and(~atom.invisible_mask, atom.proximal_mask)
            ].sum()
            proximal_volume += point_volumes[atom.proximal_mask].sum()

            invisible_area += point_areas[atom.invisible_mask].sum()
            visible_area += point_areas[~atom.invisible_mask].sum()
            proximal_visible_area += point_areas[
                np.logical_and(~atom.invisible_mask, atom.proximal_mask)
            ].sum()
            proximal_area += point_areas[atom.proximal_mask].sum()

        # Store attributes
        self.total_volume = sasa.volume
        self.proximal_volume = proximal_volume
        self.distal_volume = sasa.volume - proximal_volume
        self.invisible_volume = invisible_volume
        self.visible_volume = visible_volume
        self.proximal_visible_volume = proximal_visible_volume
        self.distal_visible_volume = visible_volume - proximal_visible_volume

        self.total_area = sasa.area
        self.proximal_area = proximal_area
        self.distal_area = sasa.area - proximal_area
        self.invisible_area = invisible_area
        self.visible_area = visible_area
        self.proximal_visible_area = proximal_visible_area
        self.distal_visible_area = visible_area - proximal_visible_area

        self._atoms = atoms

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


def cli(file: str) -> Any:
    """CLI for visible volume.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(VisibleVolume, elements, coordinates)
