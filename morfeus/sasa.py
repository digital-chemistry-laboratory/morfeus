"""Solvent accessible surface area code."""

from __future__ import annotations

from collections.abc import Iterable
import functools
import typing
from typing import Any

import numpy as np
import scipy.spatial

from morfeus.data import atomic_symbols, jmol_colors
from morfeus.geometry import Atom, Sphere
from morfeus.io import read_geometry
from morfeus.typing import Array1DFloat, Array2DFloat, ArrayLike1D, ArrayLike2D
from morfeus.utils import convert_elements, get_radii, Import, requires_dependency

if typing.TYPE_CHECKING:
    from matplotlib.colors import hex2color
    import pyvista as pv
    from pyvistaqt import BackgroundPlotter


class SASA:
    """Performs and stores results of solvent accessible surface area calculations.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        radii: VdW radii (Å)
        radii_type: Choice of vdW radii: 'bondi' or 'crc' (default)
        probe_radius: Radius of probe atom (Å)
        density: Area per point (Å²) on the vdW surface

    Attributes:
        area: Area of the solvent accessible surface.
        atom_areas: Atom areas (starting from 1)
        volume: Volume of the solvent accessible surface
    """

    area: float
    atom_areas: dict[int, float]
    volume: float
    _atoms: list[Atom]
    _density: float
    _probe_radius: float

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        radii: ArrayLike1D | None = None,
        radii_type: str = "crc",
        probe_radius: float = 1.4,
        density: float = 0.01,
    ) -> None:
        # Converting elements to atomic numbers if the are symbols
        elements = convert_elements(elements, output="numbers")
        coordinates: Array2DFloat = np.array(coordinates)

        # Getting radii if they are not supplied
        if radii is None:
            radii = get_radii(elements, radii_type=radii_type)

        # Increment the radii with the probe radius
        radii: Array1DFloat = np.array(radii)
        radii = radii + probe_radius

        # Construct list of atoms
        atoms = []
        for i, (coordinate, radius, element) in enumerate(
            zip(coordinates, radii, elements), start=1
        ):
            atom = Atom(element, coordinate, radius, i)
            atoms.append(atom)

        # Set up attributes
        self._atoms = atoms
        self._density = density
        self._probe_radius = probe_radius

        # Determine accessible and occluded points for each atom
        self._determine_accessible_points()

        # Calculate atom areas and volumes
        self._calculate()

    def _calculate(self) -> None:
        """Calculate solvent accessible surface area and volume."""
        for atom in self._atoms:
            # Get number of points of eache type
            n_accessible = len(atom.accessible_points)
            n_occluded = len(atom.occluded_points)
            n_points = len(atom.accessible_points) + len(atom.occluded_points)

            # Calculate part occluded and accessible
            ratio_occluded = n_occluded / n_points
            ratio_accessible = 1 - ratio_occluded

            # Calculate area
            area = 4 * np.pi * atom.radius**2 * ratio_accessible
            atom.area = area
            atom.point_areas = np.zeros(n_points)
            if n_accessible > 0:
                atom.point_areas[atom.accessible_mask] = atom.area / n_accessible

            # Center accessible points and normalize
            centered_points = np.array(atom.accessible_points) - atom.coordinates
            centered_points /= np.linalg.norm(centered_points, axis=1).reshape(-1, 1)

            # Add accessible points
            accessible_summed = np.sum(centered_points, axis=0)

            # Calculate volume
            volume = (4 * np.pi / 3 / n_points) * (
                atom.radius**2 * np.dot(atom.coordinates, accessible_summed)
                + atom.radius**3 * n_accessible
            )
            atom.volume = volume
            atom.point_volumes = np.zeros(n_points)
            if n_accessible > 0:
                atom.point_volumes[atom.accessible_mask] = atom.volume / n_accessible

        # Set up attributes
        self.atom_areas = {atom.index: atom.area for atom in self._atoms}
        self.area = sum([atom.area for atom in self._atoms])
        self.volume = sum([atom.volume for atom in self._atoms])

    def _determine_accessible_points(self) -> None:
        """Determine occluded and accessible points of each atom."""
        # Based on distances to all other atoms (brute force).
        for atom in self._atoms:
            # Construct sphere for atom
            sphere = Sphere(atom.coordinates, atom.radius, density=self._density)
            atom.points = sphere.points

            # Select atoms that are at a distance less than the sum of radii
            # !TODO can be vectorized
            test_atoms = []
            for test_atom in self._atoms:
                if test_atom is not atom:
                    distance = scipy.spatial.distance.euclidean(
                        atom.coordinates, test_atom.coordinates
                    )
                    radii_sum = atom.radius + test_atom.radius
                    if distance < radii_sum:
                        test_atoms.append(test_atom)

            # Select coordinates and radii for other atoms
            test_coordinates = [test_atom.coordinates for test_atom in test_atoms]
            test_radii = [test_atom.radius for test_atom in test_atoms]
            test_radii: Array1DFloat = np.array(test_radii).reshape(-1, 1)

            # Get distances to other atoms and subtract radii
            if test_coordinates:
                distances = scipy.spatial.distance.cdist(
                    test_coordinates, sphere.points
                )
                distances -= test_radii
                # Take smallest distance and perform check
                min_distances = np.min(distances, axis=0)
                atom.occluded_mask = min_distances < 0
                atom.accessible_mask = ~atom.occluded_mask
            else:
                atom.occluded_mask = np.zeros(len(atom.points), dtype=bool)
                atom.accessible_mask = np.ones(len(atom.points), dtype=bool)
            atom.occluded_points = sphere.points[atom.occluded_mask]
            atom.accessible_points = sphere.points[atom.accessible_mask]

    def print_report(self, verbose: bool = False) -> None:
        """Print report of results.

        Args:
            verbose: Whether to print atom areas
        """
        print(f"Probe radius (Å): {self._probe_radius}")
        print(f"Solvent accessible surface area (Å²): {self.area:.1f}")
        print("Volume inside solvent accessible surface (Å³): " f"{self.volume:.1f}")
        if verbose:
            print(f"{'Symbol':<10s}{'Index':<10s}{'Area (Å²)':<10s}")
            for atom, (i, area) in zip(self._atoms, self.atom_areas.items()):
                symbol = atomic_symbols[atom.element]
                print(f"{symbol:<10s}{i:<10d}{area:<10.1f}")

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
        point_color: str = "steelblue",
        opacity: float = 0.25,
        size: float = 1,
    ) -> None:
        """Draw a 3D representation.

        Draws the molecule with the solvent accessible surface area.

        Args:
            atom_scale: Scaling factor for atom size
            background_color: Background color for plot
            point_color: Color of surface points
            opacity: Point opacity
            size: Point size
        """
        # Set up plotter
        p = BackgroundPlotter()
        p.set_background(background_color)

        # Draw molecule
        for atom in self._atoms:
            color = hex2color(jmol_colors[atom.element])
            radius = atom.radius * atom_scale - self._probe_radius
            sphere = pv.Sphere(center=list(atom.coordinates), radius=radius)
            p.add_mesh(sphere, color=color, opacity=1, name=str(atom.index))

        # Draw surface points
        surface_points: Array2DFloat = np.vstack(
            [atom.accessible_points for atom in self._atoms]
        )
        p.add_points(
            surface_points, color=point_color, opacity=opacity, point_size=size
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


def cli(file: str) -> Any:
    """CLI for solvent accessible surface area.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(SASA, elements, coordinates)
