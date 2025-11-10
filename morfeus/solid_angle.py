"""Solid angle code."""

from __future__ import annotations

from collections.abc import Iterable
import functools
import typing
from typing import Any

import numpy as np

from morfeus.data import jmol_colors
from morfeus.geometry import Atom, Sphere
from morfeus.io import read_geometry
from morfeus.typing import Array1DFloat, Array2DFloat, ArrayLike1D, ArrayLike2D
from morfeus.utils import convert_elements, get_radii, Import, requires_dependency

if typing.TYPE_CHECKING:
    from matplotlib.colors import hex2color
    import pyvista as pv
    from pyvistaqt import BackgroundPlotter


class SolidAngle:
    """Calculates and stores the results of a solid angle calculation.

    As described in J. Chem. Theory Comput. 2013, 9 (12), 5734-5744 and Dalton Trans.
    2006, 33, 3991.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        metal_index: Index of central atom (1-indexed)
        radii: vdW radii (Å)
        radii_type: Type of vdW radii: 'alvarez', 'bondi', 'crc' or 'truhlar'
        density: Area per point (Å²) on the sphere surface

    Attributes:
        cone_angle: Solid cone angle (degrees)
        solid_angle: Solid angle (steradians)
        G: G parameter
    """

    G: float
    cone_angle: float
    solid_angle: float
    _atoms: list[Atom]

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        metal_index: int,
        radii: ArrayLike1D | None = None,
        radii_type: str = "crc",
        density: float = 0.001,
    ) -> None:
        if metal_index == 0:
            raise IndexError("Atom indices should not be 0 (1-indexed).")

        # Convert elements to atomic numbers if they are symbols
        elements = convert_elements(elements, output="numbers")
        coordinates: Array2DFloat = np.array(coordinates)

        # Get radii if they are not supplied
        if radii is None:
            radii = get_radii(elements, radii_type=radii_type)
        radii: Array1DFloat = np.asarray(radii)
        coordinates_metal = coordinates[metal_index - 1, :]
        coordinates -= coordinates_metal

        # Get list of atoms as Atom objects
        atoms = []
        for i, (element, coords, radius) in enumerate(
            zip(elements, coordinates, radii), start=1
        ):
            if i == metal_index:
                continue
            atom = Atom(element, coords, radius, i)
            atoms.append(atom)

        # Construct sphere and check which points on sphere are within atom cones
        sphere = Sphere(coordinates_metal, 1.0, density=density)

        mask = np.zeros(len(sphere.points), dtype=bool)
        for atom in atoms:
            atom.get_cone()
            is_inside = atom.cone.is_inside_points(sphere.points)
            mask = np.logical_or(mask, is_inside)

        # Calculate solid angle, cone angle and G
        n_occluded = sum(mask)
        n_total = len(mask)
        ratio_occluded = n_occluded / n_total
        area_occluded = ratio_occluded * sphere.area
        omega = area_occluded / sphere.radius**2
        theta = 2 * np.arccos(1 - omega / 2 / np.pi)
        G = 100 * omega / 4 / np.pi

        self.solid_angle = omega
        self.cone_angle = np.rad2deg(theta)
        self.G = G
        self._atoms = atoms
        self._sphere = sphere
        self._mask = mask

    def print_report(self) -> None:
        """Prints report of results."""
        print(f"Solid angle (sr): {self.solid_angle:.3f}")
        print(f"Cone angle (°): {self.cone_angle:.3f}")
        print(f"G: {self.G:.3f}")

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
        atom_scale: float = 1.0,
        background_color: str = "white",
        point_color: str = "steelblue",
        opacity: float = 0.25,
        size: float = 2.0,
    ) -> None:
        """Draw a 3D representation.

        Draws the molecule with the ligand shadow

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
            radius = atom.radius * atom_scale
            sphere = pv.Sphere(center=list(atom.coordinates), radius=radius)
            p.add_mesh(sphere, color=color, opacity=1, name=str(atom.index))

        coordinates: Array2DFloat = np.vstack(
            [atom.coordinates for atom in self._atoms]
        )
        radii: Array1DFloat = np.vstack([atom.radius for atom in self._atoms])
        distances = np.linalg.norm(coordinates, axis=1) + radii
        dist_max = np.max(distances)
        sphere_points = self._sphere.points[self._mask] * (dist_max + 0.1)

        # Draw surface points
        p.add_points(
            sphere_points,
            color=point_color,
            opacity=opacity,
            point_size=size,
        )

        sphere = pv.Sphere(center=[0.0, 0.0, 0.0], radius=0.5)
        p.add_mesh(sphere, color="lightsteelblue", opacity=0.75, name=str(1))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


def cli(file: str) -> Any:
    """CLI for solid angle.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(SolidAngle, elements, coordinates)
