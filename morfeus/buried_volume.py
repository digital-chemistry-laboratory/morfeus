"""Buried volume code."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
import copy
import functools
import itertools
import math
import typing
from typing import Any
import warnings

import numpy as np
import scipy.spatial

from morfeus.data import jmol_colors
from morfeus.geometry import Atom, kabsch_rotation_matrix, rotate_coordinates, Sphere
from morfeus.io import read_geometry
from morfeus.sasa import SASA
from morfeus.typing import (
    Array1DBool,
    Array1DFloat,
    Array2DFloat,
    ArrayLike1D,
    ArrayLike2D,
)
from morfeus.utils import convert_elements, get_radii, Import, requires_dependency

if typing.TYPE_CHECKING:
    from matplotlib.colors import hex2color
    import matplotlib.pyplot as plt
    import pyvista as pv
    from pyvista import BackgroundPlotter

# Quadrant and octant signs taken from
# https://en.wikipedia.org/wiki/Octant_(solid_geometry)
QUADRANT_SIGNS: dict[int, str] = {
    1: "+,+",
    2: "-,+",
    3: "-,-",
    4: "+,-",
}
"""Coventional signs for quadrants."""

OCTANT_SIGNS: dict[int, str] = {
    0: "+,+,+",
    1: "-,+,+",
    3: "+,-,+",
    2: "-,-,+",
    7: "+,+,-",
    6: "-,+,-",
    4: "+,-,-",
    5: "-,-,-",
}
"""Conventional signs for octants."""

QUADRANT_NAMES: dict[int, str] = {
    1: "NE",
    2: "NW",
    3: "SW",
    4: "SE",
}
"""Conventional names for quadrants."""

# Maps octants to quadrants
QUADRANT_OCTANT_MAP: dict[int, tuple[int, int]] = {
    1: (0, 7),
    2: (1, 6),
    3: (2, 5),
    4: (3, 4),
}
"""Map from quadrants to octants."""


class BuriedVolume:
    """Performs and stores the results of a buried volume calculation.

    Algorithm similar as to described in Organometallics 2016, 35, 2286.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        metal_index: Index of metal atom (1-indexed)
        excluded_atoms: Indices of atoms to exclude (1-indexed). Metal atom is always
            excluded and does not have to be given here.
        radii: vdW radii (Å)
        include_hs: Whether to include H atoms in the calculation
        radius: Radius of sphere (Å)
        radii_type: Type of radii to use: 'alvarez', 'bondi', 'crc' or 'truhlar'
        radii_scale: Scaling factor for radii
        density: Volume per point in the sphere (Å³)
        z_axis_atoms: Atom indices for deciding orientation of z axis (1-indexed)
        xz_plane_atoms: Atom indices for deciding orientation of xz plane (1-indexed)

    Attributes:
        buried_volume: Buried volume of sphere (Å³)
        distal_volume: Distal volume of ligand (Å³)
        fraction_buried_volume: Fraction buried volume of sphere
        free_volume: Free volume of sphere (Å³)
        octants: Results for octant analysis
        quadrants: Results for quadrant analysis
    """

    buried_volume: float
    distal_volume: float
    fraction_buried_volume: float
    free_volume: float
    molecular_volume: float
    octants: dict[str, dict[int, float]]
    quadrants: dict[str, dict[int, float]]
    _all_coordinates: Array2DFloat
    _atoms: list[Atom]
    _buried_points: Array2DFloat
    _density: float
    _excluded_atoms: set[int]
    _free_points: Array2DFloat
    _octant_limits: dict[
        int, tuple[tuple[float, float], tuple[float, float], tuple[float, float]]
    ]
    _sphere: Sphere

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        metal_index: int,
        excluded_atoms: Sequence[int] | None = None,
        radii: ArrayLike1D | None = None,
        include_hs: bool = False,
        radius: float = 3.5,
        radii_type: str = "bondi",
        radii_scale: float = 1.17,
        density: float = 0.001,
        z_axis_atoms: Sequence[int] | None = None,
        xz_plane_atoms: Sequence[int] | None = None,
    ) -> None:
        if 0 in {
            metal_index,
            *(excluded_atoms or ()),
            *(z_axis_atoms or ()),
            *(xz_plane_atoms or ()),
        }:
            raise IndexError("Atom indices should not be 0 (1-indexed).")

        # Get center and and reortient coordinate system
        coordinates: Array2DFloat = np.array(coordinates)
        center = coordinates[metal_index - 1]
        coordinates -= center

        if excluded_atoms is None:
            excluded_atoms = []
        excluded_atoms = set(excluded_atoms)

        if metal_index not in excluded_atoms:
            excluded_atoms.add(metal_index)

        if z_axis_atoms is not None and xz_plane_atoms is not None:
            z_axis_coordinates = coordinates[np.array(z_axis_atoms) - 1]
            z_point = np.mean(z_axis_coordinates, axis=0)

            xz_plane_coordinates = coordinates[np.array(xz_plane_atoms) - 1]
            xz_point = np.mean(xz_plane_coordinates, axis=0)

            v_1 = z_point - center
            v_2 = xz_point - center
            v_3: Array1DFloat = np.cross(v_2, v_1)
            real: Array2DFloat = np.vstack([v_1, v_3])
            real /= np.linalg.norm(real, axis=1).reshape(-1, 1)
            ref_1 = np.array([0.0, 0.0, -1.0])
            ref_2 = np.array([0.0, 1.0, 0.0])
            ref = np.vstack([ref_1, ref_2])
            R = kabsch_rotation_matrix(real, ref, center=False)
            coordinates = (R @ coordinates.T).T
        elif z_axis_atoms is not None:
            z_axis_coordinates = coordinates[np.array(z_axis_atoms) - 1]
            z_point = np.mean(z_axis_coordinates, axis=0)
            v_1 = z_point - center
            v_1 = v_1 / np.linalg.norm(v_1)
            coordinates = rotate_coordinates(coordinates, v_1, np.array([0, 0, -1]))

        self._z_axis_atoms = z_axis_atoms
        self._xz_plane_atoms = xz_plane_atoms

        # Save density and coordinates for steric map plotting.
        self._density = density
        self._all_coordinates = coordinates

        # Converting element ids to atomic numbers if the are symbols
        elements = convert_elements(elements, output="numbers")

        # Getting radii if they are not supplied
        if radii is None:
            radii = get_radii(elements, radii_type=radii_type, scale=radii_scale)
        radii: Array1DFloat = np.array(radii)

        # Get list of atoms as Atom objects
        atoms = []
        for i, (element, radius_, coord) in enumerate(
            zip(elements, radii, coordinates), start=1
        ):
            if i in excluded_atoms:
                continue
            elif (not include_hs) and element == 1:
                continue
            else:
                atom = Atom(element, coord, radius_, i)
                atoms.append(atom)

        # Set variables for outside access and function access.
        self._atoms = atoms
        self._excluded_atoms = set(excluded_atoms)

        # Compute buried volume
        self._compute_buried_volume(center=center, radius=radius, density=density)

    def octant_analysis(self) -> "BuriedVolume":
        """Perform octant analysis of the buried volume."""
        # Set up limits depending on the sphere radius
        lim = self._sphere.radius
        octant_limits = {
            0: ((0, lim), (0, lim), (0, lim)),
            1: ((-lim, 0), (0, lim), (0, lim)),
            3: ((0, lim), (-lim, 0), (0, lim)),
            2: ((-lim, 0), (-lim, 0), (0, lim)),
            7: ((0, lim), (0, lim), (-lim, 0)),
            6: ((-lim, 0), (0, lim), (-lim, 0)),
            4: ((0, lim), (-lim, 0), (-lim, 0)),
            5: ((-lim, 0), (-lim, 0), (-lim, 0)),
        }

        # Calculated volume for each octant.
        octant_volume = self._sphere.volume / 8

        # Do octant analysis
        percent_buried_volume = {}
        buried_volume = {}
        free_volume = {}
        for name, limits in octant_limits.items():
            buried_points = self._buried_points[
                np.logical_and.reduce(
                    [
                        self._buried_points[:, 0] > limits[0][0],
                        self._buried_points[:, 0] < limits[0][1],
                        self._buried_points[:, 1] > limits[1][0],
                        self._buried_points[:, 1] < limits[1][1],
                        self._buried_points[:, 2] > limits[2][0],
                        self._buried_points[:, 2] < limits[2][1],
                    ]
                )
            ]
            free_points = self._free_points[
                np.logical_and.reduce(
                    [
                        self._free_points[:, 0] > limits[0][0],
                        self._free_points[:, 0] < limits[0][1],
                        self._free_points[:, 1] > limits[1][0],
                        self._free_points[:, 1] < limits[1][1],
                        self._free_points[:, 2] > limits[2][0],
                        self._free_points[:, 2] < limits[2][1],
                    ]
                )
            ]
            fraction_buried = len(buried_points) / (
                len(buried_points) + len(free_points)
            )
            percent_buried_volume[name] = fraction_buried * 100
            buried_volume[name] = fraction_buried * octant_volume
            free_volume[name] = (1 - fraction_buried) * octant_volume
        self.octants = {
            "percent_buried_volume": percent_buried_volume,
            "buried_volume": buried_volume,
            "free_volume": free_volume,
        }

        # Do quadrant analysis
        percent_buried_volume = {}
        buried_volume = {}
        free_volume = {}
        for name, octants in QUADRANT_OCTANT_MAP.items():
            percent_buried_volume[name] = (
                sum(
                    [
                        self.octants["percent_buried_volume"][octant]
                        for octant in octants
                    ]
                )
                / 2
            )
            buried_volume[name] = sum(
                [self.octants["buried_volume"][octant] for octant in octants]
            )
            free_volume[name] = sum(
                [self.octants["free_volume"][octant] for octant in octants]
            )
        self.quadrants = {
            "percent_buried_volume": percent_buried_volume,
            "buried_volume": buried_volume,
            "free_volume": free_volume,
        }

        self._octant_limits = octant_limits

        return self

    def _compute_buried_volume(
        self, center: ArrayLike1D, radius: float, density: float
    ) -> None:
        """Compute buried volume."""
        center: Array1DFloat = np.array(center)
        # Construct sphere at metal center
        sphere = Sphere(
            center, radius, method="projection", density=density, filled=True
        )

        # Prune sphere points which are within vdW radius of other atoms.
        tree = scipy.spatial.cKDTree(
            sphere.points, compact_nodes=False, balanced_tree=False
        )
        mask: Array1DBool = np.zeros(len(sphere.points), dtype=bool)
        for atom in self._atoms:
            if atom.radius + sphere.radius > np.linalg.norm(atom.coordinates):
                to_prune = tree.query_ball_point(atom.coordinates, atom.radius)
                mask[to_prune] = True
        buried_points = sphere.points[mask, :]
        free_points = sphere.points[np.invert(mask), :]

        # Calculate buried_volume
        self.fraction_buried_volume = len(buried_points) / len(sphere.points)
        self.buried_volume = sphere.volume * self.fraction_buried_volume
        self.free_volume = sphere.volume - self.buried_volume
        self._sphere = sphere
        self._buried_points = buried_points
        self._free_points = free_points

    def compute_distal_volume(
        self, method: str = "sasa", octants: bool = False, sasa_density: float = 0.01
    ) -> "BuriedVolume":
        """Computes the distal volume.

        Uses either SASA or Buried volume with large radius to calculate the molecular
        volume.

        Args:
            method: Method to get total volume: 'buried_volume' or 'sasa'
            octants: Whether to compute distal volume for quadrants and octants.
                Requires method='buried_volume'
            sasa_density: Density of points on SASA surface. Ignored unless
                method='sasa'

        Returns:
            self: Self

        Raises:
            ValueError: When method is not specified correctly.
        """
        loop_coordinates: list[Array1DFloat]
        # Use SASA to calculate total volume of the molecule
        if method == "sasa":
            # Calculate total volume
            elements: list[int] = []
            loop_coordinates = []
            radii: list[float] = []
            for atom in self._atoms:
                elements.append(atom.element)
                loop_coordinates.append(atom.coordinates)
                radii.append(atom.radius)
            coordinates: Array2DFloat = np.vstack(loop_coordinates)
            sasa = SASA(
                elements,
                coordinates,
                radii=radii,
                probe_radius=0.0,
                density=sasa_density,
            )
            self.molecular_volume = sasa.volume

            # Calculate distal volume
            self.distal_volume = self.molecular_volume - self.buried_volume
        elif method == "buried_volume":
            if octants is True and self.octants is None:
                raise ValueError("Needs octant analysis.")

            # Save the values for the old buried volume calculation
            temp_bv = copy.deepcopy(self)

            # Determine sphere radius to cover the whole molecule
            loop_coordinates = []
            radii = []
            for atom in self._atoms:
                loop_coordinates.append(atom.coordinates)
                radii.append(atom.radius)
            coordinates = np.vstack(loop_coordinates)
            distances = scipy.spatial.distance.cdist(
                self._sphere.center.reshape(1, -1), coordinates
            )
            new_radius = np.max(distances + radii) + 0.5

            # Compute the distal volume
            temp_bv._compute_buried_volume(
                center=self._sphere.center,
                radius=new_radius,
                density=self._sphere.density,
            )
            self.molecular_volume = temp_bv.buried_volume
            self.distal_volume = self.molecular_volume - self.buried_volume
            if octants is True:
                temp_bv.octant_analysis()

                # Octant analysis
                distal_volume = {}
                molecular_volume = {}
                for name in self.octants["buried_volume"].keys():
                    molecular_volume[name] = temp_bv.octants["buried_volume"][name]
                    distal_volume[name] = (
                        temp_bv.octants["buried_volume"][name]
                        - self.octants["buried_volume"][name]
                    )
                self.octants["distal_volume"] = distal_volume
                self.octants["molecular_volume"] = molecular_volume

                # Quadrant analyis
                distal_volume = {}
                molecular_volume = {}
                for name, octants_ in QUADRANT_OCTANT_MAP.items():
                    distal_volume[name] = sum(
                        [self.octants["distal_volume"][octant] for octant in octants_]
                    )
                    molecular_volume[name] = sum(
                        [
                            self.octants["molecular_volume"][octant]
                            for octant in octants_
                        ]
                    )
                self.quadrants["distal_volume"] = distal_volume
                self.quadrants["molecular_volume"] = molecular_volume
        else:
            raise ValueError(f"Method {method} is not valid.")

        return self

    @requires_dependency([Import(module="matplotlib.pyplot", alias="plt")], globals())
    def plot_steric_map(  # noqa: C901
        self,
        filename: str | None = None,
        levels: float = 150,
        grid: int = 100,
        all_positive: bool = True,
        cmap: str = "viridis",
    ) -> None:
        """Plots a steric map as in the original article.

        Args:
            filename: Name of file for saving the plot.
            levels: Number of levels in the contour plot
            grid: Number of points along each axis of plotting grid
            all_positive: Whether to plot only positive values
            cmap: Matplotlib colormap for contour plot

        Raises:
            ValueError: When z-axis atoms not present
        """
        if self._z_axis_atoms is None:
            raise ValueError("Must give z-axis atoms when instantiating BuriedVolume.")
        # Set up coordinates
        atoms = self._atoms
        center: Array1DFloat = np.array(self._sphere.center)
        all_coordinates: Array2DFloat = self._all_coordinates
        coordinates: Array2DFloat = np.array([atom.coordinates for atom in atoms])

        # Translate coordinates
        all_coordinates -= center
        coordinates -= center
        center -= center

        # Get vector to midpoint of z-axis atoms
        z_axis_coordinates = all_coordinates[np.array(self._z_axis_atoms) - 1]
        point = np.mean(z_axis_coordinates, axis=0)
        vector = point - center
        vector = vector / np.linalg.norm(vector)

        # Rotate coordinate system
        coordinates = rotate_coordinates(coordinates, vector, np.array([0, 0, -1]))

        # Make grid
        r = self._sphere.radius
        x_ = np.linspace(-r, r, grid)
        y_ = np.linspace(-r, r, grid)

        # Calculate z values
        z = []
        for line in np.dstack(np.meshgrid(x_, y_)).reshape(-1, 2):
            if np.linalg.norm(line) > r:
                z.append(np.nan)
                continue
            x = line[0]
            y = line[1]
            z_list = []
            for i, atom in enumerate(atoms):
                # Check if point is within reach of atom.
                x_s = coordinates[i, 0]
                y_s = coordinates[i, 1]
                z_s = coordinates[i, 2]
                test = atom.radius**2 - (x - x_s) ** 2 - (y - y_s) ** 2
                if test >= 0:
                    z_atom = math.sqrt(test) + z_s
                    z_list.append(z_atom)
            # Take point which is furthest along z axis
            if z_list:
                z_max = max(z_list)
                # Test if point is inside the sphere. Points with positive z
                # values are included by default anyway in accordance to
                # article
                if all_positive:
                    if z_max < 0:
                        if np.linalg.norm(np.array([x, y, z_max])) >= r:
                            z_max = np.nan
                else:
                    if np.linalg.norm(np.array([x, y, z_max])) >= r:
                        z_max = np.nan
            else:
                z_max = np.nan
            z.append(z_max)

        # Create interaction surface
        z: Array2DFloat = np.array(z).reshape(len(x_), len(y_))

        # Plot surface
        fig, ax = plt.subplots()
        cf = ax.contourf(x_, y_, z, levels, cmap=cmap)
        circle = plt.Circle((0, 0), r, fill=False)
        ax.add_patch(circle)
        plt.xlabel("x (Å)")
        plt.ylabel("y (Å)")
        cf.set_clim(-r, r)
        c_bar = fig.colorbar(cf)
        c_bar.set_label("z(Å)")
        ax.set_aspect("equal", "box")

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    def print_report(self) -> None:
        """Prints a report of the buried volume."""
        print("V_bur (%):", round(self.fraction_buried_volume * 100, 1))

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
        buried_color: str = "tomato",
        free_color: str = "steelblue",
        opacity: float = 0.5,
        size: float = 5,
    ) -> None:
        """Draw a the molecule with the buried and free points.

        Args:
            atom_scale: Scaling factor for atom size
            background_color: Background color for plot
            buried_color: Color of buried points
            free_color: Color of free points
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

        # Add buried points
        p.add_points(
            self._buried_points, color=buried_color, opacity=opacity, point_size=size
        )

        # Add free points
        p.add_points(
            self._free_points, color=free_color, opacity=opacity, point_size=size
        )

        if hasattr(self, "_octant_limits"):
            for name, limits_ in self._octant_limits.items():
                limits = tuple(itertools.chain(*limits_))
                box = pv.Box(limits)
                p.add_mesh(box, style="wireframe")
                x = np.array(limits)[:2][np.argmax(np.abs(limits[:2]))]
                y = np.array(limits)[2:4][np.argmax(np.abs(limits[2:4]))]
                z = np.array(limits)[4:][np.argmax(np.abs(limits[4:]))]
                p.add_point_labels(
                    np.array([x, y, z]), [OCTANT_SIGNS[name]], text_color="black"
                )

        self._plotter = p

    @property
    def percent_buried_volume(self) -> float:
        """Deprecated attribute. Use 'fraction_buried_volume' instead."""
        warnings.warn(
            "'percent_buried_volume' is deprecated. Use 'fraction_buried_volume'.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.fraction_buried_volume

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


def cli(file: str) -> Any:
    """CLI for buried volume.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(BuriedVolume, elements, coordinates)
