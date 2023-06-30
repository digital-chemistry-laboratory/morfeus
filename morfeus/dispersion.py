"""Dispersion code."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
import functools
from os import PathLike
import typing
from typing import Any

import numpy as np
import scipy.spatial

from morfeus.calculators import D3Calculator, D4Grimme
from morfeus.data import ANGSTROM_TO_BOHR, atomic_symbols, HARTREE_TO_KCAL, jmol_colors
from morfeus.geometry import Atom
from morfeus.io import CubeParser, D3Parser, D4Parser, read_geometry, VertexParser
from morfeus.sasa import SASA
from morfeus.typing import (
    Array1DFloat,
    Array1DInt,
    Array2DFloat,
    Array2DInt,
    ArrayLike1D,
    ArrayLike2D,
)
from morfeus.utils import convert_elements, get_radii, Import, requires_dependency

if typing.TYPE_CHECKING:
    from matplotlib.colors import hex2color
    import pymeshfix
    import pyvista as pv
    from pyvistaqt import BackgroundPlotter
    import vtk


class Dispersion:
    """Calculates and stores the results for the 🍺P_int dispersion descriptor.

    The descriptor is defined in 10.1002/anie.201905439. Morfeus can compute it based on
    a surface either from vdW radii, surface vertices or the electron density.
    Dispersion can be obtained with the D3 or D4 model.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        radii: VdW radii (Å)
        radii_type: Choice of vdW radii: 'alvarez', 'bondi', 'crc', 'rahm' and 'truhlar'
        point_surface: Use point surface from vdW radii
        compute_coefficients: Whether to compute D3 coefficients with internal code
        density: Area per point (Å²) on the vdW surface
        excluded_atoms: Atoms to exclude (1-indexed). Used for substituent P_ints
        included_atoms: Atoms to include. Used for functional group P_ints

    Attributes:
        area: Area of surface (Å²)
        atom_areas: Atom indices as keys and atom areas as values (Å²)
        atom_p_int: Atom indices as keys and P_int as values (kcal¹ᐟ² mol⁻¹⸍²))
        atom_p_max: Atom indices as keys and P_max as values (kcal¹ᐟ² mol⁻¹ᐟ²)
        atom_p_min: Atom indices as keys and P_min as values( kcal¹ᐟ² mol⁻¹ᐟ²)
        p_int: P_int value for molecule (kcal¹ᐟ² mol⁻¹ᐟ²)
        p_max: Highest P value (kcal¹ᐟ² mol⁻¹ᐟ²)
        p_min: Lowest P value (kcal¹ᐟ² mol⁻¹ᐟ²)
        p_values: All P values (kcal¹ᐟ² mol⁻¹ᐟ²)
        volume: Volume of surface (Å³)

    Raises:
        Exception: When both exluded_atoms and included_atom are given
    """

    area: float
    atom_areas: dict[int, float]
    atom_p_int: dict[int, float]
    atom_p_max: dict[int, float]
    atom_p_min: dict[int, float]
    p_int: float
    p_max: float
    p_min: float
    p_values: Array1DFloat
    volume: float
    _atoms: list[Atom]
    _c_n_coefficients: dict[int, Array1DFloat]
    _density: float
    _excluded_atoms: list[int]
    _point_areas: Array1DFloat
    _point_map: Array1DInt
    _points: Array2DFloat
    _radii: Array1DFloat
    _surface: "pv.PolyData"

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        radii: ArrayLike1D | None = None,
        radii_type: str = "rahm",
        point_surface: bool = True,
        compute_coefficients: bool = True,
        density: float = 0.1,
        excluded_atoms: Sequence[int] | None = None,
        included_atoms: Sequence[int] | None = None,
    ) -> None:
        # Check that only excluded or included atoms are given
        if excluded_atoms is not None and included_atoms is not None:
            raise Exception("Give either excluded or included atoms but not both.")

        # Converting elements to atomic numbers if the are symbols
        elements = convert_elements(elements, output="numbers")
        coordinates: Array2DFloat = np.array(coordinates)

        # Set excluded atoms
        all_atoms = set(range(1, len(elements) + 1))
        if included_atoms is not None:
            included_atoms_ = set(included_atoms)
            excluded_atoms = list(all_atoms - included_atoms_)
        elif excluded_atoms is None:
            excluded_atoms = []
        else:
            excluded_atoms = list(excluded_atoms)

        self._excluded_atoms = excluded_atoms

        # Set up
        self._surface = None
        self._point_areas = None
        self._density = density

        # Getting radii if they are not supplied
        if radii is None:
            radii = get_radii(elements, radii_type=radii_type)
        radii: Array1DFloat = np.array(radii)

        self._radii = radii

        # Get vdW surface if requested
        if point_surface:
            self._surface_from_sasa(elements, coordinates)
        else:
            # Get list of atoms as Atom objects
            atoms: list[Atom] = []
            for i, (element, coord, radius) in enumerate(
                zip(elements, coordinates, radii), start=1
            ):
                atom = Atom(element, coord, radius, i)
                atoms.append(atom)
            self._atoms = atoms

        # Calculate coefficients
        if compute_coefficients:
            self.compute_coefficients(model="id3")

        # Calculatte P_int values
        if point_surface and compute_coefficients:
            self.compute_p_int()

    def _surface_from_sasa(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
    ) -> None:
        """Get surface from SASA."""
        sasa = SASA(
            elements,
            coordinates,
            radii=self._radii,
            density=self._density,
            probe_radius=0,
        )
        self._atoms = sasa._atoms
        self.area = sum(
            [
                atom.area
                for atom in self._atoms
                if atom.index not in self._excluded_atoms
            ]
        )
        self.atom_areas = sasa.atom_areas
        self.volume = sum(
            [
                atom.volume
                for atom in self._atoms
                if atom.index not in self._excluded_atoms
            ]
        )

        # Get point areas and map from point to atom
        point_areas: list[Array1DFloat] = []
        point_map = []
        for atom in self._atoms:
            n_points = len(atom.accessible_points)
            if n_points > 0:
                point_area = atom.area / n_points
            else:
                point_area = 0.0
            atom.point_areas = np.repeat(point_area, n_points)
            point_areas.extend(atom.point_areas)
            point_map.extend([atom.index] * n_points)
        self._point_areas = np.array(point_areas)
        self._point_map = np.array(point_map)

    @requires_dependency([Import(module="pyvista", alias="pv")], globals())
    def surface_from_cube(
        self,
        file: str | PathLike,
        isodensity: float = 0.001,
        method: str = "flying_edges",
    ) -> "Dispersion":
        """Adds an isodensity surface from a Gaussian cube file.

        Args:
            file: Gaussian cube file
            isodensity: Isodensity value (electrons/bohr³)
            method: Method for contouring: 'contour' or 'flying_edges

        Returns:
            self: Self
        """
        # Parse the cubefile
        parser = CubeParser(file)

        # Generate grid and fill with values
        grid = pv.UniformGrid()
        grid.dimensions = np.array(parser.X.shape)
        grid.origin = (parser.min_x, parser.min_y, parser.min_z)
        grid.spacing = (parser.step_x, parser.step_y, parser.step_z)
        grid.point_data["values"] = parser.S.flatten(order="F")
        self.grid = grid

        # Contour and process the surface
        surface = self._contour_surface(grid, method=method, isodensity=isodensity)
        self._surface = surface
        self._process_surface()

        return self

    @requires_dependency(
        [Import("pymeshfix"), Import(module="pyvista", alias="pv")], globals()
    )
    def surface_from_multiwfn(
        self, file: str | PathLike, fix_mesh: bool = True
    ) -> "Dispersion":
        """Adds surface from Multiwfn vertex file with connectivity information.

        Args:
            file: Vertex.pdb file
            fix_mesh: Whether to fix holes in the mesh with pymeshfix (recommended)

        Returns:
            self: Self
        """
        # Read the vertices and faces from the Multiwfn output file
        parser = VertexParser(file)
        vertices: Array2DFloat = np.array(parser.vertices)
        faces: Array2DInt = np.array(parser.faces)
        faces: Array2DInt = np.insert(faces, 0, values=3, axis=1)

        # Construct surface and fix it with pymeshfix
        surface = pv.PolyData(vertices, faces)
        if fix_mesh:
            meshfix = pymeshfix.MeshFix(surface)
            meshfix.repair()
            surface = meshfix.mesh

        # Process surface
        self._surface = surface
        self._process_surface()

        return self

    def _process_surface(self) -> None:
        """Extracts face center points and assigns these to atoms based on proximity."""
        # Get the area and volume
        self.area = self._surface.area
        self.volume = self._surface.volume

        # Assign face centers to atoms according to Voronoi partitioning
        coordinates: Array2DFloat = np.array([atom.coordinates for atom in self._atoms])
        points: Array2DFloat = np.array(self._surface.cell_centers().points)
        kd_tree = scipy.spatial.cKDTree(coordinates)
        _, point_regions = kd_tree.query(points, k=1)
        point_regions = point_regions + 1

        # Compute faces areas
        area_data = self._surface.compute_cell_sizes()
        areas: Array1DFloat = np.array(area_data.cell_data["Area"])

        # Assign face centers and areas to atoms
        atom_areas = {}
        for atom in self._atoms:
            atom.accessible_points = points[point_regions == atom.index]
            point_areas = areas[point_regions == atom.index]
            atom.area = np.sum(point_areas)
            atom.point_areas = point_areas
            atom_areas[atom.index] = atom.area

        # Set up attributes
        self.atom_areas = atom_areas
        self._point_areas = areas
        self._point_map = point_regions

    @requires_dependency(
        [Import(module="pyvista", alias="pv"), Import("vtk")], globals()
    )
    @staticmethod
    def _contour_surface(
        grid: "pv.Grid", method: str = "flying_edges", isodensity: float = 0.001
    ) -> "pv.PolyData":
        """Counter surface from grid.

        Args:
            grid: Electron density as PyVista Grid object
            isodensity: Isodensity value (electrons/bohr³)
            method: Method for contouring: 'contour' or 'flying_edges

        Returns:
            surface: Surface as Pyvista PolyData object
        """
        # Select method for contouring
        if method == "flying_edges":
            contour_filter = vtk.vtkFlyingEdges3D()
        elif method == "contour":
            contour_filter = vtk.vtkContourFilter()

        # Run the contour filter
        isodensity = isodensity
        contour_filter.SetInputData(grid)
        contour_filter.SetValue(0, isodensity)
        contour_filter.Update()
        surface = contour_filter.GetOutput()
        surface = pv.wrap(surface)

        return surface

    def compute_p_int(  # noqa: C901
        self, points: ArrayLike2D | None = None
    ) -> "Dispersion":
        """Compute P_int values for surface or points.

        Args:
            points: Points to compute P values for

        Returns:
            self: Self
        """
        # Set up atoms and coefficients that are part of the calculation
        atom_indices: Array1DInt = np.array(
            [
                atom.index - 1
                for atom in self._atoms
                if atom.index not in self._excluded_atoms
            ]
        )
        coordinates: Array2DFloat = np.array([atom.coordinates for atom in self._atoms])
        coordinates = coordinates[atom_indices]
        c_n_coefficients = dict(self._c_n_coefficients)
        for key, value in c_n_coefficients.items():
            c_n_coefficients[key] = np.array(value)[atom_indices] * HARTREE_TO_KCAL

        # Take surface points if none are given
        if points is None:
            points = np.vstack(
                [
                    atom.accessible_points
                    for atom in self._atoms
                    if atom.index not in self._excluded_atoms
                    and atom.accessible_points.size > 0
                ]
            )
            atomic = True
        else:
            atomic = False
            points = np.array(points)

        # Calculate p_int for each point
        dist = scipy.spatial.distance.cdist(points, coordinates) * ANGSTROM_TO_BOHR
        p = np.sum(
            [
                np.sum(np.sqrt(coefficients / (dist**order)), axis=1)
                for order, coefficients in c_n_coefficients.items()
            ],
            axis=0,
        )

        self.p_values = p

        # Take out atomic p_ints if no points are given
        if atomic is True:
            atom_p_max = {}
            atom_p_min = {}
            atom_p_int = {}
            i_start = 0
            for atom in self._atoms:
                if atom.index not in self._excluded_atoms:
                    n_points = len(atom.accessible_points)
                    if n_points > 0:
                        i_stop = i_start + n_points
                        atom_ps = p[i_start:i_stop]
                        atom.p_values = atom_ps
                        atom_p_max[atom.index] = np.max(atom_ps)
                        atom_p_min[atom.index] = np.min(atom_ps)
                        atom_p_int[atom.index] = np.sum(
                            atom_ps * atom.point_areas / atom.area
                        )
                        i_start = i_stop
                    else:
                        atom_p_max[atom.index] = 0
                        atom_p_min[atom.index] = 0
                        atom_p_int[atom.index] = 0
                        atom.p_values = np.array([])
            self.atom_p_max = atom_p_max
            self.atom_p_min = atom_p_min
            self.atom_p_int = atom_p_int

        if self._point_areas is not None:
            point_areas = self._point_areas[np.isin(self._point_map, atom_indices + 1)]
            self.p_int = np.sum(p * point_areas / self.area)

        # Calculate p_min and p_max with slight modification to Robert's
        # definitions
        self.p_min = np.min(p)
        self.p_max = np.max(p)

        # Map p_values onto surface
        if self._surface is not None:
            mapped_p = np.zeros(len(p))
            for atom in self._atoms:
                if atom.index not in self._excluded_atoms:
                    mapped_p[self._point_map == atom.index] = atom.p_values
            self._surface.cell_data["values"] = mapped_p
            self._surface = self._surface.cell_data_to_point_data()

        # Store points for later use
        self._points = points

        return self

    def compute_coefficients(
        self, model: str = "id3", order: int = 8, charge: int = 0
    ) -> "Dispersion":
        """Compute dispersion coefficients.

        Can either use internal D3 model or D4 via Grimme's dftd4 program.

        Args:
            model: Calculation model: 'id3' (default) or 'gd4'
            order: Order of the Cᴬᴬ coefficients
            charge: Molecular charge for D4 model

        Returns:
            self: Self

        Raises:
            ValueError: When model not supported
        """
        # Set up atoms and coordinates
        elements = [atom.element for atom in self._atoms]
        coordinates: Array2DFloat = np.array([atom.coordinates for atom in self._atoms])

        calculators = {
            "id3": D3Calculator,
            "gd4": D4Grimme,
        }
        calc: D3Calculator | D4Grimme
        # Calculate  D3 values with internal model
        if model == "id3":
            calc = calculators[model](elements, coordinates, order=order)
        elif model == "gd4":
            calc = calculators[model](elements, coordinates, order=order, charge=charge)
        else:
            raise ValueError(f"model={model} not supported.")
        self._c_n_coefficients = calc.c_n_coefficients

        return self

    def load_coefficients(self, file: str | PathLike, model: str) -> "Dispersion":
        """Load the C₆ and C₈ coefficients.

        Output can be read from the dftd3 and dftd4 programs by giving a file in
        combination with the corresponding model.

        Args:
            file: Output file from the dftd3 or dftd4 programs
            model: Calculation model: 'd3' or 'd4'

        Returns:
            self: Self

        Raises:
            ValueError: When model not supported
        """
        parser: D3Parser | D4Parser
        if model == "d3":
            parser = D3Parser(file)
        elif model == "d4":
            parser = D4Parser(file)
        else:
            raise ValueError(f"model={model} not supported.")
        self._c_n_coefficients = {}
        self._c_n_coefficients[6] = parser.c6_coefficients
        self._c_n_coefficients[8] = parser.c8_coefficients

        return self

    def print_report(self, verbose: bool = False) -> None:
        """Print report of results.

        Args:
            verbose: Whether to print atom P_ints
        """
        print(f"Surface area (Å²): {self.area:.1f}")
        print(f"Surface volume (Å³): {self.volume:.1f}")
        print(f"P_int (kcal¹ᐟ² mol⁻¹ᐟ²): {self.p_int:.1f}")
        if verbose:
            print(
                f"{'Symbol':<10s}{'Index':<10s}{'P_int (kcal^(1/2) mol^(-1/2))':<30s}"
            )
            for atom, (i, p_int) in zip(self._atoms, self.atom_p_int.items()):
                symbol = atomic_symbols[atom.element]
                print(f"{symbol:<10s}{i:<10d}{p_int:<10.1f}")

    def save_vtk(self, filename: str) -> "Dispersion":
        """Save surface as .vtk file.

        Args:
            filename: Name of file. Use .vtk suffix.

        Returns:
            self: Self
        """
        self._surface.save(filename)

        return self

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
        opacity: float = 1,
        display_p_int: bool = True,
        molecule_opacity: float = 1,
        atom_scale: float = 1,
    ) -> None:
        """Draw surface with mapped P_int values.

        Args:
            opacity: Surface opacity
            display_p_int: Whether to display P_int mapped onto the surface
            molecule_opacity: Molecule opacity
            atom_scale: Scale factor for atom size
        """
        # Set up plotter
        p = BackgroundPlotter()

        # Draw molecule
        for atom in self._atoms:
            color = hex2color(jmol_colors[atom.element])
            radius = atom.radius * atom_scale
            sphere = pv.Sphere(center=list(atom.coordinates), radius=radius)
            p.add_mesh(
                sphere, color=color, opacity=molecule_opacity, name=str(atom.index)
            )

        cmap: str | None
        # Set up plotting of mapped surface
        if display_p_int is True:
            color = None
            cmap = "coolwarm"
        else:
            color = "tan"
            cmap = None

        # Draw surface
        if self._surface:
            p.add_mesh(self._surface, opacity=opacity, color=color, cmap=cmap)
        else:
            point_cloud = pv.PolyData(self._points)
            point_cloud["values"] = self.p_values
            p.add_mesh(
                point_cloud,
                opacity=opacity,
                color=color,
                cmap=cmap,
                render_points_as_spheres=True,
            )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


def cli(file: str) -> Any:
    """CLI for dispersion descriptor.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    elements, coordinates = read_geometry(file)
    return functools.partial(Dispersion, elements, coordinates)
