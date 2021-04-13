"""Sterimol code."""

import math
import typing
from typing import Iterable, List, Optional, Sequence, Set, Union

import numpy as np
import scipy.spatial

from morfeus.data import jmol_colors
from morfeus.geometry import Atom, kabsch_rotation_matrix, sphere_line_intersection
from morfeus.plotting import Arrow3D
from morfeus.sasa import SASA
from morfeus.typing import Array1D, ArrayLike1D, ArrayLike2D
from morfeus.utils import convert_elements, get_radii, Import, requires_dependency

if typing.TYPE_CHECKING:
    from matplotlib.colors import hex2color
    import pyvista as pv
    from pyvistaqt import BackgroundPlotter


class Sterimol:
    """Performs and stores results of Sterimol calculation.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        dummy_index: Index of dummy atom (1-indexed)
        attached_index: Index of attached atom of substituent (1-indexed)
        radii: List of radii (Å)
        radii_type: vdW radii type: 'alvarez', 'bondi', 'crc' or 'truhlar'
        n_rot_vectors: Number of rotational vectors for determining B₁ and B₅
        excluded_atoms: Atom indices to exclude from the calculation
        calculate: Whether to calculate the Sterimol parameters directly

    Attributes:
        B_1_value: Sterimol B₁ value (Å)
        B_1: Sterimol B₁ vector (Å)
        B_5_value: Sterimol B_5 value (Å)
        B_5: Sterimol B₅ vector (Å)
        bond_length: Bond length between atom 1 and atom 2 (Å)
        L_value_uncorrected: Sterimol L value minus 0.40 (Å)
        L_value: Sterimol L value (Å)
        L: Sterimol L vector (Å)
    """

    B_1_value: float
    B_1: Array1D
    B_5_value: float
    B_5: Array1D
    bond_length: float
    L_value_uncorrected: float
    L_value: float
    L: Array1D
    _atoms: List[Atom]
    _attached_atom: Atom
    _dummy_atom: Atom
    _excluded_atoms: Set[int]
    _n_rot_vectors: int
    _origin: Array1D
    _plotter: "BackgroundPlotter"
    _points: Array1D
    _sphere_radius: float

    def __init__(
        self,
        elements: Union[Iterable[int], Iterable[str]],
        coordinates: ArrayLike2D,
        dummy_index: int,
        attached_index: int,
        radii: Optional[ArrayLike1D] = None,
        radii_type: str = "crc",
        n_rot_vectors: int = 3600,
        excluded_atoms: Optional[Sequence[int]] = None,
        calculate: bool = True,
    ) -> None:
        # Convert elements to atomic numbers if the are symbols
        elements = convert_elements(elements, output="numbers")
        coordinates = np.array(coordinates)

        if excluded_atoms is None:
            excluded_atoms = []
        excluded_atoms = np.array(excluded_atoms)

        # Get radii if they are not supplied
        if radii is None:
            radii = get_radii(elements, radii_type=radii_type)
        radii = np.array(radii)
        all_radii = radii

        # Set up coordinate array
        all_coordinates = coordinates

        # Translate coordinates so origin is at atom 2
        origin = all_coordinates[attached_index - 1]
        all_coordinates -= origin

        # Get vector pointing from atom 2 to atom 1
        vector_2_to_1 = (
            all_coordinates[attached_index - 1] - all_coordinates[dummy_index - 1]
        )
        vector_2_to_1 = vector_2_to_1 / np.linalg.norm(vector_2_to_1)

        # Get rotation quaternion that overlays vector with x-axis
        x_axis = np.array([[1.0, 0.0, 0.0]])
        R = kabsch_rotation_matrix(vector_2_to_1.reshape(1, -1), x_axis, center=False)
        all_coordinates = (R @ all_coordinates.T).T
        self._rotation_matrix = R

        # Get list of atoms as Atom objects
        atoms = []
        for i, (element, radius, coord) in enumerate(
            zip(elements, all_radii, all_coordinates), start=1
        ):
            atom = Atom(element, coord, radius, i)
            atoms.append(atom)
            if i == dummy_index:
                dummy_atom = atom
            if i == attached_index:
                attached_atom = atom

        # Set up attributes
        self._atoms = atoms
        self._excluded_atoms = set(excluded_atoms)
        self._points = np.array([])
        self._origin = origin

        self._dummy_atom = dummy_atom
        self._attached_atom = attached_atom

        self.bond_length = np.linalg.norm(vector_2_to_1)

        self._n_rot_vectors = n_rot_vectors

        if calculate:
            self.calculate()

    def set_points(self, points: Sequence[Sequence[float]], shift: bool = True) -> None:
        """Set points for calculation of Sterimol.

        Args:
            points: Points (Å)
            shift: Whether to shift the points to have origin at dummy atom.
        """
        points = np.array(points)
        if shift is True:
            points -= self._origin

    def bury(
        self,
        sphere_radius: float = 3.5,
        method: str = "delete",
        radii_scale: float = 0.5,
    ) -> None:
        """Do a Buried Sterimol calculation.

        There are three available schemes based on deletion, truncation or slicing.

        Args:
            sphere_radius: Radius of sphere (Å)
            method: Method for burying: 'delete', 'slice' or 'truncate'
            radii_scale: Scale radii for metohd='delete' calculation

        Raises:
            ValueError: When method is not specified correctly
            Exception: When using method='truncate' and sphere is too small
        """
        if method == "delete":
            # Remove all atoms outside sphere (taking vdW radii into account)
            coordinates = np.vstack([atom.coordinates for atom in self._atoms])
            radii = np.array([atom.radius for atom in self._atoms])
            distances = scipy.spatial.distance.cdist(
                self._dummy_atom.coordinates.reshape(1, -1), coordinates
            ).reshape(-1)
            distances = distances - radii * radii_scale
            excluded_atoms = set(np.array(self._atoms)[distances >= sphere_radius])
            self._excluded_atoms.update(excluded_atoms)

            # Calculate Sterimol parameters
            self.calculate()
        elif method == "truncate":
            # Calculate Sterimol parameters
            self.calculate()

            # Calculate intersection between vectors and sphere.
            atom_1_coordinates = self._dummy_atom.coordinates
            atom_2_coordinates = self._attached_atom.coordinates

            new_vectors = []
            for vector, ref_coordinates in [
                (self.L, atom_1_coordinates),
                (self.B_1, atom_2_coordinates),
                (self.B_5, atom_2_coordinates),
            ]:
                # Get intersection point
                intersection_points = sphere_line_intersection(
                    vector, atom_1_coordinates, sphere_radius
                )
                if len(intersection_points) < 1:
                    raise Exception("Sphere so small that vectors don't intersect")

                # Get vector pointing in the right direction
                trial_vectors = [
                    point - ref_coordinates for point in intersection_points
                ]
                norm_vector = vector / np.linalg.norm(vector)
                dot_products = [
                    np.dot(norm_vector, trial_vector / np.linalg.norm(trial_vector))
                    for trial_vector in trial_vectors
                ]
                new_vector = trial_vectors[int(np.argmax(dot_products))]
                new_vectors.append(new_vector)

            # Replace vectors if new ones are shorter than old ones
            if np.linalg.norm(self.L) > np.linalg.norm(new_vectors[0]):
                L = new_vectors[0]
                L_value = np.linalg.norm(L)
            if np.linalg.norm(self.B_1) > np.linalg.norm(new_vectors[1]):
                B_1 = new_vectors[1]
                B_1_value = np.linalg.norm(B_1)
            if np.linalg.norm(self.B_5) > np.linalg.norm(new_vectors[2]):
                B_5 = new_vectors[2]
                B_5_value = np.linalg.norm(B_5)

            # Set attributes
            self.L = L
            self.L_value = L_value + 0.40
            self.L_value_uncorrected = L_value

            self.B_1 = B_1
            self.B_1_value = B_1_value

            self.B_5 = B_5
            self.B_5_value = B_5_value
        elif method == "slice":
            # Remove points outside of sphere
            distances = scipy.spatial.distance.cdist(
                self._dummy_atom.coordinates.reshape(1, -1), self._points
            ).reshape(-1)
            self._points = self._points[distances <= sphere_radius]

            # Calculate Sterimol parameters
            self.calculate()
        else:
            raise ValueError(f"Method: {method} is not supported.")

        # Set attributes
        self._sphere_radius = sphere_radius

    def surface_from_radii(self, density: float = 0.01) -> None:
        """Create surface points from vdW surface.

        Args:
            density: Area per point on surface (Å²)
        """
        # Calculate vdW surface for all active atoms
        elements = []
        coordinates = []
        radii = []
        for atom in self._atoms:
            if atom.index not in self._excluded_atoms and atom is not self._dummy_atom:
                elements.append(atom.element)
                coordinates.append(atom.coordinates)
                radii.append(atom.radius)
        elements = np.array(elements)
        coordinates = np.vstack(coordinates)
        radii = radii
        sasa = SASA(elements, coordinates, radii=radii, density=density, probe_radius=0)

        # Take out points of vdW surface
        points = np.vstack(
            [
                atom.accessible_points
                for atom in sasa._atoms
                if atom.index not in self._excluded_atoms
                and atom.accessible_points.size > 0
            ]
        )
        self._points = points

    def calculate(self) -> None:
        """Calculate Sterimol parameters."""
        # Use coordinates and radii if points are not given
        if not len(self._points) > 0:
            coordinates = []
            radii = []
            for atom in self._atoms:
                if (
                    atom is not self._dummy_atom
                    and atom.index not in self._excluded_atoms
                ):
                    coordinates.append(atom.coordinates)
                    radii.append(atom.radius)
            coordinates = np.vstack(coordinates)
            radii = np.vstack(radii).reshape(-1)

        # Project coordinates onto vector between atoms 1 and 2
        vector = self._attached_atom.coordinates - self._dummy_atom.coordinates
        bond_length = np.linalg.norm(vector)
        unit_vector = vector / np.linalg.norm(vector)

        if not len(self._points) > 0:
            c_values = np.dot(unit_vector.reshape(1, -1), coordinates.T)
            projected = c_values + radii
        else:
            projected = np.dot(unit_vector.reshape(1, -1), self._points.T)

        # Get L as largest projection along the vector
        L_value = np.max(projected) + bond_length
        L = unit_vector * L_value
        L = L.reshape(-1)

        # Get rotation vectors in yz plane
        r = 1
        theta = np.linspace(0, 2 * math.pi, self._n_rot_vectors)
        x = np.zeros(len(theta))
        y = r * np.cos(theta)
        z = r * np.sin(theta)
        rot_vectors = np.column_stack((x, y, z))

        # Project coordinates onto rotation vectors
        if not len(self._points) > 0:
            c_values = np.dot(rot_vectors, coordinates.T)
            projected = c_values + radii
        else:
            projected = np.dot(rot_vectors, self._points.T)
        max_c_values = np.max(projected, axis=1)

        # Determine B1 and B5 from the smallest and largest scalar projections
        B_1_value = np.min(max_c_values)
        B_1 = rot_vectors[np.argmin(max_c_values)] * B_1_value

        B_5_value = np.max(max_c_values)
        B_5 = rot_vectors[np.argmax(max_c_values)] * B_5_value

        # Set up attributes
        self.L = L
        self.L_value = L_value + 0.40
        self.L_value_uncorrected = L_value

        self.B_1 = B_1
        self.B_1_value = B_1_value

        self.B_5 = B_5
        self.B_5_value = B_5_value

    def print_report(self, verbose: bool = False) -> None:
        """Prints the values of the Sterimol parameters.

        Args:
            verbose: Whether to print uncorrected L_value and bond length
        """
        if verbose:
            print(
                f"{'L':10s}{'B_1':10s}{'B_5':10s}" f"{'L_uncorr':10s}{'d(a1-a2)':10s}"
            )
            print(
                f"{self.L_value:<10.2f}{self.B_1_value:<10.2f}"
                f"{self.B_5_value:<10.2f}{self.L_value_uncorrected:<10.2f}"
                f"{self.bond_length:<10.2f}"
            )
        else:
            print(f"{'L':10s}{'B_1':10s}{'B_5':10s}")
            print(
                f"{self.L_value:<10.2f}{self.B_1_value:<10.2f}"
                f"{self.B_5_value:<10.2f}"
            )

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
        atom_scale: float = 0.5,
        background_color: str = "white",
        arrow_color: str = "steelblue",
    ) -> None:
        """Draw a 3D representation of the molecule with the Sterimol vectors.

        Args:
            atom_scale: Scaling factor for atom size
            background_color: Background color for plot
            arrow_color: Arrow color
        """
        # Set up plotter
        p = BackgroundPlotter()
        p.set_background(background_color)

        # Draw molecule
        for atom in self._atoms:
            color = hex2color(jmol_colors[atom.element])
            radius = atom.radius * atom_scale
            sphere = pv.Sphere(center=list(atom.coordinates), radius=radius)
            if atom.index in self._excluded_atoms:
                opacity = 0.25
            else:
                opacity = 1
            p.add_mesh(sphere, color=color, opacity=opacity, name=str(atom.index))

        # Draw sphere for Buried Sterimol
        if hasattr(self, "_sphere_radius"):
            sphere = pv.Sphere(
                center=self._dummy_atom.coordinates, radius=self._sphere_radius
            )
            p.add_mesh(sphere, opacity=0.25)

        if len(self._points) > 0:
            p.add_points(self._points, color="gray")

        # Get arrow starting points
        start_L = self._dummy_atom.coordinates
        start_B = self._attached_atom.coordinates

        # Add L arrow with label
        length = np.linalg.norm(self.L)
        direction = self.L / length
        stop_L = start_L + length * direction
        L_arrow = Arrow3D(start=start_L, direction=direction, length=length)
        p.add_mesh(L_arrow, color=arrow_color)

        # Add B_1 arrow
        length = np.linalg.norm(self.B_1)
        direction = self.B_1 / length
        stop_B_1 = start_B + length * direction
        B_1_arrow = Arrow3D(start=start_B, direction=direction, length=length)
        p.add_mesh(B_1_arrow, color=arrow_color)

        # Add B_5 arrow
        length = np.linalg.norm(self.B_5)
        direction = self.B_5 / length
        stop_B_5 = start_B + length * direction
        B_5_arrow = Arrow3D(start=start_B, direction=direction, length=length)
        p.add_mesh(B_5_arrow, color=arrow_color)

        # Add labels
        points = np.vstack([stop_L, stop_B_1, stop_B_5])
        labels = ["L", "B1", "B5"]
        p.add_point_labels(
            points,
            labels,
            text_color="black",
            font_size=30,
            bold=False,
            show_points=False,
            point_size=1,
        )

        self._plotter = p

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"
