"""Classes for performing calculations of steric descriptors of molecules."""

import copy
import itertools
import math

import numpy as np
import scipy.spatial
from scipy.io import FortranFile
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist, euclidean

from morfeus.calculators import D3Calculator, D3Grimme, D4Grimme
from morfeus.data import (AFU, AMU, ANGSTROM, ANGSTROM_TO_BOHR, BOHR,
                          BOHR_TO_ANGSTROM, DYNE, HARTREE, HARTREE_TO_EV,
                          HARTREE_TO_KCAL, C, atomic_masses, atomic_symbols,
                          jmol_colors)
from morfeus.geometry import (Angle, Atom, Bond, Cone, Dihedral,
                              InternalCoordinates, Sphere,
                              kabsch_rotation_matrix, rotate_coordinates,
                              sphere_line_intersection)
from morfeus.helpers import (Import, check_distances, convert_elements,
                             get_connectivity_matrix, get_radii,
                             requires_dependency)
from morfeus.io import CubeParser, D3Parser, D4Parser, VertexParser
from morfeus.plotting import Arrow3D, Cone3D


class Sterimol:
    """Performs and stores results of Sterimol calculation.

    Args:
        dummy_index (int): Index of dummy atom, starting at 1)
        attached_index (int): Index of attached atom of substituent, starting
                          at 1)
        coordinates (list): Coordinates (Å)
        elements (list): Elements as atomic symbols or numbers
        excluded_atoms (list): Atoms to exclude from the calculation
        n_rot_vectors (int): Number of rotational vectors for determining B_1
                             and B_5
        radii (list): List of radii (Å, optional)
        radii_type (str): Type of radii to use: 'bondi' or 'crc'
        sphere_radius (float): Radius of sphere for Buried Sterimol (Å)

    Attributes:
        B_1 (ndarray): Sterimol B_1 vector (Å)
        B_1_value (float): Sterimol B_1 value (Å)
        B_5 (ndarray): Sterimol B_5 vector (Å)
        B_5_value (float): Sterimol B_5 value (Å)
        bond_length (float): Bond length between atom 1 and atom 2 (Å)
        L (ndarray): Sterimol L vector (Å)
        L_value (float): Sterimol L value (Å)
        L_value_uncorrected (float): Sterimol L value minus 0.40 Å
    """
    def __init__(self, elements, coordinates, dummy_index, attached_index, radii=[],
                 radii_type="crc", n_rot_vectors=3600, excluded_atoms=[],
                 calculate=True):
        # Convert elements to atomic numbers if the are symbols
        elements = convert_elements(elements)

        # Get radii if they are not supplied
        if len(radii) < 1:
            radii = get_radii(elements, radii_type=radii_type)
        all_radii = np.array(radii)

        # Set up coordinate array
        all_coordinates = np.array(coordinates)

        # Translate coordinates so origin is at atom 2
        all_coordinates -= all_coordinates[attached_index - 1]

        # Get vector pointing from atom 2 to atom 1
        vector_2_to_1 = all_coordinates[attached_index - 1] \
                        - all_coordinates[dummy_index - 1]
        vector_2_to_1 = vector_2_to_1 / np.linalg.norm(vector_2_to_1)

        # Get rotation quaternion that overlays vector with x-axis
        x_axis = np.array([[1.0, 0.0, 0.0]])
        R = kabsch_rotation_matrix(vector_2_to_1.reshape(1, -1), x_axis, center=False)
        all_coordinates = (R @ all_coordinates.T).T
        self._rotation_matrix = R

        # Get list of atoms as Atom objects
        atoms = []
        for i, (element, radius, coord) in enumerate(
                zip(elements, all_radii, all_coordinates), start=1):
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

        self._dummy_atom = dummy_atom
        self._attached_atom = attached_atom

        self.L = None
        self.L_value = None
        self.L_value_uncorrected = None
        self.bond_length = np.linalg.norm(vector_2_to_1)

        self.B_1 = None
        self.B_1_value = None

        self.B_5 = None
        self.B_5_value = None

        self._n_rot_vectors = n_rot_vectors
        self._sphere_radius = None

        if calculate:
            self.calculate()

    def bury(self, sphere_radius=3.5, method="delete", radii_scale=0.5):
        """Do a Buried Sterimol calculation according to the three available
        schemes based on deletion, truncation or slicing.

        Args:
            method (str): Method for burying: 'delete', 'slice' or 'truncate'
            radii_scale (float): Scale radii for 'delete'-type calculation.
            sphere_radius (float): Radius of sphere (Å)
        """
        if method == "delete":
            # Remove all atoms outside sphere (taking vdW radii into account)
            coordinates = np.vstack([atom.coordinates for atom in self._atoms])
            radii = np.array([atom.radius for atom in self._atoms])
            distances = cdist(self._dummy_atom.coordinates.reshape(1, -1),
                              coordinates).reshape(-1)
            distances = distances - radii * radii_scale
            excluded_atoms = set(
                np.array(self._atoms)[distances >= sphere_radius])
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
            for vector, ref_coordinates in [(self.L, atom_1_coordinates), (self.B_1, atom_2_coordinates), (self.B_5, atom_2_coordinates)]:
                # Get intersection point
                intersection_points = sphere_line_intersection(vector, atom_1_coordinates, sphere_radius)
                if len(intersection_points) < 1:
                    raise Exception("Sphere so small that vectors don't intersect")

                # Get vector pointing in the right direction
                trial_vectors = [point - ref_coordinates for point in intersection_points]
                norm_vector = vector / np.linalg.norm(vector)
                dot_products = [np.dot(norm_vector, trial_vector / np.linalg.norm(trial_vector)) for trial_vector in trial_vectors]
                new_vector = trial_vectors[np.argmax(dot_products)]
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
            distances = cdist(self._dummy_atom.coordinates.reshape(1, -1), self._points).reshape(-1)
            self._points = self._points[distances <= sphere_radius]

            # Calculate Sterimol parameters
            self.calculate()
        else:
            raise Exception("Please specify valid method for burying.")

        # Set attributes
        self._sphere_radius = sphere_radius

    def surface_from_radii(self, density=0.01):
        # Calculate vdW surface for all active atoms
        elements = []
        coordinates = []
        radii = []
        for atom in self._atoms:
            if atom not in self._excluded_atoms and atom is not self._dummy_atom:
                elements.append(atom.element)
                coordinates.append(atom.coordinates)
                radii.append(atom.radius)
        elements = np.array(elements)
        coordinates = np.vstack(coordinates)
        radii = radii
        sasa = SASA(elements, coordinates, radii=radii, density=density,
                    probe_radius=0)

        # Take out points of vdW surface
        points = np.vstack([atom.accessible_points for atom in sasa._atoms
                            if atom.index not in self._excluded_atoms and
                            atom.accessible_points.size > 0])
        self._points = points

    def calculate(self):
        # Use coordinates and radii if points are not given
        if not len(self._points) > 0:
            coordinates = []
            radii = []
            for atom in self._atoms:
                if atom != self._dummy_atom and atom not in self._excluded_atoms:
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

    def print_report(self, verbose=False):
        """Prints the values of the Sterimol parameters.

        Args:
            verbose (bool): Toggles printing of uncorrected L_value and bond
                            length between atom 1 and atom 2.
        """
        if verbose:
            print(f"{'L':10s}{'B_1':10s}{'B_5':10s}"
                  f"{'L_uncorr':10s}{'d(a1-a2)':10s}")
            print(f"{self.L_value:<10.2f}{self.B_1_value:<10.2f}"
                  f"{self.B_5_value:<10.2f}{self.L_value_uncorrected:<10.2f}"
                  f"{self.bond_length:<10.2f}")
        else:
            print(f"{'L':10s}{'B_1':10s}{'B_5':10s}")
            print(f"{self.L_value:<10.2f}{self.B_1_value:<10.2f}"
                  f"{self.B_5_value:<10.2f}")

    @requires_dependency([
        Import(module="matplotlib.colors", item="hex2color"),
        Import(module="pyvista", alias="pv"),
        Import(module="pyvistaqt", item="BackgroundPlotter")
    ], globals())
    def draw_3D(self,
                atom_scale=0.5,
                background_color="white",
                arrow_color="steelblue"):
        """Draw a 3D representation of the molecule with the Sterimol vectors.
        
        Args:
            atom_scale (float): Scaling factor for atom size
            background_color (str): Background color for plot
            arrow_color (str): Arrow color
        """
        # Set up plotter
        p = BackgroundPlotter()
        p.set_background(background_color)

        # Draw molecule
        for atom in self._atoms:
            color = hex2color(jmol_colors[atom.element])
            radius = atom.radius * atom_scale
            sphere = pv.Sphere(center=list(atom.coordinates),
                               radius=radius)
            if atom in self._excluded_atoms:
                opacity = 0.25
            else:
                opacity = 1
            p.add_mesh(sphere,
                       color=color,
                       opacity=opacity,
                       name=str(atom.index))

        # Draw sphere for Buried Sterimol
        if self._sphere_radius:
            sphere = pv.Sphere(
                center=self._dummy_atom.coordinates,
                radius=self._sphere_radius
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
        p.add_point_labels(points, labels, text_color="black", font_size=30,
                           bold=False, show_points=False, point_size=1)

        self._plotter = p

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


class BuriedVolume:
    """Performs and stores the results of a buried volume calculation as
    described in Organometallics 2016, 35, 2286.

    Args:
        metal_index (int): Atom index of metal (starting from 1)
        coordinates (list): Coordinates (Å)
        density (float): Volume per point (Å**3) in the sphere
        elements (list): Elements as atomic symbols or numbers
        excluded_atoms (list): Indices of atoms to exclude from the calculation
                             (starting from 1)
        include_hs (bool): Whether to include H atoms in the calculation
        radii (list): vdW radii (Å)
        radii_scale (float): Scaling factor for radii. 1.17 from original paper.
        radii_type (str): Type of radii to use: 'bondi' (default) or 'crc'
        radius (float): Radius of sphere (Å). 3.5 from orginal paper.

    Parameters:
        buried_volume (float): Buried volume
    """

    # Quadrant and octant signs taken from
    # https://en.wikipedia.org/wiki/Octant_(solid_geometry)
    quadrant_signs = {
        1: "+,+",
        2: "-,+",
        3: "-,-",
        4: "+,-",
    }

    octant_signs = {
        0: "+,+,+",
        1: "-,+,+",
        3: "+,-,+",
        2: "-,-,+",
        7: "+,+,-",
        6: "-,+,-",
        4: "+,-,-",
        5: "-,-,-",
    }

    # Conventional names for quadrants
    quadrant_names = {
        1: "NE",
        2: "NW",
        3: "SW",
        4: "SE",
    }

    # Maps octants to quadrants
    octant_quadrant_map = {
            1: (0, 7),
            2: (1, 6),
            3: (2, 5),
            4: (3, 4)
        }

    def __init__(self, elements, coordinates, metal_index, excluded_atoms=[],
                 radii=[], include_hs=False, radius=3.5, radii_type="bondi",
                 radii_scale=1.17, density=0.001, z_axis_atoms=[],
                 xz_plane_atoms=[]):
        # Get center and and reortient coordinate system
        coordinates = np.array(coordinates)
        center = coordinates[metal_index - 1]
        coordinates -= center

        if metal_index not in excluded_atoms:
            excluded_atoms.append(metal_index)

        if len(z_axis_atoms) > 0 and len(xz_plane_atoms) > 0:
            z_axis_coordinates = coordinates[np.array(z_axis_atoms) - 1]
            z_point = np.mean(z_axis_coordinates, axis=0)

            xz_plane_coordinates = coordinates[np.array(xz_plane_atoms) - 1]
            xz_point = np.mean(xz_plane_coordinates, axis=0)

            v_1 = z_point - center
            v_2 = xz_point - center
            v_3 = np.cross(v_2, v_1)
            real = np.vstack([v_1, v_3])
            real /= np.linalg.norm(real, axis=1).reshape(-1, 1)
            ref_1 = np.array([0., 0., -1.])
            ref_2 = np.array([0., 1., 0.])
            ref = np.vstack([ref_1, ref_2])
            R = kabsch_rotation_matrix(real, ref, center=False)
            coordinates = (R @ coordinates.T).T
        elif len(z_axis_atoms) > 0:
            z_axis_coordinates = coordinates[np.array(z_axis_atoms) - 1]
            z_point = np.mean(z_axis_coordinates, axis=0)
            v_1 = z_point - center
            v_1 = v_1 / np.linalg.norm(v_1)
            coordinates = rotate_coordinates(coordinates, v_1,
                                             np.array([0, 0, -1]))

        # Save density and coordinates for steric map plotting.
        self._density = density
        self._all_coordinates = coordinates

        # Converting element ids to atomic numbers if the are symbols
        elements = convert_elements(elements)

        # Getting radii if they are not supplied
        if len(radii) == 0:
            radii = get_radii(elements, radii_type=radii_type,
                              scale=radii_scale)

        # Get list of atoms as Atom objects
        atoms= []
        for i, (element, radius_, coord) in enumerate(zip(elements, radii,
                                               coordinates), start=1):
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
        self._compute_buried_volume(center=center, radius=radius,
                                    density=density)

        # Set up attributes
        self.molecular_volume = None
        self.distal_volume = None
        self.octants = None
        self.quadrants = None
        self._octant_limits = None

        self.octants = {}
        self.quadrants = {}

    def octant_analysis(self):
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
            buried_points = self._buried_points[np.logical_and.reduce([
                self._buried_points[:, 0] > limits[0][0],
                self._buried_points[:, 0] < limits[0][1],
                self._buried_points[:, 1] > limits[1][0],
                self._buried_points[:, 1] < limits[1][1],
                self._buried_points[:, 2] > limits[2][0],
                self._buried_points[:, 2] < limits[2][1],
            ])]
            free_points = self._free_points[np.logical_and.reduce([
                self._free_points[:, 0] > limits[0][0],
                self._free_points[:, 0] < limits[0][1],
                self._free_points[:, 1] > limits[1][0],
                self._free_points[:, 1] < limits[1][1],
                self._free_points[:, 2] > limits[2][0],
                self._free_points[:, 2] < limits[2][1],
            ])]
            fraction_buried = len(buried_points) / (len(buried_points) + len(free_points))
            percent_buried_volume[name] = fraction_buried * 100
            buried_volume[name] = fraction_buried * octant_volume
            free_volume[name] = (1 - fraction_buried) * octant_volume
        self.octants["percent_buried_volume"] = percent_buried_volume
        self.octants["buried_volume"] = buried_volume
        self.octants["free_volume"] = free_volume

        # Do quadrant analysis
        percent_buried_volume = {}
        buried_volume = {}
        free_volume = {}
        for name, octants in self.octant_quadrant_map.items():
            percent_buried_volume[name] = sum([self.octants["percent_buried_volume"][octant] for octant in octants]) / 2
            buried_volume[name] = sum([self.octants["buried_volume"][octant] for octant in octants])
            free_volume[name] = sum([self.octants["free_volume"][octant] for octant in octants])
        self.quadrants["percent_buried_volume"] = percent_buried_volume
        self.quadrants["buried_volume"] = buried_volume
        self.quadrants["free_volume"] = free_volume

        self._octant_limits = octant_limits

    def _compute_buried_volume(self, center, radius, density):
        # Construct sphere at metal center
        sphere = Sphere(center, radius, method="projection", density=density,
                        filled=True)

        # Prune sphere points which are within vdW radius of other atoms.
        tree = scipy.spatial.cKDTree(sphere.points, compact_nodes=False,
                                     balanced_tree=False)
        mask = np.zeros(len(sphere.points), dtype=bool)
        for atom in self._atoms:
            if atom.radius + sphere.radius > np.linalg.norm(atom.coordinates):
                to_prune = tree.query_ball_point(atom.coordinates,
                                                 atom.radius)
                mask[to_prune] = True
        buried_points = sphere.points[mask,:]
        free_points = sphere.points[np.invert(mask),:]

        # Calculate buried_volume
        self.percent_buried_volume = len(buried_points) / len(sphere.points)
        self.buried_volume = sphere.volume * self.percent_buried_volume
        self.free_volume = sphere.volume - self.buried_volume
        self._sphere = sphere
        self._buried_points = buried_points
        self._free_points = free_points

    def compute_distal_volume(self, method="sasa", octants=False,
                              sasa_density=0.01):
        """Computes the distal volume. Uses either SASA or Buried volume with
        large radius to calculate the molecular volume.

        Args:
            method (str): Method to get total volume: 'sasa' or 'buried_volume'
            octants (bool): Whether to compute distal volume for quadrants and
                            octants. Requires method='buried_volume'
            sasa_density (float): Density of points on SASA surface. Ignored
                unless method='sasa'
        """
        # Use SASA to calculate total volume of the molecule
        if method == "sasa":
            # Calculate total volume
            elements = []
            coordinates = []
            radii = []
            for atom in self._atoms:
                elements.append(atom.element)
                coordinates.append(atom.coordinates)
                radii.append(atom.radius)
            coordinates = np.vstack(coordinates)
            sasa = SASA(elements, coordinates, radii=radii, probe_radius=0.0,
                        density=sasa_density)
            self.molecular_volume = sasa.volume

            # Calculate distal volume
            self.distal_volume = self.molecular_volume - self.buried_volume
        elif method == "buried_volume":
            # Save the values for the old buried volume calculation
            temp_bv = copy.deepcopy(self)

            # Determine sphere radius to cover the whole molecule
            coordinates = []
            radii = []
            for atom in self._atoms:
                coordinates.append(atom.coordinates)
                radii.append(atom.radius)
            coordinates = np.vstack(coordinates)
            distances = cdist(self._sphere.center.reshape(1, -1), coordinates)
            new_radius = np.max(distances + radii) + 0.5

            # Compute the distal volume
            temp_bv._compute_buried_volume(center=self._sphere.center,
                                        radius=new_radius,
                                        density=self._sphere.density)
            self.molecular_volume = temp_bv.buried_volume
            self.distal_volume = self.molecular_volume - self.buried_volume
            if octants:
                temp_bv.octant_analysis()
                # Octant analysis
                distal_volume = {}
                molecular_volume = {}
                for name in self.octants["buried_volume"].keys():
                    molecular_volume[name] = temp_bv.octants["buried_volume"][name]
                    distal_volume[name] = temp_bv.octants["buried_volume"][name] - self.octants["buried_volume"][name]
                self.octants["distal_volume"] = distal_volume
                self.octants["molecular_volume"] = molecular_volume

                # Quadrant analyis
                distal_volume = {}
                molecular_volume = {}
                for name, octants in self.octant_quadrant_map.items():
                    distal_volume[name] = sum([self.octants["distal_volume"][octant] for octant in octants])
                    molecular_volume[name] = sum([self.octants["molecular_volume"][octant] for octant in octants])
                self.quadrants["distal_volume"] = distal_volume
                self.quadrants["molecular_volume"] = molecular_volume
        else:
            raise Exception("Provide valid method.")

    @requires_dependency([Import(module="matplotlib.pyplot", alias="plt")],
                         globals())
    def plot_steric_map(self, z_axis_atoms, filename=None, levels=150, grid=100,
                        all_positive=True, cmap="viridis"):
        """Plots a steric map as in the original article.

        Args:
            all_positive (bool): Plot all positive values
            cmap (str): Colormap for contour plot
            filename (str): Name of file for saving the plot.
            grid (int): Point along each axis of plotting grid 
            levels (int): Number of levels in the contour plot
            z_axis_atoms (list): Indices of atoms for determining the
                                 orientation of the z axis (starting at 1)
        """
        # Set up coordinates
        atoms = self._atoms
        center = np.array(self._sphere.center)
        all_coordinates = self._all_coordinates
        coordinates = np.array([atom.coordinates for atom in atoms])

        # Translate coordinates
        all_coordinates -= center
        coordinates -= center
        center -= center

        # Get vector to midpoint of z-axis atoms
        z_axis_coordinates = all_coordinates[np.array(z_axis_atoms) - 1]
        point = np.mean(z_axis_coordinates, axis=0)
        vector = point - center
        vector = vector / np.linalg.norm(vector)

        #Rotate coordinate system
        coordinates = rotate_coordinates(coordinates, vector,
                                         np.array([0, 0, -1]))

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
                test = atom.radius**2 - (x - x_s)**2 - (y - y_s)**2
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
        z = np.array(z).reshape(len(x_), len(y_))

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
        ax.set_aspect('equal', 'box')

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    def print_report(self):
        """Prints a report of the buried volume for use in shell scripts"""
        print("V_bur (%):", round(self.buried_volume * 100, 1))

    @requires_dependency([
        Import(module="matplotlib.colors", item="hex2color"),
        Import(module="pyvista", alias="pv"),
        Import(module="pyvistaqt", item="BackgroundPlotter")
    ], globals())
    def draw_3D(self,
                atom_scale=1,
                background_color="white",
                buried_color="tomato",
                free_color="steelblue",
                opacity=0.05,
                size=1):
        """Draw a 3D representation of the molecule with the buried and free
        points.

        Args:
            atom_scale (float): Scaling factor for atom size
            background_color (str): Background color for plot
            buried_color (str): Color of buried points
            free_color (str): Color of free points
            opacity (float): Point opacity
            size (float): Point size
        """
        # Set up plotter
        p = BackgroundPlotter()
        p.set_background(background_color)

        # Draw molecule
        for atom in self._atoms:
            color = hex2color(jmol_colors[atom.element])
            radius = atom.radius * atom_scale
            sphere = pv.Sphere(center=list(atom.coordinates),
                               radius=radius)
            p.add_mesh(sphere, color=color, opacity=1, name=str(atom.index))

        # Add buried points
        p.add_points(self._buried_points, color=buried_color, opacity=opacity,
                     point_size=size)

        # Add free points
        p.add_points(self._free_points, color=free_color, opacity=opacity,
                     point_size=size)

        if self._octant_limits:
            for name, limits in self._octant_limits.items():
                limits = tuple(itertools.chain(*limits))
                box = pv.Box(limits)
                p.add_mesh(box, style="wireframe")
                x = np.array(limits)[:2][np.argmax(np.abs(limits[:2]))]
                y = np.array(limits)[2:4][np.argmax(np.abs(limits[2:4]))]
                z = np.array(limits)[4:][np.argmax(np.abs(limits[4:]))]
                p.add_point_labels(np.array([x, y, z]), [self.octant_signs[name]], text_color="black")

        self._plotter = p

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


class SASA:
    """Performs and stores results of solvent accessible surface area 
    calculations.

    Args:
        coordinates (list): Coordinates (Å)
        density (float): Area per point (Å**2) on the vdW surface.
        elements (list): Elements as atomic symbols or numbers.
        probe_radius (float): Radius of probe atom (Å)
        radii (list): VdW radii (Å)
        radii_type (str): Choice of vdW radii: 'bondi' or 'crc' (default)

    Attributes:
        area (float): Area of the solvent accessible surface.
        atom_areas (dict): Atom areas (starting from 1)
        atom_volumes (dict): Atom volumes (starting from 1)
        volume (float): Volume of the solvent accessible surface.
    """
    def __init__(self, elements, coordinates, radii=[], radii_type="crc",
                 probe_radius=1.4, density=0.01):
        # Converting elements to atomic numbers if the are symbols
        elements = convert_elements(elements)

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(elements, radii_type=radii_type)

        # Increment the radii with the probe radius
        radii = np.array(radii)
        radii = radii + probe_radius

        # Construct list of atoms
        atoms = []
        for i, (coordinate, radius, element) in \
                enumerate(zip(coordinates, radii, elements), start=1):
            atom = Atom(element, coordinate, radius, i)
            atoms.append(atom)

        # Determine occluded and accessible points of each atom based on
        # distances to all other atoms (brute force)
        for atom in atoms:
            # Construct sphere for atom
            sphere = Sphere(atom.coordinates, atom.radius, density=density)

            # Select atoms that are at a distance less than the sum of radii
            #!TODO can be vectorized
            test_atoms = []
            for test_atom in atoms:
                if test_atom is not atom:
                    distance = euclidean(atom.coordinates, test_atom.coordinates)
                    radii_sum = atom.radius + test_atom.radius
                    if distance < radii_sum:
                        test_atoms.append(test_atom)

            # Select coordinates and radii for other atoms
            test_coordinates = [test_atom.coordinates for
                                test_atom in test_atoms]
            test_radii = [test_atom.radius for test_atom in test_atoms]
            test_radii = np.array(test_radii).reshape(-1, 1)

            # Get distances to other atoms and subtract radii
            if test_coordinates:
                distances = cdist(test_coordinates, sphere.points)
                distances -= test_radii
                # Take smallest distance and perform check
                min_distances = np.min(distances, axis=0)

                atom.occluded_points = sphere.points[min_distances < 0]
                atom.accessible_points = sphere.points[min_distances >= 0]
            else:
                atom.accessible_points = sphere.points
                atom.occluded_points = np.empty(0)

        # Calculate atom areas and volumes
        for atom in atoms:
            # Get number of points of eache type
            n_accesible = len(atom.accessible_points)
            n_occluded = len(atom.occluded_points)
            n_points = len(atom.accessible_points) + len(atom.occluded_points)

            # Calculate part occluded and accessible
            ratio_occluded = n_occluded / n_points
            ratio_accessible = 1 - ratio_occluded

            # Calculate area
            area = 4 * np.pi * atom.radius ** 2 * ratio_accessible
            atom.area = area

            # Center accessible points around origin
            centered_points = np.array(atom.accessible_points) \
                              - atom.coordinates

            # Add accessible points
            accessible_summed = np.sum(centered_points, axis=0)

            # Calculate volume
            volume = (4 * np.pi / 3 / n_points) * (atom.radius *
                      np.dot(atom.coordinates, accessible_summed)
                      + atom.radius ** 3 * n_accesible)
            atom.volume = volume

        # Set up attributes
        self._probe_radius = probe_radius
        self.atom_areas = {atom.index: atom.area for atom in atoms}
        self.atom_volumes = {atom.index: atom.volume for atom in atoms}
        self.area = sum([atom.area for atom in atoms])
        self.volume = sum([atom.volume for atom in atoms])
        self._atoms = atoms
        self._density = density

    def print_report(self, verbose=False):
        """Print report of results

        Args:
            verbose (bool): Print atom areas
        """
        print(f"Probe radius (Å): {self._probe_radius}")
        print(f"Solvent accessible surface area (Å^2): {self.area:.1f}")
        print("Volume inside solvent accessible surface (Å^3): "
              f"{self.volume:.1f}")
        if verbose:
            print(f"{'Symbol':<10s}{'Index':<10s}{'Area (Å^2)':<10s}")
            for atom, (i, area) in zip(self._atoms, self.atom_areas.items()):
                symbol = atomic_symbols[atom.element]
                print(f"{symbol:<10s}{i:<10d}{area:<10.1f}")

    @requires_dependency([
        Import(module="matplotlib.colors", item="hex2color"),
        Import(module="pyvista", alias="pv"),
        Import(module="pyvistaqt", item="BackgroundPlotter")
    ], globals())
    def draw_3D(self, atom_scale=1, background_color="white",
                point_color="steelblue", opacity=0.25, size=1):
        """Draw a 3D representation of the molecule with the solvent accessible
        surface area

        Args:
            atom_scale (float): Scaling factor for atom size
            background_color (str): Background color for plot
            point_color (str): Color of surface points
            opacity (float): Point opacity
            size (float): Point size
        """
        # Set up plotter
        p = BackgroundPlotter()
        p.set_background(background_color)

        # Draw molecule
        for atom in self._atoms:
            color = hex2color(jmol_colors[atom.element])
            radius = atom.radius * atom_scale - self._probe_radius
            sphere = pv.Sphere(center=list(atom.coordinates),
                               radius=radius)
            p.add_mesh(sphere, color=color, opacity=1, name=str(atom.index))

        # Draw surface points
        surface_points = np.vstack([atom.accessible_points
                                    for atom in self._atoms])
        p.add_points(surface_points, color=point_color, opacity=opacity,
                     point_size=size)

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


class ConeAngle:
    """Calculates and stores the results of exact cone angle calculation as
    described in J. Comput. Chem. 2013, 34, 1189.

    Args:
        atom_1 (int): Index of central atom (starting from 1)
        coordinates (list): Coordinates (Å)
        elements (list): Elements as atomic symbols or numbers
        radii (list): vdW radii (Å)
        radii_type (str): Type of vdW radii: 'bondi' or 'crc' (default)

    Attributes:
        cone_angle (float): Exact cone angle (degrees)
        tangent_atoms (list): Atoms tangent to cone
    """
    def __init__(self, elements, coordinates, atom_1, radii=[],
                 radii_type="crc"):
        # Convert elements to atomic numbers if the are symbols
        elements = convert_elements(elements)

        # Get radii if they are not supplied
        if not radii:
            radii = get_radii(elements, radii_type=radii_type)

        # Check so that no atom is within vdW distance of atom 1
        within = check_distances(elements, coordinates, atom_1, radii=radii)
        if within:
            atom_string = ' '.join([str(i) for i in within])
            raise Exception("Atoms within vdW radius of central atom:",
                             atom_string)

        # Set up coordinate array and translate coordinates
        coordinates = np.array(coordinates)
        coordinates -= coordinates[atom_1 - 1]

        # Get list of atoms as Atom objects
        atoms = []
        for i, (element, coord, radius) in \
                enumerate(zip(elements, coordinates, radii), start=1):
            if i != atom_1:
                atom = Atom(element, coord, radius, i)
                atom.get_cone()
                atoms.append(atom)
        self._atoms = atoms

        # Search for cone over single atoms
        cone = self._search_one_cones()

        # Prune out atoms that lie in the shadow of another atom's cone
        if not cone:
            loop_atoms = list(atoms)
            remove_atoms = set()
            for cone_atom in loop_atoms:
                for test_atom in loop_atoms:
                    if cone_atom != test_atom:
                        if cone_atom.cone.is_inside(test_atom):
                            remove_atoms.add(test_atom)
            for i in remove_atoms:
                loop_atoms.remove(i)
            self._loop_atoms = loop_atoms

        # Search for cone over pairs of atoms
        if not cone:
            cone = self._search_two_cones()

        # Search for cones over triples of atoms
        if not cone:
            cone = self._search_three_cones()

        # Set attributes
        if cone:
            self._cone = cone
            self.cone_angle = math.degrees(cone.angle * 2)
            self.tangent_atoms = [atom.index for atom in cone.atoms]
        else:
            raise Exception("Cone could not be found.")

    def print_report(self):
        """Prints report of results"""
        tangent_atoms = [atom for atom in self._atoms
                         if atom.index in self.tangent_atoms]
        tangent_labels = [f'{atomic_symbols[atom.element]}{atom.index}' \
                        for atom in tangent_atoms]
        tangent_string = ' '.join(tangent_labels)
        print(f"Cone angle: {self.cone_angle:.1f}")
        print(f"No. tangent atoms: {len(tangent_atoms)}")
        print(f"Tangent to: {tangent_string}")

    def _get_upper_bound(self):
        """Calculates upper bound for apex angle

        Returns:
            upper_bound (float): Upper bound to apex angle in radians
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

    def _search_one_cones(self):
        """Searches over cones tangent to one atom

        Returns:
            max_1_cone (obj): Largest cone tangent to one atom
        """
        # Get the largest cone
        atoms = self._atoms
        alphas = np.array([atom.cone.angle for atom in atoms])
        max_1_cone = atoms[np.argmax(alphas)].cone
        self._max_1_cone = max_1_cone

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

    def _search_two_cones(self):
        """Search over cones tangent to two atoms.

        Returns:
            max_2_cone (obj): Largest cone tangent to two atoms
        """
        # Create two-atom cones
        loop_atoms = self._loop_atoms
        cones = []
        for atom_i, atom_j in itertools.combinations(loop_atoms, r=2):
            cone = self._get_two_atom_cone(atom_i, atom_j)
            cones.append(cone)

        # Select largest two-atom cone
        angles = np.array([cone.angle for cone in cones])
        max_2_cone = cones[np.argmax(angles)]
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

    def _search_three_cones(self):
        """Search over cones tangent to three atoms

        Returns:
            min_3_cone (obj): Smallest cone tangent to three atoms
        """
        # Create three-atom cones
        loop_atoms = self._loop_atoms
        cones = []
        for atom_i, atom_j, atom_k in itertools.combinations(loop_atoms, r=3):
            three_cones = self._get_three_atom_cones(atom_i, atom_j, atom_k)
            cones.extend(three_cones)

        # Get upper and lower bound to apex angle
        upper_bound = self._get_upper_bound()
        lower_bound = self._max_2_cone.angle

        # Remove cones from consideration which are outside the bounds
        remove_cones = []
        for cone in cones:
            if cone.angle - lower_bound < -1e-5 or upper_bound - cone.angle < -1e-5:
                remove_cones.append(cone)

        for cone in reversed(remove_cones):
            cones.remove(cone)

        # Keep only cones that encompass all atoms
        keep_cones = []
        for cone in cones:
            in_atoms = []
            for atom in loop_atoms:
                in_atoms.append(cone.is_inside(atom))
            if all(in_atoms):
                keep_cones.append(cone)

        # Take the smallest cone that encompasses all atoms
        cone_angles = np.array([cone.angle for cone in keep_cones])
        min_3_cone = keep_cones[np.argmin(cone_angles)]

        return min_3_cone

    @staticmethod
    def _get_two_atom_cone(atom_i, atom_j):
        """Creates a cone tangent to two atoms

        Args:
            atom_i (obj): First tangent atom 
            atom_j (obj): Second tangent atom

        Returns:
            cones (obj): Cone tangent to the two atoms
        """
        # Get the cone angle
        cone_i = atom_i.cone
        cone_j = atom_j.cone
        beta_i = cone_i.angle
        beta_j = cone_j.angle
        beta_ij = math.acos(np.dot(atom_i.cone.normal, atom_j.cone.normal))
        alpha_ij = (beta_ij + beta_i + beta_j) / 2

        # Get the cone normal
        a_ij = (1 / math.sin(beta_ij)) * \
            math.sin(0.5 * (beta_ij + beta_i - beta_j))
        b_ij = (1 / math.sin(beta_ij)) * \
            math.sin(0.5 * (beta_ij - beta_i + beta_j))
        c_ij = 0

        n = a_ij * cone_i.normal + b_ij * cone_j.normal + c_ij
        n = n / np.linalg.norm(n)

        # Create cone
        angle = alpha_ij
        normal = n
        cone = Cone(angle, [atom_i, atom_j], normal)

        return cone

    @staticmethod
    def _get_three_atom_cones(atom_i, atom_j, atom_k):
        """Creates cones tangent to three atoms

        Args:
            atom_i (obj): First tangent atom 
            atom_j (obj): Second tangent atom 
            atom_k (obj): Third tangent atom 

        Returns:
            cones (list): Cones tangent to the three atoms
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
        N = np.array([np.cross(m_j, m_k), np.cross(m_k, m_i),
                      np.cross(m_i, m_j)]).T
        P = N.T @ N
        gamma = np.dot(m_i, np.cross(m_j, m_k))

        # Set up coefficients of quadratic equation
        A = u @ P @ u
        B = v.T @ P @ v
        C = u.T @ P @ v
        D = gamma**2

        # Solve quadratic equation
        p2 = (A - B)**2 + 4 * C**2
        p1 = 2 * (A - B) * (A + B - 2 * D)
        p0 = (A + B - 2 * D)**2 - 4 * C**2

        roots = np.roots([p2, p1, p0])
        roots = np.real_if_close(roots, tol=1e10)
        roots[np.isclose(roots, 1, rtol=1e-9, atol=0.0)] = 1
        roots[np.isclose(roots, -1, rtol=1e-9, atol=0.0)] = -1

        cos_roots = [math.acos(roots[0]), 2 * np.pi - math.acos(roots[0]),
                     math.acos(roots[1]), 2 * np.pi - math.acos(roots[1])]

        # Test roots and keep only those that are physical
        angles = []
        D_tests = []
        for root in cos_roots:
            alpha = root / 2
            test = A * math.cos(alpha)**2 + B * math.sin(alpha)**2 \
                   + 2 * C * math.sin(alpha) * math.cos(alpha)
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
            a_ij = (math.cos(alpha - beta_i) - math.cos(alpha - beta_j)
                    * math.cos(beta_ij)) / math.sin(beta_ij)**2
            b_ij = (math.cos(alpha - beta_j) - math.cos(alpha - beta_i)
                    * math.cos(beta_ij)) / math.sin(beta_ij)**2
            c_ij_squared = 1 - a_ij**2 - b_ij**2 \
                           - 2 * a_ij * b_ij * math.cos(beta_ij)
            # Set c_ij_squared to 0 if negative due to numerical precision.
            if c_ij_squared < 0:
                c_ij_squared = 0
            c_ij = math.sqrt(c_ij_squared)
            p = N @ (u * math.cos(alpha) + v * math.sin(alpha)).reshape(-1)
            sign = np.sign(gamma) * np.sign(np.dot(p, np.cross(m_i, m_j)))
            if np.sign(c_ij) != sign:
                c_ij = -c_ij
            n = a_ij * m_i + b_ij * m_j + c_ij * 1 \
                / math.sin(beta_ij) * np.cross(m_i, m_j)

            # Create cone
            cone = Cone(alpha, [atom_i, atom_j, atom_k], n)
            cones.append(cone)

        return cones

    @requires_dependency([
        Import(module="matplotlib.colors", item="hex2color"),
        Import(module="pyvista", alias="pv"),
        Import(module="pyvistaqt", item="BackgroundPlotter")
    ], globals())
    def draw_3D(self, atom_scale=1, background_color="white",
                cone_color="steelblue", cone_opacity=0.75):
        """Draw a 3D representation of the molecule with the cone.

        Args:
            atom_scale (float): Scaling factor for atom size
            background_color (str): Background color for plot
            cone_color (str): Cone color
            cone_opacity (float): Cone opacity
        """
        # Set up plotter
        p = BackgroundPlotter()
        p.set_background(background_color)

        # Draw molecule
        for atom in self._atoms:
            color = hex2color(jmol_colors[atom.element])
            radius = atom.radius * atom_scale
            sphere = pv.Sphere(center=list(atom.coordinates),
                               radius=radius)
            p.add_mesh(sphere, color=color, opacity=1, name=str(atom.index))

        # Determine direction and extension of cone
        cone_angle = math.degrees(self._cone.angle)
        coordinates = np.array([atom.coordinates for atom in self._atoms])
        radii = np.array([atom.radius for atom in self._atoms])
        if cone_angle > 180:
            normal = - self._cone.normal
        else:
            normal = self._cone.normal
        projected = np.dot(normal, coordinates.T) + np.array(radii)

        max_extension = np.max(projected)
        if cone_angle > 180:
            max_extension += 1

        # Make the cone
        cone = Cone3D(center=[0, 0, 0] + (max_extension * normal) / 2,
                      direction=-normal,
                      angle=cone_angle,
                      height=max_extension,
                      capping=False,
                      resolution=100)
        p.add_mesh(cone, opacity=cone_opacity, color=cone_color)

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


class Dispersion:
    """Calculates and stores the results for the 🍺P_int dispersion descriptor.

    The descriptor is defined in Angew. Chemie Int. Ed. 2019.
    DOI: 10.1002/anie.201905439. morfeus can compute it based on a surface
    either from vdW radii, surface vertices or electron density. Dispersion can
    be obtained with the D3 or D4 model.

    Args:
        compute_coefficients (bool): Whether to compute D3 coefficients with
            internal code.
        coordinates (list): Coordinates (Å)
        density (float): Area per point (Å**2) on the vdW surface.
        elements (list): Elements as atomic symbols or numbers
        excluded_atoms (list): Atoms to exclude from the calculation. Used only
            for calculation of substituent P_ints.
        included_atoms (list): Atoms to include in the calculation. Used only
            for calculation of functional group P_ints.
        point_surface (bool): Use point surface from vdW radii.
        radii (list): VdW radii (Å)
        radii_type (str): Choice of vdW radii: 'bondi', 'crc' or 'rahm'
            (default)

    Parameters:
        area (float): Area of surface (Å^2)
        atom_areas (dict): Atom areas (Å^2, starting from 1)
        atom_p_int (dict): P_int value for atoms (kcal^(1/2) mol^(-1/2)
            starting from 1)
        atom_p_max (dict): P_max value for atoms (kcal^(1/2) mol^(-1/2)
            starting from 1)
        atom_p_min (dict): P_min value for atoms (kcal^(1/2) mol^(-1/2)
            starting from 1)
        p_int (float): P_int value for molecule (kcal^(1/2) mol^(-1/2)
        p_max (float): Highest P values (kcal^(1/2) mol^(-1/2)
        p_min (float): Lowest P values (kcal^(1/2) mol^(-1/2)
        p_values (list): All P values (kcal^(1/2) mol^(-1/2)
        volume (float): Volume of surface (Å^3)
    """
    def __init__(self, elements, coordinates, radii=[], radii_type="rahm",
                 point_surface=True, compute_coefficients=True, density=0.1,
                 excluded_atoms=[], included_atoms=[]):
        # Check that only excluded or included atoms are given
        if len(excluded_atoms) > 0 and len(included_atoms) > 0:
            raise Exception("Give either excluded or included atoms but not both.")

        # Set excluded atoms
        all_atoms = set(range(1, len(elements) + 1))
        if len(included_atoms) > 0:
            included_atoms = set(included_atoms)
            excluded_atoms = list(all_atoms - included_atoms)

        self._excluded_atoms = excluded_atoms

        # Set up
        self._surface = None
        self._density = None

        # Converting elements to atomic numbers if the are symbols
        elements = convert_elements(elements)

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(elements, radii_type=radii_type)

        # Get vdW surface if requested
        if point_surface:
            sasa = SASA(elements, coordinates, radii=radii, density=density,
                        probe_radius=0)
            self._atoms = sasa._atoms
            self.area = sum([atom.area for atom in self._atoms
                             if atom.index not in excluded_atoms])
            self.atom_areas = sasa.atom_areas
            self.volume = sum([atom.volume for atom in self._atoms
                               if atom.index not in excluded_atoms])

            # Get point areas and map from point to atom
            point_areas = []
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
            self._density = density
        else:
            # Get list of atoms as Atom objects
            atoms = []
            for i, (element, coord, radius) in \
                        enumerate(zip(elements, coordinates, radii), start=1):
                atom = Atom(element, coord, radius, i)
                atoms.append(atom)
            self._atoms = atoms

        # Calculate coefficients
        if compute_coefficients:
            self.compute_coefficients(model='id3')

        # Calculatte P_int values
        if point_surface and compute_coefficients:
            self.compute_p_int()

    @requires_dependency([Import(module="pyvista", alias="pv")], globals())
    def surface_from_cube(self, filename, isodensity=0.001,
                              method="flying_edges"):
        """Adds an isodensity surface from a Gaussian cube file.
        
        Args:
            filename (str): Name of Gaussian cube file
            isodensity (float): Isodensity value (electrons/bohr^3)
            method (str): Method for contouring: 'contour' or 'flying_edges
                          (default)
        """
        # Parse the cubefile
        parser = CubeParser(filename)

        # Generate grid and fill with values
        grid = pv.UniformGrid()
        grid.dimensions = np.array(parser.X.shape)
        grid.origin = (parser.min_x, parser.min_y, parser.min_z)
        grid.spacing = (parser.step_x, parser.step_y, parser.step_z)
        grid.point_arrays['values'] = parser.S.flatten(order='F')
        self.grid = grid

        # Contour and process the surface
        surface = self._contour_surface(grid, method=method,
                                        isodensity=isodensity)
        self._surface = surface
        self._process_surface()

    @requires_dependency(
        [Import("pymeshfix"),
         Import(module="pyvista", alias="pv")], globals())
    def surface_from_multiwfn(self, filename, fix_mesh=True):
        """Adds surface from Multiwfn vertex file with connectivity information.

        Args:
            filename (str): Name of vertex file
            fix_mesh (bool): Whether to fix holes in the mesh with pymeshfix.
        """
        # Read the vertices and faces from the Multiwfn output file
        parser = VertexParser(filename)
        vertices = np.array(parser.vertices)
        faces = np.array(parser.faces)
        faces = np.insert(faces, 0, values=3, axis=1)

        # Construct surface and fix it with pymeshfix
        surface = pv.PolyData(vertices, faces, show_edges=True)
        if fix_mesh:
            meshfix = pymeshfix.MeshFix(surface)
            meshfix.repair()
            surface = meshfix.mesh

        # Process surface
        self._surface = surface
        self._process_surface()

    def _process_surface(self):
        """Extracts face center points and assigns these to atoms based on
        proximity
        """
        # Get the area and volume
        self.area = self._surface.area
        self.volume = self._surface.volume

        # Assign face centers to atoms according to Voronoi partitioning
        coordinates = np.array([atom.coordinates for atom in self._atoms])
        points = np.array(self._surface.cell_centers().points)
        kd_tree = cKDTree(coordinates)
        _, point_regions = kd_tree.query(points, k=1)
        point_regions = point_regions + 1

        # Compute faces areas
        area_data = self._surface.compute_cell_sizes()
        areas = np.array(area_data.cell_arrays["Area"])

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

    @requires_dependency([Import(module="pyvista", alias="pv"),
                          Import("vtk")], globals())
    @staticmethod
    def _contour_surface(grid, method="flying_edges", isodensity=0.001):
        """
        Args:
            grid (obj): Electron density as PyVista grid object
            isodensity (float): Isodensity value (electrons/bohr^3)
            method (str): Method for contouring: 'contour' or 'flying_edges
                          (default)
        
        Returns:
            surface (obj): Surface as Pyvista PolyData object
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

    def compute_p_int(self, points=[]):
        """Compute P_int values for surface or points.

        Args:
            points (list): Points to compute P values for
        """
        # Set up array of points
        points = np.array(points)

        # Set up atoms and coefficients that are part of the calculation
        atom_indices = np.array([atom.index - 1 for atom in self._atoms
                        if atom.index not in self._excluded_atoms])
        coordinates = np.array([atom.coordinates for atom in self._atoms])
        coordinates = coordinates[atom_indices]
        c_n_coefficients = dict(self._c_n_coefficients)
        for key, value in c_n_coefficients.items():
            c_n_coefficients[key] = np.array(value) * HARTREE_TO_KCAL

        # Take surface points if none are given
        if points.size == 0:
            points = np.vstack([atom.accessible_points for atom in self._atoms
                                if atom.index not in self._excluded_atoms and
                                atom.accessible_points.size > 0])
            atomic = True

        # Calculate p_int for each point
        dist = cdist(points, coordinates) * ANGSTROM_TO_BOHR
        p = sum([np.sum(np.sqrt(coefficients / (dist ** order)), axis=1) for order, coefficients in c_n_coefficients.items()])
        self.p_values = p

        # Take out atomic p_ints if no points are given
        if atomic:
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
                        atom_p_int[atom.index] = np.sum(atom_ps *
                            atom.point_areas / atom.area)
                        i_start = i_stop
                    else:
                        atom_p_max[atom.index] = 0
                        atom_p_min[atom.index] = 0
                        atom_p_int[atom.index] = 0
                        atom.p_values = np.array([])
            self.atom_p_max = atom_p_max
            self.atom_p_min = atom_p_min
            self.atom_p_int = atom_p_int

        point_areas = self._point_areas[np.isin(self._point_map,
                                        atom_indices + 1)]
        self.p_int = np.sum(p * point_areas / self.area)

        # Calculate p_min and p_max with slight modification to Robert's
        # definitions
        self.p_min = np.min(p)
        self.p_max = np.max(p)

        # Map p_values onto surface
        if self._surface:
            mapped_p = np.zeros(len(p))
            for atom in self._atoms:
                if atom.index not in self._excluded_atoms:
                    mapped_p[self._point_map == atom.index] = atom.p_values
            self._surface.cell_arrays['p_int'] = mapped_p

        # Store points for later use
        self._points = points

    def compute_coefficients(self, model='id3', order=8, charge=0):
        """Compute dispersion coefficients. Can either use internal D3 model
        or D4 or D3-like model available through Grimme's dftd4 program.

        Args:
            model (str): Calculation model: 'id3'. 'gd3' or 'gd4'
        """
        # Set up atoms and coordinates
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]

        # Calculate  D3 values with internal model
        if model =="id3":
            calc = D3Calculator(elements, coordinates, order=order)
            self._c_n_coefficients = calc.c_n_coefficients
        # Calculate the D3-like values with dftd4
        if model =="gd3":
            calc = D3Grimme(elements, coordinates, order=order)
            self._c_n_coefficients = calc.c_n_coefficients
        # Calculate the D4 values with dftd4
        if model =="gd4":
            calc = D4Grimme(elements, coordinates, order=order, charge=charge)
            self._c_n_coefficients = calc.c_n_coefficients

    def load_coefficients(self, filename, model):
        """Load the C6 and C8 coefficients.

        Output can be read from the dftd3 and dftd4 programs by giving a
        filename in combination with the corresponding model.

        Args:
            filename (str): Output file from the dftd3 or dftd4 programs
            model (str): Calculation model: 'd3' or 'd4'.
        """
        if filename and model == "d3":
            # Read the data
            parser = D3Parser(filename)
            self._c_n_coefficients = {}
            self._c_n_coefficients[6] = parser.c6_coefficients
            self._c_n_coefficients[8] = parser.c8_coefficients
        elif filename and model == "d4":
            # Read the data
            parser = D4Parser(filename)
            self._c_n_coefficients = {}
            self._c_n_coefficients[6] = parser.c6_coefficients
            self._c_n_coefficients[8] = parser.c8_coefficients

    def print_report(self, verbose=False):
        """Print report of results

        Args:
            verbose (bool): Print atom P_ints
        """
        print(f"Surface area (Å^2): {self.area:.1f}")
        print(f"Surface volume (Å^3): {self.volume:.1f}")
        print(f"P_int (kcal^(1/2) mol^(-1/2): {self.p_int:.1f}")
        if verbose:
            print(f"{'Symbol':<10s}{'Index':<10s}{'P_int (kcal^(1/2) mol^(-1/2))':<30s}")
            for atom, (i, p_int) in zip(self._atoms, self.atom_p_int.items()):
                symbol = atomic_symbols[atom.element]
                print(f"{symbol:<10s}{i:<10d}{p_int:<10.1f}")

    def save_vtk(self, filename):
        """Save surface as .vtk file

        Args:
            filename (str): Name of file. Use .vtk suffix.
        """
        self._surface.save(filename)

    @requires_dependency([
        Import(module="matplotlib.colors", item="hex2color"),
        Import(module="pyvista", alias="pv"),
        Import(module="pyvistaqt", item="BackgroundPlotter")
    ], globals())
    def draw_3D(self, opacity=1, display_p_int=True, molecule_opacity=1,
                atom_scale=1):
        """Draw surface with mapped P_int values.
        
        Args:
            atom_scale (float): Scale factor for atom size
            display_p_int (bool): Display P_int mapped onto the surface or not.
            molecule_opacity (float): Molecule opacity (0-1)
            opacity (float): Surface opacity (0-1)
        """
        # Set up plotter
        p = BackgroundPlotter()

        # Draw molecule
        for atom in self._atoms:
            color = hex2color(jmol_colors[atom.element])
            radius = atom.radius * atom_scale
            sphere = pv.Sphere(center=list(atom.coordinates),
                               radius=radius)
            p.add_mesh(sphere, color=color, opacity=molecule_opacity, name=str(atom.index))

        # Set up plotting of mapped surface
        if display_p_int == True:
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
            p.add_mesh(point_cloud, opacity=opacity, color=color, cmap=cmap, render_points_as_spheres=True)

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"


class LocalForce:
    """Calculates and stores the results from local force constant 
    calculations.
    
    The method is described by Cremer in Int. J. Quantum Chem. 1998, 67, 1.
    Alternatively, the compliance matrix method can be used according to 
    J. Chem. Phys. 2010, 132, 184101.

    Args:
        coordinates (list): Coordinates (Å)
        elements (list): Elements as atomic symbols or numbers.

    Attributes:
        internal_coordinates (list): Internal coordinates
        local_force_constants (ndarray): Local force constants (mDyne/Å)
        local_frequencies (ndarray): Local mode frequencies (cm^(-1))
        n_imag (int): Number of normal modes with imaginary frequencies
    """
    def __init__(self, elements=None, coordinates=None):
        # Set up attributes
        self.n_imag = None
        self.local_force_constants = np.array([])
        self.local_frequencies = np.array([])
        self._D = np.array([])
        self._B = np.array([])
        self._force_constants = np.array([])
        self._normal_modes = np.array([])
        self._internal_coordinates = InternalCoordinates()
        self.internal_coordinates = self._internal_coordinates.internal_coordinates

        if elements is not None:
            elements = convert_elements(elements)
            masses = np.array([atomic_masses[i] for i in elements])
        else:
            masses = None
        if coordinates is not None:
            coordinates = np.array(coordinates)

        self._masses = masses
        self._elements = elements
        self._coordinates = coordinates

    def add_internal_coordinate(self, atoms):
        """Add internal coordinate composed of two (bond), three (angle) or four
        atoms (dihedral).

        Args:
            atoms (list): Atoms of internal coordinate.
        """
        self._internal_coordinates.add_internal_coordinate(atoms)

    def compute_compliance(self):
        """Compute local force constants according to the compliance matrix
        method. 
        """
        # Compute B matrix if it does not exists
        if len(self._B) == 0:
            self._B = self._internal_coordinates.get_B_matrix(self._coordinates)

        # Compute compliance matrix and get local force constants
        C = self._B @ np.linalg.pinv(self._fc_matrix) @ self._B.T
        k_s = 1 / np.diag(C) * AFU / BOHR / (DYNE / 1000) * ANGSTROM

        self.local_force_constants = k_s

    def compute_frequencies(self):
        """Compute local frequencies"""
        # Compute local frequencies
        M = np.diag(np.repeat(self._masses, 3))
        G = self._B @ np.linalg.inv(M) @ self._B.T
        frequencies = np.sqrt((np.diag(G) / AMU \
            * (self.local_force_constants / ANGSTROM * DYNE / 1000)) \
            / (4 * np.pi ** 2 * C ** 2)) / 100

        self.local_frequencies = frequencies

    def compute_local(self, project_imag=True, cutoff=1e-3):
        """Compute local force constants according to the local modes approach.
        
        Args:
            project_imag (bool): Whether to project out imaginary frequencies
            cutoff (float): Cutoff for low force constant (mDyne/Å)
        """
        # Compute D matrix from normal modes and B matrix
        if len(self._D) == 0:
            if len(self._B) == 0:
                self._B = self._internal_coordinates.get_B_matrix(self._coordinates)
            self._D = self._B @ self._normal_modes.T

        # Project out imaginary modes
        if project_imag and self.n_imag:
            # Save full D matrix and force constant vector
            self._D_full = self._D
            self._force_constants_full = self._force_constants

            # Remove force constants with imaginary force constants
            self._force_constants = self._force_constants[self.n_imag:]
            self._D = self._D[:, self.n_imag:]

        # Add cutoff to modes with low force constants
        if cutoff is not None:
            indices_small = np.where(np.array(self._force_constants)
                < cutoff)[0]
            if len(indices_small) > 0:
                self._force_constants[indices_small] \
                    = np.array(cutoff).reshape(-1)

        # Compute local mode force constants
        K = np.diag(self._force_constants)

        k_s = []
        for row in self._D:
            if not np.all(np.isclose(row, 0)):
                k = 1 / (row @ np.linalg.inv(K) @ np.conj(row).T)
                k_s.append(k)
            else:
                k_s.append(0.0)
        k_s = np.array(k_s)

        # Scale force constants due to projection of imaginary normal modes
        if project_imag and self.n_imag:
            lengths = np.linalg.norm(self._D, axis=1)
            lengths_full = np.linalg.norm(self._D_full, axis=1)
            k_s *= lengths / lengths_full

        self.local_force_constants = k_s

    def detect_bonds(self):
        """Detect bonds based on scaled sum of covalent radii"""
        if len(self._elements) > 0 and len(self._coordinates) > 0:
            self._internal_coordinates.detect_bonds(self._elements, self._coordinates)
        else:
            raise Exception("Elements or coordinates missing.")

    def get_local_force_constant(self, atoms):
        """Return the local force constant between a set of atoms.
        
        Args:
            atoms (list): Atoms of the internal coordinate
        
        Returns:
            force_constant (float): Local force constant (mDyne/Å)
        """
        coordinate = self._get_internal_coordinate(atoms)
        index = self.internal_coordinates.index(coordinate)
        if index == None:
            raise Exception(f"No internal coordinate with these atoms.")
        force_constant = self.local_force_constants[index]

        return force_constant

    def get_local_frequency(self, atoms):
        """Return the local frequency between a set of atoms.
        
        Args:
            atoms (list): Atoms in the internal coordinate
        
        Returns:
            frequency (float): Local frequency (cm^-1)
        """
        coordinate = self._get_internal_coordinate(atoms)
        index = self.internal_coordinates.index(coordinate)
        if index == None:
            raise Exception(f"No internal coordinate with these atoms.")
        frequency = self.local_frequencies[index]

        return frequency

    def load_file(self, filename, program, filetype):
        """Load data from external file
        
        Args:
            filename (str): File
            program (str): Program used to generate file: "gaussian",
                "unimovib" or "xtb"
            filetype (str): Type of file. For "gaussian": "fchk" or "log". For 
                "unimovib": "local" or "log". For "xtb": "hessian" or 
                "normal_modes"
        """
        choices = {"gaussian": {"fchk": self._parse_gaussian_fchk,
                                "log": self._parse_gaussian_log},
                   "xtb": {"hessian": self._parse_xtb_hessian,
                           "normal_modes": self._parse_xtb_normal_modes},
                   "unimovib": {"local": self._parse_unimovib_local,
                                "log": self._parse_unimovib_log,
                                "umv": self._parse_unimovib_umv,}
                   }
        choices[program][filetype](filename)

    def normal_mode_analysis(self, hessian=None, save_hessian=False):
        """Perform normal mode analysis with projection of translations
        and vibrations to get normal modes and force constants.
        
        Args:
            save_hessian (bool): Save projected Hessian for use in compliance
                matrix method.
        """
        # Set up
        coordinates = np.array(self._coordinates) * ANGSTROM_TO_BOHR
        masses = self._masses
        if not hessian:
            hessian = self._fc_matrix
        n_atoms = len(coordinates)

        # Create mass matrices
        M_minus = np.diag(np.repeat(masses, 3)**(-1/2))
        M_plus = np.diag(np.repeat(masses, 3)**(1/2))
        m_plus = np.repeat(masses, 3)**(1/2)
        m_minus = np.repeat(masses, 3)**(-1/2)

        # Mass-weight Hessian
        hessian_mw = M_minus @ hessian @ M_minus

        # Find center of mass
        com = np.sum(masses.reshape(-1, 1) * coordinates, axis=0) / np.sum(masses)

        # Shift origin to center of mass
        coordinates -= com

        # Construct translation vectors
        t_x = (np.tile(np.array([1, 0, 0]), n_atoms)).reshape(-1, 3)
        t_y = (np.tile(np.array([0, 1, 0]), n_atoms)).reshape(-1, 3)
        t_z = (np.tile(np.array([0, 0, 1]), n_atoms)).reshape(-1, 3)

        # Construct mass-weighted rotation vectors
        R_x = np.cross(coordinates, t_x).flatten() * m_plus
        R_y = np.cross(coordinates, t_y).flatten() * m_plus
        R_z = np.cross(coordinates, t_z).flatten() * m_plus

        # Mass-weight translation vectors
        T_x = t_x.flatten() * m_plus
        T_y = t_y.flatten() * m_plus
        T_z = t_z.flatten() * m_plus

        # Remove linear dependencies from translation/rotation space
        TR_vectors = np.vstack([T_x, T_y, T_z, R_x, R_y, R_z])
        Q, R = np.linalg.qr(TR_vectors.T)
        keep_indices = ~np.isclose(np.diag(R), 0, atol=1e-6, rtol=0)
        TR_vectors = Q.T[keep_indices]
        n_tr = len(TR_vectors)

        # Construct P matrix
        P = np.identity(n_atoms * 3)
        for vector in TR_vectors:
            P -= np.outer(vector, vector)

        # Project out translations and rotations
        hessian_proj = P.T @ hessian_mw @ P

        # Diagonalize
        eigenvalues, eigenvectors = np.linalg.eigh(hessian_proj)
        eigenvalues = eigenvalues[n_tr:]
        eigenvectors = eigenvectors[:, n_tr:]

        # Calculate cartesian displacements
        cart = eigenvectors.T * m_minus
        N = 1 / np.linalg.norm(cart, axis=1)
        norm_cart = cart * N.reshape(-1, 1)
        reduced_masses = N ** 2

        # Calculate frequencies and force constants
        n_imag = np.sum(eigenvalues < 0 )
        frequencies = np.sqrt(np.abs(eigenvalues) * HARTREE / BOHR ** 2 / AMU) / (2 * np.pi * C) / 100
        frequencies[:n_imag] = -frequencies[:n_imag]
        force_constants = 4 * np.pi ** 2 * (frequencies * 100) ** 2 * C ** 2 * reduced_masses * AMU / (DYNE / 1000) * ANGSTROM

        # Set up attributes
        self.n_imag = n_imag
        self._force_constants = force_constants
        self._normal_modes = norm_cart
        if save_hessian:
            self._fc_matrix = M_plus @ hessian_proj @ M_plus

    def print_report(self, angles=False, dihedrals=False, angle_units=False):
        """Print report of results.
        
        Args:
            angle_units (bool): Convert angle and dihedral force constants to
                mDyne Å rad^(-2)
            angles (bool): Whether to print angles
            dihedrals (bool): Whether to print dihedrals
        """
        # Print header
        if angle_units:
            unit = "mDyne/Å, mDyne Å rad^(-2)"
        else:
            unit = "mDyne/Å"

        string = f"{'Coordinate':30s}" + \
            f"{'Force constant ' + '(' + unit + ')':>50s}"
        if len(self.local_frequencies) > 0:
            string += f"{'Frequency (cm^-1)':>30s}"
        print(string)

        # Print results for each internal
        sorted_coordinates = sorted(self.internal_coordinates, key=lambda x: (len(x.atoms), *x.atoms))
        for coordinate in sorted_coordinates:
            # Check if internal coordinate is angle or dihedral
            if len(coordinate.atoms) == 3 and not angles:
                continue
            if len(coordinate.atoms) == 4 and not dihedrals:
                continue
            index = self.internal_coordinates.index(coordinate)
            force_constant = self.local_force_constants[index]
            if len(self.local_frequencies) > 0:
                frequency = self.local_frequencies[index]

            # Convert units for angles and dihedrals
            if len(coordinate.atoms) > 2 and angle_units:
                force_constant = force_constant * BOHR_TO_ANGSTROM ** 2

            # Print out the results
            string = f"{repr(coordinate):30s}" + f"{force_constant:50.3f}"
            if len(self.local_frequencies) > 0:
                string += f"{frequency:30.0f}"
            print(string)

    def reset_internal_coordinates(self):
        """Reset internal coordinate system"""
        self._internal_coordinates = InternalCoordinates()
        self.internal_coordinates = self._internal_coordinates.internal_coordinates

    @staticmethod
    def _get_internal_coordinate(atoms):
        # Return bond, angle or dihedral
        if len(atoms) == 2:
            return Bond(*atoms)
        elif len(atoms) == 3:
            return Angle(*atoms)
        elif len(atoms) == 4:
            return Dihedral(*atoms)

    def _parse_gaussian_fchk(self, filename):
        # Read fchk file
        with open(filename) as file:
            lines = file.readlines()

        # Set up read flags
        read_modes = False
        read_hessian = False
        read_vib_e2 = False
        read_ic = False
        read_atomic_numbers = False
        read_coordinates = False
        read_masses = False

        # Set up containers for reading data
        modes = []
        hessian = []
        vib_e2 = []
        internal_coordinates = []
        n_atoms = None
        atomic_numbers = []
        masses = []
        coordinates = []

        # Parse fchk file
        for line in lines:
            # Read normal modes
            if read_modes:
                try:
                    split_line = line.strip().split()
                    values = [float(value) for value in split_line]
                    modes.extend(values)
                except ValueError:
                    read_modes = False
            # Read cartesian force constants
            if read_hessian:
                try:
                    split_line = line.strip().split()
                    values = [float(value) for value in split_line]
                    hessian.extend(values)
                except ValueError:
                    read_hessian = False
            # Read normal mode force constants
            if read_vib_e2:
                try:
                    split_line = line.strip().split()
                    values = [float(value) for value in split_line]
                    vib_e2.extend(values)
                except ValueError:
                    read_vib_e2 = False
            # Read internal coordinates
            if read_ic:
                try:
                    split_line = line.strip().split()
                    values = [int(value) for value in split_line]
                    internal_coordinates.extend(values)
                except ValueError:
                    read_ic = False
            # Read atomic numbers
            if read_atomic_numbers:
                try:
                    split_line = line.strip().split()
                    values = [int(value) for value in split_line]
                    atomic_numbers.extend(values)
                except ValueError:
                    read_atomic_numbers = False
            # Read atomic masses
            if read_masses:
                try:
                    split_line = line.strip().split()
                    values = [float(value) for value in split_line]
                    masses.extend(values)
                except ValueError:
                    read_masses= False
            # Read coordinates
            if read_coordinates:
                try:
                    split_line = line.strip().split()
                    values = [float(value) for value in split_line]
                    coordinates.extend(values)
                except ValueError:
                    read_coordinates = False
            # Read number of atoms
            if "Number of atoms" in line:
                n_atoms = int(line.strip().split()[4])
            # Read number of normal modes
            if "Number of Normal Modes" in line:
                n_modes = int(line.strip().split()[5])
            # Read number of internal coordinates
            if "Redundant internal coordinates" in line:
                n_redundant = int(line.strip().split()[5])
            # Detect when to read data
            if "Vib-Modes" in line:
                read_modes = True
            if "Vib-AtMass" in line:
                read_masses = True
            if "Cartesian Force Constants" in line:
                read_hessian = True
            if "Vib-E2 " in line:
                read_vib_e2 = True
            if "Redundant internal coordinate indices" in line:
                read_ic = True
            if "Atomic numbers" in line:
                read_atomic_numbers = True
            if "Current cartesian coordinates " in line:
                read_coordinates = True

        # Take out normal mode force constants
        force_constants = np.array(vib_e2[n_modes * 2:n_modes * 3])

        # Construct force constant matrix from lower triangular matrix
        fc_matrix = np.zeros((n_atoms * 3, n_atoms * 3))
        fc_matrix[np.tril_indices_from(fc_matrix)] = hessian
        fc_matrix = np.triu(fc_matrix.T, 1) + fc_matrix

        # Take out the internal coordinates
        internal_coordinates = np.array(internal_coordinates)
        for i, coordinate in enumerate(np.split(internal_coordinates,
                                                n_redundant)):
            if all(coordinate >= 0): # Sort out linear bends
                atoms = [i for i in coordinate if i != 0]
                self.add_internal_coordinate(atoms)

        # Convert coordinates to right Ångström
        coordinates = np.array(coordinates).reshape(-1, 3) * BOHR_TO_ANGSTROM

        # Set up attributes
        self._fc_matrix = fc_matrix
        self._normal_modes = np.array(modes).reshape(n_modes, n_atoms * 3)
        self._force_constants = force_constants
        self._elements = convert_elements(atomic_numbers)
        self._coordinates = coordinates
        self._masses = np.array(masses)

    def _parse_gaussian_log(self, filename):
        # Read the log file
        with open(filename) as file:
            lines = file.readlines()

        # Set up read flags
        read_b_atoms = False
        read_b_vectors = False
        read_fc_matrix = False
        read_ifc_matrix = False
        read_input_orientation = False
        read_standard_orientation = False
        read_internal = False
        read_internal_modes = False
        read_hp_modes = False
        read_masses = False

        # Set up containers for reading data
        B_atom_map = {}
        B_vectors = {}
        normal_modes = []
        internal_modes = []
        force_constants = []
        masses = []
        fc_matrix = np.array([])
        ifc_matrix = np.array([])
        input_coordinates = []
        standard_coordinates = []
        n_imag = None
        n_atoms = None
        internal_indices = {}
        atomic_numbers = []

        # Parse through log file content
        counter = 0
        internal_names = []
        internal_vector = []
        coordinates = []
        for line in lines:
            # Read internal coordinate definitions
            if read_internal:
                if counter > 1:
                    if " --------------------------------------------------------------------------------" in line:
                        read_internal = False
                        n_internals = len(internal_indices.items())
                    else:
                        split_line = line.strip().split()
                        name = split_line[1]
                        internal_names.append(name)
                        atoms = split_line[2][2:].replace(")", "").split(",")
                        atoms = [int(atom) for atom in atoms]
                        internal_indices[frozenset(atoms)] = counter - 2
                counter += 1
            # Read Cartesian force constant matrix
            if read_fc_matrix:
                if  " FormGI is forming" in line:
                    read_fc_matrix = False
                    fc_matrix = np.triu(fc_matrix.T, 1) + fc_matrix
                elif "          " in line:
                    column_indices = [int(value) for value
                                      in line.strip().split()]
                elif "D" in line:
                    split_line = line.strip().split()
                    row_index = int(split_line[0]) - 1
                    values = split_line[1:]
                    values = [float(value.replace("D", "E")) for value
                              in values]
                    for i, value in enumerate(values):
                        column_index = column_indices[i] - 1
                        fc_matrix[row_index, column_index] = value
            # Read internal force constant matrix
            if read_ifc_matrix:
                if  "Leave Link  716" in line:
                    read_ifc_matrix = False
                    ifc_matrix = np.triu(ifc_matrix.T, 1) + ifc_matrix
                elif "          " in line:
                    column_indices = [int(value) for value
                                      in line.strip().split()]
                elif "D" in line:
                    split_line = line.strip().split()
                    row_index = int(split_line[0]) - 1
                    values = split_line[1:]
                    values = [float(value.replace("D", "E")) for value
                              in values]
                    for i, value in enumerate(values):
                        column_index = column_indices[i] - 1
                        ifc_matrix[row_index, column_index] = value
            # Read atoms for creation of B matrix
            if read_b_atoms:
                if " B Matrix in FormBX:" in line:
                    read_b_atoms = False
                else:
                    if counter == 0:
                        atoms = [int(value) for value in line.strip().split()]
                        for atom in atoms:
                            B_atom_map[atom] = []
                    if counter > 0 and counter < 5:
                        values = [int(value) for value
                                  in line.strip().split()[1:]]
                        for atom, value in zip(atoms, values):
                            B_atom_map[atom].append(value)
                    counter += 1
                    if counter == 5:
                        counter = 0
            # Read values of B matrix
            if read_b_vectors:
                if  " IB Matrix in Red2BG:" in line or "Iteration" in line \
                        or " G Matrix:" in line:
                    read_b_vectors = False
                else:
                    if counter == 0:
                        atoms = [int(value) for value in line.strip().split()]
                        for atom in atoms:
                            B_vectors[atom] = []
                    if counter > 0 and counter < 13:
                        values = [float(value.replace("D", "E")) for value
                                  in line.strip().split()[1:]]
                        for atom, value in zip(atoms, values):
                            B_vectors[atom].append(value)
                    counter += 1
                    if counter == 13:
                        counter = 0
            # Read atomic coordinates in input orientation
            if read_input_orientation:
                if counter > 3:
                    if "---------------------------------------------------------------------" in line:
                        read_input_orientation = False
                    else:
                        strip_line = line.strip().split()
                        values = [float(value) for value in strip_line[3:]]
                        input_coordinates.append(values)
                        atomic_numbers.append(int(strip_line[1]))
                counter += 1
            # Read atomic coordinates in standard orientation
            if read_standard_orientation:
                if counter > 3:
                    if "---------------------------------------------------------------------" in line:
                        read_standard_orientation = False
                    else:
                        strip_line = line.strip().split()
                        values = [float(value) for value in strip_line[3:]]
                        standard_coordinates.append(values)
                counter += 1
            # Read decomposition of normal modes in internal coordinates
            if read_internal_modes:
                if counter > 3:
                    if "--------------------------------------------------------------------------------" in line:
                        read_internal_modes = False
                        internal_modes.append(internal_vector)
                    else:
                        value = float(line.strip().split()[3])
                        internal_vector.append(value)
                counter += 1
            # Read high-precision normal modes
            if read_hp_modes:
                if counter < n_atoms * 3:
                    strip_line = line.strip().split()
                    values = [float(value) for value in strip_line[3:]]
                    coordinates.append(values)
                if counter == n_atoms * 3:
                    coordinates = np.array(coordinates)
                    normal_modes.append(coordinates)
                    read_hp_modes = False
                counter += 1
            # Read atomic masses
            if read_masses:
                if "Molecular mass: " in line:
                    read_masses = False
                elif "and mass" in line:
                    masses.append(float(line.strip().split()[8]))
            # Read number of atoms
            if " NAtoms=" in line and not n_atoms :
                n_atoms = int(line.strip().split()[1])
            # Read normal mode force constants
            if "Frc consts" in line:
                split_line = line.strip().split()
                values = [float(value) for value in split_line[3:]]
                force_constants.extend(values)
            # Read number of imaginary frequencies
            if "imaginary frequencies (negative Signs)" in line:
                n_imag = int(line.strip().split()[1])
            # Detect when to read data
            if "Name  Definition              Value          Derivative Info." in line:
                read_internal = True
                counter = 1
                internal_names = []
            if "- Thermochemistry -" in line:
                read_masses = True
            if " IB Matrix in FormBX:" in line:
                read_b_atoms = True
                counter = 0
            if " B Matrix in FormBX:" in line:
                read_b_vectors = True
                counter = 0
            if " Force constants in Cartesian coordinates: " in line:
                read_fc_matrix = True
                fc_matrix = np.zeros((3 * n_atoms, 3 * n_atoms))
                counter = 0
            if " Force constants in internal coordinates: " in line:
                read_ifc_matrix = True
                ifc_matrix = np.zeros((n_internals, n_internals))
                counter = 0
            if "Input orientation: " in line:
                read_input_orientation = True
                counter = 0
            if "Standard orientation: " in line:
                read_standard_orientation = True
                counter = 0
            if "Normal Mode" in line:
                read_internal_modes = True
                counter = 1
                internal_vector = []
            if " Coord Atom Element:" in line:
                read_hp_modes = True
                coordinates = []
                counter = 0

        # Process internal coordinates
        if len(internal_indices) > 0:
            for name, indices in zip(internal_names, internal_indices):
                if name[0] == "R" and len(indices) == 2:
                    self._internal_coordinates.add_internal_coordinate(indices)
                if name[0] == "A" and len(indices) == 3:
                    self._internal_coordinates.add_internal_coordinate(indices)
                if name[0] == "D" and len(indices) == 4:
                    self._internal_coordinates.add_internal_coordinate(indices)

        # Construct the B matrix from atoms and vectors
        if B_vectors:
            n_cartesian = n_atoms * 3
            B = np.zeros((n_internals, n_cartesian))
            for i in range(n_internals):
                for j, atom in enumerate(B_atom_map[i + 1]):
                    if atom:
                        B[i][(atom - 1) * 3:(atom - 1) * 3 + 3] \
                            = B_vectors[i + 1][j * 3 :j * 3 + 3]
            B_inv = np.linalg.pinv(B)
            # Detect whether the internal coordinate system is redundant
            if B.shape[0] == len(force_constants):
                self._redundant = False
            else:
                self._redundant = True
        else:
            B = np.array([])
            B_inv = np.array([])

        # Detect whether the input coordinates have been rotated. If so, rotate
        # B matrix and its inverse.
        input_coordinates = np.array(input_coordinates).reshape(-1, 3)
        standard_coordinates = np.array(standard_coordinates).reshape(-1, 3)

        if not np.array_equal(input_coordinates, standard_coordinates) \
                and standard_coordinates.size > 0:
            rotation_i_to_s = kabsch_rotation_matrix(input_coordinates,
                                                     standard_coordinates)
            if len(B) > 0:
                B = (rotation_i_to_s @ B.reshape(-1, 3).T).T.reshape(
                        n_internals, -1)
                B_inv = np.linalg.pinv(B)

        # Set up attributes
        self._B = B
        self._B_inv = B_inv
        self._fc_matrix = fc_matrix
        if fc_matrix.size == 0 and ifc_matrix.size > 0:
            self._fc_matrix = B.T @ ifc_matrix @ B
        self._ifc_matrix = ifc_matrix
        self._force_constants = np.array(force_constants)
        if len(normal_modes) > 0:
            self._normal_modes = np.hstack(normal_modes).T
        else:
            self._normal_modes = np.array([])
        if len(internal_modes) > 0:
            self._D = np.vstack(internal_modes).T
        self.n_imag = n_imag
        self._masses = np.array(masses)
        self._input_coordinates = input_coordinates
        self._standard_coordinates = standard_coordinates
        self._coordinates = input_coordinates
        self._elements = convert_elements(atomic_numbers)

    def _parse_unimovib_local(self, filename):
        # Read file
        with open(filename) as file:
            lines = file.readlines()

        # Set up flags for reading
        read_masses = False
        read_atomic_numbers = False
        read_coordinates = False
        read_hessian = False
        read_normal_modes = False

        # Set up containers for data
        masses = []
        atomic_numbers = []
        coordinates = []
        hessian = []
        normal_modes = []

        # Parse file
        for line in lines:
            # Read atomic masses
            if read_masses:
                if " $ZA  $END" in line:
                    read_masses = False
                else:
                    split_line = line.strip().split()
                    masses.extend([float(value.replace("D", "E")) for value in split_line])
            # Read atomic numbers
            if read_atomic_numbers:
                if " $XYZ  $END" in line:
                    read_atomic_numbers = False
                else:
                    split_line = line.strip().split()
                    atomic_numbers.extend([int(float(value.replace("D", "E"))) for value in split_line])
            # Read coordinates
            if read_coordinates:
                if " $FFX  $END" in line:
                    read_coordinates = False
                else:
                    split_line = line.strip().split()
                    coordinates.extend([float(value.replace("D", "E")) for value in split_line])
            # Read Hessian
            if read_hessian:
                if " $NMMODE  $END" in line:
                    read_hessian = False
                else:
                    split_line = line.strip().split()
                    hessian.extend([float(value.replace("D", "E")) for value in split_line])
            # Read normal modes
            if read_normal_modes:
                if " $APT" in line:
                    read_normal_modes = False
                else:
                    split_line = line.strip().split()
                    normal_modes.extend([float(value.replace("D", "E")) for value in split_line])
            # Read number of atoms and normal modes
            if " $CONTRL" in line:
                split_line = line.strip().split()
                n_atoms = int(split_line[1].split("=")[1])
                n_modes = int(split_line[2].split("=")[1])
            # Detect when to read data
            if " $AMASS  $END" in line:
                read_masses = True
            if " $ZA  $END" in line:
                read_atomic_numbers = True
            if " $XYZ  $END" in line:
                read_coordinates = True
            if " $FFX  $END" in line:
                read_hessian = True
            if " $NMMODE  $END" in line:
                read_normal_modes = True

        # Convert data to right format
        coordinates = np.array(coordinates).reshape(-1, 3) * BOHR_TO_ANGSTROM
        hessian = np.array(hessian).reshape(n_atoms * 3, n_atoms * 3)
        normal_modes = np.array(normal_modes).reshape(-1, n_atoms * 3)[:n_modes]
        elements = convert_elements(atomic_numbers)

        # Set up attributes
        self._normal_modes = normal_modes
        self._masses = np.array(masses)
        self._fc_matrix = hessian
        self._coordinates = coordinates
        self._elements = elements

    def _parse_unimovib_log(self, filename):
        # Read file
        with open(filename) as file:
            lines = file.readlines()

        # Set up read flags
        read_coordinates = False
        read_vibrations = False
        read_normal_modes = False

        # Set up data containers
        atomic_numbers = []
        masses = []
        coordinates = []
        force_constants = []
        normal_modes = []
        frequencies = []

        # Parse file
        normal_modes_chunk = [[]]
        counter = 0
        n_modes_chunk = 0
        for line in lines:
            # Read coordinates, atomic numbers and atomic masses
            if read_coordinates:
                if counter > 0:
                    if " ------------------------------------------------------------------------------------------" in line:
                        read_coordinates = False
                    else:
                        split_line = line.strip().split()
                        atomic_numbers.append(int(split_line[2]))
                        coordinates.extend([float(value) for value in split_line[3:6]])
                        masses.append(float(split_line[6]))
                counter += 1
            # Read normal modes and force constants
            if read_vibrations:
                if " Results of translations and rotations:" in line:
                    read_vibrations = False
                if read_normal_modes:
                    split_line = line.strip().split()
                    if len(split_line) == 0:
                        read_normal_modes = False
                        normal_modes.extend(normal_modes_chunk)
                    else:
                        values = [float(value) for value in split_line[2:]]
                        for i in range(n_modes_chunk):
                            normal_modes_chunk[i].append(values[i * 3: i * 3 + 3])
                elif "Irreps" in line:
                    split_line = line.strip().split()
                    n_modes_chunk = len(split_line) - 1
                    normal_modes_chunk = [[] for i in range(n_modes_chunk)]
                elif "Force constants" in line:
                    split_line = line.strip().split()
                    force_constants.extend([float(value) for value in split_line[2:]])
                elif "Frequencies" in line:
                    split_line = line.strip().split()
                    frequencies.extend([float(value) for value in split_line[2:]])
                elif "        Atom  ZA" in line:
                    read_normal_modes = True
            # Detect when to read data
            if "No.   Atom    ZA" in line:
                read_coordinates = True
                counter = 0
            if "Results of vibrations:" in line:
                read_vibrations = True

        # Convert quantities to right format
        n_atoms = len(masses)
        coordinates = np.array(coordinates).reshape(-1, 3)
        normal_modes = np.array(normal_modes).reshape(-1, n_atoms * 3)
        elements = convert_elements(atomic_numbers)
        force_constants = np.array(force_constants)
        frequencies = np.array(frequencies)

        # Detect imaginary modes
        self.n_imag = np.sum(frequencies < 0)

        # Set up attributes
        self._elements = elements
        self._coordinates = coordinates
        self._masses = masses
        self._force_constants = force_constants
        self._normal_modes = normal_modes

    def _parse_unimovib_umv(self, filename):
        # Read file
        with open(filename) as file:
            lines= file.readlines()

        # Set up flags for reading
        read_n_atoms = False
        read_masses = False
        read_atomic_numbers = False
        read_coordinates = False
        read_hessian = False

        # Set up containers for data
        n_atoms = None
        masses = []
        atomic_numbers = []
        coordinates = []
        hessian = []

        # Parse file
        for line in lines:
            # Read number of atoms
            if read_n_atoms:
                if "AMASS" in line:
                    read_n_atoms = False
                else:
                    n_atoms = int(line.strip())
            # Read atomic masses
            if read_masses:
                if "ZA" in line:
                    read_masses = False
                else:
                    split_line = line.strip().split()
                    masses.extend([float(value.replace("D", "E")) for value in split_line])
            # Read atomic numbers
            if read_atomic_numbers:
                if "XYZ" in line:
                    read_atomic_numbers = False
                else:
                    split_line = line.strip().split()
                    atomic_numbers.extend([int(float(value.replace("D", "E"))) for value in split_line])
            # Read coordinates
            if read_coordinates:
                if "FFX" in line:
                    read_coordinates = False
                else:
                    split_line = line.strip().split()
                    coordinates.extend([float(value.replace("D", "E")) for value in split_line])
            # Read Hessian
            if read_hessian:
                if "APT" in line:
                    read_hessian = False
                else:
                    split_line = line.strip().split()
                    hessian.extend([float(value.replace("D", "E")) for value in split_line])
            # Read number of atoms and normal modes
            if "NATM" in line:
                read_n_atoms = True
            # Detect when to read data
            if "AMASS" in line:
                read_masses = True
            if "ZA" in line:
                read_atomic_numbers = True
            if "XYZ" in line:
                read_coordinates = True
            if "FFX" in line:
                read_hessian = True

        # Convert data to right format
        coordinates = np.array(coordinates).reshape(-1, 3) * BOHR_TO_ANGSTROM
        hessian = np.array(hessian).reshape(n_atoms * 3, n_atoms * 3)
        elements = convert_elements(atomic_numbers)

        # Set up attributes
        self._masses = np.array(masses)
        self._fc_matrix = hessian
        self._coordinates = coordinates
        self._elements = elements

    def _parse_xtb_hessian(self, filename):
        # Read hessian file
        with open(filename) as file:
            lines = file.readlines()

        # Parse file
        hessian = []
        for line in lines:
            try:
                hessian.extend([float(value) for value in line.strip().split()])
            except ValueError:
                pass
        # Set up force constant matrix
        dimension = int(np.sqrt(len(hessian)))
        self._fc_matrix = np.array(hessian).reshape(dimension, dimension)

    def _parse_xtb_normal_modes(self, filename):
        # Read in data from xtb normal modes file
        with FortranFile(filename) as file:
            n_modes = file.read_ints()[0]
            frequencies = file.read_reals()
            reduced_masses = file.read_reals()
            xtb_normal_modes = file.read_reals().reshape(n_modes, - 1)

        # Find the number of rotations and translations
        n_tr_1  = np.sum(frequencies == 0)
        n_tr_2 = np.sum(reduced_masses == 0)
        n_tr_3 = np.sum(np.sum(xtb_normal_modes, axis=1) == 0)
        if n_tr_1 == n_tr_2 == n_tr_3:
            frequencies = frequencies[n_tr_1:]
            reduced_masses = reduced_masses[n_tr_1:]
            xtb_normal_modes = xtb_normal_modes[n_tr_1:]
        else:
            raise Exception("Can't determine the number of translations and rotations.")

        # Un-mass-weight normal modes, calculate reduced masses and normalize
        m_inverse = 1 / np.repeat(self._masses, 3)
        cart_disp = xtb_normal_modes * np.sqrt(m_inverse)
        N = 1 / np.linalg.norm(cart_disp, axis=1)
        norm_cart_disp = cart_disp * N.reshape(-1, 1)
        reduced_masses = N ** 2

        # Calculate force constants
        force_constants = 4 * np.pi ** 2 * (frequencies * 100) ** 2 * C ** 2 * reduced_masses * AMU
        force_constants = force_constants / (DYNE / 1000) * ANGSTROM

        # Detect imaginary modes
        n_imag = np.sum(frequencies < 0)

        self._normal_modes = norm_cart_disp
        self._force_constants = force_constants
        self.n_imag = n_imag

    def __repr__(self):
        n_internal = len(self.internal_coordinates)
        return f"{self.__class__.__name__}({n_internal!r} internal coordinates)"


class Pyramidalization:
    """Calculates and stores results of pyramidalization and alpha angle as 
    described in Struct. Chem. 1991, 2, 107 and pyramidalization according to
    bond angles as in J. Comput. Chem. 2012, 33 (27), 2173–2179.

    Args:
        atom_index (int): Index of pyramidalized atom (1-indexed)
        coordinates (list): Coordinates (Å)
        elements (list): Elements as atomic symbols or numbers. 
        excluded_atoms (list): Indices of atoms to exclude. 
        method (str): Method for detecting neighbors. Either 'connectivity' or 
                      'distance'. Ignored if neighbor_indices is given.
        neighbor_indices (list): Indices of neighbors to pyramidalized atom.
        radii (list): Covalent radii used to determine connectivity (Å)
        radii_type (str): Choice of covalent radii: 'pyykko' (default)

    Attributes:
        alpha (float): Average alpha angle (deg)
        alphas (float): Alpha angles for all permutations of neighbors (deg)
        P (float): Pyramidalization according to Radhakrishnan
        P_angle (float): Pyramidalization according to Gavrish
    """
    def __init__(self, coordinates, atom_index, neighbor_indices=[], elements=[], radii=[], radii_type="pyykko", excluded_atoms=[],
                 method="distance"):
        coordinates = np.array(coordinates)
        atom_coordinates = coordinates[atom_index - 1]
        excluded_atoms = np.array(excluded_atoms, dtype=int)

        # Get 3 closest neighbors
        if len(neighbor_indices) > 0:
            if len(neighbor_indices) != 3:
                raise Exception(f"Only {len(neighbor_indices)} neighbors.")
            neighbors = np.array(neighbor_indices) - 1
        elif method == "distance":
            # Generate mask for excluded atoms
            mask = np.zeros(len(coordinates), dtype=bool)
            mask[excluded_atoms - 1] = True
            mask[atom_index - 1] = True

            # Get three closest atoms not in the excluded atoms
            distances = cdist(atom_coordinates.reshape(1, -1), coordinates).reshape(-1)
            distances[mask] = np.inf
            neighbors = np.argsort(distances)[:3]
        elif method == "connectivity":
            # Construct connectivity matrix and get closest neighbors.
            if not (len(elements) > 0 or len(radii) > 0):
                raise Exception("Connectivity requires elements or radii.")
            connectivity_matrix = get_connectivity_matrix(elements, coordinates, radii_type="pyykko")
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
        alphas = []
        vectors = []
        cos_alphas = []
        thetas = []
        for v_1, v_2, v_3 in itertools.permutations([a, b, c], 3):
            # Calculate cos_alpha
            normal = np.cross(v_1, v_2)
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


@requires_dependency([Import("xtb"), Import("xtb.interface")], globals())
class XTB:
    """Calculates electronic properties with the xtb program using the
    xtb-python package.

    Args:
        charge (int): Molecular charge
        coordinates (list): Coordinates (Å)
        electronic_temperature (float): Electronic temperature (K)
        elements (list): Elements as atomic symbols or numbers. 
        n_unpaired (int): Number of unpaired electrons
        solvent (str): Solvent. See xtb-python documentation.
        version (str): Version of xtb to use. Currently works with '1' or '2'.
    """
    _ipea_corrections = {"1": 5.700, "2": 4.846}

    def __init__(self, elements, coordinates, version="2", charge=0,
                 n_unpaired=None, solvent=None, electronic_temperature=None):
        # Converting elements to atomic numbers if the are symbols
        self._elements = np.array(convert_elements(elements, output="numbers"))

        # Store settings
        self._coordinates = np.array(coordinates)
        self._version = version
        self._charge = charge
        self._solvent = solvent
        self._n_unpaired = n_unpaired
        self._electronic_temperature = electronic_temperature

        # Set up results dictionary
        self._results = {-1: None,
                         0: None,
                         1: None,
                         }

        self._params = {"1": xtb.interface.Param.GFN1xTB,
                        "2": xtb.interface.Param.GFN2xTB,
                        }

    def get_bond_order(self, i, j):
        """Calculate bond order.
        
        Args:
            i (int): Index of atom 1 (1-indexed)
            j (int): Index of atom 2 (1-indexed)
        
        Returns:
            bond_order (float): Bond order
        """
        bo_matrix = self.get_bond_orders()
        bond_order = bo_matrix[i - 1, j - 1]

        return bond_order

    def get_bond_orders(self, charge_state=0):
        self._check_results(charge_state)
        return self._results[charge_state].get_bond_orders()

    def get_charges(self, charge_state=0):
        self._check_results(charge_state)
        return self._results[charge_state].get_charges()

    def get_dipole(self):
        """Calculate dipole vector (a.u.).
        
        Returns:
            dipole (ndarray): Dipole vector. 
        """
        self._check_results(0)
        dipole = self._results[0].get_dipole()

        return dipole

    def get_ea(self, corrected=False):
        """Calculate electron affinity.

        Args:
            corrected (bool): Apply correction term.
        
        Returns:
            ea (float): Electron affinity (eV)
        """
        # Calculate energies
        energy_neutral = self._get_energy(0)
        energy_anion = self._get_energy(-1)

        # Calculate electron affinity
        ea = (energy_neutral - energy_anion) * HARTREE_TO_EV
        if corrected:
            ea -= self._ipea_corrections[self._version]

        return ea

    def get_fukui(self, variety):
        """Calculate Fukui coefficients
        
        Args:
            variety (str): Type of Fukui coefficient: 'nucleophilicity', 
                'electrophilicity', 'radical', 'dual', 'local_nucleophilicity'
                or 'local_electrophilicity'.
        
        Returns
            fukui (ndarray): Array of atom Fukui coefficients.
        """
        if variety == "nucleophilicity":
            fukui = self.get_charges(0) - self.get_charges(1)
        elif variety == "electrophilicity":
            fukui = self.get_charges(-1) - self.get_charges(0)
        elif variety == "radical":
            fukui = (self.get_charges(-1) - self.get_charges(0)) / 2
        elif variety == "dual":
            fukui = 2 * self.get_charges(0) - self.get_charges(1) - \
                self.get_charges(-1)
        elif variety == "local_nucleophilicity":
            fukui = self.get_fukui("nucleophilicity")
        elif variety == "local_electrophilicity":
            chem_pot = -(self.get_ip() + self.get_ea()) / 2
            hardness = self.get_ip() - self.get_ea()
            fukui = -(chem_pot / hardness) * self.get_fukui("radical") + \
                1 / 2 * (chem_pot / hardness) ** 2 * self.get_fukui("dual")
        else:
            raise Exception("Choose variety.")

        return fukui

    def get_global_descriptor(self, variety, corrected=False):
        """Calculate global reactivity descriptors.
        
        Args:
            corrected (bool): Apply correction term
            variety (str): Type of descriptor: 'electrophilicity',
                'nucleophilicity', 'electrofugality' or 'nucleofugality'.
        
        Returns:
            descriptor (float): Global reactivity descriptor (eV).
        """
        if variety == "electrophilicity":
            descriptor = (self.get_ip(corrected=corrected) + \
                 self.get_ea(corrected=corrected)) ** 2 / \
                     (8 * (self.get_ip(corrected=corrected) \
                         - self.get_ea(corrected=corrected)))
        elif variety == "nucleophilicity":
            descriptor = -self.get_ip(corrected=corrected)
        elif variety == "electrofugality":
            descriptor = (3 * self.get_ip(corrected=corrected) - \
                self.get_ea(corrected=corrected)) ** 2 / \
                    (8 * (self.get_ip(corrected=corrected) - \
                        self.get_ea(corrected=corrected)))
        elif variety == "nucleofugality":
            descriptor = (self.get_ip(corrected=corrected) - \
                3 * self.get_ea(corrected=corrected)) ** 2 / \
                    (8 * (self.get_ip(corrected=corrected) - \
                        self.get_ea(corrected=corrected)))

        return descriptor

    def get_homo(self):
        """Calculate HOMO energy.
        
        Returns:
            homo_energy (float): HOMO energy (a.u.)
        """
        eigenvalues = self._get_eigenvalues()
        occupations = self._get_occupations()
        homo_index = int(occupations.sum().round(0) / 2) - 1
        homo_energy = eigenvalues[homo_index]

        return homo_energy

    def get_ip(self, corrected=False):
        """Calculate ionization potential.
        
        Args:
            corrected (bool): Apply correction term

        Returns:
            ip (float): Ionization potential (eV)
        """
        # Calculate energies
        energy_neutral = self._get_energy(0)
        energy_cation = self._get_energy(1)

        # Calculate ionization potential
        ip = (energy_cation - energy_neutral) * HARTREE_TO_EV
        if corrected:
            ip -= self._ipea_corrections[self._version]

        return ip

    def get_lumo(self):
        """Calculate LUMO energy.
        
        Returns:
            lumo_energy (float): LUMO energy (a.u.)
        """
        eigenvalues = self._get_eigenvalues()
        occupations = self._get_occupations()
        homo_index = int(occupations.sum().round(0) / 2) - 1
        lumo_index = homo_index + 1
        lumo_energy = eigenvalues[lumo_index]

        return lumo_energy

    def _check_results(self, charge_state):
        if self._results[charge_state] is None:
            self._sp(charge_state)

    def _get_eigenvalues(self):
        self._check_results(0)
        return self._results[0].get_orbital_eigenvalues()

    def _get_energy(self, charge_state=0):
        self._check_results(charge_state)
        return self._results[charge_state].get_energy()

    def _get_occupations(self):
        self._check_results(0)
        return self._results[0].get_orbital_occupations()

    def _sp(self, charge_state=0):
        # Set up calculator
        calc = xtb.interface.Calculator(self._params[self._version],
                                        self._elements,
                                        self._coordinates * ANGSTROM_TO_BOHR,
                                        charge=self._charge + charge_state,
                                        uhf=self._n_unpaired)
        calc.set_verbosity(xtb.libxtb.VERBOSITY_MUTED)

        # Set solvent
        if self._solvent:
            solvent = xtb.utils.get_solvent(self._solvent)
            if solvent is None:
                raise Exception(f"Solvent '{self._solvent}' not recognized")
            calc.set_solvent(solvent)

        # Set electronic temperature
        if self._electronic_temperature:
            calc.set_electronic_temperature(self._electronic_temperature)

        # Do singlepoint calculation and store the result
        res = calc.singlepoint()
        self._results[charge_state] = res
