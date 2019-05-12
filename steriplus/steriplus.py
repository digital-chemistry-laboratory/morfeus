"""This module contains classes for performing calculations of steric
descriptors of molecules.

Classes:
    Atom: Atom class.
    BuriedVolume: Calculates buried volumes
    ConeAngle: Calculates exact cone angles.
    SASA: Calculates solvent accessible surface area.
    Sterimol: Calculates Sterimol parameters
"""

import math
import itertools
import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

import scipy.linalg
import scipy.spatial
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

from steriplus.data import atomic_numbers, atomic_symbols, \
    bondi_radii, crc_radii, jmol_colors
from steriplus.geometry import rotate_coordinates, Atom, Sphere, Cone
from steriplus.helpers import convert_elements, check_distances, get_radii
from steriplus.io import read_gjf, read_xyz
from steriplus.plotting import MoleculeScene

class Sterimol:
    """Performs and stores results of Sterimol calculation.

    Args:
        atom_1 (int): Index of atom 1 (dummy atom, starting at 1)
        atom_2 (int): Index of atom 2 (connected atom of substituent, starting
                      at 1)
        coordinates (list): Coordinates (Å)
        elements (list): Elements as atomic symbols or numbers
        n_rot_vectors (int): Number of rotational vectors for determining B_1
                             and B_5
        radii (list): List of radii (Å, optional)
        radii_type (str): Type of radii to use: 'bondi' or 'crc'

    Attributes:
        atom_1 (int): Index of atom 1 (dummy atom, starting at 1)
        atom_2 (int): Index of atom 2 (connected atom of substituent, starting
                      at 1)
        B_1 (ndarray): Sterimol B_1 vector (Å)
        B_1_value (float): Sterimol B_1 value (Å)
        B_5 (ndarray): Sterimol B_5 vector (Å)
        B_5_value (float): Sterimol B_5 value (Å)
        bond_length (float): Bond length between atom 1 and atom 2 (Å)
        L (ndarray): Sterimol L vector (Å)
        L_value (float): Sterimol L value (Å)
        L_value_uncorrected (float): Sterimol L value minus 0.40 Å
        vector (ndarray): Vector between atom 1 and atom 2 (Å)
    """
    def __init__(self, elements, coordinates, atom_1, atom_2, radii=[],
                 radii_type="crc", n_rot_vectors=3600):
        # Convert elements to atomic numbers if the are symbols
        element_ids = convert_elements(elements)

        # Get radii if they are not supplied
        if not radii:
            radii = get_radii(element_ids, radii_type=radii_type)
        radii = np.array(radii)

        # Set up coordinate array
        all_coordinates = np.array(coordinates)

        # Translate coordinates so origin is at atom 2
        all_coordinates -= all_coordinates[atom_2 - 1]

        # Get vector pointing from atom 2 to atom 1
        vector_2_to_1 = all_coordinates[atom_2 - 1] \
                        - all_coordinates[atom_1 - 1]
        vector_2_to_1 = vector_2_to_1 / np.linalg.norm(vector_2_to_1)

        # Get rotation quaternion that overlays vector with x-axis
        x_axis = np.array([1, 0, 0])
        all_coordinates = rotate_coordinates(all_coordinates, vector_2_to_1,
                                             x_axis)

        # Get list of atoms as Atom objects
        atoms = []
        for i, (element, radius, coord) in enumerate(
                zip(elements, radii, all_coordinates), start=1):
            atom = Atom(element, coord, radius, i)
            atoms.append(atom)

        coordinates = np.delete(all_coordinates, atom_1 - 1, axis=0)
        radii = np.delete(radii, atom_1 - 1)

        # Project coordinates onto vector between atoms 1 and 2
        vector = all_coordinates[atom_2 - 1] - all_coordinates[atom_1 - 1]
        bond_length = np.linalg.norm(vector)
        unit_vector = vector / np.linalg.norm(vector)

        c_values = np.dot(unit_vector.reshape(1, -1), coordinates.T)
        projected = c_values + radii

        # Get L as largest projection along the vector
        L_value = np.max(projected) + bond_length
        L = vector * L_value

        # Get rotation vectors in yz plane
        r = 1
        theta = np.linspace(0, 2 * math.pi, n_rot_vectors)
        x = np.zeros(len(theta))
        y = r * np.cos(theta)
        z = r * np.sin(theta)
        rot_vectors = np.column_stack((x, y, z))
    
        # Project coordinates onto rotation vectors
        c_values = np.dot(rot_vectors, coordinates.T)
        projected = c_values + radii
        max_c_values = np.max(projected, axis=1)
    
        # Determine B1 and B5 from the smallest and largest scalar projections
        B_1_value = np.min(max_c_values)
        B_1 = rot_vectors[np.argmin(max_c_values)] * B_1_value
    
        B_5_value = np.max(max_c_values)
        B_5 = rot_vectors[np.argmax(max_c_values)] * B_5_value

        # Set up attributes
        self._atoms = atoms
        self.vector = vector.reshape(-1)
        
        self.atom_1 = atom_1
        self.atom_2 = atom_2
        
        self.L = L.reshape(-1)
        self.L_value = L_value + 0.40
        self.L_value_uncorrected = L_value
        self.bond_length = bond_length
        
        self.B_1 = B_1
        self.B_1_value = B_1_value

        self.B_5 = B_5
        self.B_5_value = B_5_value

    def draw_3D(self):
        """Draw a 3D representation of the molecule with the Sterimol vectors"""
        # Set up lists for drawing
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]
        radii = [atom.radius for atom in self._atoms]
        indices = [atom.index for atom in self._atoms]
        
        # Draw molecule scene
        scene = MoleculeScene(elements, coordinates, radii, indices)

        # Set the size of the atoms
        scene.set_scale(0.3)

        # Draw L vector
        L_start = coordinates[self.atom_1 - 1]
        L_stop = self.L
        L_length = self.L_value 
        scene.add_arrow(L_start, L_stop, L_length, "L")

        # Draw B_1 vector
        B_1_start = coordinates[self.atom_2 - 1]
        B_1_stop = self.B_1
        B_1_length = self.B_1_value
        scene.add_arrow(B_1_start, B_1_stop, B_1_length, "B_1")

        # Draw B_5 vector
        B_5_start = coordinates[self.atom_2 - 1]
        B_5_stop = self.B_5
        B_5_length = self.B_5_value
        scene.add_arrow(B_5_start, B_5_stop, B_5_length, "B_5")       

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

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"

class BuriedVolume:
    """Performs and stores the results of a buried volume calculation.

    Args:
        atom_1 (int)                :   Atom index of metal (starting from 1)
        coordinates (list)          :   List of coordinates in Å
        density (float)             :   Density of points on the vdW surface in
                                        Å^-2
        element_ids (list)          :   List of element_ids, either as atomic
                                        numbers or symbols
        exclude_list (list)         :   List of atom indices to exclude from the
                                        calculation (one-indexed)
        include_hs (bool)           :   Flag to include H atoms when calculating
                                        the buried volume or not

        radii (list)                :   List of radii (optional)
        radii_scale (float)         :   Scaling factor for radii. Default set
                                        from original paper.
        radii_type (str)            :   Type of radii to use, either "crc" or
                                        "bondi"

    Parameters:
        atoms (list)                :   List of atoms objects for the ligand
                                        (excluding the exclude list)
        atom_coordinates (ndarray)  :   Array of atom coordinates
        buried_points (ndarray)     :   Array of points buried by the ligand
        buried_volume (float)       :   Buried volume
        free_points (ndarray)       :   Array of free points
        sphere (object)             :   Sphere object at the metal center

    """
    def __init__(self, element_ids, coordinates, atom_1, exclude_list=[],
                 radii=[], include_hs=False, radius=3.5, radii_type="bondi",
                 radii_scale=1.17, density=0.001):
        # Get the coordinates for the central atom
        center = np.array(coordinates[atom_1 - 1])

        # Construct sphere at metal center
        s = Sphere(center, radius, method="projection", density=density,
                   filled=True)

        # Save all coordinates before deleting some atoms
        self._density = density
        self._excluded_atoms = exclude_list
        self._all_coordinates = np.array(coordinates)

        # Make copies of lists to avoid deleting originals
        coordinates = list(coordinates)
        element_ids = list(element_ids)

        # Remove atoms from the exclude list when assessing buried volume
        if exclude_list:
            for entry in sorted(exclude_list, reverse=True):
                del element_ids[entry - 1]
                del coordinates[entry - 1]

        # Converting element ids to atomic numbers if the are symbols
        element_ids = convert_element_ids(element_ids)

        # Exclude Hs if option set
        exclude_H_list = []
        for i, element in enumerate(element_ids):
            if element == 1:
                exclude_H_list.append(i)
        for i in reversed(exclude_H_list):
            del element_ids[i]
            del coordinates[i]

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(element_ids, radii_type=radii_type,
                              scale=radii_scale)

        # Setting up coordinate array and center array
        atom_coordinates = np.array(coordinates)

        # Get list of atoms as Atom objects
        atom_list = []
        for i, (element_id, radius, coord) in enumerate(zip(element_ids, radii,
                                             atom_coordinates), start=1):
            coord = np.array(coord)
            atom = Atom(element_id, radius, coord, i)
            atom_list.append(atom)

        # Prune sphere points which are within vdW radius of other atoms.
        tree = scipy.spatial.cKDTree(s.points, compact_nodes=False,
                                     balanced_tree=False)
        mask = np.zeros(len(s.points), dtype=bool)
        for atom in atom_list:
            if atom.radius + s.radius > np.linalg.norm(atom.coordinates):
                to_prune = tree.query_ball_point(atom.coordinates,
                                                 atom.radius)
                mask[to_prune] = True
        buried_points = s.points[mask,:]
        free_points = s.points[np.invert(mask),:]

        # Calculate buried_volume
        self.buried_volume = len(buried_points) / len(s.points)

        # Set variables for outside access and function access.
        self._atoms = atom_list
        self._atom_coordinates = atom_coordinates
        self._sphere = s
        self._buried_points = buried_points
        self._free_points = free_points

    def draw_3D(self, full_density=False):
        """Draw a 3D representation of the molecule with buried and free 
        points
        
        Args:
            full_density (bool): Requests drawing of all points (slow).
        """
        # Set up lists for drawing
        elements = []
        coordinates = []
        radii = []
        indices = []
        for atom in self._atoms:
            elements.append(atom.element_id)
            coordinates.append(atom.coordinates)
            radii.append(atom.radius)
            indices.append(atom.index)
        
        coordinates = np.vstack(coordinates)
        
        # Draw molecule scene
        scene = MoleculeScene(elements, coordinates, radii, indices)

        # Take out reasonable amount of points
        if full_density:
            step = 1
        else:
            step = math.ceil(1 / self._density / 25)  

        buried_points = np.array(self._buried_points)
        np.random.shuffle(buried_points)
        buried_points = buried_points[::step, :]

        free_points = np.array(self._free_points)
        np.random.shuffle(free_points)
        free_points = free_points[::step, :]

        # Add points
        scene.add_points(buried_points, color='#1f77b4')
        scene.add_points(free_points, color='#ff7f0e')

    def print_report(self):
        """Prints a report of the buried volume for use in shell scripts"""
        print("V_bur (%):", round(self.buried_volume * 100, 1))

    def plot_steric_map(self, z_axis_atoms, filename=None, levels=150, grid=100, all_positive=True, cmap="viridis"):
        """Plots a steric map as in the original article.

        Args:
            all_positive (bool)     :   Plot all positive values, even if they
                                        are outside the sphere.
            cmap (str)              :   Colormap for contour plot
            filename (str)          :   Name of file for saving the plot.
            grid (int)              :   Dimension of plotting grid
            levels (int)            :   Number of levels in the contour plot
            z_axis_atoms (list)     :   Atom indices for determining the
                                        orientation of the z axis (one-indexed)
        """
        # Set up coordinates
        atoms = self._atoms
        center = np.array(self._sphere.center)
        all_coordinates = np.array(self._all_coordinates)
        coordinates = np.array(self._atom_coordinates)

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
                thingy = atom.radius**2 - (x - x_s)**2 - (y - y_s)**2
                if thingy >= 0:
                    z_atom = math.sqrt(thingy) + z_s
                    z_list.append(z_atom)
            # Take point which is furthest along z axis
            if z_list:
                z_max = max(z_list)
                # Test if point is inside the sphere. Points with positive z
                # values are included by default anyway in accordance to article
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
        circle = plt.Circle((0,0), r, fill=False)
        ax.add_patch(circle)
        plt.xlabel("X (Å)")
        plt.ylabel("Y (Å)")
        cf.set_clim(-r, r)
        c_bar = fig.colorbar(cf)
        c_bar.set_label("Z(Å)")
        ax.set_aspect('equal', 'box')

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

class SASA:
    """Performs and stores results of solvent accessible surface area 
    calculations.

    Args:
        coordinates (list):     List of coordinates in Å.
        density (float):        Density of points on the vdW atomic
                                spheres in Å^-2
        element_ids (list):     List of element identifiers, either as symbols
                                or numbers.
        probe_radius (float):   Radius of probe atom in Å.
        radii (list):           List of radii (optional)
        radii_type (str):       Choice of vdW radii, either "crc" or "bondi"

    Attributes:
        area (float):           Area of the solvent accessible surface.
        atom_areas (dict):      Dictionary of atom areas (1-indexed)
        atom_volumes (dict):    Dictionary of atom volumes (1-indexed)
        volume (float):         Volume of the solvent accessible surface.
    """
    def __init__(self, element_ids, coordinates, radii=[], radii_type="crc", 
                 probe_radius=1.4, density=0.01):
        # Converting element ids to atomic numbers if the are symbols
        element_ids = convert_element_ids(element_ids)

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(element_ids, radii_type=radii_type)
        
        # Increment the radii with the probe radius
        orig_radii = np.array(radii)
        radii = orig_radii + probe_radius
        
        # Center coordinate system at geometric center
        coordinates = np.array(coordinates)
        center = np.mean(coordinates, axis=0)
        coordinates -= center
        
        # Construct list of atoms
        atoms = []
        orig_atoms = []
        for i, (coordinate, radius, orig_radius, element_id) in \
                enumerate(zip(coordinates, radii, orig_radii, element_ids), start=1):
            sasa_atom = SASAAtom(element_id, radius, coordinate, i)
            orig_atom = Atom(element_id, orig_radius, coordinate, i)
            atoms.append(sasa_atom)
            orig_atoms.append(orig_atom)
        
        # Determine occluded and accessible points of each atom based on
        # distances to all other atoms (brute force)
        for atom in atoms:
            # Construct sphere for atom
            sphere = Sphere(atom.coordinates, atom.radius, density=density)

            # Select atoms to test against
            test_atoms = [test_atom for test_atom in atoms 
                if test_atom is not atom and np.linalg.norm(atom.coordinates - test_atom.coordinates) < atom.radius + test_atom.radius]

            # Select coordinates and radii for other atoms
            test_coordinates = [test_atom.coordinates for test_atom in test_atoms]
            test_radii = [test_atom.radius for test_atom in test_atoms]
            test_radii = np.array(test_radii).reshape(-1, 1)

            # Get distances to other atoms and subtract radii
            distances = cdist(test_coordinates, sphere.points)
            distances -= test_radii

            # Take smallest distance and perform check
            min_distances = np.min(distances, axis=0)

            atom.occluded_points = sphere.points[min_distances <= 0]
            atom.accessible_points = sphere.points[min_distances > 0]

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
        self._orig_atoms = orig_atoms
        self._density = density

    def draw_3D(self, full_density=False, highlight=[]):
        """Draw a 3D representation of the molecule with the cone
        
        Args:
            full_density (bool) : Requests drawing of all points (slow)
            highlight (list)    : List of atom indices for highlighting surface
        """
        # Set up lists for drawing
        elements = []
        coordinates = []
        radii = []
        indices = []
        for atom in self._orig_atoms:
            elements.append(atom.element_id)
            coordinates.append(atom.coordinates)
            radii.append(atom.radius)
            indices.append(atom.index)
        
        coordinates = np.vstack(coordinates)
        
        # Draw molecule scene
        scene = MoleculeScene(elements, coordinates, radii, indices)

        # Take out a reasonable amount of points
        for atom in self._atoms:
            points = np.array(atom.accessible_points)
            if full_density:
                step = 1
            else:
                step = math.ceil(1 / self._density / 25)
            np.random.shuffle(points)
            points = points[::step, :]
        
            # Add points
            if atom.index in highlight:
                color = '#ff7f0e'
            else:
                color = '#1f77b4'
            if np.any(points):
                scene.add_points(points, color=color)
       
    def print_report(self, verbose=False):
        """Print report of results

        Args:
            verbose(bool):  Prints areas per atom if set to True.
        """
        print(f"Probe radius (Å): {self._probe_radius}")
        print(f"Solvent accessible surface area (Å^2): {self.area:.1f}")
        print("Volume inside solvent accessible surface (Å^3): "
              f"{self.volume:.1f}")
        if verbose:
            print(f"{'Atom':<10s}{'Index':<10s}{'Area (Å^2)':<10s}")
            for atom, (i, area) in zip(self._atoms, self.atom_areas.items()):
                symbol = atomic_symbols[atom.element_id]
                print(f"{symbol:<10s}{i:<10d}{area:<10.1f}")

class ConeAngle:
    """Calculates and stores the results of exact cone angle calculation as
    described in J. Comput. Chem. 2013, 34, 1189.

    Args:
        atom_1 (int)            :   Index of central atom (starting from 1)
        coordinates (list)      :   List of atom coordinates (Å)
        element_ids (list)      :   Elements as atomic numbers or symbols
        radii (list)            :   vdW radii (Å) (optional)
        radii_type (str)        :   Type of radii, "crc" or "bondi"

    Attributes:
        atoms (list)            :   List of atoms objects
        cone (object)           :   Smallest cone encompassing all atoms
        cone_angle (float)      :   Exact cone angle (degrees)
        element_ids (list)      :   Element ids as atomic numbers or symbols
        radii (list)            :   vdW radii (Å) (optional)
        tangent_atoms (list)    :   Atoms tangent to cone
    """
    def __init__(self, element_ids, coordinates, atom_1, radii=[], radii_type="crc"):
        # Converting element ids to atomic numbers if the are symbols
        element_ids = convert_element_ids(element_ids)

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(element_ids, radii_type=radii_type)

        # Check distances to atom_1 so that it is not within their vdW radii
        within = check_distances(element_ids, coordinates, atom_1, radii=radii)
        if within:
            atom_string = ' '.join([str(i) for i in within])
            raise Exception("Atoms within vdW radius of atom 1:", atom_string)

        # Setting up coordinate array and translate coordinates
        atom_coordinates = np.array(coordinates)
        atom_coordinates -= coordinates[atom_1 - 1]

        # Get list of atoms as Atom objects
        atoms = []
        for i, (coord, radius, element_id) in enumerate(zip(atom_coordinates, radii, element_ids), start=1):
            if i != atom_1:
                coord = np.array(coord)
                atom = ConeAngleAtom(coord, radius, i, element_id)
                atoms.append(atom)
                atom.get_cone()
        self._atoms = atoms

        # Search for cone over single atoms
        cone = self._search_one_cones()
        if cone:
            self.cone = cone
            self.cone_angle = math.degrees(cone.angle * 2)
            self.tangent_atoms = [atom for atom in cone.atoms]
        else:
            # Prune out atoms that lie in the shadow of another atom's cone
            loop_list = list(atoms)
            remove_set = set()
            for cone_atom in loop_list:
                for test_atom in loop_list:
                    if cone_atom != test_atom:
                        if cone_atom.cone.is_inside(test_atom):
                            remove_set.add(test_atom)
            for i in remove_set:
                loop_list.remove(i)
            self._loop_list = loop_list

        # Search for cone over pairs of atoms
        if not cone:
            cone = self._search_two_cones()
            if cone:
                self.cone = cone
                self.cone_angle = math.degrees(cone.angle * 2)
                self.tangent_atoms = [atom for atom in cone.atoms]

        # Search for cones over triples of atoms
        if not cone:
            cone = self._search_three_cones()
            if cone:
                self.cone = cone
                self.cone_angle = math.degrees(cone.angle * 2)
                self.tangent_atoms = [atom for atom in cone.atoms]
    
    def draw_3D(self):
        """Draw a 3D representation of the molecule with the cone"""
        # Set up lists for drawing
        elements = []
        coordinates = []
        radii = []
        indices = []
        for atom in self._atoms:
            elements.append(atom.element_id)
            coordinates.append(atom.coordinates)
            radii.append(atom.radius)
            indices.append(atom.index)
        
        coordinates = np.vstack(coordinates)
        
        # Draw molecule scene
        scene = MoleculeScene(elements, coordinates, radii, indices)

        # Determine direction and extension of cone
        if self.cone_angle > 180:
            normal = - self.cone.normal 
        else:
            normal = self.cone.normal
        projected = np.dot(normal, coordinates.T) + np.array(radii)

        max_extension = np.max(projected)
        if self.cone_angle > 180:
            max_extension += 1
        
        # Add cone
        scene.add_cone([0, 0, 0], normal, self.cone.angle, max_extension)

    def _search_one_cones(self):
        """Searches over cones tangent to one atom

        Returns:
            max_1_cone (object)     :   Largest cone tangent to one atom
        """
        # Get the largest cone
        atoms = self._atoms
        alphas = np.array([atom.cone.angle for atom in atoms])
        max_1_cone = atoms[np.argmax(alphas)].cone
        self._max_1_cone = max_1_cone

        # Check if all atoms are contained in cone. If yes, return cone,
        # otherwise, return None.
        in_list = []
        test_list = [atom for atom in atoms if atom not in max_1_cone.atoms]
        for atom in test_list:
            in_list.append(max_1_cone.is_inside(atom))
        if all(in_list):
            return max_1_cone
        else:
            return None

    def _search_two_cones(self):
        """Search over cones tangent to two atoms

        Returns:
            max_2_cone (object)     :   Largest cone tangent to two atoms
        """
        # Create two-atom cones
        loop_list = self._loop_list
        cone_list = []
        for atom_i, atom_j in itertools.combinations(loop_list, r=2):
            cone = self._get_two_atom_cone(atom_i, atom_j)
            cone_list.append(cone)

        # Select largest two-atom cone
        angles = np.array([cone.angle for cone in cone_list])
        max_2_cone = cone_list[np.argmax(angles)]
        self._max_2_cone = max_2_cone

        # Check if all atoms are contained in cone. If yes, return cone,
        # otherwise, return None
        in_list = []
        for atom in loop_list:
            in_list.append(max_2_cone.is_inside(atom))

        if all(in_list):
            return max_2_cone
        else:
            return None

    def _search_three_cones(self):
        """Search over cones tangent to three atoms

        Returns:
            min_3_cone (object)     :   Smallest cone tangent to three atoms
                                        encompassing all atoms
        """
        # Create three-atom cones
        loop_list = self._loop_list
        cone_list = []
        for atom_i, atom_j, atom_k in itertools.combinations(loop_list, r=3):
            cones = self._get_three_atom_cones(atom_i, atom_j, atom_k)
            cone_list.extend(cones)

        # Get upper and lower bound to apex angle
        upper_bound = self._get_upper_bound()
        lower_bound = self._max_2_cone.angle

        # Remove cones from consideration which are outside the bounds
        remove_list = []
        for cone in cone_list:
            if cone.angle - lower_bound < -1e-5 or upper_bound - cone.angle < -1e-5:
                remove_list.append(cone)

        for cone in reversed(remove_list):
            cone_list.remove(cone)

        # Keep cones that encompass all atoms
        keep_list = []
        for cone in cone_list:
            in_list = []
            for atom in loop_list:
                in_list.append(cone.is_inside(atom))
            if all(in_list):
                keep_list.append(cone)

        # Take the smallest cone that encompasses all atoms
        cone_angles = np.array([cone.angle for cone in keep_list])
        min_3_cone = keep_list[np.argmin(cone_angles)]

        return min_3_cone

    def _get_upper_bound(self):
        """Calculates upper bound for apex angle

        Returns:
            upper_bound (float)     :   Upper bound to apex angle in radians
        """
        loop_list = self._loop_list

        # Calculate unit vector to centroid
        coordinates = np.array([atom.coordinates for atom in self._atoms])
        centroid_vector = np.mean(coordinates, axis=0)
        centroid_unit_vector = centroid_vector / np.linalg.norm(centroid_vector)

        # Getting sums of angle to centroid and vertex angle.
        angle_sum_list = []
        for atom in self._atoms:
            cone = atom.cone
            cos_angle = np.dot(centroid_unit_vector, cone.normal)
            angle = math.acos(cos_angle)
            angle_sum = cone.angle + angle
            angle_sum_list.append(angle_sum)

        # Select upper bound as the maximum angle
        upper_bound = max(angle_sum_list)

        return upper_bound

    @staticmethod
    def _get_two_atom_cone(atom_i, atom_j):
        """Creates a cone tangent to two atoms

        Args:
            atom_i (object) :   First tangent atom object
            atom_j (object) :   Second tangent atom object

        Returns:
            cones (object)  :   Cone tangent to the two atoms
        """
        # Get the cone angle
        cone_i = atom_i.cone
        cone_j = atom_j.cone
        beta_i = cone_i.angle
        beta_j = cone_j.angle
        beta_ij = math.acos(np.dot(atom_i.cone.normal, atom_j.cone.normal))
        alpha_ij = (beta_ij + beta_i + beta_j) / 2

        # Get the cone normal
        a_ij = (1 / math.sin(beta_ij)) * math.sin(0.5 * (beta_ij + beta_i - beta_j))
        b_ij = (1 / math.sin(beta_ij)) * math.sin(0.5 * (beta_ij - beta_i + beta_j))
        c_ij = 0

        n = a_ij * cone_i.normal + b_ij * cone_j.normal + c_ij
        n = n / np.linalg.norm(n)

        # Create cone
        angle = alpha_ij
        normal = n
        cone = ConeAngleCone(angle, [atom_i, atom_j], normal)

        return cone

    @staticmethod
    def _get_three_atom_cones(atom_i, atom_j, atom_k):
        """Creates a list of cones tangent to three atoms

        Args:
            atom_i (object)     :   First tangent atom object
            atom_j (object)     :   Second tangent atom object
            atom_k (object)     :   Third tangent atom object

        Returns:
            cone_list (list)    :   List of cones tangent to the three atoms

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
        N = np.array([np.cross(m_j, m_k), np.cross(m_k, m_i), np.cross(m_i, m_j)]).T
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
        cos_roots = [math.acos(roots[0]), 2 * np.pi - math.acos(roots[0]), math.acos(roots[1]), 2 * np.pi - math.acos(roots[1])]

        # Test roots and keep only those that are physical
        angle_list = []
        test_list = []
        for root in cos_roots:
            alpha = root / 2
            test = A * math.cos(alpha)**2 + B * math.sin(alpha)**2 + 2 * C * math.sin(alpha) * math.cos(alpha)
            test_D = abs(test - D)
            angle_list.append(alpha)
            test_list.append(test_D)
        angles = np.array(angle_list)                
        tests = np.array(test_list)
        physical_angles = angles[np.argsort(tests)][:2]

        # Create cones for physical angles
        cone_list = []
        for alpha in physical_angles:
            # Calculate normal vector
            a_ij = (math.cos(alpha - beta_i) - math.cos(alpha - beta_j) * math.cos(beta_ij)) / math.sin(beta_ij)**2
            b_ij = (math.cos(alpha - beta_j) - math.cos(alpha - beta_i) * math.cos(beta_ij)) / math.sin(beta_ij)**2
            c_ij_squared = 1 - a_ij**2 - b_ij**2 - 2 * a_ij * b_ij * math.cos(beta_ij)
            # Set c_ij_squared to 0 if negative due to numerical precision.
            if c_ij_squared < 0:
                c_ij_squared = 0
            c_ij = math.sqrt(c_ij_squared)
            p = N @ (u * math.cos(alpha) + v * math.sin(alpha)).reshape(-1)
            sign = np.sign(gamma) * np.sign(np.dot(p, np.cross(m_i, m_j)))
            if np.sign(c_ij) != sign:
                c_ij = -c_ij
            n = a_ij * m_i + b_ij * m_j + c_ij * 1 / math.sin(beta_ij) * np.cross(m_i, m_j)

            # Create cone
            cone = ConeAngleCone(alpha, [atom_i, atom_j, atom_k], n)
            cone_list.append(cone)

        return cone_list

    def print_report(self):
        """Prints report of results"""
        tangent_list = [f'{atomic_symbols[atom.element_id]}{atom.index}' \
                        for atom in self.tangent_atoms]
        tangent_string = ' '.join(tangent_list)
        print(f"Cone angle: {self.cone_angle:.1f}")
        print(f"No. tangent atoms: {len(self.tangent_atoms)}")
        print(f"Tangent to: {tangent_string}")