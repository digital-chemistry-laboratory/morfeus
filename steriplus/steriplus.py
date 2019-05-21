"""Classes for performing calculations of steric descriptors of molecules.

Classes:
    BuriedVolume: Calculates buried volumes
    ConeAngle: Calculates exact cone angles.
    SASA: Calculates solvent accessible surface area.
    Sterimol: Calculates Sterimol parameters
"""
import math
import itertools

import matplotlib.pyplot as plt
import numpy as np

import scipy.spatial
from scipy.spatial.distance import cdist

from steriplus.data import atomic_symbols
from steriplus.geometry import Atom, Cone, rotate_coordinates, Sphere
from steriplus.helpers import check_distances, convert_elements, get_radii
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
        B_1 (ndarray): Sterimol B_1 vector (Å)
        B_1_value (float): Sterimol B_1 value (Å)
        B_5 (ndarray): Sterimol B_5 vector (Å)
        B_5_value (float): Sterimol B_5 value (Å)
        bond_length (float): Bond length between atom 1 and atom 2 (Å)
        L (ndarray): Sterimol L vector (Å)
        L_value (float): Sterimol L value (Å)
        L_value_uncorrected (float): Sterimol L value minus 0.40 Å
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
        
        self._atom_1 = atom_1
        self._atom_2 = atom_2
        
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
        L_start = coordinates[self._atom_1 - 1]
        L_stop = self.L
        L_length = self.L_value 
        scene.add_arrow(L_start, L_stop, L_length, "L")

        # Draw B_1 vector
        B_1_start = coordinates[self._atom_2 - 1]
        B_1_stop = self.B_1
        B_1_length = self.B_1_value
        scene.add_arrow(B_1_start, B_1_stop, B_1_length, "B_1")

        # Draw B_5 vector
        B_5_start = coordinates[self._atom_2 - 1]
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
    """Performs and stores the results of a buried volume calculation as
    described in Organometallics 2016, 35, 2286.

    Args:
        atom_1 (int): Atom index of metal (starting from 1)
        coordinates (list): Coordinates (Å)
        density (float): Volume per point (Å**3) in the sphere
        elements (list): Elements as atomic symbols or numbers
        exclude_list (list): Indices of atoms to exclude from the calculation
                             (starting from 1)
        include_hs (bool): Whether to include H atoms in the calculation
        radii (list): vdW radii (Å)
        radii_scale (float): Scaling factor for radii. 1.17 from original paper.
        radii_type (str): Type of radii to use: 'bondi' (default) or 'crc'
        radius (float): Radius of sphere (Å). 3.5 from orginal paper.

    Parameters:
        buried_volume (float): Buried volume
    """
    def __init__(self, elements, coordinates, atom_1, exclude_list=[],
                 radii=[], include_hs=False, radius=3.5, radii_type="bondi",
                 radii_scale=1.17, density=0.001):
        # Get the coordinates for the central atom
        center = np.array(coordinates[atom_1 - 1])

        # Construct sphere at metal center
        sphere = Sphere(center, radius, method="projection", density=density,
                   filled=True)

        # Save density and coordinates for steric map plotting.
        self._density = density
        self._all_coordinates = np.array(coordinates)

        # Converting element ids to atomic numbers if the are symbols
        elements = convert_elements(elements)

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(elements, radii_type=radii_type,
                              scale=radii_scale)

        # Get list of atoms as Atom objects
        atoms= []
        for i, (element, radius, coord) in enumerate(zip(elements, radii,
                                               coordinates), start=1):
            if i not in exclude_list and element != 1:
                atom = Atom(element, coord, radius, i)
                atoms.append(atom)

        # Prune sphere points which are within vdW radius of other atoms.
        tree = scipy.spatial.cKDTree(sphere.points, compact_nodes=False,
                                     balanced_tree=False)
        mask = np.zeros(len(sphere.points), dtype=bool)
        for atom in atoms:
            if atom.radius + sphere.radius > np.linalg.norm(atom.coordinates):
                to_prune = tree.query_ball_point(atom.coordinates,
                                                 atom.radius)
                mask[to_prune] = True
        buried_points = sphere.points[mask,:]
        free_points = sphere.points[np.invert(mask),:]

        # Calculate buried_volume
        self.buried_volume = len(buried_points) / len(sphere.points)

        # Set variables for outside access and function access.
        self._atoms = atoms
        self._sphere = sphere
        self._buried_points = buried_points
        self._free_points = free_points

    def draw_3D(self, full_density=False):
        """Draw a 3D representation of the molecule with buried and free 
        points
        
        Args:
            full_density (bool): Requests drawing of all points (slow).
        """
        # Set up lists for drawing
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]
        radii = [atom.radius for atom in self._atoms]
        indices = [atom.index for atom in self._atoms]
        
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
        
        # Center coordinate system at geometric center
        coordinates = np.array(coordinates)
        center = np.mean(coordinates, axis=0)
        coordinates -= center
        
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
            test_atoms = []
            for test_atom in atoms:
                if test_atom is not atom:
                    distance = np.linalg.norm(atom.coordinates 
                                              - test_atom.coordinates)
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

                atom.occluded_points = sphere.points[min_distances <= 0]
                atom.accessible_points = sphere.points[min_distances > 0]
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

    def draw_3D(self, full_density=False, highlight=[]):
        """Draw a 3D representation of the molecule with the solvent accessible
        surface areas as dots.
        
        Args:
            full_density (bool): Requests drawing of all points (slow)
            highlight (list): Atom indices for highlighting surface
        """
        # Set up lists for drawing
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]
        radii = [atom.radius for atom in self._atoms]
        indices = [atom.index for atom in self._atoms]
       
        coordinates = np.vstack(coordinates)
        radii = np.array(radii) - self._probe_radius
        
        # Draw molecule scene
        scene = MoleculeScene(elements, coordinates, radii, indices)

        # Take out a reasonable amount of points
        for atom in self._atoms:
            points = atom.accessible_points
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

    def draw_3D(self):
        """Draw a 3D representation of the molecule with the cone"""
        # Set up lists for drawing
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]
        radii = [atom.radius for atom in self._atoms]
        indices = [atom.index for atom in self._atoms]
        
        coordinates = np.vstack(coordinates)
        
        # Draw molecule scene
        scene = MoleculeScene(elements, coordinates, radii, indices)

        # Determine direction and extension of cone
        if self.cone_angle > 180:
            normal = - self._cone.normal 
        else:
            normal = self._cone.normal
        projected = np.dot(normal, coordinates.T) + np.array(radii)

        max_extension = np.max(projected)
        if self.cone_angle > 180:
            max_extension += 1
        
        # Add cone
        scene.add_cone([0, 0, 0], normal, self._cone.angle, max_extension)

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

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"