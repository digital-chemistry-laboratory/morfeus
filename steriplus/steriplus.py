"""Classes for performing calculations of steric descriptors of molecules.

Classes:
    BuriedVolume: Calculates buried volumes.
    ConeAngle: Calculates exact cone angles.
    LocalForce: Calculates local force constants.
    SASA: Calculates solvent accessible surface area.
    Sterimol: Calculates Sterimol parameters.
"""
import math
import itertools

# Matplotlib is required for plotting steric maps for buried volumes
try:
    import matplotlib.pyplot as plt
    from matplotlib.colors import hex2color
except ImportError:
    _has_matplotlib = False
else:
    _has_matplotlib = True
_warning_matplotlib = "Install matplotlib to use this function."

import numpy as np
from subprocess import Popen, DEVNULL, PIPE

# VTK and PyVista are required for 3D visualization and for generating surfaces
# for dispersion calculations.
try:
    import vtk
    import pyvista as pv
    import pymeshfix
    from steriplus.plotting import Arrow_3D, Cone_3D
except ImportError:
    _has_vtk = False
else:
    _has_vtk = True
_warning_vtk = "Install pyvista and vtk to use this function."

import scipy.spatial
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist, euclidean

from steriplus.data import atomic_symbols, au_to_kcal, angstrom_to_bohr
from steriplus.data import jmol_colors
from steriplus.geometry import Atom, Cone, rotate_coordinates, Sphere
from steriplus.geometry import kabsch_rotation_matrix
from steriplus.helpers import check_distances, convert_elements, get_radii
from steriplus.helpers import D3Calculator, conditional
from steriplus.io import CubeParser, D3Parser, D4Parser, VertexParser

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
        elements = convert_elements(elements)

        # Get radii if they are not supplied
        if not radii:
            radii = get_radii(elements, radii_type=radii_type)
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
        L = unit_vector * L_value

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

    @conditional(_has_vtk, _warning_vtk)
    @conditional(_has_matplotlib, _warning_matplotlib) 
    def draw_3D(self, atom_scale=0.5, background_color="white",
                arrow_color="steelblue"):
        """Draw a 3D representation of the molecule with the Sterimol vectors.
        
        Args:
            atom_scale (float): Scaling factor for atom size
            background_color (str): Background color for plot
            arrow_color (str): Arrow color
        """
        # Set up plotter
        p = pv.BackgroundPlotter()
        p.set_background(background_color)

        # Set up lists for drawing
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]
        radii = [atom.radius for atom in self._atoms]
        colors = [hex2color(jmol_colors[i]) for i in elements]

        # Draw molecule
        for coordinate, radius, color in zip(coordinates, radii, colors):
            sphere = pv.Sphere(center=list(coordinate),
                               radius=radius * atom_scale)
            p.add_mesh(sphere, color=color, opacity=1)
        
        # Get arrow starting points
        start_L = self._atoms[self._atom_1 - 1].coordinates
        start_B = self._atoms[self._atom_2 - 1].coordinates

        # Add L arrow with label 
        length = np.linalg.norm(self.L)
        direction = self.L / length
        stop_L = start_L + length * direction 
        L_arrow = Arrow_3D(start=start_L, direction=direction, length=length)
        p.add_mesh(L_arrow, color=arrow_color)

        # Add B_1 arrow
        length = np.linalg.norm(self.B_1)
        direction = self.B_1 / length
        stop_B_1 = start_B + length * direction 
        B_1_arrow = Arrow_3D(start=start_B, direction=direction, length=length)
        p.add_mesh(B_1_arrow, color=arrow_color)

        # Add B_5 arrow
        length = np.linalg.norm(self.B_5)
        direction = self.B_5 / length
        stop_B_5 = start_B + length * direction 
        B_5_arrow = Arrow_3D(start=start_B, direction=direction, length=length)
        p.add_mesh(B_5_arrow, color=arrow_color)

        # Add labels
        points = np.vstack([stop_L, stop_B_1, stop_B_5])
        labels = ["L", "B1", "B5"]
        p.add_point_labels(points, labels, text_color="black", font_size=30,
                           bold=False, show_points=False, point_size=1)

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
    
    @conditional(_has_matplotlib, _warning_matplotlib)
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

    @conditional(_has_vtk, _warning_vtk)
    @conditional(_has_matplotlib, _warning_matplotlib) 
    def draw_3D(self, atom_scale=1, background_color="white",
                buried_color="tomato", free_color="steelblue", opacity=0.05,
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
        p = pv.BackgroundPlotter()
        p.set_background(background_color)

        # Set up lists for drawing
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]
        radii = [atom.radius for atom in self._atoms]
        colors = [hex2color(jmol_colors[i]) for i in elements]

        # Draw molecule
        for coordinate, radius, color in zip(coordinates, radii, colors):
            sphere = pv.Sphere(center=list(coordinate),
                               radius=radius * atom_scale)
            p.add_mesh(sphere, color=color, opacity=1)
        
        # Add buried points
        p.add_points(self._buried_points, color=buried_color, opacity=opacity)

        # Add free points
        p.add_points(self._free_points, color=free_color, opacity=opacity,
                     size=size)

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

    @conditional(_has_vtk, _warning_vtk)
    @conditional(_has_matplotlib, _warning_matplotlib) 
    def draw_3D(self, atom_scale=1, background_color="white",
                point_color="steelblue", opacity=1, size=1):
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
        p = pv.BackgroundPlotter()
        p.set_background(background_color)

        # Set up lists for drawing
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]
        radii = np.array([atom.radius for atom in self._atoms]) - \
            self._probe_radius
        colors = [hex2color(jmol_colors[i]) for i in elements]

        # Draw molecule
        for coordinate, radius, color in zip(coordinates, radii, colors):
            sphere = pv.Sphere(center=list(coordinate),
                               radius=radius * atom_scale)
            p.add_mesh(sphere, color=color, opacity=1)
        
        # Draw surface points
        surface_points = np.vstack([atom.accessible_points 
                                    for atom in self._atoms])
        p.add_points(surface_points, color=point_color, opacity=opacity,
                     size=size)
        
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

    @conditional(_has_vtk, _warning_vtk)
    @conditional(_has_matplotlib, _warning_matplotlib) 
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
        p = pv.BackgroundPlotter()
        p.set_background(background_color)

        # Set up lists for drawing
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]
        radii = [atom.radius for atom in self._atoms]
        colors = [hex2color(jmol_colors[i]) for i in elements]

        # Draw molecule
        for coordinate, radius, color in zip(coordinates, radii, colors):
            sphere = pv.Sphere(center=list(coordinate),
                               radius=radius * atom_scale)
            p.add_mesh(sphere, color=color, opacity=1)
        
        # Determine direction and extension of cone
        cone_angle = math.degrees(self._cone.angle)
        coordinates = np.array(coordinates)
        if cone_angle > 180:
            normal = - self._cone.normal 
        else:
            normal = self._cone.normal
        projected = np.dot(normal, coordinates.T) + np.array(radii)

        max_extension = np.max(projected)
        if cone_angle > 180:
            max_extension += 1
        
        # Make the cone
        cone = Cone_3D(center=[0, 0, 0] + (max_extension * normal) / 2,
                       direction=-normal, angle=cone_angle, 
                       height=max_extension, capping=False, resolution=100)
        p.add_mesh(cone, opacity=cone_opacity, color=cone_color)

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"

class Dispersion:
    """Calculates and stores the results for the P_int dispersion descriptor.

    The descriptor is defined in Angew. Chemie Int. Ed. 2019.
    DOI: 10.1002/anie.201905439. Steriplus can calculate it based on a surface
    either from vdW radii, surface vertices or electron density. Dispersion can
    be obtained with the D3 or D4 model.

    Args:
        calculate_coefficents (bool): Whether to calculate D3 coefficients with
            internal code.
        coordinates (list): Coordinates (Å)
        density (float): Area per point (Å**2) on the vdW surface.
        elements (list): Elements as atomic symbols or numbers
        excluded_atoms (list): Atoms to exclude from the calculation. Used only
            for calculation of substituent P_ints.
        point_surface (bool): Use point surface from vdW radii.
        radii (list): VdW radii (Å)
        radii_type (str): Choice of vdW radii: 'bondi', 'crc' or 'rahm'
            (default)

    Parameters:
        area (float): Area of surface (Å^2)
        atom_areas (dict): Atom areas (Å^2, starting from 1)
        atom_p_ints (dict): P_int value for atoms (kcal^(1/2) Bohr^(-1/2)
            starting from 1)
        p_int (float): P_int value for molecule (kcal^(1/2) Bohr^(-1/2)
        p_max (float): Mean of 10 highest P values (kcal^(1/2) Bohr^(-1/2)
        p_min (float): Mean of 10 lowest P values (kcal^(1/2) Bohr^(-1/2)
        p_values (list): All P values (kcal^(1/2) Bohr^(-1/2)
        volume (float): Volume of surface (Å^3)
    """
    def __init__(self, elements, coordinates, radii=[], radii_type="rahm",
                 point_surface=True, calculate_coefficients=True, density=0.1,
                 excluded_atoms=[]):
        # Set up
        self._surface = None
        self._excluded_atoms = excluded_atoms
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
        if calculate_coefficients:
            self.get_coefficients(model='id3')
        
        # Calculatte P_int values
        if point_surface and calculate_coefficients:
            self.calculate_p_int()

    @conditional(_has_vtk, _warning_vtk)
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
    
    @conditional(_has_vtk, _warning_vtk)
    def surface_from_multiwfn(self, filename):
        """Adds surface from Multiwfn vertex file with connectivity information.

        Args:
            filename (str): Name of vertex file
        """
        # Read the vertices and faces from the Multiwfn output file
        parser = VertexParser(filename)
        vertices = np.array(parser.vertices)
        faces = np.array(parser.faces)
        faces = np.insert(faces, 0, values=3, axis=1)
        
        # Construct surface and fix it with pymeshfix
        surface = pv.PolyData(vertices, faces, show_edges=True)
        meshfix = pymeshfix.MeshFix(surface)
        meshfix.repair()
        surface = meshfix.mesh

        # Process surface
        self._surface = surface
        self._process_surface()

    @conditional(_has_vtk, _warning_vtk)
    def surface_from_radii(self, step_size=0.313, smoothing=76,
                           radii_scale=1.138, method="flying_edges"):
        """Construct surface by smoothening vdW surface from atomic radii.
        Method described by Tom Goddard.
        https://www.cgl.ucsf.edu/chimera/data/surface-oct2013/surface.html        

        Args:
            step_size (float): Step size of sampling grid for surface
                               construction (Å)
            smoothing (int): Iterations of VTK smoothing filter
            radii_scale (float): Scaling factor for vdW radii
            method (str): Method for contouring: 'contour' or 'flying_edges
                          (default)
        """
        # Set up coordinates and radii
        coordinates = np.array([atom.coordinates for atom in self._atoms])
        radii = np.array([atom.radius for atom in self._atoms]) * radii_scale
        radii_x3 = np.tile(radii.reshape(-1, 1), (1, 3))

        # Extract boundaries of grid
        min_x, min_y, min_z = np.min(coordinates - radii_x3, axis=0) - step_size * 5
        max_x, max_y, max_z = np.max(coordinates + radii_x3, axis=0) + step_size * 5

        # Construct grid and extract points
        x = np.arange(min_x, max_x + step_size, step_size)
        y = np.arange(min_y, max_y + step_size, step_size)
        z = np.arange(min_z, max_z + step_size, step_size)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        points = np.stack((X.flatten(order="F"), Y.flatten(order="F"),
                           Z.flatten(order="F")), axis=1)

        # Take out distances to vdW surface from to each grid point
        dists = cdist(points, coordinates) - radii
        min_dists = np.min(dists, axis=1)
        min_dists = min_dists

        # Construct grid and add points
        grid = pv.UniformGrid(X, Y, Z)
        grid.dimensions = np.array((x.shape, y.shape, z.shape))
        grid.origin = (min_x, min_y, min_z)
        grid.spacing = (step_size, step_size, step_size)
        grid.point_arrays["values"] = min_dists

        # Countour the surface
        surface = self._contour_surface(grid, method=method, isodensity=0)
        self._surface = surface.smooth(smoothing)
        self._process_surface()
    
    @conditional(_has_vtk, _warning_vtk)
    def _process_surface(self):
        """Extracts face center points and assigns these to atoms based on
        proximity
        """
         # Get the area and volume
        self.area = self._surface.area
        self.volume = self._surface.volume

        # Assign face centers to atoms according to Voronoi partitioning
        coordinates = np.array([atom.coordinates for atom in self._atoms])
        points = self._surface.cell_centers().points
        kd_tree = cKDTree(coordinates)
        _, point_regions = kd_tree.query(points, k=1)
        point_regions = point_regions + 1
        
        # Compute faces areas 
        area_data = self._surface.compute_cell_sizes()
        areas = area_data.cell_arrays["Area"]
        
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
    
    @conditional(_has_vtk, _warning_vtk)
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

    def calculate_p_int(self, points=[]):
        """Calculate P_int values for surface or points.

        Args:
            points (list): Points to calculate P values for
        """
        # Set up array of points
        points = np.array(points)

        # Set up atoms and coefficients that are part of the calculation
        atom_indices = [atom.index - 1 for atom in self._atoms 
                        if atom.index not in self._excluded_atoms]
        coordinates = np.array([atom.coordinates for atom in self._atoms])
        coordinates = coordinates[atom_indices]
        c6_coefficients = np.array(self._c6_coefficients)
        c6_coefficients = c6_coefficients[atom_indices] * au_to_kcal
        c8_coefficients = np.array(self._c8_coefficients)
        c8_coefficients = c8_coefficients[atom_indices] * au_to_kcal

        # Take surface points if none are given
        if points.size == 0:
            points = np.vstack([atom.accessible_points for atom in self._atoms 
                                if atom.index not in self._excluded_atoms and
                                atom.accessible_points.size > 0])
            atomic = True

        # Calculate p_int for each point
        dist = cdist(points, coordinates) * angstrom_to_bohr 
        p = np.sum(np.sqrt(c6_coefficients/(dist**6)), axis=1) \
            + np.sum(np.sqrt(c8_coefficients/(dist**8)), axis=1)

        self.p_values = p
    
        # Take out atomic p_ints if no points are given
        if atomic:
            atom_p_ints = {}
            i_start = 0
            for atom in self._atoms:
                if atom.index not in self._excluded_atoms:
                    n_points = len(atom.accessible_points)
                    if n_points > 0:
                        i_stop = i_start + n_points
                        atom_ps = p[i_start:i_stop]
                        atom.p_values = atom_ps
                        atom_p_ints[atom.index] = np.sum(atom_ps * 
                            atom.point_areas / atom.area)
                        i_start = i_stop
                    else:
                        atom_p_ints[atom.index] = 0
                        atom.p_values = np.array([])
            self.atom_p_ints = atom_p_ints

        self.p_int = np.sum(p * self._point_areas / self.area)

        # Calculate p_min and p_max with slight modification to Robert's 
        # definitions
        p_sorted = np.sort(p)
        #self.p_min = np.median(p_sorted[:100]) # Robert's definition
        self.p_min = np.mean(p_sorted[:10])
        self.p_max = np.mean(p_sorted[-10:])

        # Map p_values onto surface
        if self._surface:
            mapped_p = np.zeros(len(p))
            for atom in self._atoms:
                if atom.index not in self._excluded_atoms:
                    mapped_p[self._point_map == atom.index] = atom.p_values
            self._surface.cell_arrays['p_int'] = mapped_p

    def get_coefficients(self, filename=None, model='id3'):
        """Get the C6 and C8 coefficients.

        The default model is the internal D3 calculator. Output can be read from
        the dftd3 and dftd4 programs by giving a filename in combination with
        the corresponding 'model' keyword argument.

        Args:
            filename (str): Output file from the dftd3 or dftd4 programs
            model (str): Calculation model: 'id3' (default), 'd3' or 'd4'.
        """
        if not filename and model =="id3":
            # Set up atoms and coordinates
            elements = [atom.element for atom in self._atoms]
            coordinates = [atom.coordinates for atom in self._atoms]

            # Calculate the D3 values
            calc = D3Calculator(elements, coordinates)
            self._c6_coefficients = calc.c6_coefficients
            self._c8_coefficients = calc.c8_coefficients
        elif filename and model == "d3":
            # Read the data
            parser = D3Parser(filename)
            self._c6_coefficients = parser.c6_coefficients
            self._c8_coefficients = parser.c8_coefficients
        elif filename and model == "d4":
            # Read the data
            parser = D4Parser(filename)
            self._c6_coefficients = parser.c6_coefficients
            self._c8_coefficients = parser.c8_coefficients            

    def print_report(self, verbose=False):
        """Print report of results

        Args:
            verbose (bool): Print atom P_ints
        """
        print(f"Surface area (Å^2): {self.area:.1f}")
        print(f"Surface volume (Å^3): {self.volume:.1f}")
        print(f"P_int (kcal^(1/2) Bohr^(-1/2): {self.p_int:.1f}")
        if verbose:
            print(f"{'Symbol':<10s}{'Index':<10s}{'P_int (kcal^(1/2) Bohr^(-1/2))':<30s}")
            for atom, (i, p_int) in zip(self._atoms, self.atom_p_ints.items()):
                symbol = atomic_symbols[atom.element]
                print(f"{symbol:<10s}{i:<10d}{p_int:<10.1f}")        

    @conditional(_has_vtk, _warning_vtk)
    @conditional(_has_matplotlib, _warning_matplotlib)
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
        p = pv.BackgroundPlotter()

        # Set up lists for drawing
        elements = [atom.element for atom in self._atoms]
        coordinates = [atom.coordinates for atom in self._atoms]
        radii = [atom.radius for atom in self._atoms]
        colors = [hex2color(jmol_colors[i]) for i in elements]

        # Draw molecule
        for coordinate, radius, color in zip(coordinates, radii, colors):
            sphere = pv.Sphere(center=coordinate, radius=radius * atom_scale)
            p.add_mesh(sphere, color=color, opacity=molecule_opacity)
        
        # Set up plotting of mapped surface
        if display_p_int == True:
            if self._surface:
                self._surface.set_active_scalar('p_int')
            color = None
            cmap = "coolwarm"
        else:
            color = "tan"
            cmap = None

        # Draw surface
        if self._surface:
            p.add_mesh(self._surface, opacity=opacity, color=color, cmap=cmap)

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"

class LocalForce:
    """Calculates and stores the results from local force constant 
    calculations.
    
    The method is described by Cremer in Int. J. Quantum Chem. 1998, 67, 1.
    Alternatively, the compliance matrix method can be used according to 
    J. Chem. Phys. 2010, 132, 184101.

    Args:
        log_file (str): Name of Gaussian log file
        fchk_file (str): Name of Gaussian fchk file
        pes_file (stry): Name of Gaussian PES file
        project_imag (bool): Whether to project out imaginary frequencies
            when using the local modes method
        cutoff (float): Cutoff for low force constant when using the local
            modes method
        method (str): Method for calculating the local force constants:
             'compliance' or 'local' (default)

    Attributes:
        internal_names (list): Gaussian names of internal indices
        internal_indices (dict): Mapping of internal coordinates to indices
            for force constants and frequencies
        local_force_constants (ndarray): Local mode force constants
            (mDyne/Å)
        local_frequencies (ndarray): Local mode frequencies (cm^(-1))
        n_imag (int): Number of normal modes with imaginary frequencies
    """
    def __init__(self, log_file, fchk_file=None, pes_file=None, 
        project_imag=True, cutoff=1e-3, method="local"):
        # Set up attributes
        self._log_file = log_file
        self.n_imag = None
        self.local_force_constants = np.array([])
        self.local_frequencies = np.array([])
        self.internal_indices = {}
        self.internal_names = []
        self._project_imag = project_imag
        self._cutoff = cutoff
        self._rotation_matrix = None

        # Read in data from Gaussian log file
        self._parse_log_file(log_file)

        # Read in high-precision normal modes force constants from fchk file
        if fchk_file:
            self._parse_fchk_file(fchk_file)
        
        # Read in high-precision force constants from PES file.
        if pes_file:
            self._parse_pes_file(pes_file)
        
        # Compute local mode properties from local modes or compliance matrix
        if method == "local":
            self._compute_local()
        if method == "compliance":
            self._compute_compliance()

    def get_local_force_constant(self, atoms):
        """Return the local force constant between a set of atoms.
        
        Args:
            atoms (list): Atoms in the internal coordinate
        
        Returns:
            force_constant (float): Local force constant (mDyne/Å)
        """
        atoms = frozenset(atoms)
        index = self.internal_indices.get(atoms)
        if index == None:
            raise Exception(f"No bond internal coordinate with these atoms.")
        force_constant = self.local_force_constants[index]
        
        return force_constant

    def get_local_frequency(self, atoms):
        """Return the local frequency between a set of atoms.
        
        Args:
            atoms (list): Atoms in the internal coordinate
        
        Returns:
            frequency (float): Local frequency (cm^-1)
        """
        atoms = frozenset(atoms)
        index = self.internal_indices.get(atoms)
        if index == None:
            raise Exception(f"No bond internal coordinate with these atoms.")
        frequency = self.local_frequencies[index]
        
        return frequency

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
        
        print(f"{'Atom_1':>10s}{'Atom_2':>10s}{'Atom_3':>10s}{'Atom_4':>10s}"
            f"{'Force constant' + '(' + unit + ')':>50s}"
            f"{'Frequency (cm^-1)':>40s}")
        
        # Print results for each internal
        for atoms, index in self.internal_indices.items():
            # Check if internal coordinate is angle or dihedral
            if len(atoms) == 3 and not angles:
                continue
            if len(atoms) == 4 and not dihedrals:
                continue
            force_constant = self.local_force_constants[index]
            frequency = self.local_frequencies[index]

            # Convert units for angles and dihedrals
            if len(atoms) > 2 and angle_units:
                force_constant *= 0.52918**2
            
            # Print out the results
            atom_iter = iter(sorted(atoms))
            atom_string = ""
            for i in range(4):
                atom = next(atom_iter, "")
                atom_string += f"{atom:>10}"
            print(atom_string + f"{force_constant:50.3f}" \
                + f"{frequency:40.0f}")

    def _add_cutoff(self):
        # Cut off low force constants at given value
        indices_small = np.where(np.array(self._force_constants)
            < self._cutoff)[0]
        if len(indices_small) > 0:
            self._force_constants[indices_small] \
                = np.array(self._cutoff).reshape(-1)

    def _compute_compliance(self):
        # Compute force constant matrix and get local force constants
        C = self._B @ np.linalg.pinv(self._fc_matrix) @ self._B.T
        k_s = 1 / np.diag(C) * 15.569141 #TODO update constants
        self.local_force_constants = k_s
        
    def _compute_local(self):
        # Compute D matrix from normal modes and B matrix
        if self._B.size > 0  and self._normal_modes.size > 0:
            self._D = self._B @ self._normal_modes.T
        elif len(self._D) < 1:
            raise Exception("Need both B matrix and normal modes,"
                            " or read local modes.")
        
        # Project out imaginary modes
        if self._project_imag:
            self._prune_imaginary()
        
        # Add cutoff to modes with low force constants
        if self._cutoff is not None:
            self._add_cutoff()
                
        # Compute local mode force constants
        K = np.diag(self._force_constants)
        
        k_s = []
        for row in self._D:
            k = 1 / (row @ np.linalg.inv(K) @ np.conj(row).T)
            k_s.append(k)
        k_s = np.array(k_s)

        # Scale force constants due to projection of imaginary normal modes
        if self._project_imag and self.n_imag:
            lengths = np.linalg.norm(self._D, axis=1)
            lengths_full = np.linalg.norm(self._D_full, axis=1)
            k_s *= lengths / lengths_full
        
        # Compute local mode frequencies TODO add these constants somewhere
        G = self._D @ self._D.T
        frequencies = np.sqrt((np.diag(G) * (1.66053904e-27)**-1 \
            * (k_s * (1e-10)**-1 * 1e-8)) / (4 * np.pi**2 * 299792458**2)) / 100
        
        self.local_force_constants = k_s
        self.local_frequencies = frequencies

    def _prune_imaginary(self):
        # Save full D matrix and force constant vector
        self._D_full = self._D
        self._force_constants_full = self._force_constants

        # Remove force constants with imaginary force constants
        self._force_constants = self._force_constants[self.n_imag:]
        self._D = self._D[:, self.n_imag:]

    def _parse_fchk_file(self, file):
        # Read fchk file
        lines = open(file).readlines()

        # Set up read flags
        read_modes = False
        read_hessian = False
        read_vib_e2 = False
        read_ic = False

        # Set up containers for reading data
        modes = []
        hessian = []
        vib_e2 = []
        internal_coordinates = []
        n_atoms = None

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
            # Read number of atoms
            if "Number of atoms" in line:
                n_atoms = int(line.strip().split()[4])
            # Read number of normal modes
            if "Number of Normal Modes" in line:
                n_modes = int(line.strip().split()[5])
            # Read number of internal coordinates
            if "Redundant internal coordinates" in line:
                n_redundant = int(line.strip().split()[5])        
            if "Vib-Modes" in line:
                read_modes = True
            if "Cartesian Force Constants" in line:
                read_hessian = True
            if "Vib-E2 " in line:
                read_vib_e2 = True
            if "Redundant internal coordinate indices" in line:
                read_ic = True

        # Take out normal mode force constants
        force_constants = np.array(vib_e2[n_modes * 2:n_modes * 3])

        # Construct force constant matrix from lower triangle       
        fc_matrix = np.zeros((n_atoms * 3, n_atoms * 3))
        fc_matrix[np.tril_indices_from(fc_matrix)] = hessian
        fc_matrix = np.triu(fc_matrix.T, 1) + fc_matrix

        # Take out the internal coordinates
        internal_coordinates = np.array(internal_coordinates)
        internal_indices = {}
        for i, coordinate in enumerate(np.split(internal_coordinates,
                                                n_redundant)):
            pruned_coordinate = [i for i in coordinate if i != 0]
            internal_indices[frozenset(pruned_coordinate)] = i
        
        # Set up attributes
        self._fc_matrix = fc_matrix
        self._normal_modes = np.array(modes).reshape(n_modes, n_atoms * 3)
        self._force_constants = force_constants
        self.internal_indices = internal_indices

    def _parse_log_file(self, file):
        # Read the log file
        lines = open(file).readlines()
        
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
        read_atomic_masses = False

        # Set up containers for reading data
        B_atom_map = {}
        B_vectors = {}
        normal_modes = []
        internal_modes = []
        force_constants = []
        atomic_masses = []
        fc_matrix = np.array([])
        ifc_matrix = np.array([])
        input_coordinates = []
        standard_coordinates = []
        n_imag = None
        n_atoms = None
        internal_indices = {}

        # Parse through log file content
        for line in lines:
            # Read in internal coordinate definitions
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
            # Read in Cartesian force constant matrix
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
            # Read in internal force constant matrix 
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
            # Read in atoms for creation of B matrix       
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
            # Read in values of B matrix
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
            if read_atomic_masses:
                if "Molecular mass: " in line:
                    read_atomic_masses = False
                elif "and mass" in line:
                    atomic_masses.append(float(line.strip().split()[8]))
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
            if "Name  Definition              Value          Derivative Info." in line:
                read_internal = True
                counter = 1
                internal_names = []
            if "- Thermochemistry -" in line:
                read_atomic_masses = True
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
        
        # Construct the B matrix from atoms and vectors
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
        
        # Detect whether the input coordinates have been rotated. If so, rotate
        # B matrix and its inverse.
        input_coordinates = np.array(input_coordinates).reshape(-1, 3)
        standard_coordinates = np.array(standard_coordinates).reshape(-1, 3)
        
        if not np.array_equal(input_coordinates, standard_coordinates) \
                and standard_coordinates.size > 0:
            rotation_i_to_s = kabsch_rotation_matrix(input_coordinates,
                                                     standard_coordinates)
            B = (rotation_i_to_s @ B.reshape(-1, 3).T).T.reshape(n_internals,
                                                                 -1)
            B_inv = np.linalg.pinv(B)
            self._rotation_matrix = rotation_i_to_s
        
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
            self._D = np.array(internal_modes).T
        self.n_imag = n_imag
        self._atomic_masses = np.array(atomic_masses)
        self.internal_indices = internal_indices
        self.internal_names = internal_names
        self._input_coordinates = input_coordinates
        self._standard_coordinate = standard_coordinates
    
    def _parse_pes_file(self, file):
        # Read PES file
        lines = open(file).readlines()
        
        # Set up read flags
        read_hessian = False

        # Set up containers for reading data
        hessian = []

        for line in lines:
            # Read Cartesian force constant matrix
            if read_hessian:
                try:
                    split_line = line.strip().split()
                    values = [float(value.replace("D", "E")) for value
                              in split_line]
                    hessian.extend(values)
                except ValueError:
                    read_hessian = False
            # Read number of atoms
            if " NAtoms= " in line:
                n_atoms = int(line.strip().split()[1]) 
            if "Cartesian Force constants:" in line:
                read_hessian = True
                hessian = []

        # Construct force constant matrix from lower triangle
        fc_matrix = np.zeros((n_atoms * 3, n_atoms * 3))
        fc_matrix[np.tril_indices_from(fc_matrix)] = hessian
        fc_matrix = np.triu(fc_matrix.T, 1) + fc_matrix
        fc_matrix = (self._rotation_matrix @ fc_matrix.T.reshape(-1, 3, 3) @ 
                     self._rotation_matrix).reshape(12, 12).T[:3]

        # Set up attributes
        self._fc_matrix = fc_matrix

    def __repr__(self):
        return f"{self.__class__.__name__}({self._log_file!r})"    