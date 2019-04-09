import math
from typing import NamedTuple
import time
import itertools

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

import scipy.linalg
import scipy.spatial
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

from steriplus.data import atomic_numbers, bondi_radii, crc_radii, jmol_colors
from steriplus.io import read_gjf, read_xyz
from steriplus.plotting import ax_3D, coordinate_axes, set_axes_equal
from steriplus.geometry import rotate_coordinates, Sphere, Cone, ConeAngleCone, ConeAngleAtom, SASAAtom

class Sterimol:
    """Performs and stores results of Sterimol calculation.

    Args:
        atom_1 (int)                :   Index of atom 1 (dummy atom)
        atom_2 (int)                :   Index of atom 2 (connected atom of
                                        substituent)
        coordinates (list)          :   List of coordinates in Å
        density (float)             :   Density of points on the vdW surface in
                                        Å^-2
        element_ids (list)          :   List of element_ids, either as atomic
                                        numbers or symbols
        n_rot_vectors (int)         :   Number of rotational vectors for
                                        determining B_1 and B_5
        radii (list)                :   List of radii (optional)
        radii_type (str)            :   Type of radii to use, either "crc" or
                                        "bondi"

    Attributes:
        hull_area (float)           :   Area of the convex hull enclosing the
                                        vdW surface in Å^2
        atom_1 (int)                :   Index of atom 1 (dummy atom)
        atom_2 (int)                :   Index of atom 2 (connected atom of
                                        substituent)
        atoms (list)                :   List of atoms objects
        B_1 (ndarray)               :   Sterimol B_1 vector
        B_1_value (float)           :   Sterimol B_1 value in Å
        B_5 (ndarray)               :   Sterimol B_5 vector
        B_5_value (float)           :   Sterimol B_5 value in Å
        bond length (float)         :   Bond length between atom 1 and atom 2 in
                                        Å
        coordinates (ndarray)       :   List of coordinates on vdW surface
        density (float)             :   Density of points on the vdW surface in
                                        Å^-2
        element_ids (list)          :   List of element_ids, either as atomic
                                        numbers
        L (ndarray)                 :   Sterimol L vector
        L_value (float)             :   Sterimol L value in Å
        L_value_uncorrected (float) :   Sterimol L value minus 0.40 Å
        radii (list)                :   List of radii
        vector (ndarray)            :   Vector between atom 1 and atom 2
        hull_volume (float)         :   Hull volume containing the vdW surface
                                        in Å^3
    """
    def __init__(self, element_ids, coordinates, atom_1, atom_2, radii=[], radii_type="crc", density=0.01, n_rot_vectors=360):
        # Converting element ids to atomic numbers if the are symbols
        if type(element_ids[0]) == str:
            element_ids = [atomic_numbers[element_id] for element_id in element_ids]
        self.element_ids = element_ids

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(element_ids, radii_type=radii_type)
        self.radii = radii

        # Setting up coordinate array
        atom_coordinates = np.array(coordinates)

        # Translate coordinates so origin is at atom 1
        atom_coordinates -= atom_coordinates[atom_1 - 1]

        # Get vector pointing from atom 2 to atom 1
        vector_2_to_1 = atom_coordinates[atom_2 - 1] - atom_coordinates[atom_1 - 1]
        vector_2_to_1 = vector_2_to_1 / np.linalg.norm(vector_2_to_1)

        # Get rotation quaternion that overlays vector with x-asis
        x_axis = np.array([1, 0, 0])
        real = np.dot(vector_2_to_1, x_axis).reshape(1) + 1

        #  Handle case of antiparallel vectors
        if real < 1e-6:
            w = np.cross(vector_2_to_1, np.array([0, 0, 1]))
            if np.linalg.norm(w) < 1e-6:
                w = np.cross(vector_2_to_1, np.array([1, 0, 0]))
        else:
            w = np.cross(vector_2_to_1, x_axis)

        q = np.concatenate((w, real))
        q = q / np.linalg.norm(q)

        # Rotate atomic coordinates
        rot = Rotation.from_quat(q)
        atom_coordinates = rot.apply(atom_coordinates)

        # Get list of atoms as Atom objects
        atom_list = []
        for element_id, radius, coord in zip(element_ids, radii, atom_coordinates):
            coord = np.array(coord).reshape(1,-1)
            atom = Atom(element_id, radius, coord)
            atom_list.append(atom)

        # Get list of spheres as Sphere objects
        sphere_list = []
        for i, atom in enumerate(atom_list, start=1):
            if not i == atom_1:
                sphere = Sphere(atom.coordinates, atom.radius, density=density)
                sphere_list.append(sphere)

        # Construct array of coordinates on the vdW surface of the spheres.
        coordinate_list = []
        for sphere in sphere_list:
            coordinate_list.append(sphere.points)
        coordinates = np.vstack(coordinate_list)

        # Prune coordinates that are not on the molecular vdW surface
        # Hull algorithm also removes more points that are note important
        # for determining Sterimol parameters.
        hull = ConvexHull(coordinates)
        coordinates = coordinates[hull.vertices]
        self.hull_area = hull.area
        self.hull_volume = hull.volume

        # Get vector from atom 1 to atom 2
        vector = atom_list[atom_2 - 1].coordinates - atom_list[atom_1 - 1].coordinates
        vector = vector / np.linalg.norm(vector)

        # Project coordinates onto vector
        c_values = vector.dot(coordinates.T)

        # Get L as largest projection along the vector
        L_value = np.max(c_values)
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
        max_c_values = np.max(c_values, axis=1)

        # Determine B1 and B5 from the smallest and largest scalar projections
        B_1_value = np.min(max_c_values)
        B_1 = rot_vectors[np.argmin(max_c_values)] * B_1_value

        B_5_value = np.max(max_c_values)
        B_5 = rot_vectors[np.argmax(max_c_values)] * B_5_value

        # Set up attributes
        self.atoms = atom_list
        self.coordinates = coordinates
        self.density = density
        self.vector = vector.reshape(-1)

        self.atom_1 = atom_1
        self.atom_2 = atom_2

        self.L = L.reshape(-1)
        self.L_value = L_value + 0.40
        self.L_value_uncorrected = L_value
        self.bond_length = np.linalg.norm(atom_list[atom_2 - 1].coordinates - atom_list[atom_1 - 1].coordinates)

        self.B_1 = B_1
        self.B_1_value = B_1_value

        self.B_5 = B_5
        self.B_5_value = B_5_value

    def print_report(self, verbose=False):
        """Prints the values of the Sterimol parameters.

        Args:
            verbose (bool) :    Toggles printing of uncorrected L_value and bond
                                length between atom 1 and atom 2.
        """
        if verbose:
            print(f"{'L':10s}{'B_1':10s}{'B_5':10s}{'L_uncorr':10s}{'d(a1-a2)':10s}")
            print(f"{self.L_value:<10.2f}{self.B_1_value:<10.2f}{self.B_5_value:<10.2f}{self.L_value_uncorrected:<10.2f}{self.bond_length:<10.2f}")
        else:
            print(f"{'L':10s}{'B_1':10s}{'B_5':10s}")
            print(f"{self.L_value:<10.2f}{self.B_1_value:<10.2f}{self.B_5_value:<10.2f}")

    def plot_3D(self, density=1, filename=None):
        """Plots a 3D representation of the molecule, vdW surface and Sterimol
        vectors.

        Args:
            density (float) :   Factor for adjusting the density of points in the
                                plot.
            filename (str)  :   Filename for saving the plot to file.
        """
        with ax_3D() as ax:
            ax.set_aspect('equal')
            # Setting up
            atom_1 = self.atom_1
            atom_2 = self.atom_2
            B_1 = self.B_1
            B_5 = self.B_5
            L = self.L
            atom_list = self.atoms
            vector = self.vector
            coordinates = self.coordinates

            # Set the density of points for plotting with input factor
            n_points = self.hull_area / self.density / density
            step = round(n_points / len(coordinates))

            # Set up dictionary for coloring atomic centers
            element_list = [atom.element_id for atom in atom_list]
            color_dict = {element_id: jmol_colors[element_id] for element_id in set(element_list)}

            # Plot the coordinates on the vdW surface
            ax.scatter(coordinates[::step,0], coordinates[::step,1], coordinates[::step,2], alpha=0.1)

            # Plot the atomic centers
            for atom in atom_list:
                x, y, z = atom.coordinates[:,0], atom.coordinates[:,1], atom.coordinates[:,2]
                ax.scatter(x, y, z, color=color_dict[atom.element_id], s=50, edgecolors="k",)

            # Plot the B_1 and B_5 vectors
            x, y, z = atom_list[atom_2 - 1].coordinates[:,0], atom_list[atom_2 - 1].coordinates[:,1], atom_list[atom_2 - 1].coordinates[:,2]
            ax.quiver(x, y, z, *vector)
            ax.quiver(x, y, z, *B_1)
            ax.quiver(x, y, z, *B_5)
            B_1_pos = B_1 + atom_list[atom_2 - 1].coordinates
            ax.text(*B_1_pos[0], "B1")
            B_5_pos = B_5 + atom_list[atom_2 - 1].coordinates
            ax.text(*B_5_pos[0], "B5")

            # Plot the L vector
            x, y, z = atom_list[atom_1 - 1].coordinates[:,0], atom_list[atom_1 - 1].coordinates[:,1], atom_list[atom_1 - 1].coordinates[:,2]
            ax.quiver(x, y, z, *L)
            L_pos = L + atom_list[atom_1 - 1].coordinates
            ax.text(*L_pos[0], "L")

            if filename:
                plt.savefig(filename)
            else:
                plt.show()

    def plot_2D(self, plane="yz", density=1, filename=None):
        """Plots a 2D representation of the molecule, vdW surface and Sterimol
        vectors.

        Args:
            density (float) :   Factor for adjusting the density of points in the
                                plot
            filename (str)  :   Filename for saving the plot to file
            plane (str)     :   Plotting plane ("xy", "xz" or "yz")
        """
        plt.figure(figsize=(15,10))

        # Set equal aspects for axes
        plt.axes().set_aspect('equal', 'datalim')

        atom_1 = self.atom_1
        atom_2 = self.atom_2
        atom_list = self.atoms

        # Select coordinates for plotting
        if plane == "yz":
            i_1 = 1
            i_2 = 2
        if plane == "xy":
            i_1 = 0
            i_2 = 1
        if plane == "xz":
            i_1 = 0
            i_2 = 2

        # Get coordinates according to plotting plane
        coordinates = self.coordinates[:, [i_1, i_2]]
        B_1 = self.B_1[[i_1, i_2]]
        B_5 = self.B_5[[i_1, i_2]]
        L = self.L[[i_1, i_2]]
        coordinates = self.coordinates[:, [i_1, i_2]]
        atom_coordinates_list = [atom.coordinates[:, [i_1, i_2]] for atom in atom_list]
        atom_2_coordinates = atom_list[atom_2 - 1].coordinates[:,[i_1, i_2]]
        atom_1_coordinates = atom_list[atom_1 - 1].coordinates[:,[i_1, i_2]]
        atom_coordinates = np.vstack(atom_coordinates_list)

        # Set the density of points for plotting with input factor
        n_points = self.hull_area / self.density / density
        step = max(round(n_points / len(coordinates)), 1)

        # Set up dictionary for coloring atomic centers
        element_list = [atom.element_id for atom in atom_list]
        color_dict = {element_id: jmol_colors[element_id] for element_id in set(element_list)}

        # Plot the coordinates on the vdW surface
        plt.scatter(coordinates[::step,0], coordinates[::step,1], alpha=0.1)

        # Plot the atomic centers
        for i, atom in enumerate(atom_list, start=0):
            plt.scatter(atom_coordinates[i,0], atom_coordinates[i,1], color=color_dict[atom.element_id], s=50, edgecolors="k",)

        # Plot the B_1 and B_5 vectors
        coord_1, coord_2 = atom_2_coordinates[:,0], atom_2_coordinates[:,1]
        plt.arrow(coord_1[0], coord_2[0], B_1[0], B_1[1], fc="k", ec="k", head_width=0.1, head_length=0.2, length_includes_head=True)
        plt.arrow(coord_1[0], coord_2[0], B_5[0], B_5[1], fc="k", ec="k", head_width=0.1, head_length=0.2, length_includes_head=True)
        B_1_pos = B_1 + atom_2_coordinates
        plt.text(B_1_pos[:,0][0], B_1_pos[:,1][0], "B_1")
        B_5_pos = B_5 + atom_2_coordinates
        plt.text(B_5_pos[:,0][0], B_5_pos[:,1][0], "B_5")

        # Plot the L vector in the case of xy or xz plotting plane
        if plane != "yz":
            coord_1, coord_2 = atom_1_coordinates[:,0], atom_2_coordinates[:,1]
            plt.arrow(coord_1[0], coord_2[0], L[0], L[1], fc="k", ec="k", head_width=0.1, head_length=0.2, length_includes_head=True)
            L_pos = L + atom_1_coordinates
            plt.text(L_pos[:,0][0], L_pos[:,1][0], "L")
        
        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self.element_ids)!r} atoms)"

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
        self.density = density
        self.excluded_atoms = exclude_list
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
        if type(element_ids[0]) == str:
            element_ids = [atomic_numbers[element_id] for element_id in
                           element_ids]

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
        for element_id, radius, coord in zip(element_ids, radii,
                                             atom_coordinates):
            coord = np.array(coord).reshape(1,-1)
            atom = Atom(element_id, radius, coord)
            atom_list.append(atom)

        # Prune sphere points which are within vdW radius of other atoms.
        tree = scipy.spatial.cKDTree(s.points, compact_nodes=False,
                                     balanced_tree=False)
        mask = np.zeros(len(s.points), dtype=bool)
        for atom in atom_list:
            if atom.radius + s.radius > np.linalg.norm(atom.coordinates):
                to_prune = tree.query_ball_point(atom.coordinates,
                                                 atom.radius)[0]
                mask[to_prune] = True
        buried_points = s.points[mask,:]
        free_points = s.points[np.invert(mask),:]

        # Calculate buried_volume
        self.buried_volume = len(buried_points) / len(s.points)

        # Set variables for outside access and function access.
        self.atoms = atom_list
        self.atom_coordinates = atom_coordinates
        self.sphere = s
        self.buried_points = buried_points
        self.free_points = free_points

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
        atoms = self.atoms
        center = np.array(self.sphere.center)
        all_coordinates = np.array(self._all_coordinates)
        coordinates = np.array(self.atom_coordinates)

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
        r = self.sphere.radius
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

    def plot_3D(self, filename=None, density=1):
        """Plots a 3D view of the sphere volume that is free
        Args:
            density (float)    :    Parameter that regulates the density of points of the
                                    sphere
            filename (str)     :    Name of file for saving the plot.
        """
        buried_points = np.array(self.buried_points)
        np.random.shuffle(buried_points)
        free_points = np.array(self.free_points)
        np.random.shuffle(free_points)
        atom_list = self.atoms

        # Set the density of points for plotting with input factor
        n_points = self.sphere.volume / self.density / density * 100
        step = round(n_points / (len(free_points) + len(buried_points)))
        
        # Set up dictionary for coloring atomic centers
        element_list = [atom.element_id for atom in atom_list]
        color_dict = {element_id: jmol_colors[element_id] for element_id in set(element_list)}

        with ax_3D() as ax:
            ax.scatter(buried_points[::step,0], buried_points[::step,1], buried_points[::step,2], s=1, c="r")
            ax.scatter(free_points[::step,0], free_points[::step,1], free_points[::step,2], s=1, c="b")

            # Plot the atomic centers
            for atom in atom_list:
                x, y, z = atom.coordinates[:,0], atom.coordinates[:,1], atom.coordinates[:,2]
                ax.scatter(x, y, z, color=color_dict[atom.element_id], s=50, edgecolors="k",)

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

class SASA:
    def __init__(self, element_ids, coordinates, radii=[], radii_type="crc", 
                 probe_radius=1.4, density=0.01):
        # Converting element ids to atomic numbers if the are symbols
        if type(element_ids[0]) == str:
            element_ids = [atomic_numbers[element_id] for
                           element_id in element_ids]
        self._element_ids = element_ids

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(element_ids, radii_type=radii_type)
        
        # Increment the radii with the probe radius
        radii = np.array(radii)
        radii = radii + probe_radius
        
        # Center coordinate system at geometric center
        coordinates = np.array(coordinates)
        center = np.mean(coordinates, axis=0)
        coordinates -= center
        
        # Construct list of atoms
        atoms = []
        for i, (coordinate, radius, element_id) in enumerate(zip(coordinates, radii, element_ids), start=1):
            atom = SASAAtom(element_id, radius, coordinate, i)
            atoms.append(atom)
        
        # Determine occluded and accesible points of each atom based on
        # distances to all other atoms (brute force)
        for atom in atoms:
            # Construct sphere for atom
            sphere = Sphere(atom.coordinates, atom.radius, density=density)
            
            # Select coordinates and radii for other atoms
            other_coordinates = [test_atom.coordinates for test_atom 
                                 in atoms if test_atom is not atom]
            other_radii = [test_atom.radius for test_atom 
                           in atoms if test_atom is not atom]
            other_radii = np.array(other_radii).reshape(-1, 1)

            # Get distances to other atoms and subtract radii
            distances = cdist(other_coordinates, sphere.points)
            distances -= other_radii

            # Take smallest distance and perform check
            min_distances = np.min(distances, axis=0)
            atom.occluded_points = sphere.points[min_distances < 0]
            atom.accessible_points = sphere.points[min_distances > 0]
        
        # Calculate atom areas and volumes
        for atom in atoms:
            # Get number of points of eache type
            n_accesible = len(atom.accessible_points)
            n_occluded = len(atom.occluded_points)
            n_points = len(atom.accessible_points) + len(atom.occluded_points)

            # Calculate part occluded and accessible
            ratio_occluded = n_occluded / n_points
            ratio_accesible = 1 - ratio_occluded
            
            # Calculate area
            area = 4 * np.pi * atom.radius ** 2 * ratio_accesible
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
        self.probe_radius = probe_radius 
        self.atom_areas = {atom.index: atom.area for atom in atoms}
        self.atom_volumes = {atom.index: atom.volume for atom in atoms}
        self.area = sum([atom.area for atom in atoms])
        self.volume = sum([atom.volume for atom in atoms])
        self._atoms = atoms
        self._density = density
    
    def plot_3D(self, highlight=[]):

        # Set up dictionary for coloring atomic centers
        element_list = [atom.element_id for atom in self._atoms]
        color_dict = {element_id: jmol_colors[element_id] 
                      for element_id in set(element_list)}
        
        with ax_3D() as ax:
            for atom in self._atoms:
                step = int(max(round(1 / self._density / 20), 1))
                accessible_points = np.array(atom.accessible_points)
                np.random.shuffle(accessible_points)
                x = accessible_points[::step, 0]
                y = accessible_points[::step, 1]
                z = accessible_points[::step, 2]
                if atom.index in highlight:
                    ax.scatter(x, y, z, c=color_dict[atom.element_id], 
                               alpha=1, edgecolors="k", linewidths=0.5)
                else:
                    ax.scatter(x, y, z, c=color_dict[atom.element_id],
                               alpha=0.3, edgecolor="k", linewidths=0.2)
    
    def print_report(self, verbose=False):
        print(f"Solvent-accessible surface area (Å^2): {self.area:.1f}")
        print("Volume inside solvent-accessible surface (Å^3): ",  
              f"{self.volume:.1f}")
        if verbose:
            print("Areas per atom (Å^2):")
            for i, area in self.atom_areas.items():
                print(f"{i:5d}: {area:5.1f}")

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
        # Make copies to avoid deletion
        coordinates = list(coordinates)
        element_ids = list(element_ids)
        
        # Removing central atom
        center_coordinates = np.array(coordinates[atom_1 - 1])
        del coordinates[atom_1 - 1]
        del element_ids[atom_1 - 1]

        # Converting element ids to atomic numbers if the are symbols
        if type(element_ids[0]) == str:
            element_ids = [atomic_numbers[element_id] for element_id in element_ids]
        self.element_ids = element_ids

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(element_ids, radii_type=radii_type)
        self.radii = radii

        # Setting up coordinate array
        atom_coordinates = np.array(coordinates)

        # Translate coordinates so origin is at atom 1
        atom_coordinates -= center_coordinates

        # Get list of atoms as Atom objects
        atoms = []
        for i, (coord, radius, element_id) in enumerate(zip(atom_coordinates, radii, element_ids), start=1):
            coord = np.array(coord)
            if np.linalg.norm(coord) < radius:
                raise Exception(f"Center within vdW surface of atom {i}")
            atom = ConeAngleAtom(coord, radius, i, element_id)
            atoms.append(atom)
            atom.get_cone()

        self.atoms = atoms

        # Search for cone over single atoms
        cone = self._search_one_cones()
        if cone:
            self.cone = cone
            self.cone_angle = math.degrees(cone.angle * 2)
            self.tangent_atoms = [atom.index for atom in cone.atoms]
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
                self.tangent_atoms = [atom.index for atom in cone.atoms]

        # Search for cones over triples of atoms
        if not cone:
            cone = self._search_three_cones()
            if cone:
                self.cone = cone
                self.cone_angle = math.degrees(cone.angle * 2)
                self.tangent_atoms = [atom.index for atom in cone.atoms]

    def _search_one_cones(self):
        """Searches over cones tangent to one atom

        Returns:
            max_1_cone (object)     :   Largest cone tangent to one atom
        """
        # Get the largest cone
        atoms = self.atoms
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
        test_list = [atom for atom in loop_list if atom not in max_2_cone.atoms]
        for atom in test_list:
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
            if cone.angle < lower_bound or cone.angle > upper_bound:
                remove_list.append(cone)

        for cone in reversed(remove_list):
            cone_list.remove(cone)

        # Keep cones that encompass all atoms
        keep_list = []
        for cone in cone_list:
            in_list = []
            test_list = [atom for atom in loop_list]
            for atom in test_list:
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
        coordinates = np.array([atom.coordinates for atom in self.atoms])
        centroid_vector = np.mean(coordinates, axis=0)
        centroid_unit_vector = centroid_vector / np.linalg.norm(centroid_vector)

        # Getting sums of angle to centroid and vertex angle.
        angle_sum_list = []
        for atom in loop_list:
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
        a_ij = (1 / math.sin(beta_ij)) * (0.5 * (beta_ij + beta_i - beta_j))
        b_ij = (1 / math.sin(beta_ij)) * (0.5 * (beta_ij - beta_i + beta_j))
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
        physical_angles = angles[np.argsort(tests)]

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
        print(f"Cone angle: {self.cone_angle:.1f}")
        print(f"No. tangent atoms: {len(self.tangent_atoms)}")

    def plot_3D(self, height=5, plot_cone=True):
        """Plot 3D representation of points on convex hull together with cone

        Args:
            height (float)      :   Height of the cone in Å
            plot_cone (bool)    :   Whether to plot the cone or just the points
        """
        # Set up dictionary for coloring atomic centers
        element_list = [atom.element_id for atom in self.atoms]
        color_dict = {element_id: jmol_colors[element_id] for element_id in set(element_list)}

        # Construct cone and set direction
        angle = np.degrees(self.cone.angle)
        if angle * 2 > 180:
            angle = 180 - angle
            normal = -self.cone.normal
        else:
            angle = angle
            normal = self.cone.normal

        cone = Cone(angle, height, normal, spacing=0.1)

        with ax_3D() as ax:
            ax.set_aspect('equal')
            # Plot vdW surface of atoms
            for atom in self.atoms:
                sphere = Sphere(atom.coordinates, atom.radius, density=0.5)
                cvx = ConvexHull(sphere.points)
                x, y, z = sphere.points.T
                tri = Triangulation(x, y, triangles=cvx.simplices)
                ax.plot_trisurf(tri, z, color=color_dict[atom.element_id])
            # Plot cone
            if plot_cone:
                ax.plot_trisurf(cone.points[:,0], cone.points[:,1], cone.points[:,2], alpha=0.1)
                ax.quiver(0, 0, 0, *normal, color="r")
            set_axes_equal(ax)
        
        plt.show()

class Atom(NamedTuple):
    """Atom class"""
    element_id: int
    radius: float
    coordinates: list

def get_radii(element_id_list, radii_type="crc", scale=1):
    """Gets radii for list of element ids

    Args:
        element_id_list (list)  :   List of element ids as atomic numbers or
                                    symbols
        radii_type (str)        :   Type of vdW radius, "bondi" or "crc"

    Returns:
        radii_list (list)       :   List of atomic radii
    """
    # Set up dictionary of atomic radii for all elements present in the list
    element_set = set(element_id_list)
    if radii_type == "bondi":
        radii_dict = {element_id: bondi_radii.get(element_id) for element_id in element_set}
    elif radii_type == "crc":
        radii_dict = {element_id: crc_radii.get(element_id) for element_id in element_set}

    # Get radii for elements in the element id list. Set 2 as default if does not exist
    radii_list = [radii_dict.get(element_id) * scale if radii_dict.get(element_id, 2.0) else 2.0 * scale for element_id in element_id_list ]

    return radii_list
