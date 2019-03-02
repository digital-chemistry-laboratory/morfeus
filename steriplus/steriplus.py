import math
from typing import NamedTuple

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import scipy.spatial
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation

from steriplus.data import atomic_numbers, bondi_radii, crc_radii, jmol_colors
from steriplus.io import read_gjf, read_xyz
from steriplus.plotting import ax_3D

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
        area (float)                :   Area of vdW surface
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
    """
    def __init__(self, element_ids, coordinates, atom_1, atom_2, radii=[], radii_type="crc", density=0.005, n_rot_vectors=360):
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

        #  Handle case of parallel vectors
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
        self.area = hull.area

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
            n_points = self.area / self.density / density
            step = round(n_points / len(coordinates))

            # Set up dictionary for coloring atomic centers
            element_list = [atom.element_id for atom in atom_list]
            color_dict = {element_id: jmol_colors[element_id] for element_id in set(element_list)}

            # Plot the coordinates on the vdW surface
            ax.scatter(coordinates[::step,0], coordinates[::step,1], coordinates[::step,2], alpha=0.1)

            # Plot the atomic centers
            for i, atom in enumerate(atom_list, start=1):
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
        plt.figure()

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
        coodinates = self.coordinates[:, [i_1, i_2]]
        B_1 = self.B_1[[i_1, i_2]]
        B_5 = self.B_5[[i_1, i_2]]
        L = self.L[[i_1, i_2]]
        vector = self.vector[[i_1, i_2]]
        coordinates = self.coordinates[:, [i_1, i_2]]
        atom_coordinates_list = [atom.coordinates[:, [i_1, i_2]] for atom in atom_list]
        atom_2_coordinates = atom_list[atom_2 - 1].coordinates[:,[i_1, i_2]]
        atom_1_coordinates = atom_list[atom_1 - 1].coordinates[:,[i_1, i_2]]
        atom_coordinates = np.vstack(atom_coordinates_list)

        # Set the density of points for plotting with input factor
        n_points = self.area / self.density / density
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

class Atom(NamedTuple):
    """Atom class"""
    element_id: int
    radius: float
    coordinates: list

class Sphere:
    """Sphere class for creating and holding points on vdW surface.

    Args:
        center (list)      :    Center of sphere
        density (float)    :    Density of points in Å^-2
        radius (float)     :    Radius in Å

    Attributes:
        area (float)            :   Area of sphere in Å^2
        center (list)           :   Center of sphere
        circumference (float)   :   Circumference in Å
        points (ndarray)        :   Points on vdW surface of sphere
        radius (list)           :   Radius in Å
    """

    def __init__(self, center, radius, density=0.005):
        self.center = center
        self.radius = radius
        self.circumference = math.pi * radius * 2
        self.area = 4 * radius**2 * math.pi
        self.points = self.get_points_polar(density=density)

    def get_points_polar(self, density):
        """Calculates points on the vdW surface of the sphere.

        Args:
            density (float)     :   Density of points on the vdW surface in Å^-2

        Returns:
            points (ndarray)    :   Points on vdW surface of sphere
        """
        # Calculate number of points
        n = round((self.area / density / 2)**(1 / 2))

        # Get the range of theta and phi
        theta = np.linspace(0, math.pi, n)
        phi = np.linspace(0, 2 * math.pi, 2 * n)

        # Combine together all the possible combinations of theta and phi
        combied_theta_phi = np.dstack(np.meshgrid(theta, phi)).reshape(-1, 2)

        # Get the Cartesian coordinates
        theta = combied_theta_phi[:,0]
        phi = combied_theta_phi[:,1]
        points = self.get_cartesian_coordinates(self.radius, theta, phi)
        points = points + self.center

        return points

    @staticmethod
    def get_cartesian_coordinates(r, theta, phi):
        """Converts polar to Cartesian coordinates.

        Args:
            r (float)           :   Radius in Å
            theta (ndarray)     :   Array of theta angles in radians
            phi (ndarray)       :   Array of phi angles in radians

        Returns:
            points (ndarray)    :   Array of xyz points
        """
        # Calculate x, y and z coordinates
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        # Stack coordinates as columns
        points = np.column_stack((x, y, z))

        return points

    def __repr__(self):
        return f"{self.__class__.__name__}(center: {self.center}, radius: {self.radius})"

def get_radii(element_id_list, radii_type="crc"):
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

    # Get radii for elements in the element id list
    radii_list = [radii_dict[element_id] for element_id in element_id_list]

    return radii_list
