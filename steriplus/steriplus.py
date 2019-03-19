import math
from typing import NamedTuple

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import scipy.spatial
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation

from steriplus.data import atomic_numbers, bondi_radii, crc_radii, jmol_colors
from steriplus.io import read_gjf, read_xyz, create_rdkit_mol
from steriplus.plotting import ax_3D, coordinate_axes, set_axes_equal
from steriplus.geometry import rotate_coordinates
try:
    from rdkit.Chem import rdFreeSASA
except:
    rdFreeSASA = None
    print("SASA calculations not available on Windows.")

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

        set_axes_equal(ax)

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self.element_ids)!r} atoms)"

class BuriedVolume:
    """Performs and stores the results of a buried volume calculation.

    Args:
        center (list)               :   List of x, y, z coordinates in Å for the
                                        metal center
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
    def __init__(self, element_ids, coordinates, center, exclude_list=[],
                 radii=[], include_hs=False, radius=3.5, radii_type="bondi",
                 radii_scale=1.17, density=0.001):

        # Construct sphere at metal center
        s = Sphere(center, radius, method="projection", density=density,
                   filled=True)

        # Save all coordinates before deleting some atoms
        self.all_coordinates = np.array(coordinates)
        self.density = density

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
            if element == 0:
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
        center = np.array(center)

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
        print("%V_bur:", round(self.buried_volume * 100, 1))

    def plot_steric_map(self, z_axis_atoms, filename=None, levels=30, grid=100, all_positive=True, cmap="jet"):
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
        buried_points = self.buried_points
        center = np.array(self.sphere.center)
        coordinates = self.all_coordinates
        atom_coordinates = self.atom_coordinates

        z_axis_atom_coordinates = coordinates[np.array(z_axis_atoms) - 1]
        point = np.mean(z_axis_atom_coordinates)
        vector = point - center
        vector = vector / np.linalg.norm(vector)

        # Translate coordinates
        atom_coordinates -= center
        center -= center

        #Rotate coordinate system
        atoms = self.atoms
        atom_coordinates = rotate_coordinates(atom_coordinates, vector, np.array([0, 0, -1]))

        # Determine z_values
        r = self.sphere.radius
        x_ = np.linspace(-r, r, grid)
        y_ = np.linspace(-r, r, grid)

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
                x_s = atom_coordinates[i, 0]
                y_s = atom_coordinates[i, 1]
                z_s = atom_coordinates[i, 2]
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
                        if np.linalg.norm(np.array([x, y, z_max])) > r:
                            z_max = np.nan
                else:
                        if np.linalg.norm(np.array([x, y, z_max])) > r:
                            z_max = np.nan
            else:
                z_max = np.nan
            z.append(z_max)

        # Create interaction surface
        z = np.array(z).reshape(len(x_), len(y_))
        self.interaction_surface_x = x_
        self.interaction_surface_y = y_
        self.interaction_surface_z = z

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
        buried_points = self.buried_points
        free_points = self.free_points
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
            for i, atom in enumerate(atom_list, start=1):
                x, y, z = atom.coordinates[:,0], atom.coordinates[:,1], atom.coordinates[:,2]
                ax.scatter(x, y, z, color=color_dict[atom.element_id], s=50, edgecolors="k",)

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

class SASA:
    """Calculates the solvent accesible surface area (SASA) with RDKit from a
    list of element ids and coordinates.
    def __init__(self, , center, exclude_list=[],
                 radii=[], include_hs=False, radius=3.5, radii_type="bondi",
                 radii_scale=1.17, density=0.001):
    Args:
        mol                     :    rdkit mol object
        atom_index_list (list)  :    List of atom indices (starting at 1)
    """
    def __init__(self, element_ids, coordinates, radii=None, radii_type="crc"):
        # Check if rdFreeSASA is loaded
    #    if not rdFreeSASA:
    #        raise Exception("rdFreeSASA not loaded")

        # Converting element ids to atomic numbers if the are symbols
        if type(element_ids[0]) == str:
            element_ids = [atomic_numbers[element_id] for element_id in element_ids]
        self.element_ids = element_ids

        # Getting radii if they are not supplied
        if not radii:
            radii = get_radii(element_ids, radii_type=radii_type)
        self.radii = radii

        # Create rdkit mol object from the element_ids and coordinates
        self.mol = create_rdkit_mol(element_ids, coordinates)

        # Calculate total area
        self.total_area = rdFreeSASA.CalcSASA(mol, radii=radii)

        atom_areas = {}
        # Calculate area for each atom
        for atom in mol.GetAtoms():
            atom.SetProp("SASAClassName", "Polar")
            atom.SetProp("SASAClass", "2")
            q_atom = rdFreeSASA.MakeFreeSasaPolarAtomQuery()
            atom_area = rdFreeSASA.CalcSASA(mol, radii=radii, query=q_atom)
            atom_areas[atom.GetIdx()] = atom_area
            atom.ClearProp("SASAClassName")
            atom.ClearProp("SASAClass")

        self.atom_areas = atom_areas

    def print_report(self):
        pass

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

    def __init__(self, center, radius, density=0.005, method="fibonacci", filled=False):
        self.center = center
        self.radius = radius
        self.circumference = math.pi * radius * 2
        self.area = 4 * radius**2 * math.pi
        self.volume = 4 * radius**3 * math.pi / 3

        if filled:
            self.points = self.get_points_projected(density=density, filled=True)

        if method == "polar":
            self.points = self.get_points_polar(density=density)
        elif method =="projection":
            self.points = self.get_points_projected(density=density)
        elif method == "fibonacci":
            self.points = self.get_points_fibonacci(density=density)

    def get_points_fibonacci(self, density):
        rnd = 1
        n = round((self.area / density))
        point_list = []
        offset = 2.0 / n
        increment = math.pi * (3.0 - math.sqrt(5.0));

        for i in range(n):
            y = ((i * offset) - 1) + (offset / 2);
            r = math.sqrt(1 - pow(y, 2))

            phi = ((i + rnd) % n) * increment

            x = math.cos(phi) * r
            z = math.sin(phi) * r

            point_list.append([x,y,z])

        points = np.array(point_list)
        points = points / np.linalg.norm(points, axis=1).reshape(-1,1) * self.radius
        points = points + self.center

        return points

    def get_points_projected(self, density, filled=True):
        if not filled:
            n = round((self.area / density * 6 / math.pi)**(1 / 3))
        else:
            n = round((self.volume / density * 6 / math.pi)**(1 / 3))
        r = self.radius
        x = np.linspace(-r, r, n)
        y = np.linspace(-r, r, n)
        z = np.linspace(-r, r, n)
        points = np.stack(np.meshgrid(x, y, z), -1).reshape(-1, 3)
        numpoints = len(points)
        lengths = np.linalg.norm(points, axis=1)
        points = points[lengths <= r]
        if not filled:
            points = points / np.linalg.norm(points, axis=1).reshape(-1,1) * r
        points = points + self.center
        return points

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
