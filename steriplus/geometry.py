"""This module contains help classes and funcitons related to geometry

Classes:
    ConeAngleAtom: Atom used in cone angle calculations.
    ConeAngleCone: Cone used in cone angle calculations.
    Cone: Cone class.
    Sphere: Sphere class.

Functions:
    rotatate_coordinates: Rotate coordinates given a vector and an axis
"""

import numpy as np
from scipy.spatial.transform import Rotation
import math

class ConeAngleCone:
    """Class for supporting cone angle calculations.

    Args:
        angle (float)       :   Cone angle in radians
        atoms (atoms)       :   List of 1-index atoms that are tangent to cone
        normal (ndarray)    :   Normal vector of cone

    Attributes:
        angle (float)       :   Cone angle in radians
        atoms (atoms)       :   List of 1-index atoms that are tangent to cone
        normal (ndarray)    :   Normal vector of cone
    """
    def __init__(self, angle, atoms, normal):
        self.angle = angle
        self.atoms = atoms
        self.normal = normal

    def is_inside(self, atom):
        """Tests if atom lies inside the cone
        Args:
            atom        :   Atom to test.
        Returns:
            (bool)      :   True if inside, False if outside.
        """
        # Get vertex angle of atom
        beta = atom.cone.angle

        # Calculate angle between cone normal vector and normal vector to atom
        cos = np.dot(atom.cone.normal, self.normal)

        # Take into account numerical problems that sometimes gives a value
        # somewhat above 1
        if 1 - cos > 0 and 1 - cos < 1e-5:
            cos = 1
        angle = math.acos(cos)

        # Check if atom lies inside cone, within numerical reason
        diff = self.angle - (beta + angle)
        if diff > -1e-5:
            return True
        else:
            return False

class ConeAngleAtom:
    """Class for supporting cone angle calculations.

    Args:
        coordinates (ndarray)   :   Atomic coordinates in Å
        element_id (int)        :   Atomic number (0-indexed)
        index (int)             :   Atom index starting from 1
        radius (float)          :   vdW radius of atom in Å

    Attributes:
        cone (object)           :   Cone tangent to atom
        coordinates (ndarray)   :   Atomic coordinates in Å
        index (int)             :   Atom index starting from 1
        radius (float)          :   vdW radius of atom in Å

    """
    def __init__(self, coordinates, radius, index, element_id):
        self.coordinates = coordinates
        self.radius = radius
        self.cone = None
        self.index = index
        self.element_id = element_id

    def get_cone(self):
        """Constructs cone tangent to atom."""
        vector = self.coordinates
        n = vector / np.linalg.norm(vector)
        sin_alpha = self.radius / np.linalg.norm(vector)
        alpha = math.asin(sin_alpha)
        cone = ConeAngleCone(alpha, [self], n)
        self.cone = cone

class SASAAtom:
    """Class for supporting SASA calculations.

    Args:
        coordinates (ndarray)   :
        element_id (int)        :
        index (int)             :
        radius (float)          :

    Attributes:
        accessible_points (list)    :   Points accessible to solvent 
        area (float)                :   Solvent-accessible surface area (Å^2)
        coordinates (ndarray)       :   Coordinates (Å)
        element_id (int)            :   Element id as atomic number
        index (int)                 :   Atomic index (1-indexed)
        occluded_points (list)      :   Points occluded by other atoms
        radius (float)              :   Radius (Å)
        volume (float)              :   Volume inside solvent-accessible surface
                                        area (Å^3)
    """
    def __init__(self, element_id, radius, coordinates, index):
        self.coordinates = coordinates
        self.radius = radius
        self.index = index
        self.element_id = element_id
        self.occluded_points = None
        self.accessible_points = None
        self.area = None
        self.volume = None

class Cone:
    def __init__(self, angle, height, direction, spacing=0.1):
        circle_list = []
        dists = np.arange(0, height, spacing)
        theta = angle * math.pi / 180
        point_list = []
        for dist in dists:
            s = np.linspace(0, math.pi * 2, 50)
            t = np.linspace(0, math.pi * 2, 50)
            x = dist * np.tan(theta) * np.cos(t)
            y = dist * np.tan(theta) * np.sin(t)
            z = np.full((50,1), dist)
            points = np.column_stack([x, y, z])
            point_list.append(points)
        points = np.vstack(point_list)
        self.angle = angle
        self.center = [0, 0, 0]
        self.height = height
        self.points = points
        self.rotate(direction)
        self.direction = direction / np.linalg.norm(direction)

    def rotate(self, vector):
        vector = vector /np.linalg.norm(vector)

         # Get rotation quaternion that overlays vector with x-asis
        x_axis = vector
        vector = np.array([0, 0, 1])
        real = np.dot(vector, x_axis).reshape(1) + 1

        #  Handle case of antiparallel vectors
        if real < 1e-6:
            w = np.cross(vector, np.array([0, 0, 1]))
            if np.linalg.norm(w) < 1e-6:
                w = np.cross(vector, np.array([1, 0, 0]))
        else:
            w = np.cross(vector, x_axis)

        q = np.concatenate((w, real))
        q = q / np.linalg.norm(q)

        # Rotate atomic coordinates
        rot = Rotation.from_quat(q)
        self.points = rot.apply(self.points)

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
        offset = 2.0 / n
        increment = math.pi * (3.0 - math.sqrt(5.0));

        i = np.arange(n)
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - np.square(y))
        phi = np.mod((i + rnd), n) * increment
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        points = np.column_stack((x, y, z))
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

def rotate_coordinates(coordinates, vector, axis):
    """Rotates coordinates with the rotation that aligns vector with axis.

    Args:
        axis (ndarray)                  :   Vector to align vector to
        coordinates (ndarray)           :   Coordinates to rotate
        vector (ndarray)                :   Vector to align

    Returns:
        rotated_coordinates (ndarray)   :   Rotated coordinates.
    """
    # Get rotation quaternion that overlays vector with x-asis
    real = np.dot(vector, axis).reshape(1) + 1
    #  Handle case of antiparallel vectors
    if real < 1e-6:
        w = np.cross(vector, np.array([0, 0, 1]))
        if np.linalg.norm(w) < 1e-6:
            w = np.cross(vector, np.array([1, 0, 0]))
    else:
        w = np.cross(vector, axis)

    q = np.concatenate((w, real))
    q = q / np.linalg.norm(q)

    # Rotate atomic coordinates
    rot = Rotation.from_quat(q)
    rotated_coordinates = rot.apply(coordinates)
    return rotated_coordinates
