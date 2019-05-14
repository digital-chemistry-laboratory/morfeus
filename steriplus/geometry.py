"""Help classes and functions related to geometry

Classes:
    Atom: Atom class.
    Cone: Cone used in cone angle calculations.
    Sphere: Sphere for representing atoms

Functions:
    rotatate_coordinates: Rotate coordinates given a vector and an axis
"""
import math

import numpy as np
from scipy.spatial.transform import Rotation

class Atom:
    """Atom common for Steriplus calculations.

     Args:
        coordinates (ndarray): Coordinates (Å)
        element (int): Atomic number (starting from 1) 
        index (int): Atom index (starting from 1)
        radius (float): vdW radius (Å)

    Attributes:
        accessible_points (list): Points accessible to solvent 
        area (float): Solvent-accessible surface area (Å**2)
        cone (object): Cone tangent to atom
        coordinates (ndarray): Coordinates (Å)
        element (int): Atomic number (starting from 1)
        index (int): Atom index (starting from 1)
        occluded_points (list): Points occluded by other atoms
        radius (float): vdW radius (Å)
        volume (float): Volume inside solvent-accessible surface area (Å**3)   
    """
    def __init__(self, element, coordinates, radius, index):
        # Set up initial attributes
        self.coordinates = coordinates
        self.element = element
        self.index = index
        self.radius = radius

        # Set up settable attributes
        self.accessible_points = None
        self.area = None
        self.cone = None
        self.occluded_points = None
        self.volume = None
    
    def get_cone(self):
        """Construct cone for atom"""
        vector = self.coordinates
        normal = vector / np.linalg.norm(vector)
        sin_alpha = self.radius / np.linalg.norm(vector)
        alpha = math.asin(sin_alpha)
        cone = Cone(alpha, [self], normal)
        self.cone = cone

    def __repr__(self):
        return f"{self.__class__.__name__}({self.index!r})"

class Cone:
    """Cone used in cone angle calculations.

    Args:
        angle (float): Cone angle in radians
        atoms (list): List of 1-index atoms that are tangent to cone
        normal (ndarray): Normal vector of cone

    Attributes:
        angle (float): Cone angle in radians
        atoms (list): List of 1-index atoms that are tangent to cone
        normal (ndarray): Normal vector of cone
    """
    def __init__(self, angle, atoms, normal):
        self.angle = angle
        self.atoms = atoms
        self.normal = normal

    def is_inside(self, atom):
        """Tests if atom lies inside the cone

        Args:
            atom (object)   :   Atom to test.

        Returns:
            (bool)          :   True if inside, False if outside.
        """
        # Get vertex angle of atom
        beta = atom.cone.angle

        # Calculate angle between cone normal vector and unit vector to atom
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

    def __repr__(self):
        atoms = ', '.join([str(atom.index) for atom in self.atoms])
        return f"{self.__class__.__name__}(Tangent atoms: {atoms})"

class Sphere:
    """Sphere class for creating and holding points on vdW surface.

    Args:
        center (list): Coordinates for center (Å)
        density (float): Area per point (Å**2) for empty sphere and volume per
                         point (Å**3) for filled sphere.
        filled (bool): Whether a sphere with internal points should be
                       constructed (works only with projection)
        method (str): Method for generating points:
                     'fibonacci' (default), 'polar' or 'projection'
        radius (float): Radius (Å)

    Attributes:
        area (float): Area of sphere (Å**2)
        center (list): Coordinates for center (Å)
        circumference (float): Circumference (Å)
        points (ndarray): Points on sphere surface
        radius (float): Radius (Å)
    """
    def __init__(self, center, radius, density=0.005, method="fibonacci",
                 filled=False):
        self.center = center
        self.radius = radius
        self.circumference = math.pi * radius * 2
        self.area = 4 * radius**2 * math.pi
        self.volume = 4 * radius**3 * math.pi / 3

        if method == "polar":
            self.points = self._get_points_polar(density=density)
        elif method =="projection":
            self.points = self._get_points_projected(density=density,
                                                     filled=filled)
        elif method == "fibonacci":
            self.points = self._get_points_fibonacci(density=density)

    @staticmethod
    def _get_cartesian_coordinates(r, theta, phi):
        """Converts polar to Cartesian coordinates.

        Args:
            phi (ndarray): Phi angles (rad)
            r (float): Radius (Å)
            theta (ndarray): Theta angles (rad)

        Returns:
            points (ndarray): Cartesian points (Å)
        """
        # Calculate x, y and z coordinates
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        # Stack coordinates as columns
        points = np.column_stack((x, y, z))

        return points

    def _get_points_fibonacci(self, density):
        """Construct points on sphere surface by the Fibonacci golden spiral
        method. Method vectorized from 
        https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere

        Args:
            density (float): Area per point (Å**2)
        
        Returns:
            points (ndarray): Array of surface points.
        """
        # Generate points on unit sphere
        rnd = 1
        n = round((self.area / density))
        offset = 2.0 / n
        increment = math.pi * (3.0 - math.sqrt(5.0))

        i = np.arange(n)
        y = ((i * offset) - 1) + (offset / 2)
        r = np.sqrt(1 - np.square(y))
        phi = np.mod((i + rnd), n) * increment
        x = np.cos(phi) * r
        z = np.sin(phi) * r

        # Generate points and adjust with radius and center
        points = np.column_stack((x, y, z))
        points = points * self.radius
        points = points + self.center

        return points

    def _get_points_polar(self, density):
        """Construct points on sphere by polar coordinate method.

        Args:
            density (float): Area per point (Å**2) for empty sphere and volume
                             per point (Å**3) for filled sphere.
      
        Returns:
            points (ndarray): Array of surface points.
        """
        # Calculate number of points
        n = round((self.area / density / 2)**(1 / 2))

        # Set up points along theta and phi
        theta = np.linspace(0, math.pi, n)
        phi = np.linspace(0, 2 * math.pi, 2 * n)

        # Combine together all the possible combinations of theta and phi
        combined_theta_phi = np.dstack(np.meshgrid(theta, phi)).reshape(-1, 2)
        theta = combined_theta_phi[:,0]
        phi = combined_theta_phi[:,1]

        # Get the Cartesian coordinates
        points = self._get_cartesian_coordinates(self.radius, theta, phi)

        # Adjust to sphere center
        points = points + self.center

        return points

    def _get_points_projected(self, density, filled=False):
        """Construct points on sphere surface by projection.

        Args:
            density (float): Area per point (Å**2) for empty sphere and volume
                             per point (Å**3) for filled sphere.
            filled (bool): Whether to generate internal points
        
        Returns:
            points (ndarray): Array of surface points.
        """
        # Calculate number of points from density of empty or filled sphere.
        if filled:
            n = round((self.volume / density * 6 / math.pi)**(1 / 3))
        else:
            n = round((self.area / density * 6 / math.pi)**(1 / 3))

        # Generate points in box
        r = self.radius
        x = np.linspace(-r, r, n)
        y = np.linspace(-r, r, n)
        z = np.linspace(-r, r, n)
        points = np.stack(np.meshgrid(x, y, z), -1).reshape(-1, 3)

        # Remove points outside of sphere
        lengths = np.linalg.norm(points, axis=1)
        points = points[lengths <= r]
        
        # Project points onto sphere surface if sphere is empty
        if not filled:
            points = points / np.linalg.norm(points, axis=1).reshape(-1,1) * r
        
        # Adjust with sphere center
        points = points + self.center

        return points

    def __repr__(self):
        return (f"{self.__class__.__name__}(center: {self.center}, ",
                "radius: {self.radius})")

def rotate_coordinates(coordinates, vector, axis):
    """Rotates coordinates by the rotation that aligns vector with axis.

    Args:
        axis (ndarray): Reference vector to align input vector to
        coordinates (ndarray): Coordinates to rotate
        vector (ndarray): Vector to align

    Returns:
        rotated_coordinates (ndarray): Rotated coordinates.
    """
    # Get real part of quaternion
    real = np.dot(vector, axis).reshape(1) + 1

    # Get imaginary dimensions and handle case of antiparallel vectors
    if real < 1e-6:
        w = np.cross(vector, np.array([0, 0, 1]))
        if np.linalg.norm(w) < 1e-6:
            w = np.cross(vector, np.array([1, 0, 0]))
    else:
        w = np.cross(vector, axis)

    # Form quaternion and normalize
    q = np.concatenate((w, real))
    q = q / np.linalg.norm(q)

    # Rotate atomic coordinates
    rotation = Rotation.from_quat(q)
    rotated_coordinates = rotation.apply(coordinates)

    return rotated_coordinates
