"""Help classes and functions related to geometry."""

import math

import numpy as np
from scipy.spatial.transform import Rotation
import scipy.spatial

from morfeus.helpers import convert_elements, get_connectivity_matrix
from morfeus.data import ANGSTROM_TO_BOHR, cov_radii_pyykko

class Atom:
    """Atom common for morfeus calculations.

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
        self.coordination_number = None
        self.occluded_points = None
        self.point_areas = None
        self.p_values = None
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
        self.density = density

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
        n = int(round((self.area / density)))
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
        n = int(round((self.area / density / 2)**(1 / 2)))

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
            n = int(round((self.volume / density * 6 / math.pi)**(1 / 3)))
        else:
            n = int(round((self.area / density * 6 / math.pi)**(1 / 3)))

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
        return (f"{self.__class__.__name__}(center: {self.center}, "
                f"radius: {self.radius})")

class InternalCoordinates:
    def __init__(self):
        self.internal_coordinates = []
    
    def add_bond(self, atom_1, atom_2):
        bond = Bond(atom_1, atom_2)
        if not bond in self.internal_coordinates:
            self.internal_coordinates.append(bond)
    
    def add_angle(self, atom_1, atom_2, atom_3):
        angle = Angle(atom_1, atom_2, atom_3)
        if not angle in self.internal_coordinates:
            self.internal_coordinates.append(angle)
    
    def add_dihedral(self, atom_1, atom_2, atom_3, atom_4):
        dihedral = Dihedral(atom_1, atom_2, atom_3, atom_4)
        if not dihedral in self.internal_coordinates:
            self.internal_coordinates.append(dihedral)
    
    def add_internal_coordinate(self, atoms):
        if len(atoms) == 2:
            self.add_bond(*atoms)
        elif len(atoms) == 3:
            self.add_angle(*atoms)
        elif len(atoms) == 4:
            self.add_dihedral(*atoms)
    
    def detect_bonds(self, elements, coordinates):
        # Detect bonds based on covalent radii
        connectivity_matrix = get_connectivity_matrix(elements, coordinates, radii_type="pyykko")
        indices = np.where(connectivity_matrix)

        bonds = set()
        for i, j in zip(*indices):
            bond = frozenset([i + 1, j + 1])
            bonds.add(bond)
        
        # Add each bond as an internal coordinate
        for bond in bonds:
            i, j = sorted(bond)
            self.add_bond(i, j)
    
    def get_B_matrix(self, coordinates):
        b_vectors = []
        for internal_coordinate in self.internal_coordinates:
            b_vectors.append(internal_coordinate.get_b_vector(coordinates))
        B_matrix = np.vstack(b_vectors)
        
        return B_matrix
        
    def __repr__(self):
        return (f"{self.__class__.__name__}({len(self.internal_coordinates)} coordinates)")        

class Bond:
    def __init__(self, atom_1, atom_2):
        self.i = atom_1
        self.j = atom_2
        self.atoms = [atom_1, atom_2]
    
    def get_b_vector(self, coordinates):
        i, j = self.i - 1, self.j -1
        v = (coordinates[i] - coordinates[j]) * ANGSTROM_TO_BOHR
        r = np.linalg.norm(v)
        grad  = np.array([v / r, -v / r])
        
        b_vector = np.zeros(coordinates.size)
        b_vector[i * 3:i * 3 + 3] = grad[0]
        b_vector[j * 3:j * 3 + 3] = grad[1]
        
        return b_vector

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.i == other.i or self.i == other.j) and \
                   (self.j == other.i or self.j == other.j)
        else:
            return False
    
    def __hash__(self):
        return hash((self.i, self.j))

    def __repr__(self):
        atoms = ", ".join(str(i) for i in sorted(self.atoms))
        return (f"{self.__class__.__name__}({atoms})")

class Angle:
    def __init__(self, atom_1, atom_2, atom_3):
        self.i = atom_1
        self.j = atom_2
        self.k = atom_3
        self.atoms = [atom_1, atom_2, atom_3]
    
    def get_b_vector(self, coordinates):
        i, j, k = self.i - 1, self.j - 1, self.k - 1
        v_1 = (coordinates[i] - coordinates[j]) * ANGSTROM_TO_BOHR
        v_2 = (coordinates[k] - coordinates[j]) * ANGSTROM_TO_BOHR
        dot_product = np.dot(v_1, v_2) / (np.linalg.norm(v_1) * np.linalg.norm(v_2))
        if dot_product < -1:
            dot_product = -1
        elif dot_product > 1:
            dot_product = 1
        phi = np.arccos(dot_product)
        if abs(phi) > np.pi - 1e-6:
            grad = [
                (np.pi - phi) / (2 * np.linalg.norm(v_1) ** 2) * v_1,
                (1 / np.linalg.norm(v_1) - 1 / np.linalg.norm(v_2)) * (np.pi - phi) / (2 * np.linalg.norm(v_1)) * v_1,
                (np.pi - phi) / (2 * np.linalg.norm(v_2) ** 2) * v_2,
            ]
        else:
            grad = [
                1 / np.tan(phi) * v_1 / np.linalg.norm(v_1) ** 2
                - v_2 / (np.linalg.norm(v_1) * np.linalg.norm(v_2) * np.sin(phi)),
                (v_1 + v_2) / (np.linalg.norm(v_1) * np.linalg.norm(v_2) * np.sin(phi))
                - 1 / np.tan(phi) * (v_1 / np.linalg.norm(v_1) ** 2 + v_2 / np.linalg.norm(v_2) ** 2),
                1 / np.tan(phi) * v_2 / np.linalg.norm(v_2) ** 2
                - v_1 / (np.linalg.norm(v_1) * np.linalg.norm(v_2) * np.sin(phi)),
            ]
        
        b_vector = np.zeros(coordinates.size)
        b_vector[i * 3:i * 3 + 3] = grad[0]
        b_vector[j * 3:j * 3 + 3] = grad[1]
        b_vector[k * 3:k * 3 + 3] = grad[2]
        
        return b_vector

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.i == other.i or self.i == other.k) and \
                   (self.k == other.i or self.k == other.k) and \
                    self.j == other.j
        else:
            return False
    
    def __hash__(self):
        return hash((self.i, self.j, self.k))

    def __repr__(self):
        atoms = ", ".join(str(i) for i in sorted(self.atoms))
        return (f"{self.__class__.__name__}({atoms})")

class Dihedral:
    def __init__(self, atom_1, atom_2, atom_3, atom_4):
        self.i = atom_1
        self.j = atom_2
        self.k = atom_3
        self.l = atom_4
        self.atoms = [atom_1, atom_2, atom_3, atom_4]

    def get_b_vector(self, coordinates):       
        i, j, k, l = self.i - 1, self.j - 1, self.k - 1, self.l - 1
        v_1 = (coordinates[i] - coordinates[j]) * ANGSTROM_TO_BOHR
        v_2 = (coordinates[l] - coordinates[k]) * ANGSTROM_TO_BOHR
        w = (coordinates[k] - coordinates[j]) * ANGSTROM_TO_BOHR
        ew = w / np.linalg.norm(w)
        a_1 = v_1 - np.dot(v_1, ew) * ew
        a_2 = v_2 - np.dot(v_2, ew) * ew
        sgn = np.sign(np.linalg.det(np.array([v_2, v_1, w])))
        sgn = sgn or 1
        dot_product = np.dot(a_1, a_2) / (np.linalg.norm(a_1) * np.linalg.norm(a_2))
        if dot_product < -1:
            dot_product = -1
        elif dot_product > 1:
            dot_product = 1
        phi = np.arccos(dot_product) * sgn

        if abs(phi) > np.pi - 1e-6:
            g = np.cross(w, a_1)
            g = g / np.linalg.norm(g)
            A = np.dot(v_1, ew) / np.linalg.norm(w)
            B = np.dot(v_2, ew) / np.linalg.norm(w)
            grad = [
                g / (np.linalg.norm(g) * np.linalg.norm(a_1)),
                -((1 - A) / np.linalg.norm(a_1) - B / np.linalg.norm(a_2)) * g,
                -((1 + B) / np.linalg.norm(a_2) + A / np.linalg.norm(a_1)) * g,
                g / (np.linalg.norm(g) * np.linalg.norm(a_2)),
            ]
        elif abs(phi) < 1e-6:
            g = np.cross(w, a_1)
            g = g / np.linalg.norm(g)
            A = np.dot(v_1, ew) / np.linalg.norm(w)
            B = np.dot(v_2, ew) / np.linalg.norm(w)
            grad = [
                g / (np.linalg.norm(g) * np.linalg.norm(a_1)),
                -((1 - A) / np.linalg.norm(a_1) + B / np.linalg.norm(a_2)) * g,
                ((1 + B) / np.linalg.norm(a_2) - A / np.linalg.norm(a_1)) * g,
                -g / (np.linalg.norm(g) * np.linalg.norm(a_2)),
            ]
        else:
            A = np.dot(v_1, ew) / np.linalg.norm(w)
            B = np.dot(v_2, ew) / np.linalg.norm(w)
            grad = [
                1 / np.tan(phi) * a_1 / np.linalg.norm(a_1) ** 2
                - a_2 / (np.linalg.norm(a_1) * np.linalg.norm(a_2) * np.sin(phi)),
                ((1 - A) * a_2 - B * a_1) / (np.linalg.norm(a_1) * np.linalg.norm(a_2) * np.sin(phi))
                - 1
                / np.tan(phi)
                * ((1 - A) * a_1 / np.linalg.norm(a_1) ** 2 - B * a_2 / np.linalg.norm(a_2) ** 2),
                ((1 + B) * a_1 + A * a_2) / (np.linalg.norm(a_1) * np.linalg.norm(a_2) * np.sin(phi))
                - 1
                / np.tan(phi)
                * ((1 + B) * a_2 / np.linalg.norm(a_2) ** 2 + A * a_1 / np.linalg.norm(a_1) ** 2),
                1 / np.tan(phi) * a_2 / np.linalg.norm(a_2) ** 2
                - a_1 / (np.linalg.norm(a_1) * np.linalg.norm(a_2) * np.sin(phi)),
            ]
        b_vector = np.zeros(coordinates.size)
        b_vector[i * 3:i * 3 + 3] = grad[0]
        b_vector[j * 3:j * 3 + 3] = grad[1]
        b_vector[k * 3:k * 3 + 3] = grad[2]
        b_vector[l * 3:l * 3 + 3] = grad[3]
        
        return b_vector        

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return all([self.i == other.i, self.j == other.j,
                        self.k == other.k, self.l == other.l]) or \
                   all([self.i == other.l, self.j == other.k,
                        self.k == other.j, self.l == other.i])
        else:
            return False

    def __hash__(self):
        return hash((self.i, self.j, self.k, self.l))

    def __repr__(self):
        atoms = ", ".join(str(i) for i in sorted(self.atoms))
        return (f"{self.__class__.__name__}({atoms})")

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

def kabsch_rotation_matrix(P, Q, center=True):
    """Construct the rotation matrix that overlays the points in P with the
    points in Q with minimum RMSD.
    https://en.wikipedia.org/wiki/Kabsch_algorithm

    Args:
        P (ndarray): Coordinates to be rotated
        Q (ndarray): Reference coordinates

    Returns:
        R (ndarray): Rotation matrix
    """
    # Calculate centroids and center coordinates
    if center:
        centroid_P = np.mean(P, axis=0)
        centroid_Q = np.mean(Q, axis=0)
        P -= centroid_P
        Q -= centroid_Q

    # Calculate cross-covariance matrix
    H = P.T @ Q

    # Perform SVD
    U, S, V_T = np.linalg.svd(H)

    # Determine correction for right-hand coordinate system
    d = np.sign(np.linalg.det(V_T.T @ U.T))

    # Obtain rotation matrix
    R = V_T.T @ np.array([[1, 0, 0], [0, 1, 0], [0, 0, d]]) @ U.T
    
    return R

def sphere_line_intersection(vector, center, radius):
    """Get points of intersection between line and sphere. Follows the
    procedure outlined here: http://paulbourke.net/geometry/circlesphere/.

    Args:
        vector (ndarray): Vector giving direction of line
        center (ndarray): Center of sphere
        radius (float): Radius of sphere

    Returns:
        intersection_points (list): Intersection points
    """
    # Set up points
    p_1 = vector
    p_2 = vector * 2

    # Set up terms for quadratic equation
    a = (p_2[0] - p_1[0]) ** 2 + (p_2[1] - p_1[1]) ** 2 \
        + (p_2[2] - p_1[2]) ** 2
    b = 2 * ((p_2[0] - p_1[0]) * (p_1[0] - center[0]) \
        + (p_2[1] - p_1[1]) * (p_1[1] - center[1]) \
        + (p_2[2] - p_1[2]) * (p_1[2] - center[2]))
    c = center[0] ** 2 + center[1] ** 2 + center[2] ** 2 \
        + p_1[0] ** 2 + p_1[1] ** 2 + p_1[2] ** 2 \
        - 2 * (center[0] * p_1[0] + center[1] * p_1[1] + center[2] * p_1[2]) \
        - radius ** 2

    # Determine value within the square root and select cases
    within_sqrt = b ** 2 - 4 * a * c
    if within_sqrt < 0:
        us = []
    elif within_sqrt == 0:
        us = [-b / (2 * a)]
    elif within_sqrt > 0:
        us = [(-b + math.sqrt(within_sqrt)) / (2 * a), (-b - math.sqrt(within_sqrt)) / (2 * a)]
    
    # Calculate intersection points.
    intersection_points = [p_1 + u * (p_2 - p_1) for u in us]
    
    return intersection_points