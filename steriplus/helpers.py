"""Help functions

Functions:
    check_distances: Controls distances for steric clashes.
    convert_elements: Converts atomic symbols to atomic numbers.
    get_radii: Returns radii for list of elements.
"""
import numpy as np
import pickle
import pkg_resources
from scipy.spatial.distance import cdist

from steriplus.data import atomic_numbers, atomic_symbols
from steriplus.data import radii_bondi, radii_crc, radii_rahm
from steriplus.data import cov_radii_pyykko, r2_r4
from steriplus.geometry import Atom

class D3Calculator:
    """Calculates C6(AA) and C8(AA) coefficients based on the procedure in 
    J. Chem. Phys. 2010, 132, 154104.

    Args:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)

    Attributes:
        c6_coefficients (list): C6(AA) coefficients (au)
        c8_coefficients (list): C8(AA) coefficients (au)
    """
    def __init__(self, elements, coordinates):
        # Convert elements to atomic numbers
        elements = convert_elements(elements)

        # Load the covalent radii
        radii = get_radii(elements, radii_type='pyykko')

        # Set up the atoms objects
        atoms = []
        for i, (element, coordinate, radius) in enumerate(zip(elements,
                coordinates, radii), start=1):
            atom = Atom(element, coordinate, radius, i)
            atoms.append(atom)
        self._atoms = atoms
            
        # Calculate the coordination numbers according to Grimme's recipe.
        k_1 = 16
        k_2 = 4 / 3
        k_3 = 4
        for cn_atom in atoms:
            # Get coordinates and radii of all other atoms
            coordinates = np.array([atom.coordinates for 
                                    atom in atoms if atom is not cn_atom])
            radii = np.array([atom.radius for
                              atom in atoms if atom is not cn_atom])

            # Calculate distances
            dists = cdist(np.array(cn_atom.coordinates).reshape(-1, 3),
                          coordinates.reshape(-1, 3))
            
            # Calculate coordination number
            coordination_number = np.sum(1 / (1 + np.exp(-k_1 * 
                (k_2 * (cn_atom.radius + radii) / dists - 1))))
            cn_atom.coordination_number = coordination_number
        
        # Load the reference data
        data_file = pkg_resources.resource_filename('steriplus',
            '../data/c6_reference_data.pickle')
        with open(data_file, "rb") as file:
            c6_reference_data = pickle.load(file)
        
        # Calculate the C6 coefficients
        c6_coefficients = []
        c8_coefficients = []
        for atom in atoms:
            # Get the reference data for atom
            reference_data = c6_reference_data[atom.element]
            cn = atom.coordination_number
            
            # Take out the coordination numbers and c6(aa) values
            c6_ref = reference_data[:,0]
            cn_1 = reference_data[:,1]
            cn_2 = reference_data[:,2]

            # Calculate c6  and c8 according to the recipe
            r = (cn - cn_1)**2 + (cn - cn_2)**2
            L = np.exp(-k_3 * r)
            W = np.sum(L)
            Z = np.sum(c6_ref * L)
            c6 = Z / W
            c8 = 3 * c6 * r2_r4[atom.element]**2
            c6_coefficients.append(c6)
            c8_coefficients.append(c8)

        # Set up attributes
        self._atoms = atoms
        self.c6_coefficients = c6_coefficients
        self.c8_coefficients = c8_coefficients

        def __repr__(self):
            return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"

def check_distances(elements, coordinates, check_atom, radii=[], check_radius=0,
                    exclude_list=[], epsilon=0, radii_type="crc"):
    """
    Args:
        elements (list): Elements as atomic symbols or numbers
        check_atom (int): Index of atom to check against (starting from 1)
        check_radius (float): Radius to use for check_atom (Å)
        coordinates (list): Coordinates (Å)
        epsilon (float): Numeric term to adjust distance check (Å). Positive
                         values add to the radii to make the test more
                         discriminatory.
        exclude_list (list): Indices of atoms to exclude from check (starting at
                             1)
        radii (list): vdW radii (Å)
        radii_type (str): Type of radii to use: 'bondi' or 'crc' (default)
    
    Returns:
        within_list (list): List of atoms within vdW distance of check atom.
                            Returns none if list is empty.
    """
    # Convert elements to atomic numbers if the are symbols
    elements = convert_elements(elements)

    # Get radii if they are not supplied
    if not radii:
        radii = get_radii(elements, radii_type=radii_type)
    radii = np.array(radii)

    atom_coordinates = np.array(coordinates)
    check_coordinates = np.array(coordinates[check_atom - 1]).reshape(-1, 3)

    # Calculate distances between check atom and all atoms
    distances = cdist(atom_coordinates, check_coordinates) - \
        radii.reshape(-1, 1) - check_radius - epsilon
    distances = distances.reshape(-1)
    
    # Determine atoms which are within a vdW distance from the check atom
    within_distance = list(np.argwhere(distances < 0).reshape(-1))
    
    # Remove check atom and atoms in the exclude list
    within_distance.remove(check_atom - 1)
    within_distance = [i for i in within_distance if i not in exclude_list]

    # Return atoms which are within vdW distance from check atom
    if within_distance:
        within_list = [i + 1 for i in within_distance]
        return within_list
    else:
        return None

def convert_elements(elements, output='numbers'):
    """Converts elements to atomic symbols or numbers.

    Args:
        elements (list): Elements as atomic symbols or numbers
        output (str): Output format: 'numbers' (default) or 'symbols'.
    
    Returns:
        elements (list): Converted elements
    """
    try:
        elements = [element.capitalize() for element in elements]
    except AttributeError:
        pass

    if output == "numbers":
        try:
            elements = [atomic_numbers[element] for 
                                  element in elements]
        except KeyError:
            pass
    if output == "symbols":
        try:
            elements = [atomic_symbols[element] for 
                                  element in elements]
        except KeyError:
            pass

    return elements

def get_radii(elements, radii_type="crc", scale=1):
    """Gets radii from element identifiers

    Args:
        elements (list): Elements as atomic symbols or numbers
        radii_type (str): Type of radius: 'bondi', 'crc' (default), 'rahm' or
                          'pyykko'

    Returns:
        radii (list): vdW radii (Å)
    """
    elements = convert_elements(elements)

    # Set up dictionary of radii types
    radii_choice = {'bondi': radii_bondi,
                    'crc': radii_crc,
                    'rahm': radii_rahm,
                    'pyykko': cov_radii_pyykko}
    
    # Get the radii. Replace with 2.0 if it the radius doesn't exist.
    radii = [radii_choice[radii_type].get(element, 2.0) * scale for element in elements]

    return radii
