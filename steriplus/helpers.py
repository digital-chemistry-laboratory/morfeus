"""Help functions

Functions:
    check_distances: Controls distances for steric clashes.
    convert_elements: Converts atomic symbols to atomic numbers.
    get_radii: Returns radii for list of elements.
"""
import numpy as np
from scipy.spatial.distance import cdist

from steriplus.data import atomic_numbers, bondi_radii, crc_radii

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
    distances = cdist(atom_coordinates, check_coordinates) - radii.reshape(-1, 1) - check_radius - epsilon
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

def convert_elements(elements):
    """Converts elements to atomic numbers.

    Args:
        elements (list): Elements as atomic symbols or numbers
    
    Returns:
        converted_elements (list): Elements as atomic numbers
    """
    if type(elements[0]) == str:
        cap_elements = [element.capitalize() for element in elements]
        converted_elements = [atomic_numbers[element] for 
                              element in cap_elements]
    else:
        converted_elements = elements

    return converted_elements

def get_radii(elements, radii_type="crc", scale=1):
    """Gets radii from element identifiers

    Args:
        elements (list): Elements as atomic symbols or numbers
        radii_type (str): Type of vdW radius, 'bondi' or 'crc' (default)

    Returns:
        radii (list): vdW radii (Å)
    """
    elements = convert_elements(elements)

    # Set up dictionary of atomic radii for all elements present in the list
    element_set = set(elements)
    if radii_type == "bondi":
        radii_dict = {element: bondi_radii.get(element) for element in element_set}
    elif radii_type == "crc":
        radii_dict = {element: crc_radii.get(element) for element in element_set}

    # Get radii for elements in the element id list. Set 2 as default if does not exist
    radii = [radii_dict.get(element) * scale if radii_dict.get(element, 2.0) else 2.0 * scale for element in elements]

    return radii
