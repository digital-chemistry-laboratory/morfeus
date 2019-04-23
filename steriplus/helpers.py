"""Module containing helper functions

Functions:
    check_distances: Controls distances for steric clashes.
    convert_elements: Converts elements to atomic numbers.
    get_radii: Helper function to get radii from list of elements.
"""

from steriplus.data import atomic_numbers, atomic_symbols, \
    bondi_radii, crc_radii, jmol_colors
import numpy as np
from scipy.spatial.distance import cdist

def check_distances(element_ids, coordinates, dummy_atom, exclude_list=[], dummy_radius=0, epsilon=0, radii=[], radii_type="crc"):
    # Converting element ids to atomic numbers if the are symbols
    if type(element_ids[0]) == str:
        element_ids = [element.capitalize() for element in element_ids]
        element_ids = [atomic_numbers[element_id] for element_id in element_ids]
    
    # Getting radii if they are not supplied
    if not radii:
        radii = get_radii(element_ids, radii_type=radii_type)
    
    radii = np.array(radii)
    atom_coordinates = np.array(coordinates)
    dummy_coordinates = np.array(coordinates[dummy_atom - 1]).reshape(-1, 3)

    distances = cdist(atom_coordinates, dummy_coordinates) - radii.reshape(-1, 1) - dummy_radius - epsilon
    distances = distances.reshape(-1)

    within_distance = list(np.argwhere(distances < 0).reshape(-1))
    within_distance.remove(dummy_atom - 1)
    within_distance = [i for i in within_distance if i not in exclude_list]
    
    if any(within_distance):
        within_list = [i + 1 for i in within_distance]
        return within_list
    else:
        return None

def convert_element_ids(element_ids):
    """Converts elements to list of atomic numbers.

    Args:
        element_ids (list): List of element identifiers
    """
    if type(element_ids[0]) == str:
        element_ids = [element.capitalize() for element in element_ids]
        element_ids = [atomic_numbers[element_id] for
                       element_id in element_ids]
    return element_ids

def get_radii(element_ids, radii_type="crc", scale=1):
    """Gets radii for list of element ids

    Args:
        element_id_list (list)  :   List of element ids as atomic numbers or
                                    symbols
        radii_type (str)        :   Type of vdW radius, "bondi" or "crc"

    Returns:
        radii_list (list)       :   List of atomic radii
    """
    element_ids = convert_element_ids(element_ids)

    # Set up dictionary of atomic radii for all elements present in the list
    element_set = set(element_ids)
    if radii_type == "bondi":
        radii_dict = {element_id: bondi_radii.get(element_id) for element_id in element_set}
    elif radii_type == "crc":
        radii_dict = {element_id: crc_radii.get(element_id) for element_id in element_set}

    # Get radii for elements in the element id list. Set 2 as default if does not exist
    radii = [radii_dict.get(element_id) * scale if radii_dict.get(element_id, 2.0) else 2.0 * scale for element_id in element_ids]

    return radii
