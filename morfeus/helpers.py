"""Help functions."""

import warnings

import numpy as np
from scipy.spatial.distance import cdist
import scipy.spatial

from morfeus.data import atomic_numbers, atomic_symbols
from morfeus.data import (radii_alvarez, radii_bondi, radii_crc, radii_rahm, 
    radii_truhlar, cov_radii_pyykko)

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

def conditional(condition, warning=None):
    """Decorator factory to control if functions are presented.

    Args:
        condition (bool): Whether to import function or not  
        warning (warning): Warning to print if function not imported
    
    Returns:
        noop_decorator (obj): Decorator that passes the original function
        neutered_function (obj): Decorator that passes nothing and prints warning.

    """
    def noop_decorator(function):
        """Returns function unchanged."""
        return function

    def neutered_function(function):
        """Returns function that just prints warning"""
        def neutered(*args, **kw):
            if warning:
                warnings.warn(warning)
            return
        return neutered

    return noop_decorator if condition else neutered_function

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
        radii_type (str): Type of radius: 'bondi', 'crc' (default), 'rahm',
                          'truhlar' or 'pyykko'

    Returns:
        radii (list): vdW radii (Å)
    """
    elements = convert_elements(elements)

    # Set up dictionary of radii types
    radii_choice = {'alvarez': radii_alvarez,
                    'bondi': radii_bondi,
                    'crc': radii_crc,
                    'rahm': radii_rahm,
                    'pyykko': cov_radii_pyykko,
                    'truhlar': radii_truhlar}
    
    # Get the radii. Replace with 2.0 if it the radius doesn't exist.
    radii = [radii_choice[radii_type].get(element, 2.0) * scale for element in elements]

    return radii

def get_connectivity_matrix(elements, coordinates, radii=[], radii_type="pyykko"):
    elements = convert_elements(elements)

    if len(radii) < 1:
        radii = get_radii(elements, radii_type="pyykko")
    distance_matrix = scipy.spatial.distance_matrix(coordinates, coordinates)
    radii_matrix = np.add.outer(radii, radii) * 1.2
    connectivity_matrix = (distance_matrix < radii_matrix)  - np.identity(len(elements))

    return connectivity_matrix

