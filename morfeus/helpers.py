"""Help functions."""

import shutil
from importlib import import_module

import numpy as np
import scipy.spatial
from scipy.spatial.distance import cdist

from morfeus.data import (atomic_numbers, atomic_symbols, cov_radii_pyykko,
                          radii_alvarez, radii_bondi, radii_crc, radii_rahm,
                          radii_truhlar)


def check_distances(elements,
                    coordinates,
                    check_atom,
                    radii=None,
                    check_radius=0,
                    exclude_list=None,
                    epsilon=0,
                    radii_type="crc"):
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
    if radii is None:
        radii = get_radii(elements, radii_type=radii_type)
    else:
        radii = np.array(radii)

    if exclude_list is None:
        exclude_list
    else:
        exclude_list = list(exclude_list)

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


def requires_executable(executables):
    """Decorator factory to control optional executables.

    Args:
        executables (list): Names of executables

    Returns:
        decorator (function): Either 'noop_decorator' returns the original
            function or 'error_decorator' which raises an Exception and lists
            absent executables.
    """
    def noop_decorator(function):
        """Returns function unchanged."""
        return function

    def error_decorator(function):
        """Raises error"""
        def error(*args, **kwargs):
            error_msg = "Required executables not found in path:"
            for exe_error in exe_errors:
                error_msg += f" {exe_error}"
            raise Exception(error_msg)
        return error

    # Try to find excetubles in path
    exe_errors = []
    for executable in executables:
        if not shutil.which(executable):
            exe_errors.append(executable)

    return error_decorator if len(exe_errors) > 0 else noop_decorator


def requires_dependency(modules, _globals):
    """Decorator factory to control optional dependencies.

    Args:
        modules (list): Modules as (module, name) tuple pairs
        _globals (dict): Global symbol table from calling module.

    Returns:
        decorator (function): Either 'noop_decorator' returns the original
            function or 'error_decorator' which raises an ImportError and lists
            absent dependencies.
    """
    def noop_decorator(function):
        """Returns function unchanged."""
        return function

    def error_decorator(function):
        """Raises error"""
        def error(*args, **kwargs):
            error_msg = "Install extra requirements to use this function:"
            for e in import_errors:
                error_msg += f" {e.name}"
            raise ImportError(error_msg)
        return error

    # Try to import dependencies
    import_errors = []
    for module, name in modules:
        try:
            _globals[name] = import_module(module)
        except ImportError as import_error:
            import_errors.append(import_error)

    return error_decorator if len(import_errors) > 0 else noop_decorator


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
    radii = [radii_choice[radii_type].get(element, 2.0) * scale for element
        in elements]

    return radii

def get_connectivity_matrix(elements, coordinates, radii=[],
    radii_type="pyykko"):
    elements = convert_elements(elements)

    if len(radii) < 1:
        radii = get_radii(elements, radii_type="pyykko")
    distance_matrix = scipy.spatial.distance_matrix(coordinates, coordinates)
    radii_matrix = np.add.outer(radii, radii) * 1.2
    connectivity_matrix = (distance_matrix < radii_matrix) \
        - np.identity(len(elements))

    return connectivity_matrix
