"""Help functions."""

from dataclasses import dataclass
from importlib import import_module
from numbers import Integral
import shutil
from typing import (
    Any,
    Callable,
    cast,
    Iterable,
    List,
    Literal,
    Optional,
    overload,
    Sequence,
    Union,
)

import numpy as np
import scipy.spatial

from morfeus.data import (
    atomic_numbers,
    atomic_symbols,
    cov_radii_pyykko,
    radii_alvarez,
    radii_bondi,
    radii_crc,
    radii_rahm,
    radii_truhlar,
)
from morfeus.typing import Array2D, ArrayLike1D, ArrayLike2D


# TODO change exclude_list, within_list, return empty list instead of None
def check_distances(
    elements: Union[Iterable[int], Iterable[str]],
    coordinates: ArrayLike2D,
    check_atom: int,
    radii: Optional[ArrayLike1D] = None,
    check_radius: float = 0,
    exclude_list: Optional[Sequence[int]] = None,
    epsilon: float = 0,
    radii_type: str = "crc",
) -> Optional[List[int]]:
    """Check which atoms are within clashing vdW radii distances.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        check_atom: Index of atom to check against (1-indexed)
        radii: vdW radii (Å)
        check_radius: Radius to use for check_atom (Å)
        exclude_list: Atom indices to exclude (1-indexed)
        epsilon: Numeric term add to the radii (Å)
        radii_type: Radii type: 'bondi' or 'crc'

    Returns:
        within_list: Atom indices within vdW distance of check atom.
            Returns none if list is empty.
    """
    # Convert elements to atomic numbers if the are symbols
    elements = convert_elements(elements, output="numbers")

    # Get radii if they are not supplied
    if radii is None:
        radii = get_radii(elements, radii_type=radii_type)
    radii = np.array(radii)

    if exclude_list is None:
        exclude_list = []
    else:
        exclude_list = list(exclude_list)

    coordinates = np.array(coordinates)
    atom_coordinates = np.array(coordinates)
    check_coordinates = np.array(coordinates[check_atom - 1]).reshape(-1, 3)

    # Calculate distances between check atom and all atoms
    distances = (
        scipy.spatial.distance.cdist(atom_coordinates, check_coordinates)
        - radii.reshape(-1, 1)
        - check_radius
        - epsilon
    )
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


def requires_executable(executables: Sequence[str]) -> Callable[..., Callable]:
    """Decorator factory to control optional executables.

    Args:
        executables: Names of executables

    Returns:
        decorator: Either 'noop_decorator' that returns the original function or
            'error_decorator' that raises an OSError and lists absent executables.
    """

    def noop_decorator(function: Callable[..., Any]) -> Callable[..., Any]:
        """Returns function unchanged."""
        return function

    def error_decorator(function: Callable[..., Any]) -> Callable[..., Any]:
        """Raises error."""

        def error(*args, **kwargs) -> OSError:
            error_msg = "Required executables not found in path:"
            for exe_error in exe_errors:
                error_msg += f" {exe_error}"
            raise OSError(error_msg)

        return error

    # Try to find excetubles in path
    exe_errors = []
    for executable in executables:
        if not shutil.which(executable):
            exe_errors.append(executable)

    return error_decorator if len(exe_errors) > 0 else noop_decorator


@dataclass
class Import:
    """Class for handling optional dependency imports."""

    module: str
    item: Optional[str] = None
    alias: Optional[str] = None


def requires_dependency(  # noqa: C901
    imports: Sequence[Import], _globals: dict
) -> Callable[..., Callable]:
    """Decorator factory to control optional dependencies.

    Args:
        imports: Imports
        _globals: Global symbol table from calling module.

    Returns:
        decorator: Either 'noop_decorator' that returns the original function or
            'error_decorator' that raises an ImportError and lists absent dependencies.
    """

    def noop_decorator(function: Callable[..., Any]) -> Callable[..., Any]:
        """Returns function unchanged."""
        return function

    def error_decorator(function: Callable[..., Any]) -> Callable[..., Any]:
        """Raises error."""

        def error(*args, **kwargs) -> ImportError:
            error_msg = "Install extra requirements to use this function:"
            for e in import_errors:
                error_msg += f" {e.name}"
            raise ImportError(error_msg)

        return error

    import_errors = []
    for imp in imports:
        # Import module
        try:
            module = import_module(imp.module)

            # Try to import item as attribute
            if imp.item is not None:
                try:
                    item = getattr(module, imp.item)
                except AttributeError:
                    item = import_module(f"{imp.module}.{imp.item}")
                name = imp.item
            else:
                item = module
                name = imp.module

            # Convert item name to alias
            if imp.alias is not None:
                name = imp.alias

            _globals[name] = item
        except ImportError as import_error:
            import_errors.append(import_error)

    return error_decorator if len(import_errors) > 0 else noop_decorator


@overload
def convert_elements(
    elements: Union[Iterable[int], Iterable[str]], output: Literal["numbers"]
) -> List[int]:
    ...


@overload
def convert_elements(
    elements: Union[Iterable[int], Iterable[str]], output: Literal["symbols"]
) -> List[str]:
    ...


def convert_elements(
    elements: Union[Iterable[int], Iterable[str]], output: str = "numbers"
) -> Union[List[int], List[str]]:
    """Converts elements to atomic symbols or numbers.

    Args:
        elements: Elements as atomic symbols or numbers
        output: Output format: 'numbers' (default) or 'symbols'.

    Returns:
        elements: Converted elements

    Raises:
        TypeError: When input type not supported
        ValueError: When output not supported
    """
    if output not in ["numbers", "symbols"]:
        raise ValueError(f"ouput={output} not supported. Use 'numbers' or 'symbols'")

    if all(isinstance(element, str) for element in elements):
        elements = cast(List[str], elements)
        if output == "numbers":
            elements = [atomic_numbers[element.capitalize()] for element in elements]
        return elements
    elif all(isinstance(element, Integral) for element in elements):
        elements = cast(List[int], elements)
        if output == "symbols":
            elements = [atomic_symbols[element] for element in elements]
        return elements
    else:
        raise TypeError("elements must be all integers or all strings.")


def get_radii(
    elements: Union[Iterable[int], Iterable[str]],
    radii_type: str = "crc",
    scale: float = 1,
) -> List[float]:
    """Gets radii from element identifiers.

    Args:
        elements: Elements as atomic symbols or numbers
        radii_type: Radii type: 'bondi', 'crc', 'pyykko' 'rahm' or 'truhlar'
        scale: Scaling factor

    Returns:
        radii: Radii (Å)
    """
    elements = convert_elements(elements, output="numbers")

    # Set up dictionary of radii types
    radii_choice = {
        "alvarez": radii_alvarez,
        "bondi": radii_bondi,
        "crc": radii_crc,
        "rahm": radii_rahm,
        "pyykko": cov_radii_pyykko,
        "truhlar": radii_truhlar,
    }

    # Get the radii. Replace with 2.0 if it the radius doesn't exist.
    radii = [radii_choice[radii_type].get(element, 2.0) * scale for element in elements]

    return radii


def get_connectivity_matrix(
    coordinates: ArrayLike2D,
    elements: Optional[Union[Iterable[int], Iterable[str]]] = None,
    radii: Optional[ArrayLike1D] = None,
    radii_type: str = "pyykko",
    scale_factor: float = 1.2,
) -> Array2D:
    """Get connectivity matrix from covalent radii.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        radii: Radii (Å)
        radii_type: Radii type: 'pyykko'
        scale_factor: Factor for scaling covalent radii

    Returns:
        connectivity_matrix: Connectivity matrix

    Raises:
        RuntimeError: When neither elements nor radii given
    """
    coordinates = np.array(coordinates)
    n_atoms = len(coordinates)
    if radii is None:
        if elements is None:
            raise RuntimeError("Either elements or radii needed.")
        elements = convert_elements(elements, output="numbers")
        radii = get_radii(elements, radii_type=radii_type)
    radii = np.array(radii)
    distance_matrix = scipy.spatial.distance_matrix(coordinates, coordinates)
    radii_matrix = np.add.outer(radii, radii) * scale_factor
    connectivity_matrix: np.ndarray = (distance_matrix < radii_matrix) - np.identity(
        n_atoms
    ).astype(int)

    return connectivity_matrix