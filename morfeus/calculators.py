"""Internal calculators."""

import pickle
import pkg_resources

import numpy as np
from scipy.spatial.distance import cdist

from morfeus.helpers import convert_elements, get_radii
from morfeus.geometry import Atom
from morfeus.data import r2_r4

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
        data_file = pkg_resources.resource_filename('morfeus',
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
