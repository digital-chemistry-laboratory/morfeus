"""Internal calculators."""

import pickle
import pkg_resources

import numpy as np
from scipy.spatial.distance import cdist

from morfeus.helpers import convert_elements, get_radii, conditional
from morfeus.geometry import Atom
from morfeus.data import r2_r4, HARTREE, BOHR, EV, ANGSTROM

# Matplotlib is required for plotting steric maps for buried volumes
try:
    from dftd4.calculators import D4_model, D3_model
    from dftd4.utils import extrapolate_c_n_coeff
    import ase.io
except ImportError:
    _has_dftd4 = False
else:
    _has_dftd4 = True
_warning_dftd4 = "Install dftd4 and ase python packages for this function."


@conditional(_has_dftd4, _warning_dftd4)
class D3Grimme:
    """Calculates Cᴬᴬ dispersion coefficients with the dftd4 program and a
    D3-like method.

    Args:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
        order (int): Maximum order for the CN coefficients.

    Attributes:
        c_n_coefficients (dict): Cᴬᴬ coefficients (a.u.)
        polarizabilities (ndarray): Atomic polarizabilities (a.u.)
    """    
    def __init__(self, elements, coordinates, order=6):
        #if order %

        # Convert elements to atomic numbers
        elements = convert_elements(elements)

        # Do calculation 
        atoms = ase.Atoms(numbers=elements, positions=coordinates)
        calc = D3_model(energy=False, forces=False)
        polarizabilities = calc.get_property("polarizibilities", atoms=atoms) * (ANGSTROM / BOHR) ** 3
        c6_coefficients_all = calc.get_property('c6_coefficients', atoms=atoms) * EV / HARTREE * (ANGSTROM / BOHR) ** 6
        c6_coefficients = np.diag(c6_coefficients_all)
        c_n_coefficients = {}
        c_n_coefficients[6] = c6_coefficients

        # Extrapolate
        for i in range(8, order + 1, 2):
            c_n_coefficients[i] = np.array([extrapolate_c_n_coeff(c6, element, element, i) for c6, element in zip(c6_coefficients, elements)])

        # Store attributes
        self.polarizabilities = polarizabilities
        self.c_n_coefficients = c_n_coefficients
        self._atoms = atoms 

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"

@conditional(_has_dftd4, _warning_dftd4)
class D4Grimme:
    """Calculates Cᴬᴬ dispersion coefficients with the D4 method and the dftd4 program.
    
    Args:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
        order (int): Maximum order for the CN coefficients.

    Attributes:
        c_n_coefficients (dict): Cᴬᴬ coefficients (a.u.)
        charges (ndarray): Atomic charges.
        polarizabilities (ndarray): Atomic polarizabilities (a.u.)
    """        
    def __init__(self, elements, coordinates, order=8, charge=0):
        # Convert elements to atomic numbers
        elements = convert_elements(elements)
        
        # Set up atoms object
        charges = np.zeros(len(elements))
        charges[0] = charge
        atoms = atoms = ase.Atoms(numbers=elements, positions=coordinates, charges=charges)

        # Do calculation
        calc = D4_model(energy=False, forces=False)
        polarizabilities = calc.get_property("polarizibilities", atoms=atoms) * (ANGSTROM / BOHR) ** 3
        charges = calc.get_property("charges", atoms=atoms)
        c6_coefficients_all = calc.get_property('c6_coefficients', atoms=atoms) * EV / HARTREE * (ANGSTROM / BOHR) ** 6
        c6_coefficients = np.diag(c6_coefficients_all)
        c_n_coefficients = {}
        c_n_coefficients[6] = c6_coefficients

        # Extrapolate
        for i in range(8, order + 1, 2):
            c_n_coefficients[i] = np.array([extrapolate_c_n_coeff(c6, element, element, i) for c6, element in zip(c6_coefficients, elements)])

        # Store attributes
        self.polarizabilities = polarizabilities
        self.charges = charges        
        self.c_n_coefficients = c_n_coefficients
        self._atoms = atoms 

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"

class D3Calculator:
    """Calculates Cᴬᴬ dispersion coefficients with the D3 method based on
    the procedure in J. Chem. Phys. 2010, 132, 154104.

    Args:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
        order (int): Maximum order for the CN coefficients.

    Attributes:
        c_n_coefficients (dict): Cᴬᴬ coefficients (a.u.)
        coordination_numbers (ndarray): Atomic coordination numbers.
        polarizabilities (ndarray): Atomic polarizabilities (a.u.)
    """
    def __init__(self, elements, coordinates, order=8):
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
        
        # Calculate the C_N coefficients
        c_n_coefficients = {i: []  for i in range(6, order + 1, 2)}
        for atom in atoms:
            # Get the reference data for atom
            reference_data = c6_reference_data[atom.element]
            cn = atom.coordination_number
            
            # Take out the coordination numbers and c6(aa) values
            c_6_ref = reference_data[:,0]
            cn_1 = reference_data[:,1]
            cn_2 = reference_data[:,2]

            # Calculate c6  and c8 according to the recipe
            r = (cn - cn_1)**2 + (cn - cn_2)**2
            L = np.exp(-k_3 * r)
            W = np.sum(L)
            Z = np.sum(c_6_ref * L)
            c_6 = Z / W

            temp_coefficients = {i: None for i in range(6, order + 1, 2)}
            for i in range(6, order + 1, 2):
                if i == 6:
                    temp_coefficients[i] = c_6
                if i == 8:
                    c_8 = 3 * c_6 * r2_r4[atom.element]**2
                    temp_coefficients[i] = c_8
                elif i == 10:
                    c_10 = 49.0 / 40.0 * c_8 ** 2 / c_6
                    temp_coefficients[i] = c_10
                elif i > 10:
                    c_n = temp_coefficients[i - 6] * (temp_coefficients[i - 2] / temp_coefficients[i - 4]) ** 3 
                    temp_coefficients[i] = c_n
            for key, value in temp_coefficients.items():
                c_n_coefficients[key].append(value)

        # Set up attributes
        coordination_numbers = np.array([atom.coordination_number for atom in self._atoms])

        self._atoms = atoms
        self.c_n_coefficients = c_n_coefficients
        self.coordination_numbers = coordination_numbers

    def __repr__(self):
        return f"{self.__class__.__name__}({len(self._atoms)!r} atoms)"
