"""Interface to quantum-chemical programs."""

import numpy as np

from morfeus.data import ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM
from morfeus.helpers import conditional, convert_elements

try:
    import qcelemental as qcel
    import qcengine as qcng
    _HAS_QCNG = True
except ImportError:
    _HAS_QCNG = False
_MSG_QCNG = "Install QCEngine to use this function."


@conditional(_HAS_QCNG, _MSG_QCNG)    
def optimize_qc_engine(
        elements, coordinates, charge=None, multiplicity=None,
        connectivity_matrix=None, program="xtb", model=None, keywords=None,
        local_options=None, procedure="berny", return_trajectory=False):
    """Optimize molecule with QCEngine.

    Args:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
        charge (int): Molecular charge
        multiplicity (int): Multiplicity
        connectivity_matrix (ndarray): Connectivity matrix
        program (str): QCEngine program
        model (dict): QCEngine model
        keywords (dict): QCEngine keywords
        local_options (dict): QCEngine local options
        procedure (str): QCEngine procedure
        return_trajectory (bool): Return coordinates for all steps

    Returns:
        opt_coordinates (ndarray): Conformer coordinates (Å)
        energies (ndarray): Energies for all steps (a.u.)
    """
    # Set defaults
    if model is None:
        model = {"method": "GFN2-xTB"}
    if keywords is None:
        keywords = {}
    if local_options is None:
        local_options = {}
    
    # Create molecule object
    molecule = _generate_qcel_molecule(elements, coordinates, charge,
        multiplicity, connectivity_matrix)
    
    # Create optimization input
    opt_input = {
        "keywords": {
            "program": program
        },
        "input_specification": {
            "driver": "gradient",
            "model": model,
            "keywords": keywords,
        },
        "initial_molecule": molecule
    }
    
    # Perform optimization
    opt = qcng.compute_procedure(opt_input, procedure=procedure,
        local_options=local_options)
    if not opt.success:
        raise Exception(opt.error.error_message)

    # Take out results
    energies = np.array(opt.energies)
    if return_trajectory:
        opt_coordinates = np.array([result.molecule.geometry 
            for result in opt.trajectory])
    else:
        opt_coordinates = opt.final_molecule.geometry
    opt_coordinates *= BOHR_TO_ANGSTROM
    
    return opt_coordinates, energies


@conditional(_HAS_QCNG, _MSG_QCNG)    
def sp_qc_engine(
        elements, coordinates, charge=None, multiplicity=None,
        connectivity_matrix=None, program="xtb", model=None, keywords=None,
        local_options=None):
    """Single-point calculation with QCEngine.

    Args:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
        charge (int): Molecular charge
        multiplicity (int); Molecular multiplicity
        connectivity_matrix (ndarray): Connectivity matrix
        program (str): QCEngine program
        model (dict): QCEngine model
        keywords (dict): QCEngine keywords
        local_options (dict): QCEngine local options

    Returns:
        energy (float): Energy (a.u.)
    """
    # Set defaults
    if model is None:
        model = {"method": "GFN2-xTB"}
    if keywords is None:
        keywords = {}
    if local_options is None:
        local_options = {}
    
    # Crate molecule object
    molecule = _generate_qcel_molecule(elements, coordinates, charge,
        multiplicity, connectivity_matrix)
    
    # Create sp input
    sp_input = qcel.models.AtomicInput(
        molecule=molecule,
        driver="energy",
        model=model,
        keywords=keywords,
    )    

    # Perform sp calculation
    sp = qcng.compute(sp_input, program=program, local_options=local_options)
    if not sp.success:
        raise Exception(sp.error.error_message)

    # Take out results
    energy = sp.return_result

    return energy


@conditional(_HAS_QCNG, _MSG_QCNG)    
def _generate_qcel_molecule(elements, coordinates, charge=None,
    multiplicity=None, connectivity_matrix=None):
    """Generate QCElemental molecule object.

    Args:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
        charge (int): Molecular charge
        multiplicity (int): Molecular multiplicity
        connectivity_matrix (ndarray): Connectivity matrix
    
    Returns:
        molecule (obj): QCElemental molecule object.
    """
    # Generate bond order list from connectivity matrix
    if connectivity_matrix is not None:
        bos = []
        i, j = np.tril_indices_from(connectivity_matrix)
        for k, l in zip(i, j):
            if k != l:
                bo = int(connectivity_matrix[k, l])
                if bo != 0:
                    bos.append((k, l, bo))
    else:
        bos = None
    
    # Create molecule object
    elements = np.array(convert_elements(elements, output="symbols"))
    coordinates = np.array(coordinates) * ANGSTROM_TO_BOHR
    molecule = qcel.models.Molecule(symbols=elements,
        geometry=coordinates, molecular_charge=charge, connectivity=bos,
        molecular_multiplicity=multiplicity)

    return molecule
