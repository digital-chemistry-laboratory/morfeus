import numpy as np

from morfeus.data import (HARTREE, K_B, KCAL_TO_HARTREE, HARTREE_TO_KCAL,
    HARTREE_TO_KJ, KJ_TO_HARTREE)
from morfeus.helpers import conditional, convert_elements
from morfeus.io import get_xyz_string
from morfeus.qc import optimize_qc_engine, sp_qc_engine, _generate_qcel_molecule

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _has_rdkit = True
except ImportError:
    _has_rdkit = False
_warning_rdkit = "Install RDKit to use this function."

try:
    import spyrmsd
    from spyrmsd import rmsd
    _has_spyrmsd = True
except ImportError:
    _has_spyrmsd = False
_warning_spyrmsd = "Install spyrmsd to use this function."    

class Conformer:
    """Conformer with coordinates, energy and properties.

    Args:
        coordinates (list): Coordinates (Å)
        energy (float): Energy (a.u.)
        properties (dict): Conformers properties.
    
    Attributes:
        coordinates (list): Coordinates (Å)
        energy (float): Energy (a.u.)
        properties (dict): Conformers properties.    
    """
    def __init__(self, coordinates, energy=None, properties=None):
        if properties is None:
            properties = {}
        self.coordinates = np.array(coordinates)
        self.energy = energy
        self.properties = properties

    def __repr__(self):
        n_atoms = self.coordinates.shape[0]
        return f"{self.__class__.__name__}({n_atoms!r} atoms)"

class ConformerEnsemble:
    """Conformer ensemble object that supports sorting, pruning, optimization
    and single-point calculations.

    Args:
        conformer_coordinates (list): Conformer coordinates (Å)
        connectivity_matrix (ndarray): Connectivity matrix
        elements (list): Elements as atomic symbols or numbers
        energies (list): Energy (a.u.)
        properties (dict): Conformers properties.
        charge (int): Molecular charge
        multiplicity (int): Molecular multiplicity.

    Attributes:
        charge (int): Molecular charge
        conformers (list): Conformers
        connectivity_matrix (ndarray): Connectivity matrix
        elements (list): Elements as atomic symbols or numbers
        multiplicity (int): Molecular multiplicity.
    """
    def __init__(self, elements, conformer_coordinates=None, energies=None,
        connectivity_matrix=None, properties=None, charge=None,
        multiplicity=None):
        if conformer_coordinates is None:
            conformer_coordinates = []
            
        # Store conformers
        self.conformers = []
        self._add_conformers(conformer_coordinates, energies, properties)
        
        # Set up attributes
        self.elements = np.array(convert_elements(elements))
        self.charge = charge
        self.multiplicity = multiplicity
        self.connectivity_matrix = connectivity_matrix
        
    def add_conformers(self, coordinates, properties=None):
        """Add conformer to ensemble.

        Args:
            coordinates (list): Conformer coordinates (Å)
            properties (dict): Conformer properties.
        """
        self._add_conformers(coordinates, properties)
    
    def _add_conformers(self, conformer_coordinates, energies=None,
        properties=None):
        if properties is  None:
            properties = {}
        if energies is  None:
            energies = [None] * len(conformer_coordinates)
        for coordinates, energy in zip(conformer_coordinates, energies):
            conformer = Conformer(coordinates, energy)
            self.conformers.append(conformer)
        for key, value in properties.items():
            self.set_properties(key, value)
    
    def boltzmann_statistic(self, property_name, temperature=298.15,
        statistic="avg"):
        """Calculate Boltzmann staistic of property over ensemble.

        Args:
            property_name (str): Name of property
            statistic (str): Boltzmann statistic: 'avg', 'var' or 'std'
            temperature (float): Temperature (K)
        
        Returns:
            statistic (float): Boltzmann statistic
        """
        energies = self.get_energies()
        properties = self.get_properties()[property_name]
        statistic = boltzmann_statistic(properties, energies,
            temperature=temperature, statistic=statistic)
        
        return statistic
    
    def boltzmann_average_dT(self, property_name, temperature=298.15):
        """Calculate temperature derivative of Boltzmann average of property.
        
        Args:
            property_name (str): Name of property
            temperature (float): Temperature (K)
        
        Returns:
            derivative (float): Derivative of Boltzmann average.
        """
        energies = self.get_energies()
        properties = self.get_properties()[property_name]
        derivative = boltzmann_average_dT(properties, energies,
            temperature=temperature)
        
        return derivative
    
    def boltzmann_weights(self, temperature=298.15):
        """Calculate Boltzmann weights for ensemble.

        Args:
            temperature (float): Temperature (K)

        Returns:
            weights (ndarray): Conformer weights (normalized to unity)
        """
        energies = self.get_energies()
        weights = boltzmann_weights(energies, temperature)
        
        return weights

    def get_properties(self):
        """Get conformer properties.
        
        Returns:
            properties (dict): Conformer properties.
        """
        properties = {}
        for conformer in self.conformers:
            for key, value in conformer.properties.items():
                properties.setdefault(key, []).append(value)
        
        for key, value in properties.items():
            properties[key] = np.array(value)
        
        return properties
    
    def set_properties(self, property_name, values):
        """Set conformer properties.

        Args:
            property_name (str): Name of property
            values (list): Property values
        """
        for conformer, value in zip(self.conformers, values):
            conformer.properties[property_name] = value
    
    def get_energies(self, unit="hartree", relative=False):
        """Get conformer energies.
        Args:
            relative (bool): Return energies relative to lowest-energy
                conformer.
            unit (str): Unit of returned energies: 'hartree' (default),
                'kcal/mol' or 'kJ/mol'.

        Returns:
            energies (ndarray): Conformer energies.
        """
        energies = np.array([conformer.energy for conformer in
            self.conformers])

        # Adjust unit
        if unit.lower() == "hartree":
            pass
        if unit.lower() == "kcal/mol":
            energies *= HARTREE_TO_KCAL
        elif unit.lower() == "kj/mol":
            energies *= HARTREE_TO_KJ
        
        if relative:
            energies -= np.min(energies)
        
        return energies

    def set_energies(self, energies):
        """Set conformer energies.

        Args:
            energies (list): Conformer energies (a.u.)
        """
        for conformer, energy in zip(self.conformers, energies):
            conformer.energy = energy
    
    def get_coordinates(self):
        """Get conformer coordinates.

        Returns:
            conformer_coordinates (ndarray): Conformer coordinates (Å)
        """
        conformer_coordinates = np.array([conformer.coordinates for
            conformer in self.conformers])
        
        return conformer_coordinates

    def __repr__(self):
        n_conformers = len(self.conformers)
        return f"{self.__class__.__name__}({n_conformers!r} conformers)"


    @conditional(_has_spyrmsd, _warning_spyrmsd)    
    def get_rmsd(self, i, j, include_hs=False, symmetry=True):
        """Get RSMD between two conformers.
        
        Args:
            i (int): Index of conformer 1
            include_hs (bool): Include H atoms in RMSD calculation.
            j (int): Index of conformer 2
            symmetry (bool): Consider symmetry (requires connectivity matrix)

        Returns:
            rmsd (float): RSMD (Å)
        """
        # Construct mask for H atoms
        if not include_hs:
            mask = self.elements != 1
        else:
            mask = np.ones(len(self.elements), dtype=np.bool)
        
        # Pick out conformer
        conformer_1 = self.conformers[i - 1]
        conformer_2 = self.conformers[j - 1]
        
        # Calculate RMSD. Use spyrmsd with no. heavy atoms < 3 due to unstable
        # performance of spyrmsd.
        if np.sum(self.elements != 1) < 3:
            molecule_1 = _generate_qcel_molecule(self.elements[mask], conformer_1.coordinates[mask])
            molecule_2 = _generate_qcel_molecule(self.elements[mask], conformer_2.coordinates[mask])
            rmsd = molecule_1.align(molecule_2)[1]["rmsd"]
        else:            
            if symmetry:              
                rmsd = spyrmsd.rmsd.symmrmsd(
                    conformer_1.coordinates[mask],
                    conformer_2.coordinates[mask],
                    self.elements[mask],
                    self.elements[mask],
                    self.connectivity_matrix[mask,:][:,mask],
                    self.connectivity_matrix[mask,:][:,mask],
                    center=True, minimize=True
                )
            else:
                rmsd = spyrmsd.rmsd.rmsd(
                    conformer_1.coordinates[mask],
                    conformer_2.coordinates[mask],
                    self.elements[mask],
                    self.elements[mask],
                    center=True,
                    minimize=True
                )
        
        return rmsd

    @conditional(_has_spyrmsd, _warning_spyrmsd) 
    def prune_rmsd(self, prune_rmsd=0.35, include_hs=False, symmetry=True):
        """Prune conformers based on RMSD.
        
        Args:
            include_hs (bool): Include H atoms in RMSD calculation
            prune_rmsd (float): Threshold for RSMD pruning (Å)
            symmetry (bool): Consider symmetry (requires connectivity matrix)
        """
        # Select conformers to keep
        keep_list = []
        for i in range(len(self.conformers)):
            keep = True
            for j in keep_list:
                rmsd = self.get_rmsd(i + 1, j + 1, include_hs=include_hs,
                    symmetry=symmetry)
                if rmsd < prune_rmsd:
                    keep = False
                    break
            if keep:
                keep_list.append(i)

        # Update conformer list
        self.conformers = [conformer for i, conformer in
            enumerate(self.conformers) if i in keep_list]
    
    def prune_energy(self, threshold=3.0, unit="kcal/mol"):
        """Prune conformers based on energy compared to minimum energy
        conformer.

        Args:
            threshold (float): Energy threshold for pruning. 
            unit (str): Unit for input energy threshold 'hartree', 'kcal/mol'
                (default) or 'kJ/mol.
        """
        # Convert threshold to Hartree.
        if unit.lower() == "hartree":
            pass
        if unit.lower() == "kcal/mol":
            threshold *= KCAL_TO_HARTREE
        elif unit.lower() == "kj/mol":
            threshold *= KJ_TO_HARTREE
        
        # Prune conformers.
        energies = self.get_energies()
        energies -= np.min(energies)
        remove_list = np.where(energies > threshold)[0]
        
        for i in reversed(remove_list):
            del self.conformers[i]    
            
    def sort(self):
        """Sort conformers based on energy."""
        self.conformers = sorted(self.conformers, key=lambda x: x.energy)
    
    def sp_qc_engine(self, ids=None, program="xtb", model=None, keywords=None,
        local_options=None):
        """Calculate conformer energies with QCEngine interface.

        Args:
            ids (list): Conformer indices to optimize. If None, all are 
                optimized. 1-index.
            keywords (dict): QCEngine keywords
            local_options (dict): QCEngine local options
            model (dict): QCEngine model
            procedure (str): QCEngine procedure
            program (str): QCEngine program
        """
        # Set defaults
        if model is None:
            model = {"method": "GFN2-xTB"}
        if keywords is None:
            keywords = {}
        if local_options is None:
            local_options = {}
        if ids is not None:
            conformers = [self.conformers[i] for i in ids]
        else:
            conformers = self.conformers

        for conformer in conformers:
            energy = sp_qc_engine(
                self.elements,
                conformer.coordinates,
                charge=self.charge,
                multiplicity=self.multiplicity,
                connectivity_matrix=self.connectivity_matrix,
                program=program,
                model=model,
                keywords=keywords,
                local_options=local_options,
            )
            conformer.energy = energy

    def optimize_qc_engine(self, ids=None, program="xtb", model=None,
        keywords=None, local_options=None, procedure="berny"):
        """Optimize conformers with QCEngine interface.

        Args:
            ids (list): Conformer indices to optimize. If None, all are 
                optimized. 1-index.
            keywords (dict): QCEngine keywords
            local_options (dict): QCEngine local options
            model (dict): QCEngine model
            procedure (str): QCEngine procedure
            program (str): QCEngine program
        """
        # Set defaults
        if model is None:
            model = {"method": "GFN2-xTB"}
        if keywords is None:
            keywords = {}
        if local_options is None:
            local_options = {}
        if ids is not None:
            conformers = [self.conformers[i] for i in ids]
        else:
            conformers = self.conformers
        
        for conformer in conformers:
            opt_coordinates, energies = optimize_qc_engine(
                self.elements,
                conformer.coordinates,
                charge=self.charge,
                multiplicity=self.multiplicity,
                connectivity_matrix=self.connectivity_matrix,
                program=program,
                model=model,
                keywords=keywords,
                local_options=local_options,
                procedure=procedure,
            )
            conformer.coordinates = opt_coordinates
            conformer.energy = energies[-1]
        
    def write_conformers(self, filename, unit="kcal/mol", relative=True,
        separate=False):
        """Write conformers to xyz file.

        Args:
            filename (str): Filename
            unit (str): Output unit for energies in xyz comment field
            relative (bool): Give energies relative to lowest energy conformer
            separate (bool): Write conformers to separate xyz files.
        """
        # Retrieve symbols and energies
        symbols = convert_elements(self.elements, output="symbols")
        energies = self.get_energies(unit=unit, relative=relative)

        # Get xyz strings
        write_strings = []
        for conformer, energy in zip(self.conformers, energies):
            if energy is not None:
                energy = energy
            else:
                energy = 0.0
            xyz_string = get_xyz_string(symbols, conformer.coordinates,
                comment=f"{energy:.5f}")
            write_strings.append(xyz_string)
        
        # Write to file
        if not separate:
            with open(filename, "w") as file:
                write_string = "".join(write_strings)
                file.write(write_string)
        else:
            for i, write_string in enumerate(write_strings, start=1):
                conf_filename = filename.split(".")[0] + f"_{i}.xyz"
                with open(conf_filename, "w") as file:
                    file.write(write_string)                


def boltzmann_weights(energies, temperature=298.15):
    """Compute Boltzmann weights.
    
    Args:
        energies (list): Conformer energies (a.u.)
        temperature (float): Temperature (K)
    
    Returns:
        weights (ndarray): Conformer weights (normalized to unity)
    """
    energies = np.array(energies)
    energies -= energies.min()
    terms = np.exp(-energies / (K_B / HARTREE * temperature))
    weights = terms / np.sum(terms)
    
    return weights


def boltzmann_statistic(properties, energies, temperature=298.15,
    statistic="avg"):
    """Compute Boltzmann statistic.
    
    Args:
        energies (list): Conformer energies (a.u.)
        properties (list): Conformer properties
        statistic (str): Statistic to compute: 'avg', 'var' or 'std'
        temperature (float): Temperature (K)
    
    Returns:
        result (float): Boltzmann statistic
    """ 
    properties = np.array(properties)
    
    # Get conformer weights
    weights = boltzmann_weights(energies, temperature)

    # Compute Boltzmann weighted statistic
    if statistic == "avg":
        result = np.average(properties, weights=weights)
    elif statistic == "var":
        avg = np.average(properties, weights=weights)
        result = np.sum(weights * (properties - avg) ** 2)
    elif statistic == "std":
        avg = np.average(properties, weights=weights)
        var = np.sum(weights * (properties - avg) ** 2)
        result = np.sqrt(var)
    
    return result


def boltzmann_average_dT(properties, energies, temperature=298.15):
    """Return the derivative of the Boltzmann average.

    Args:
        energies (list): Conformer energies (a.u.)
        properties (list): Conformer properties
        temperature (float): Temperature (K)

    Returns:
        derivative (float): Derivative of Boltzmann average.
    """
    energies = np.array(energies)

    # Calculate Boltzmann averaged properties
    avg_prop_en = boltzmann_statistic(properties * energies * HARTREE,
        energies, temperature=temperature, statistic="avg")
    avg_en = boltzmann_statistic(energies * HARTREE, energies,
        temperature=temperature, statistic="avg")
    avg_prop = boltzmann_statistic(properties, energies,
        temperature=temperature, statistic="avg")

    # Calculate derivative
    derivative = (avg_prop_en - avg_en * avg_prop) / (K_B * temperature ** 2)

    return derivative


@conditional(_has_rdkit, _warning_rdkit)    
def generate_conformers(smiles, n_confs=None, version=2, optimize=None,
    small_rings=True, macrocycles=True, random_seed=None, prune_rmsd=0.35,
    n_threads=1):
    """Generates conformers for an RDKit mol object. Recipe based on 
    J. Chem. Inf. Modeling 2012, 52, 1146.

    Args:
        macrocycles (bool): Impose macrocycle torsion angle preferences.
        n_confs (int): Number of conformers to generate. If None, a 
            reasonable number will be set depending on the number of
            rotatable bonds.
        n_threads (int): Number of threads for force field optimization.
        optimize (str): Force field used for conformer optimization: 'MMFF',
            'UFF' or None
        prune_rmsd (float): Pruning RMSD threshold (Å).
        random_seed (int): Random seed for conformer generation.
        small_rings (bool): Impose small ring torsion angle preferences.
        smiles (str): SMILES string of molecule
        version (int): Version of the experimental torsion-angle preferences.

    Returns:
        conformer_coordinates (ndarray): Coordinates for all conformers (Å)
        connectivity_matrix (ndarray): Connectivity matrix
        elements (list): Atomic symbols
        energies (list): Conformer energies (a.u.) 
    """
    # Generate mol object
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # If n_confs is not set, set number of conformers based on number of
    # rotatable bonds 
    if not n_confs:
        n_rot_bonds = AllChem.CalcNumRotatableBonds(mol)

        if n_rot_bonds <= 7:
            n_confs = 50
        elif n_rot_bonds >= 8 and n_rot_bonds <= 12:
            n_confs = 200
        else:
            n_confs = 300
    
    # Generate conformers
    if prune_rmsd is not None:
        rdkit_prune_rmsd = prune_rmsd
    else:
        rdkit_prune_rmsd = -1
    if random_seed is not None:
        rdkit_random_seed = random_seed
    else:
        rdkit_random_seed = -1
    AllChem.EmbedMultipleConfs(
        mol, numConfs=n_confs, randomSeed=rdkit_random_seed,
        useSmallRingTorsions=small_rings, useMacrocycleTorsions=macrocycles,
        ETversion=version, pruneRmsThresh=rdkit_prune_rmsd
    )
    
    # Optimize with force fields
    results = None
    if optimize is not None:
        if optimize == "MMFF":
            results = AllChem.MMFFOptimizeMoleculeConfs(mol,
                numThreads=n_threads)
        if optimize == "UFF":
            results = AllChem.UFFOptimizeMoleculeConfs(mol,
                numThreads=n_threads)
    
    # Set energies
    if results is not None:
        energies = np.array([i[1] for i in results]) * KCAL_TO_HARTREE
    else:
        energies = None
    
    # Take out elements, coordinates and connectivity matrix
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]

    conformer_coordinates = []
    for conformer in mol.GetConformers():
        coordinates = conformer.GetPositions()
        conformer_coordinates.append(coordinates)
    conformer_coordinates = np.array(conformer_coordinates)

    connectivity_matrix = Chem.GetAdjacencyMatrix(mol)

    return elements, conformer_coordinates, energies, connectivity_matrix

