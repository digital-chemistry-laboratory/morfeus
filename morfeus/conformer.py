"""Conformer tools."""

import numbers
import subprocess
import tempfile
from collections import Counter
from copy import copy, deepcopy
from pathlib import Path

import numpy as np

from morfeus.data import (HARTREE, HARTREE_TO_KCAL, HARTREE_TO_KJ, K_B,
                          KCAL_TO_HARTREE, KJ_TO_HARTREE)
from morfeus.helpers import conditional, convert_elements
from morfeus.io import get_xyz_string
from morfeus.qc import (_generate_qcel_molecule, optimize_qc_engine,
                        sp_qc_engine)

# Optional dependencies
try:
    from openbabel import openbabel as ob
    _HAS_OPENBABEL = True
except ImportError:
    _HAS_OPENBABEL = False
_MSG_OPENBABEL = "Install OpenBabel to use this function."

try:
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False
_MSG_RDKIT = "Install RDKit to use this function."

try:
    import spyrmsd.rmsd
    _HAS_SPYRMSD = True
except ImportError:
    _HAS_SPYRMSD = False
_MSG_SPYRMSD = "Install spyrmsd to use this function."


def boltzmann_average_dT(properties, energies, temperature=298.15):
    """Return the derivative of the Boltzmann average.

    Args:
        properties (list): Conformer properties
        energies (list): Conformer energies (a.u.)
        temperature (float): Temperature (K)

    Returns:
        derivative (float): Derivative of Boltzmann average.
    """
    energies = np.array(energies)

    # Calculate Boltzmann averaged properties
    avg_prop_en = boltzmann_statistic(properties * energies * HARTREE,
                                      energies,
                                      temperature=temperature,
                                      statistic="avg")
    avg_en = boltzmann_statistic(energies * HARTREE,
                                 energies,
                                 temperature=temperature,
                                 statistic="avg")
    avg_prop = boltzmann_statistic(properties,
                                   energies,
                                   temperature=temperature,
                                   statistic="avg")

    # Calculate derivative
    derivative = (avg_prop_en - avg_en * avg_prop) / (K_B * temperature ** 2)

    return derivative


def boltzmann_statistic(properties, energies, temperature=298.15,
                        statistic="avg"):
    """Compute Boltzmann statistic.

    Args:
        properties (list): Conformer properties
        energies (list): Conformer energies (a.u.)
        temperature (float): Temperature (K)
        statistic (str): Statistic to compute: 'avg', 'var' or 'std'

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


class Conformer:
    """Conformer with coordinates, energy and properties.

    Args:
        elemennts (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
        energy (float): Energy (a.u.)
        degeneracy (int): Degeneracy
        properties (dict): Conformers properties.
        cip_label (tuple): Tuple of CIP labels for all atoms: 'R', 'S'
            or ''

    Attributes:
        coordinates (list): Coordinates (Å)
        degeneracy (int): Degeneracy
        elemennts (list): Elements as atomic symbols or numbers
        energy (float): Energy (a.u.)
        properties (dict): Conformers properties.
        cip_label (tuple): Tuple of CIP labels for all atoms: 'R', 'S'
            or ''
    """
    def __init__(self, elements, coordinates, energy=None, degeneracy=1,
                 properties=None, cip_label=None):
        if properties is None:
            properties = {}
        self.coordinates = np.array(coordinates)
        self.elements = np.array(elements)
        self.energy = energy
        self.properties = properties
        self.degeneracy = degeneracy
        self.cip_label = cip_label

    def __repr__(self):
        n_atoms = self.coordinates.shape[0]
        return f"{self.__class__.__name__}({n_atoms!r} atoms)"


class ConformerEnsemble:
    """Conformer ensemble object that supports sorting, pruning, optimization
    and single-point calculations.

    Args:
        elements (list): Elements as atomic symbols or numbers
        conformer_coordinates (list): Conformer coordinates (Å)
        energies (list): Energy (a.u.)
        connectivity_matrix (ndarray): Connectivity matrix
        degeneracies (list): Degeneracies
        properties (dict): Conformers properties.
        charge (int): Molecular charge
        multiplicity (int): Molecular multiplicity.
        ref_cip_label (tuple): Tuple of CIP labels for all atoms: 'R', 'S'
            or ''

    Attributes:
        elements (list): Elements as atomic symbols or numbers
        charge (int): Molecular charge
        conformers (list): Conformers
        connectivity_matrix (ndarray): Connectivity matrix
        elements (list): Elements as atomic symbols or numbers
        multiplicity (int): Molecular multiplicity.
        ref_cip_label (tuple): Tuple of CIP labels for all atoms: 'R', 'S'
            or ''
    """
    def __init__(
            self, elements, conformer_coordinates=None, energies=None,
            connectivity_matrix=None, degeneracies=None, properties=None,
            charge=None, multiplicity=None, ref_cip_label=None):
        if conformer_coordinates is None:
            conformer_coordinates = []

        # Store conformers
        self.elements = np.array(convert_elements(elements))
        self.conformers = []
        self._add_conformers(conformer_coordinates, energies, properties,
                             degeneracies)

        # Set up attributes
        self.charge = charge
        self.multiplicity = multiplicity
        self.connectivity_matrix = np.array(connectivity_matrix)
        self._mol = None
        self.ref_cip_label = ref_cip_label

    def add_conformers(self, coordinates, energies=None, degeneracies=None,
                       properties=None,):
        """Add conformer to ensemble.

        Args:
            coordinates (list): Conformer coordinates (Å)
            energies (list): Energies (a.u.)
            degeneracies (list): Degeneracies
            properties (dict): Conformer properties.
        """
        self._add_conformers(coordinates, energies, degeneracies, properties)

    def add_inverted(self):
        """Add inverted images of all conformers.

        Scrambles stereochemistry and leads to redundant conformers so use with
        care and prune based on RMSD as postprocessing.
        """
        conformers = self.conformers
        conformer_coordinates = []
        energies = []
        for conformer in conformers:
            coordinates = conformer.coordinates * np.array([-1, -1, -1])
            conformer_coordinates.append(coordinates)
            energies.append(conformer.energy)
        self.add_conformers(conformer_coordinates, energies)

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
        derivative = boltzmann_average_dT(properties,
                                          energies,
                                          temperature=temperature)

        return derivative

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
        statistic = boltzmann_statistic(
            properties, energies, temperature=temperature, statistic=statistic)

        return statistic

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

    def detect_enantiomers(self, thres=0.001):
        """Detect enantiomers in ensemble.

        Args:
            thres (float): RMSD threshold for detecting enantiomers in terms
                of coordinates.

        Returns:
            enantiomers (dict): Mapping of enantiomer with conformer id as keys
                and enantiomer ids as values.
        """
        # Add set of inverted conformers
        n_conformers = self.n_conformers
        self.add_inverted()

        # Map conformers to enantiomers
        enantiomers = {i: set() for i in range(n_conformers)}
        rmsds = self.get_rmsd()
        for i in range(n_conformers):
            # Do test that inverted conformer should have 0 RMSD to be
            # enantiomer
            test = np.where(rmsds[i, n_conformers:] < thres)[0]
            if len(test) > 0:
                if i not in test:
                    enantiomers[i].update(test)
                    for j in test:
                        enantiomers[j].add(i)
        enantiomers = {key: list(value) for key, value in enantiomers.items()}

        return enantiomers

    @classmethod
    def from_rdkit(cls, *args, **kwargs):
        """Generate conformer ensemble from RDKit.

        See the documentation for the function generate_conformers_rdkit for
        more information.

        Args:
            *args: Positional arguments
            **kwargs Keyword arguments
        """
        # Run RDKit conformer search and generate ensemble.
        elements, conformer_coordinates, energies, connectivity_matrix, mol = \
            generate_conformers_rdkit(*args, **kwargs)
        ce = cls(elements, conformer_coordinates, energies,
                 connectivity_matrix)
        ce._mol = mol

        # Set reference CIP label if enantiomerically pure.
        cip_labels = ce.get_cip_labels()
        if len(set(cip_labels)) == 1:
            ce.ref_cip_label = cip_labels[0]

        return ce

    @conditional(_HAS_RDKIT, _MSG_RDKIT)
    def get_cip_labels(self):
        """Generate tuples of CIP labels for conformer.

        Returns:
            cip_label (list): Tuples of CIP labels for each conformer.
        """
        # Update RDKit Mol object with current conformers.
        self._reset_mol()

        # Generate CIP labels with new RDKit
        cip_labels = []
        mol = self._mol
        for i in range(mol.GetNumConformers()):
            Chem.AssignStereochemistryFrom3D(mol, i)
            labels = []
            for atom in mol.GetAtoms():
                if atom.HasProp("_CIPCode"):
                    label = atom.GetProp("_CIPCode")
                else:
                    label = ""
                labels.append(label)
            cip_labels.append(tuple(labels))

        return cip_labels

    def get_coordinates(self):
        """Get conformers coordinates.

        Returns:
            conformer_coordinates (ndarray): Conformer coordinates (Å)
        """
        conformer_coordinates = np.array(
            [conformer.coordinates for conformer in self.conformers])

        return conformer_coordinates

    def get_degeneracies(self):
        """Get degeneracies.

        Returns:
            degeneriacies (ndarray): Degeneracies
        """
        degeneracies = np.array(
            [conformer.degeneracy for conformer in self.conformers])

        return degeneracies

    def get_energies(self):
        """Get energies.

        Returns:
            energies (ndarray): Energy (a.u.)
        """
        energies = np.array(
            [conformer.energy for conformer in self.conformers])
        return energies

    def get_properties(self):
        """Get conformer properties

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

    def get_relative_energies(self, unit="kcal/mol", relative=True):
        """Get conformer energies with choice of units and reference value.

        Args:
            unit (str): Unit of returned energies: 'hartree' (default),
                'kcal/mol' or 'kJ/mol'.
            relative (bool): Return energies relative to lowest-energy
                conformer.

        Returns:
            energies (ndarray): Conformer energies.
        """
        energies = self.get_energies()

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

    @property
    def n_conformers(self):
        "n_conformers (int): Number of conformers"
        n_conformers = len(self.conformers)

        return n_conformers

    def prune_enantiomers(self, keep="original", ref_label=None):
        """Prune out conformers so that only one enantiomer is present in the
        ensemble

        Args:
            keep (str): Which enantiomer to keep: 'original' (default), 'most
                common' or 'specified'. Choice of 'original' requires that the
                ref_cip_label attribute is set. Choice of 'specified' requires
                ref_label to be given.
            ref_label (tuple): Reference CIP labels for all atoms.
        """
        cip_labels = self.get_cip_labels()

        # Set up reference label
        if keep == "original":
            ref_label = self.ref_cip_label
            if ref_label is None:
                keep = "most common"
        elif keep == "specified":
            pass
        if keep == "most common":
            counter = Counter(cip_labels)
            ref_label = counter.most_common(n=1)[0][0]

        # Prune conformers
        to_keep = []
        for i, cip_label in enumerate(cip_labels):
            if cip_label == ref_label:
                to_keep.append(i)
        self.conformers = [self.conformers[i] for i in to_keep]

    def set_coordinates(self, conformer_coordinates):
        """Set conformer coordinates.

        Args:
            conformer_coordinates (ndarray): Conformer coordinates (Å)
        """
        if len(conformer_coordinates) != self.n_conformers:
            msg = f"Number of coordinates ({len(conformer_coordinates)}) " \
                f"!= number of conformers ({self.n_conformers})."
            raise ValueError(msg)
        for conformer, coordinates in zip(self.conformers,
                                          conformer_coordinates):
            conformer.coordinates = coordinates

    def set_degeneracies(self, degeneracies):
        """Set degeneracies.

        Args:
            degeneriacies (ndarray): Degeneracies
        """
        if len(degeneracies) != self.n_conformers:
            msg = f"Number of degeneracies ({len(degeneracies)}) " \
                f"!= number of conformers ({self.n_conformers})."
            raise ValueError(msg)
        for conformer, degeneracy in zip(self.conformers, degeneracies):
            conformer.degeneracy = degeneracy

    def set_energies(self, energies):
        """Set energies.

        Args:
            energies (ndarray): Energy (a.u.)
        """
        if len(energies) != self.n_conformers:
            msg = f"Number of energies ({len(energies)}) != number of " \
                f"conformers ({self.n_conformers})."
            raise ValueError(msg)
        for conformer, energy in zip(self.conformers, energies):
            conformer.energy = energy

    def set_properties(self, key, values):
        """Set conformer properties.

        Args:
            key (str): Name of property
            values (list): Property values
        """
        for conformer, value in zip(self.conformers, values):
            conformer.properties[key] = value

    def get_rmsd(self, i_s=None, j_s=None, include_hs=False, symmetry=True,
                 method="openbabel"):
        """Get RSMD between two conformers.

        For very small systems 'openbabel' or 'spyrmsd' work well. For larger
        systems a significant speed-up is attained with 'obrms-batch' or
        'obrms-iter'.

        Args:
            i_s (list): Indices of conformer 1
            j_s (list): Indices of conformer 2
            include_hs (bool): Include H atoms in RMSD calculation. Ignored for
                'obrms-iter' and 'obrms-batch' which only use heavy atoms.
            symmetry (bool): Consider symmetry (requires connectivity matrix).
                Ignored for 'obrms-iter' and 'obrms-batch' which always use
                symmetry.
            method (str): RMSD calculation method: 'obrms-batch', 'obrms-iter',
                openbabel' (default) or 'spyrmsd'.

        Returns:
            rmsds (ndarray): RSMDs (Å)
        """
        if i_s is None:
            i_s = np.arange(1, len(self.conformers) + 1)
        else:
            i_s = np.array(i_s)
        if j_s is None:
            j_s = np.arange(1, len(self.conformers) + 1)
        else:
            j_s = np.array(j_s)

        if method == "obrms-batch":
            rmsds = self._get_rmsd_obrms_batch(i_s, j_s)
        if method == "obrms-iter":
            rmsds = self._get_rmsd_obrms_iter(i_s, j_s)
        if method == "openbabel":
            rmsds = self._get_rmsd_openbabel(i_s, j_s, include_hs, symmetry)
        elif method == "spyrmsd":
            rmsds = self._get_rmsd_spyrmsd(i_s, j_s, include_hs, symmetry)
        return rmsds

    def optimize_qc_engine(
            self, ids=None, program=None, model=None, keywords=None,
            local_options=None, procedure="berny"):
        """Optimize conformers with QCEngine interface.

        Args:
            ids (list): Conformer indices to optimize. If None, all are
                optimized. 1-index.
            program (str): QCEngine program
            model (dict): QCEngine model
            keywords (dict): QCEngine keywords
            local_options (dict): QCEngine local options
            procedure (str): QCEngine procedure
        """
        # Set defaults
        if model is None:
            model = {}
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

    def prune_rmsd(self, thres=0.35, include_hs=False, symmetry=True,
                   method="openbabel"):
        """Prune conformers based on RMSD.

        Args:
            thres (float): Threshold for RSMD pruning (Å)
            include_hs (bool): Include H atoms in RMSD calculation. Ignored for
                'obrms-iter' and 'obrms-batch' which only use heavy atoms.
            symmetry (bool): Consider symmetry (requires connectivity matrix).
                Ignored for 'obrms-iter' and 'obrms-batch' which always use
                symmetry.
            method (str): RMSD calculation method: 'obrms-batch', 'obrms-iter',
                openbabel' (default) or 'spyrmsd'
        """
        # Select conformers to keep
        candidates = np.arange(len(self.conformers))
        keep_list = []

        # Prune conformers iteratively if possible to reduce number of RMSD
        # evaluations. For large systems obrms might be faster in batch.
        if method in ("openbabel", "spyrmsd"):
            while len(candidates) > 0:
                keeper = candidates[0]
                keep_list.append(keeper)
                rmsd = self.get_rmsd(
                    [keeper + 1],
                    candidates + 1, include_hs=include_hs, symmetry=symmetry,
                    method=method)
                candidates = candidates[rmsd[0] > thres]
        elif method == "obrms-batch":
            rmsds = self.get_rmsd(include_hs=include_hs,
                                  symmetry=symmetry, method=method)
            working_array = rmsds
            while len(working_array) > 0:
                keeper = candidates[0]
                keep_list.append(keeper)
                rmsd = working_array[0]
                mask = rmsd > thres
                candidates = candidates[mask]
                working_array = working_array[mask, :][:, mask]
        elif method == "obrms-iter":
            while len(candidates) > 0:
                keeper = candidates[0]
                keep_list.append(keeper)
                rmsd = self.get_rmsd(
                    [keeper + 1],
                    candidates + 1, include_hs=include_hs, symmetry=symmetry,
                    method=method)
                candidates = candidates[rmsd[0] > thres]

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
        elif unit.lower() == "kJ/mol":
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
            program (str): QCEngine program
            model (dict): QCEngine model
            keywords (dict): QCEngine keywords
            local_options (dict): QCEngine local options
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

    def write_conformers(self, file, ids=None, unit="kcal/mol", relative=True,
                         separate=False):
        """Write conformers to xyz file.

        Args:
            file (str): Filename or path object. Needs filename if
                `separate=True`
            ids (list): Conformer indices (1-indexed)
            unit (str): Output unit for energies in xyz comment field
            relative (bool): Give energies relative to lowest energy conformer
            separate (bool): Write conformers to separate xyz files.
        """
        if ids is None:
            ids = np.arange(len(self.conformers))
        else:
            ids = np.array(ids) - 1

        # Retrieve symbols and energies
        symbols = convert_elements(self.elements, output="symbols")
        conformers = [conformer for i, conformer in
                      enumerate(self.conformers) if i in ids]
        energies = self.get_relative_energies(
            unit=unit, relative=relative)[ids]

        # Get xyz strings
        write_strings = []
        for conformer, energy in zip(conformers, energies):
            if energy is None:
                energy = 0.0
            xyz_string = get_xyz_string(symbols, conformer.coordinates,
                                        comment=f"{energy:.5f}")
            write_strings.append(xyz_string)

        # Write to file
        if not separate:
            with open(file, "w") as f:
                write_string = "".join(write_strings) + "\n"
                f.write(write_string)
        else:
            for i, write_string in zip(ids, write_strings):
                conf_filename = file.split(".")[0] + f"_{i + 1}.xyz"
                with open(conf_filename, "w") as f:
                    f.write(write_string)

    def _add_conformers(self, conformer_coordinates, energies=None,
                        properties=None, degeneracies=None):
        n_conformers = len(conformer_coordinates)
        if properties is None:
            properties = {}
        if energies is None:
            energies = [np.nan] * n_conformers
        if degeneracies is None:
            degeneracies = np.ones(n_conformers)
        for coordinates, energy, degeneracy in zip(conformer_coordinates,
                                                   energies, degeneracies):
            conformer = Conformer(
                self.elements, coordinates, energy, degeneracy)
            self.conformers.append(conformer)
        for key, value in properties.items():
            self.set_properties(key, value)

    @conditional(_HAS_OPENBABEL, _MSG_OPENBABEL)
    def _get_rmsd_obrms_batch(self, i_s, j_s):
        """Calculate RMSD with obrms in batch mode, first calculating the
        matrix of all pairwise RMSDs and then taking out those of interest."""
        with tempfile.NamedTemporaryFile(suffix=".xyz") as ref_file:
            p_ref = Path(ref_file.name)
            self.write_conformers(p_ref, unit="hartree", relative=False)
            process = subprocess.run(
                f"obrms {p_ref.as_posix()} --cross "
                "--minimize".split(" "),
                capture_output=True)
        rmsds = np.genfromtxt(process.stdout.splitlines(),
                              delimiter=",")[:, 1:]
        rmsds = rmsds[i_s - 1, :][:, j_s - 1]
        return rmsds

    @conditional(_HAS_OPENBABEL, _MSG_OPENBABEL)
    def _get_rmsd_obrms_iter(self, i_s, j_s):
        """Calculate RMSD with obrms in iterative row-wise mode for heavy atoms
        and witout symmetry."""
        rmsds = []
        for i in i_s:
            with tempfile.NamedTemporaryFile(suffix=".xyz") as ref_file, \
                    tempfile.NamedTemporaryFile(suffix=".xyz") as \
                    test_file:
                p_ref = Path(ref_file.name)
                p_test = Path(test_file.name)
                self.write_conformers(p_ref, ids=j_s)
                self.write_conformers(p_test, ids=[i])
                process = subprocess.run(
                    f"obrms {p_ref.as_posix()} "
                    f"{p_test.as_posix()} --minimize".split(" "),
                    capture_output=True)
            row_rmsds = np.genfromtxt(process.stdout.splitlines(),
                                      usecols=(-1))
            rmsds.append(row_rmsds)
        rmsds = np.vstack(rmsds)
        return rmsds

    @conditional(_HAS_OPENBABEL, _MSG_OPENBABEL)
    def _get_rmsd_openbabel(self, i_s, j_s, include_hs, symmetry):
        """Calculate RMSD row-wise with openbabel python interface."""
        rmsds = []
        for i in i_s:
            conformer_1 = self.conformers[i - 1]
            ob_mol_1 = _get_ob_mol(self.elements, conformer_1.coordinates,
                                   self.connectivity_matrix)
            align = ob.OBAlign(include_hs, symmetry)
            align.SetRefMol(ob_mol_1)
            rmsds_row = []
            for j in j_s:
                conformer_2 = self.conformers[j - 1]
                ob_mol_2 = _get_ob_mol(self.elements, conformer_2.coordinates,
                                       self.connectivity_matrix)
                align.SetTargetMol(ob_mol_2)
                align.Align()
                rmsd = align.GetRMSD()
                rmsds_row.append(rmsd)
            rmsds.append(rmsds_row)
        rmsds = np.array(rmsds)
        return rmsds

    @conditional(_HAS_SPYRMSD, _MSG_SPYRMSD)
    def _get_rmsd_spyrmsd(self, i_s, j_s, include_hs, symmetry):
        """Calculate RMSD row-wise with spyrmsd"""
        # Construct mask for H atoms
        if not include_hs:
            mask = self.elements != 1
        else:
            mask = np.ones(len(self.elements), dtype=np.bool)

        # Calculate RMSD. Use qcel with no. heavy atoms < 3 due to unstable
        # performance of spyrmsd.
        if np.sum(self.elements != 1) < 3:
            rmsds = []
            for i in i_s:
                ref_coordinates = self.conformers[i - 1].coordinates
                ref_molecule = _generate_qcel_molecule(self.elements[mask],
                                                       ref_coordinates[mask])
                row_rmsds = []

                for j in j_s:
                    test_coordinates = self.conformers[j - 1].coordinates
                    test_molecule = _generate_qcel_molecule(
                        self.elements[mask],
                        test_coordinates[mask])
                    rmsd = ref_molecule.align(test_molecule)[1]["rmsd"]
                    row_rmsds.append(rmsd)
            rmsds.append(np.array(row_rmsds))
            rmsds = np.vstack(rmsds)
        else:
            rmsds = []
            for i in i_s:
                ref_coordinates = self.conformers[i - 1].coordinates[mask]
                if symmetry:
                    test_coordinates = [
                        self.conformers[j - 1].coordinates[mask] for j in j_s]
                    row_rmsds = spyrmsd.rmsd.symmrmsd(
                        ref_coordinates,
                        test_coordinates,
                        self.elements[mask],
                        self.elements[mask],
                        self.connectivity_matrix[mask, :][:, mask],
                        self.connectivity_matrix[mask, :][:, mask],
                        center=True,
                        minimize=True)
                else:
                    row_rmsds = []
                    for j in j_s:
                        test_coordinates = self.conformers[j -
                                                           1].coordinates[mask]
                        rmsd = spyrmsd.rmsd.rmsd(ref_coordinates,
                                                 test_coordinates,
                                                 self.elements[mask],
                                                 self.elements[mask],
                                                 center=True,
                                                 minimize=True)
                        row_rmsds.append(rmsd)
                rmsds.append(np.array(row_rmsds))
            rmsds = np.vstack(rmsds)
        return rmsds

    @conditional(_HAS_RDKIT, _MSG_RDKIT)
    def _reset_mol(self):
        self._mol.RemoveAllConformers()
        conformer_coordinates = self.get_coordinates()
        _add_conformers_to_mol(self._mol, conformer_coordinates)

    def __copy__(self):
        # Generate copy where conformers and mol object are shared with old
        # ensemble.
        cls = type(self)
        ce = cls(
            self.elements,
            connectivity_matrix=self.connectivity_matrix,
            charge=self.charge,
            multiplicity=self.multiplicity,
            ref_cip_label=self.ref_cip_label,
        )
        ce._mol = self._mol
        ce.conformers = self.conformers
        return ce

    def __deepcopy__(self, memo):
        # Generate copy where conformers and mol object are new
        cls = type(self)
        ce = cls(
            self.elements,
            conformer_coordinates=self.get_coordinates(),
            energies=self.get_energies(),
            properties=self.get_properties(),
            degeneracies=self.get_degeneracies(),
            connectivity_matrix=self.connectivity_matrix,
            charge=self.charge,
            multiplicity=self.multiplicity,
            ref_cip_label=self.ref_cip_label,
        )
        ce._mol = deepcopy(ce._mol, memo)
        return ce

    def __delitem__(self, index):
        del self.conformers[index]

    def __getitem__(self, index):
        cls = type(self)
        if isinstance(index, slice):
            # Generate copy of ensemble with selected conformers.
            ce = copy(self)
            ce.conformers = ce.conformers[index]
            return ce
        elif isinstance(index, numbers.Integral):
            return self.conformers[index]
        else:
            msg = f'{cls.__name__} indices must be integers'
            raise TypeError(msg)

    def __len__(self):
        return len(self.conformers)

    def __repr__(self):
        n_conformers = len(self.conformers)
        return f"{self.__class__.__name__}({n_conformers!r} conformers)"


@conditional(_HAS_RDKIT, _MSG_RDKIT)
def generate_conformers_rdkit(
        smiles, n_confs=None, optimize=None, version=2, small_rings=True,
        macrocycles=True, random_seed=None, rmsd_thres=0.35, n_threads=1,
        add_rotamers=False):
    """Generates conformers for an RDKit mol object. Recipe based on
    J. Chem. Inf. Modeling 2012, 52, 1146.

    Args:
        smiles (str): SMILES string of molecule
        n_confs (int): Number of conformers to generate. If None, a
            reasonable number will be set depending on the number of
            rotatable bonds.
        optimize (str): Force field used for conformer optimization: 'MMFF',
            'UFF' or None
        version (int): Version of the experimental torsion-angle preferences.
        small_rings (bool): Impose small ring torsion angle preferences.
        macrocycles (bool): Impose macrocycle torsion angle preferences.
        random_seed (int): Random seed for conformer generation.
        rmsd_thres (float): Pruning RMSD threshold (Å).
        n_threads (int): Number of threads.
        add_rotamers (bool): Whether to add H rotamers after RDKit step.

    Returns:
        elements (list): Atomic symbols
        conformer_coordinates (ndarray): Coordinates for all conformers (Å)
        energies (list): Conformer energies (a.u.)
        connectivity_matrix (ndarray): Connectivity matrix with bond orders
        mol (obj): RDKit Mol object. Only returned when `return_mol=True`.
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
    if rmsd_thres is not None:
        rdkit_prune_rmsd = rmsd_thres
    else:
        rdkit_prune_rmsd = -1
    if random_seed is not None:
        rdkit_random_seed = random_seed
    else:
        rdkit_random_seed = -1
    AllChem.EmbedMultipleConfs(
        mol, numConfs=n_confs, randomSeed=rdkit_random_seed,
        useSmallRingTorsions=small_rings, useMacrocycleTorsions=macrocycles,
        ETversion=version, pruneRmsThresh=rdkit_prune_rmsd,
        numThreads=n_threads,
    )

    if add_rotamers:
        _add_rotamers_to_mol(mol)

    # Optimize with force fields
    results = None
    if optimize is not None:
        if optimize == "MMFF94":
            results = AllChem.MMFFOptimizeMoleculeConfs(
                mol, mmffVariant="MMFF94", numThreads=n_threads)
        elif optimize == "MMFF94s":
            results = AllChem.MMFFOptimizeMoleculeConfs(
                mol, mmffVariant="MMFF94s", numThreads=n_threads)
        if optimize == "UFF":
            results = AllChem.UFFOptimizeMoleculeConfs(
                mol, numThreads=n_threads)

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

    connectivity_matrix = Chem.GetAdjacencyMatrix(mol, useBO=True)

    return elements, conformer_coordinates, energies, connectivity_matrix, mol


@conditional(_HAS_RDKIT, _MSG_RDKIT)
def _add_conformers_to_mol(mol, conformer_coordinates):
    """Add conformers to RDKit Mol object.

    Args:
        conformer_coordinates (list): Conformer coordinates (Å)
    """
    conformer_coordinates = np.array(conformer_coordinates)
    if len(conformer_coordinates.shape) == 2:
        conformer_coordinates.reshape(-1, conformer_coordinates.shape[0], 3)

    for coordinates in conformer_coordinates:
        conformer = Chem.Conformer()
        for i, coord in enumerate(coordinates):
            point = rdkit.Geometry.Point3D(*coord)
            conformer.SetAtomPosition(i, point)
        mol.AddConformer(conformer, assignId=True)


@conditional(_HAS_RDKIT, _MSG_RDKIT)
def _add_rotamers_to_mol(mol):
    """Add simple CH₃, NH₂ etc. rotamers to conformers of RDKit Mol object

    Args:
        mol (obj): RDKit Mol object
    """
    dihedrals = _get_dihedrals(mol)
    for conformer in mol.GetConformers():
        ref_id = conformer.GetId()
        for dihedral in dihedrals:
            conformer = mol.GetConformer(ref_id)
            for increment in [-120, 120]:
                new_id = mol.AddConformer(conformer, assignId=True)
                conformer = mol.GetConformer(new_id)
                angle = AllChem.GetDihedralDeg(conformer, *dihedral)
                AllChem.SetDihedralDeg(conformer, *dihedral, angle + increment)


@conditional(_HAS_RDKIT, _MSG_RDKIT)
def _get_dihedrals(mol):
    """Return dihedral which ends in saturated atom of type CH₃, NH₂ etc.

    Args:
        mol (obj): Molecule as RDKit Mol object.

    Returns:
        dihedral_indices:
            (0-indexed)
    """
    # Set up SMARTS query and match bonds to satured XH3, XH2, or XH1.
    smarts = "[Av4H3,Av3H2,Av2H1]-[*!D1]"
    query = Chem.MolFromSmarts(smarts)
    matches = mol.GetSubstructMatches(query)

    # Take out dihedral indices for match.
    dihedral_indices = []
    for match in matches:
        begin_atom = mol.GetAtomWithIdx(match[0])
        end_atom = mol.GetAtomWithIdx(match[1])
        j = begin_atom.GetIdx()
        k = end_atom.GetIdx()
        i = [atom for atom in begin_atom.GetNeighbors() if atom.GetIdx()
             is not end_atom.GetIdx()][0].GetIdx()
        l = [atom for atom in end_atom.GetNeighbors() if atom.GetIdx()
             is not begin_atom.GetIdx()][0].GetIdx()

        # Order indices so that the end atom is the saturated AHX group.
        dihedral_indices.append((l, k, j, i))

    return dihedral_indices


@conditional(_HAS_OPENBABEL, _MSG_OPENBABEL)
def _get_ob_mol(elements, coordinates, connectivity_matrix):
    """Generate OpenBabel OBMol object.

    Args:
        elements (list): Elements as atomic symbols or numbers.
        coordinates (list): Coordinates (Å)
        connectivity_matrix (ndarray): Connectivity matrix with bond orders

    Returns:
        mol (obj): OpenBabel OBMol object
    """
    elements = convert_elements(elements)

    # Construct OBMol object based on atoms and bonds
    mol = ob.OBMol()
    for element, coordinate in zip(elements, coordinates):
        a = mol.NewAtom()
        a.SetAtomicNum(int(element))
        a.SetVector(*coordinate)
    for i, j in zip(*np.tril_indices_from(connectivity_matrix)):
        if i != j:
            bo = connectivity_matrix[i, j]
            if bo != 0:
                mol.AddBond(int(i + 1), int(j + 1), int(bo))
    return mol
