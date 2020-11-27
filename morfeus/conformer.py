"""Conformer tools."""

import numbers
import subprocess
import tempfile
import warnings
from collections import Counter
from copy import copy, deepcopy
from pathlib import Path

import numpy as np
from pkg_resources import parse_version

from morfeus.data import (HARTREE, HARTREE_TO_KCAL, HARTREE_TO_KJ, K_B,
                          KCAL_TO_HARTREE, KJ_TO_HARTREE)
from morfeus.helpers import (Import, convert_elements, requires_dependency,
                             requires_executable)
from morfeus.io import write_xyz
from morfeus.qc import optimize_qc_engine, sp_qc_engine


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
        formal_charges (list): Atomic formal charges
        ref_cip_label (tuple): Tuple of CIP labels for all atoms: 'R', 'S'
            or ''

    Attributes:
        elements (ndarray): Elements as atomic symbols or numbers
        charge (int): Molecular charge
        conformers (list): Conformers
        connectivity_matrix (ndarray): Connectivity matrix
        formal_charges (ndarray): Atomic formal charges.
        mol (obj): RDKit mol object.
        multiplicity (int): Molecular multiplicity.
        ref_cip_label (tuple): Tuple of CIP labels for all atoms: 'R', 'S'
            or ''
    """

    def __init__(
            self, elements, conformer_coordinates=None, energies=None,
            connectivity_matrix=None, degeneracies=None, properties=None,
            charge=None, multiplicity=None, formal_charges=None,
            ref_cip_label=None):
        if conformer_coordinates is None:
            conformer_coordinates = []
        if formal_charges is None:
            formal_charges = np.zeros(len(elements))
        if charge is None:
            charge = np.sum(formal_charges)
        else:
            if charge != np.sum(formal_charges):
                msg = f"Charge ({charge}) is different from sum of partial " \
                    f"charges ({np.sum(formal_charges)})"
                raise Exception(msg)
        if multiplicity is None:
            multiplicity = 1

        # Store conformers
        self.elements = np.array(convert_elements(elements))
        self.conformers = []
        self._add_conformers(conformer_coordinates, energies, properties,
                             degeneracies)

        # Set up attributes
        self.charge = charge
        self.formal_charges = np.array(formal_charges)
        self.multiplicity = multiplicity
        self.connectivity_matrix = np.array(connectivity_matrix)
        self.mol = None
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

    def align_conformers(self):
        """Align conformers with RDKit."""
        self.update_mol()
        AllChem.AlignMolConformers(self.mol)

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

    def detect_enantiomers(self, thres=0.01, method="rdkit", include_hs=False):
        """Detect enantiomers in ensemble.

        Args:
            thres (float): RMSD threshold for detecting enantiomers in terms
                of coordinates.
            method (str): RMSD calculation method: 'obrms-batch', 'obrms-iter',
                openbabel' (default) or 'spyrmsd'.

        Returns:
            enantiomers (dict): Mapping of enantiomer with conformer id as keys
                and enantiomer ids as values.
        """
        # Add set of inverted conformers
        n_conformers = self.n_conformers
        self.add_inverted()

        # Map conformers to enantiomers
        enantiomers = {i: set() for i in range(n_conformers)}
        rmsds = self.get_rmsd(method=method,
                              include_hs=include_hs,
                              symmetry=False)
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

        # Reset number of conformers
        self.conformers = self.conformers[:n_conformers]

        return enantiomers

    @classmethod
    def from_rdkit(cls, *args, **kwargs):
        """Generate conformer ensemble from RDKit.

        See the documentation for the function generate_conformers_rdkit for
        more information.

        Args:
            *args: Positional arguments
            **kwargs Keyword arguments

        Returns:
            ce (obj): Conformer ensemble object.
        """
        # Run RDKit conformer search and generate ensemble.
        elements, conformer_coordinates, energies, connectivity_matrix, \
            charges, mol = generate_conformers_rdkit(*args, **kwargs)
        ce = cls(elements, conformer_coordinates, energies=energies,
                 connectivity_matrix=connectivity_matrix,
                 formal_charges=charges)
        ce.mol = mol
        ce.set_multiplicity_from_mol()

        # Set reference CIP label if enantiomerically pure.
        cip_labels = ce.get_cip_labels()
        if len(set(cip_labels)) == 1:
            ce.ref_cip_label = cip_labels[0]

        return ce

    @classmethod
    def from_openbabel_ga(
            cls, *args, generate_rdkit_mol=False, update_charges=True,
            update_connectivity=True, update_multiplicity=True, **kwargs):
        """Generate conformer ensemble from OpenBabel with GA method.

        See the documentation for the function generate_conformers_openbabel_ga
        for more information.

        Args:
            *args: Positional arguments for generate_conformers_openbabel_ga
            generate_rdkit_mol (bool): Generate RDKit mol object for ensemble
            update_charges (bool): Update formal charges from generated RDKit
                Mol object. Only used if generate_rdkit_mol is True.
            update_connectivity (bool): Update connectivity from generated
                RDKit Mol object. Only used if generate_rdkit_mol is True.
            update_multiplicity (bool): Update multiplicity from generated
                RDKit Mol object. Only used if generate_rdkit_mol is True.
            **kwargs: Keyword arguments for generate_conformers_openbabel_ga

        Returns:
            ce (obj): Conformer ensemble.
        """
        # Run Openbabel conformer search and generate ensemble.
        (elements, conformer_coordinates, connectivity_matrix, charges,
         ob_mol) = generate_conformers_openbabel_ga(*args, **kwargs)
        ce = cls(elements, conformer_coordinates,
                 connectivity_matrix=connectivity_matrix,
                 formal_charges=charges)

        # Generate RDKit mol object and CIP labels
        if generate_rdkit_mol:
            ce.generate_mol(update_charges=update_charges,
                            update_connectivity=update_connectivity,
                            update_multiplicity=update_multiplicity)

            cip_labels = ce.get_cip_labels()
            if len(set(cip_labels)) == 1:
                ce.ref_cip_label = cip_labels[0]

        return ce

    @classmethod
    def from_openbabel_ff(
            cls, *args, generate_rdkit_mol=False, update_charges=True,
            update_connectivity=True, update_multiplicity=True, **kwargs):
        """Generate conformer ensemble from OpenBabel with FF method.

        See the documentation for the function generate_conformers_openbabel_ff
        for more information.

        Args:
            *args: Positional arguments for generate_conformers_openbabel_ga
            generate_rdkit_mol (bool): Generate RDKit mol object for ensemble
            update_charges (bool): Update formal charges from generated RDKit
                Mol object. Only used if generate_rdkit_mol is True.
            update_connectivity (bool): Update connectivity from generated
                RDKit Mol object. Only used if generate_rdkit_mol is True.
            update_multiplicity (bool): Update multiplicity from generated
                RDKit Mol object. Only used if generate_rdkit_mol is True.
            **kwargs: Keyword arguments for generate_conformers_openbabel_ga

        Returns:
            ce (obj): Conformer ensemble.
        """
        # Run Openbabel conformer search and generate ensemble.
        (elements, conformer_coordinates, connectivity_matrix, charges,
         ob_mol) = generate_conformers_openbabel_ff(*args, **kwargs)
        ce = cls(elements, conformer_coordinates,
                 connectivity_matrix=connectivity_matrix,
                 formal_charges=charges)

        # Generate RDKit mol object and CIP labels
        if generate_rdkit_mol:
            ce.generate_mol(update_charges=update_charges,
                            update_connectivity=update_connectivity,
                            update_multiplicity=update_multiplicity)

            cip_labels = ce.get_cip_labels()
            if len(set(cip_labels)) == 1:
                ce.ref_cip_label = cip_labels[0]

        return ce

    @requires_dependency([Import(module="rdkit", item="Chem")], globals())
    def generate_mol(self,
                     update_charges=True,
                     update_connectivity=True,
                     update_multiplicity=True):
        """Generate RDKit Mol object"""
        mol = _get_rdkit_mol(self.elements, self.get_coordinates(),
                             self.connectivity_matrix, self.formal_charges)
        if update_charges:
            self.formal_charges = np.array([
                atom.GetFormalCharge() for atom in mol.GetAtoms()
            ])
        if update_connectivity:
            self.connectivity_matrix = np.array(
                Chem.GetAdjacencyMatrix(mol, useBO=True))
        if update_multiplicity:
            self.set_multiplicity_from_mol()

        self.mol = mol

    @requires_dependency([Import(module="rdkit", item="Chem")], globals())
    def get_cip_labels(self):
        """Generate tuples of CIP labels for conformer.

        Returns:
            cip_label (list): Tuples of CIP labels for each conformer.
        """
        # Update RDKit Mol object with current conformers.
        self.update_mol()

        # Generate CIP labels with new RDKit
        cip_labels = []
        mol = self.mol
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

    def get_rmsd(self, i_s=None, j_s=None, include_hs=False, symmetry=False,
                 method="rdkit"):
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
        elif method == "obrms-iter":
            rmsds = self._get_rmsd_obrms_iter(i_s, j_s)
        elif method == "openbabel":
            rmsds = self._get_rmsd_openbabel(i_s, j_s, include_hs, symmetry)
        elif method == "rdkit":
            rmsds = self._get_rmsd_rdkit(i_s, j_s, include_hs)
        elif method == "spyrmsd":
            rmsds = self._get_rmsd_spyrmsd(i_s, j_s, include_hs, symmetry)
        return rmsds

    @property
    def n_conformers(self):
        "n_conformers (int): Number of conformers"
        n_conformers = len(self.conformers)

        return n_conformers

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
        if np.any(self.formal_charges != 0) and (program.lower() == "rdkit"):
            raise Exception("QCEngine using RDKit does not work with charges.")

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

    def condense_enantiomeric(self, thres=None):
        cip_labels = self.get_cip_labels()
        if "".join(["".join(labels) for labels in cip_labels]) != "":
            raise Exception("Chiral molecule. Not safe to condense.")

        enantiomers = self.detect_enantiomers(thres=0.01)
        to_keep = []
        to_remove = []
        for key, value in enantiomers.items():
            if key not in to_remove:
                conformer = self.conformers[key]
                conformer.degeneracy = len(value) + 1
                to_keep.append(conformer)
                to_remove.extend(value)

        self.conformers = to_keep

    def prune_enantiomers(self, keep="original", ref_label=None):
        """Prune out conformers so that only one enantiomer is present in the
        ensemble.

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

    def prune_rmsd(self, thres=0.35, include_hs=False, symmetry=False,
                   method="rdkit"):
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
        if method in ("obrms-iter", "openbabel", "spyrmsd", "rdkit"):
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

        # Update conformer list
        self.conformers = [conformer for i, conformer in
                           enumerate(self.conformers) if i in keep_list]

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

    @requires_dependency([Import(module="rdkit.Chem", item="Descriptors")],
                         globals())
    def set_multiplicity_from_mol(self):
        """Sets multiplicity based on unpaired electrons in Mol object."""
        num_radical = Descriptors.NumRadicalElectrons(self.mol)

        # Assume maximum spin pairing
        if (num_radical % 2) == 0:
            multiplicity = 1
        else:
            multiplicity = 2

        self.multiplicity = multiplicity

    def set_properties(self, key, values):
        """Set conformer properties.

        Args:
            key (str): Name of property
            values (list): Property values
        """
        for conformer, value in zip(self.conformers, values):
            conformer.properties[key] = value

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
        if np.any(self.formal_charges != 0) and (program.lower() == "rdkit"):
            raise Exception("QCEngine using RDKit does not work with charges.")

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

    def write_xyz(self, file, ids=None, unit="kcal/mol", relative=True,
                  separate=False, n_decimals=3):
        """Write conformers to xyz file.

        Args:
            file (str): Filename or path object. Needs filename if
                `separate=True`
            ids (list): Conformer indices (1-indexed)
            unit (str): Output unit for energies in xyz comment field
            relative (bool): Give energies relative to lowest energy conformer
            separate (bool): Write conformers to separate xyz files.
            n_decimals (int): Number of decimals for energies.
        """
        if ids is None:
            ids = np.arange(len(self.conformers))
        else:
            ids = np.array(ids) - 1

        # Retrieve symbols, coordinates and energies
        symbols = convert_elements(self.elements, output="symbols")
        conformer_coordinates = self.get_coordinates()[ids]
        energies = self.get_relative_energies(
            unit=unit, relative=relative)[ids].round(n_decimals)

        # Write conformers
        if separate:
            for i, coordinates, energy in zip(ids, conformer_coordinates,
                                              energies):
                conf_filename = file.split(".")[0] + f"_{i + 1}.xyz"
                write_xyz(conf_filename,
                          symbols,
                          conformer_coordinates,
                          comments=energies)
        else:
            write_xyz(file, symbols, conformer_coordinates, comments=energies)

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

    @requires_executable(["obrms"])
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
        result = np.genfromtxt(process.stdout.splitlines(),
                               delimiter=",")
        # Reshape array to 2D if 1D
        if len(result.shape) == 1:
            result = result.reshape(1, -1)
        rmsds = result[:, 1:]

        rmsds = rmsds[i_s - 1, :][:, j_s - 1]

        return rmsds

    @requires_executable(["obrms"])
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

    @requires_dependency([Import(module="openbabel.openbabel", alias="ob")],
                         globals())
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

    @requires_dependency([
        Import(module="rdkit", item="Chem"),
        Import(module="rdkit.Chem", item="AllChem")
    ], globals())
    def _get_rmsd_rdkit(self, i_s, j_s, include_hs):
        """Calculate RMSD row-wise with RDKit."""
        # Update mol object from conformers and make copy.
        self.update_mol()
        mol = Chem.Mol(self.mol)

        # Construct atom list for rmsd: heavy atoms or all atoms
        if include_hs:
            atom_ids = [atom.GetIdx() for atom in mol.GetAtoms()]
        else:
            atom_ids = [atom.GetIdx() for atom in mol.GetAtoms()
                        if atom.GetAtomicNum() != 1]

        # Calculated RMSD row-wise with RDKit
        rmsds = []
        conformers = mol.GetConformers()
        for i in i_s:
            ref_conformer = conformers[i - 1]
            mol.RemoveAllConformers()
            mol.AddConformer(ref_conformer)
            for j in j_s:
                conformer = conformers[j - 1]
                mol.AddConformer(conformer)
            rmsds_row = []
            AllChem.AlignMolConformers(
                mol, atomIds=atom_ids, RMSlist=rmsds_row)
            rmsds.append(rmsds_row)
        rmsds = np.array(rmsds)

        return rmsds

    @requires_dependency(
        [Import("spyrmsd"), Import("spyrmsd.rmsd")], globals())
    def _get_rmsd_spyrmsd(self, i_s, j_s, include_hs, symmetry):
        """Calculate RMSD row-wise with spyrmsd"""
        # Construct mask for H atoms
        if not include_hs:
            mask = self.elements != 1
        else:
            mask = np.ones(len(self.elements), dtype=np.bool)

        # Calculate RMSD row-wise with spyrmsd
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

    def update_mol(self):
        """Update Mol object with conformers"""
        self.mol.RemoveAllConformers()
        conformer_coordinates = self.get_coordinates()
        _add_conformers_to_mol(self.mol, conformer_coordinates)

    def __copy__(self):
        # Generate copy where conformers and mol object are shared with old
        # ensemble.
        cls = type(self)
        ce = cls(
            self.elements,
            connectivity_matrix=self.connectivity_matrix,
            charge=self.charge,
            multiplicity=self.multiplicity,
            formal_charges=self.formal_charges,
            ref_cip_label=self.ref_cip_label,
        )
        ce.mol = self.mol
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
            formal_charges=self.formal_charges,
            ref_cip_label=self.ref_cip_label,
        )
        ce.mol = deepcopy(ce.mol, memo)
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


@requires_dependency([
    Import(module="openbabel.openbabel", alias="ob"),
    Import(module="openbabel.pybel", alias="pybel"),
    Import("openbabel")
], globals())
def generate_conformers_openbabel_ff(mol,
                                     num_conformers=30,
                                     ff="MMFF94",
                                     method="systematic",
                                     rings=False):
    """Generates conformers based on the force field algorithm in OpenBabel.

    Follows the recipe of the command line script obabel --conformer:
    https://github.com/openbabel/openbabel/blob/master/src/ops/conformer.cpp
    If an OBMol object with 3D coordinates is given, the conformer search will
    start from that structure.

    Args:
        mol (str or OBMol): Molecule either as SMILES string or OBMol object.
        num_conformers (int): Maximum number of conformers
        ff (str): Force field supported by OpenBabel
        method (str): 'fast', random', 'systematic' (default) or 'weighted'
        rings (bool): Sample ring torsions.

    """
    # Create 3D structure if not given
    try:
        py_mol = pybel.readstring("smi", mol)
        py_mol.make3D()
        ob_mol = py_mol.OBMol
    except TypeError:
        ob_mol = mol
        py_mol = pybel.Molecule(ob_mol)

    # Set up FF
    ff = ob.OBForceField_FindForceField(ff)
    ff.EnableCutOff(True)
    ff.SetVDWCutOff(10.0)
    ff.SetElectrostaticCutOff(20.0)
    ff.SetUpdateFrequency(10)

    # Do conformer search
    if method == "systematic":
        ff.SystematicRotorSearch(10, rings)
    elif method == "fast":
        ff.FastRotorSearch(True)
    elif method == "random":
        ff.RandomRotorSearch(num_conformers, 10, rings)
    elif method == "weighted":
        ff.WeightedRotorSearch(num_conformers, 10, rings)
    ff.GetConformers(ob_mol)

    # Extract information
    (elements, conformer_coordinates, connectivity_matrix,
     charges) = _extract_from_ob_mol(ob_mol)

    return (elements, conformer_coordinates, connectivity_matrix, charges,
            ob_mol)


@requires_dependency([
    Import(module="openbabel.openbabel", alias="ob"),
    Import(module="openbabel.pybel", alias="pybel"),
    Import("openbabel")
], globals())
def generate_conformers_openbabel_ga(mol,
                                     num_conformers=None,
                                     num_children=None,
                                     mutability=None,
                                     convergence=None,
                                     score="rmsd",
                                     filter_method="steric",
                                     cutoff=0.8,
                                     vdw_factor=0.5,
                                     check_hydrogens=True):
    """Generates conformers based on the genetic algorithm in OpenBabel.

    Follows the recipe of the command line script obabel --conformer:
    https://github.com/openbabel/openbabel/blob/master/src/ops/conformer.cpp
    If an OBMol object with 3D coordinates is given, the conformer search will
    start from that structure.

    Args:
        mol (str or OBMol): Molecule either as SMILES string or OBMol object.
        num_conformers (int): Maximum number of conformers
        num_children (int): Number of children to generate for each parent
        mutability (float): Mutation frequency
        convergence (int): Number of identical generations before convergence
            is reached
        score (str): Scoring function: 'rmsd' (default), 'min_rmsd', 'energy',
            'min_energy'
        filter_method (str): Filtering algorithm: 'steric' (default)
        cutoff (float): Absolute distance in Ånström below which atoms are
            considered to clash.
        vdw_factor (float): Scale factor applied to van der Waals radii for
            detecting clashes."
        check_hydrogens (bool): Detect clashes with hydrogen atoms
    """
    # Create 3D structure if not given
    try:
        py_mol = pybel.readstring("smi", mol)
        py_mol.make3D()
        ob_mol = py_mol.OBMol
    except TypeError:
        ob_mol = mol
        py_mol = pybel.Molecule(ob_mol)

    # Create search object and set parameters
    conf_search = ob.OBConformerSearch()
    conf_search.Setup(ob_mol)
    if num_conformers is not None:
        conf_search.SetNumConformers(num_conformers)
    if num_children is not None:
        conf_search.SetNumChildren(num_children)
    if convergence is not None:
        conf_search.SetConvergence(convergence)
    if mutability is not None:
        conf_search.SetMutability(mutability)
    if score is not None:
        # Scorers don't work with earlier versions of openbabel
        if not (parse_version(openbabel.__version__) > parse_version("3.1.0")):
            warnings.warn(
                "Scorer only works with openbabel version > 3.1.0. "
                "Proceeding without scorer.")
        else:
            if score == "rmsd":
                scorer = ob.OBRMSDConformerScore()
            if score == "min_rmsd":
                scorer = ob.OBMinimizingRMSDConformerScore()
            if score == "energy":
                scorer = ob.OBEnergyConformerScore()
            if score == "min_energy":
                scorer = ob.OBMinimizingEnergyConformerScore()
            conf_search.SetScore(scorer)
    if filter_method == "steric":
        # Filters don't work with earlier versions of openbabel
        if not (parse_version(openbabel.__version__) > parse_version("3.1.0")):
            warnings.warn(
                "Filter only works with openbabel version > 3.1.0. "
                "Proceeding without filter.")
        else:
            ob_filter = ob.OBStericConformerFilter(
                cutoff, vdw_factor, check_hydrogens)
            conf_search.SetFilter(ob_filter)

    # Do conformational search
    conf_search.Search()
    conf_search.GetConformers(ob_mol)

    # Extract information
    (elements, conformer_coordinates, connectivity_matrix,
     charges) = _extract_from_ob_mol(ob_mol)

    return (elements, conformer_coordinates, connectivity_matrix, charges,
            ob_mol)


@requires_dependency([
    Import(module="rdkit", item="Chem"),
    Import(module="rdkit.Chem", item="AllChem")
], globals())
def generate_conformers_rdkit(mol,
                              n_confs=None,
                              optimize=None,
                              version=2,
                              small_rings=True,
                              macrocycles=True,
                              random_seed=None,
                              rmsd_thres=0.35,
                              n_threads=1):
    """Generates conformers for an RDKit mol object. Recipe based on
    J. Chem. Inf. Modeling 2012, 52, 1146.

    Args:
        mol (str or Mol): Molecule either as SMILES string or RDKit Mol object.
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

    Returns:
        elements (list): Atomic symbols
        conformer_coordinates (ndarray): Coordinates for all conformers (Å)
        energies (list): Conformer energies (a.u.)
        connectivity_matrix (ndarray): Connectivity matrix with bond orders
        mol (obj): RDKit Mol object. Only returned when `return_mol=True`.
    """
    if optimize not in ("MMFF94", "MMFF94s", "UFF", None):
        raise Exception(f"Force field {optimize} not found. Choose one of"
                        "MMFF94, MMFF94s, UFF.")
    # Generate mol object
    try:
        mol = Chem.MolFromSmiles(mol)
    except TypeError:
        pass
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

    # Extract information from mol
    (elements, conformer_coordinates, connectivity_matrix,
     charges) = _extract_from_mol(mol)

    return (elements, conformer_coordinates, energies, connectivity_matrix,
            charges, mol)


@requires_dependency([
    Import(module="rdkit"),
    Import(module="rdkit", item="Chem"),
], globals())
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


@requires_dependency([Import(module="rdkit", item="Chem")], globals())
def _extract_from_mol(mol):
    """Extract information from RDKit Mol object with conformers."""
    # Take out elements, coordinates and connectivity matrix
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    charges = [atom.GetFormalCharge() for atom in mol.GetAtoms()]

    conformer_coordinates = []
    for conformer in mol.GetConformers():
        coordinates = conformer.GetPositions()
        conformer_coordinates.append(coordinates)
    conformer_coordinates = np.array(conformer_coordinates)

    connectivity_matrix = Chem.GetAdjacencyMatrix(mol, useBO=True)

    return elements, conformer_coordinates, connectivity_matrix, charges


@requires_dependency([
    Import(module="openbabel.openbabel", alias="ob"),
    Import(module="openbabel.pybel", alias="pybel"),
    Import("openbabel")
], globals())
def _extract_from_ob_mol(ob_mol):
    """Extract information from Openbabel OBMol object with conformers."""
    py_mol = pybel.Molecule(ob_mol)
    elements = np.array([atom.atomicnum for atom in py_mol.atoms])
    charges = np.array([atom.formalcharge for atom in py_mol.atoms])

    n_atoms = len(py_mol.atoms)
    connectivity_matrix = np.zeros((n_atoms, n_atoms))
    for bond in ob.OBMolBondIter(ob_mol):
        i = bond.GetBeginAtomIdx() - 1
        j = bond.GetEndAtomIdx() - 1
        bo = bond.GetBondOrder()
        connectivity_matrix[i, j] = bo
        connectivity_matrix[j, i] = bo

    # Retrieve conformer coordinates
    conformer_coordinates = []
    for i in range(ob_mol.NumConformers()):
        ob_mol.SetConformer(i)
        coordinates = np.array([atom.coords for atom in py_mol.atoms])
        conformer_coordinates.append(coordinates)
    conformer_coordinates = np.array(conformer_coordinates)

    return elements, conformer_coordinates, connectivity_matrix, charges


@requires_dependency([Import(module="openbabel.openbabel", alias="ob")],
                     globals())
def _get_ob_mol(elements, coordinates, connectivity_matrix, charges=None):
    """Generate OpenBabel OBMol object.

    Args:
        elements (list): Elements as atomic symbols or numbers.
        coordinates (list): Coordinates (Å)
        connectivity_matrix (ndarray): Connectivity matrix with bond orders

    Returns:
        mol (obj): OpenBabel OBMol object
    """
    elements = convert_elements(elements)
    if charges is None:
        charges = np.zeros(len(elements))

    mol = ob.OBMol()

    # Add atoms
    for element, coordinate, charge in zip(elements, coordinates, charges):
        a = mol.NewAtom()
        a.SetAtomicNum(int(element))
        a.SetVector(*coordinate)
        a.SetFormalCharge(int(charge))

    # Add bonds
    for i, j in zip(*np.tril_indices_from(connectivity_matrix)):
        if i != j:
            bo = connectivity_matrix[i, j]
            if bo != 0:
                mol.AddBond(int(i + 1), int(j + 1), int(bo))
    return mol


@requires_dependency([Import(module="rdkit", item="Chem")], globals())
def _get_rdkit_mol(elements, coordinates, connectivity_matrix, charges=None):
    _RDKIT_BOND_TYPES = {
        1.0: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2.0: Chem.BondType.DOUBLE,
        3.0: Chem.BondType.TRIPLE,
        4.0: Chem.BondType.QUADRUPLE,
        5.0: Chem.BondType.QUINTUPLE,
        6.0: Chem.BondType.HEXTUPLE,
    }
    if charges is None:
        charges = np.zeros(len(elements))
    elements = convert_elements(elements, output="symbols")

    mol = Chem.RWMol()

    # Add atoms
    for element, charge in zip(elements, charges):
        atom = Chem.Atom(element)
        atom.SetFormalCharge(int(charge))
        mol.AddAtom(atom)

    # Add bonds
    for i, j in zip(*np.tril_indices_from(connectivity_matrix)):
        if i != j:
            bo = connectivity_matrix[i, j]
            if bo != 0:
                bond_type = _RDKIT_BOND_TYPES[float(bo)]
                mol.AddBond(int(i), int(j), bond_type)

    # Add conformers
    _add_conformers_to_mol(mol, coordinates)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)

    return mol
