"""Conformer tools."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from copy import copy, deepcopy
import functools
import numbers
from os import PathLike
from pathlib import Path
import subprocess
import tempfile
import typing
from typing import Any, cast, Type, TypeVar
import warnings

import numpy as np
from packaging.version import parse

from morfeus.data import (
    HARTREE,
    HARTREE_TO_KCAL,
    HARTREE_TO_KJ,
    K_B,
    KCAL_TO_HARTREE,
    KJ_TO_HARTREE,
)
from morfeus.io import CrestParser, write_xyz
from morfeus.qc import optimize_qc_engine, sp_qc_engine
from morfeus.typing import (
    Array1DFloat,
    Array1DInt,
    Array2DFloat,
    Array2DInt,
    Array3DFloat,
    ArrayLike1D,
    ArrayLike2D,
    ArrayLike3D,
)
from morfeus.utils import (
    convert_elements,
    Import,
    requires_dependency,
    requires_executable,
)

if typing.TYPE_CHECKING:
    import openbabel
    import openbabel.openbabel as ob
    import openbabel.pybel as pybel
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdDistGeom
    import spyrmsd


def boltzmann_average_dT(
    properties: ArrayLike1D, energies: ArrayLike1D, temperature: float = 298.15
) -> float:
    """Return the derivative of the Boltzmann average.

    Args:
        properties: Conformer properties
        energies: Conformer energies (a.u.)
        temperature: Temperature (K)

    Returns:
        derivative: Derivative of Boltzmann average.
    """
    energies: Array1DFloat = np.array(energies)
    properties: Array1DFloat = np.array(properties)

    # Calculate Boltzmann averaged properties
    avg_prop_en = boltzmann_statistic(
        properties * energies * HARTREE,
        energies,
        temperature=temperature,
        statistic="avg",
    )
    avg_en = boltzmann_statistic(
        energies * HARTREE, energies, temperature=temperature, statistic="avg"
    )
    avg_prop = boltzmann_statistic(
        properties, energies, temperature=temperature, statistic="avg"
    )

    # Calculate derivative
    derivative = (avg_prop_en - avg_en * avg_prop) / (K_B * temperature**2)

    return derivative


def boltzmann_statistic(
    properties: ArrayLike1D,
    energies: ArrayLike1D,
    temperature: float = 298.15,
    statistic: str = "avg",
) -> float:
    """Compute Boltzmann statistic.

    Args:
        properties: Conformer properties
        energies: Conformer energies (a.u.)
        temperature: Temperature (K)
        statistic: Statistic to compute: 'avg', 'var' or 'std'

    Returns:
        result: Boltzmann statistic

    Raises:
        ValueError: When statistic not specified correctly
    """
    properties: Array1DFloat = np.array(properties)

    # Get conformer weights
    weights = boltzmann_weights(energies, temperature)

    # Compute Boltzmann weighted statistic
    result: float
    if statistic == "avg":
        result = np.average(properties, weights=weights, axis=0)
    elif statistic == "var":
        avg = boltzmann_statistic(
            properties, energies, temperature=temperature, statistic="avg"
        )
        result = np.average((properties - avg) ** 2, weights=weights, axis=0)
    elif statistic == "std":
        avg = boltzmann_statistic(
            properties, energies, temperature=temperature, statistic="avg"
        )
        var = boltzmann_statistic(
            properties, energies, temperature=temperature, statistic="var"
        )
        result = np.sqrt(var)
    else:
        raise ValueError("Choose between: 'avg', 'var' and 'std'")

    return result


def boltzmann_weights(
    energies: ArrayLike1D, temperature: float = 298.15
) -> Array1DFloat:
    """Compute Boltzmann weights.

    Args:
        energies: Conformer energies (a.u.)
        temperature: Temperature (K)

    Returns:
        weights: Conformer weights (normalized to unity)
    """
    energies: Array1DFloat = np.array(energies)
    energies -= energies.min()
    terms = np.exp(-energies / (K_B / HARTREE * temperature))
    weights: Array1DFloat = terms / np.sum(terms)

    return weights


class Conformer:
    """Conformer with coordinates, energy and properties.

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        energy: Energy (a.u.)
        degeneracy: Degeneracy
        properties: Conformers properties.
        cip_label: Tuple of CIP labels for all atoms: 'R', 'S' or ''

    Attributes:
        cip_label: Tuple of CIP labels for all atoms: 'R', 'S' or ''
        coordinates: Coordinates (Å)
        degeneracy: Degeneracy
        elements: Elements as atomic symbols or numbers
        energy: Energy (a.u.)
        properties: Conformers properties.
    """

    cip_label: tuple[str, ...] | None
    coordinates: Array2DFloat
    degeneracy: int
    elements: Array1DFloat
    energy: float | None
    properties: dict[str, float]

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        coordinates: ArrayLike2D,
        energy: float | None = None,
        degeneracy: int = 1,
        properties: dict[str, float] | None = None,
        cip_label: tuple[str] | None = None,
    ) -> None:
        if properties is None:
            properties = {}
        self.coordinates = np.array(coordinates)
        self.elements = np.array(elements)
        self.energy = energy
        self.properties = properties
        self.degeneracy = degeneracy
        self.cip_label = cip_label

    def __repr__(self) -> str:
        n_atoms = self.coordinates.shape[0]
        return f"{self.__class__.__name__}({n_atoms!r} atoms)"


T = TypeVar("T", bound="ConformerEnsemble")


class ConformerEnsemble:
    """Conformer ensemble object.

    Supports sorting, pruning, optimization and single-point calculations.

    Args:
        elements: Elements as atomic symbols or numbers
        conformer_coordinates: Conformer coordinates (Å)
        energies: Energy (a.u.)
        connectivity_matrix: Connectivity matrix
        degeneracies: Degeneracies
        properties: Conformers properties.
        charge: Molecular charge
        multiplicity: Molecular multiplicity.
        formal_charges: Atomic formal charges
        ref_cip_label: Tuple of CIP labels for all atoms: 'R', 'S' or ''

    Attributes:
        charge: Molecular charge
        conformers: Conformers
        connectivity_matrix: Connectivity matrix
        elements: Elements as atomic symbols or numbers
        formal_charges: Atomic formal charges.
        mol: RDKit mol object.
        multiplicity: Molecular multiplicity.
        ref_cip_label: Tuple of CIP labels for all atoms: 'R', 'S' or ''
    """

    charge: int
    confomers: list[Conformer]
    connectivity_matrix: Array2DInt
    elements: Array1DInt
    formal_charges: Array1DInt
    mol: Chem.Mol
    multiplicity: int
    ref_cip_label: tuple[str, ...] | None

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str],
        conformer_coordinates: ArrayLike3D | None = None,
        energies: ArrayLike1D | None = None,
        connectivity_matrix: ArrayLike2D | None = None,
        degeneracies: ArrayLike1D | None = None,
        properties: Mapping[str, ArrayLike1D] | None = None,
        charge: int | None = None,
        multiplicity: int | None = None,
        formal_charges: ArrayLike1D | None = None,
        ref_cip_label: tuple[str, ...] | None = None,
    ) -> None:
        elements: Array1DInt = np.array(convert_elements(elements, output="numbers"))
        self.elements = elements

        if conformer_coordinates is None:
            conformer_coordinates = []
        if formal_charges is None:
            formal_charges = np.zeros(len(elements), dtype=int)
        if charge is None:
            charge = int(np.sum(formal_charges))
        else:
            if charge != np.sum(formal_charges):
                msg = (
                    f"Charge ({charge}) is different from sum of partial "
                    f"charges ({np.sum(formal_charges)})"
                )
                raise Exception(msg)
        if multiplicity is None:
            multiplicity = 1

        # Store conformers
        self.conformers: list[Conformer] = []
        self._add_conformers(conformer_coordinates, energies, degeneracies, properties)

        # Set up attributes
        self.charge = charge
        self.formal_charges = np.array(formal_charges)
        self.multiplicity = multiplicity
        if connectivity_matrix is not None:
            self.connectivity_matrix = np.array(connectivity_matrix)
        else:
            self.connectivity_matrix = None
        self.mol = None
        self.ref_cip_label = ref_cip_label

    def add_conformers(
        self,
        conformer_coordinates: ArrayLike3D,
        energies: ArrayLike1D | None = None,
        degeneracies: ArrayLike1D | None = None,
        properties: Mapping[str, ArrayLike1D] | None = None,
    ) -> ConformerEnsemble:
        """Add conformer to ensemble.

        Args:
            conformer_coordinates: Conformer coordinates (Å)
            energies: Energies (a.u.)
            degeneracies: Degeneracies
            properties: Conformer properties.

        Returns:
            self: Self
        """
        self._add_conformers(conformer_coordinates, energies, degeneracies, properties)

        return self

    def add_inverted(self) -> ConformerEnsemble:
        """Add inverted images of all conformers.

        Scrambles stereochemistry and leads to redundant conformers so use with
        care and prune based on RMSD as postprocessing.

        Returns:
            self: Self
        """
        conformers = self.conformers
        conformer_coordinates: list[Array2DFloat] = []
        energies: list[float | None] = []
        for conformer in conformers:
            coordinates = conformer.coordinates * np.array([-1, -1, -1])
            conformer_coordinates.append(coordinates)
            energies.append(conformer.energy)

        if not all([isinstance(energy, float) for energy in energies]):
            self._add_conformers(conformer_coordinates)
        else:
            energies_ = cast("list[float]", energies)
            self.add_conformers(conformer_coordinates, energies_)

        return self

    @requires_dependency([Import(module="rdkit.Chem", item="AllChem")], globals())
    def align_conformers(self) -> ConformerEnsemble:
        """Align conformers with RDKit."""
        self.update_mol()
        AllChem.AlignMolConformers(self.mol)

        return self

    def boltzmann_average_dT(
        self, property_name: str, temperature: float = 298.15
    ) -> float:
        """Calculate temperature derivative of Boltzmann average of property.

        Args:
            property_name: Name of property
            temperature: Temperature (K)

        Returns:
            derivative: Derivative of Boltzmann average
        """
        energies = self.get_energies()
        properties = self.get_properties()[property_name]
        derivative = boltzmann_average_dT(properties, energies, temperature=temperature)

        return derivative

    def boltzmann_statistic(
        self, property_name: str, temperature: float = 298.15, statistic: str = "avg"
    ) -> float:
        """Calculate Boltzmann staistic of property over ensemble.

        Args:
            property_name: Name of property
            statistic: Boltzmann statistic: 'avg', 'var' or 'std'
            temperature: Temperature (K)

        Returns:
            statistic: Boltzmann statistic
        """
        energies = self.get_energies()
        properties = self.get_properties()[property_name]
        statistic = boltzmann_statistic(
            properties, energies, temperature=temperature, statistic=statistic
        )

        return statistic

    def boltzmann_weights(self, temperature: float = 298.15) -> Array1DFloat:
        """Calculate Boltzmann weights for ensemble.

        Args:
            temperature: Temperature (K)

        Returns:
            weights: Conformer weights (normalized to unity)
        """
        energies = self.get_energies()
        weights = boltzmann_weights(energies, temperature)

        return weights

    def detect_enantiomers(
        self, thres: float = 0.01, method: str = "rdkit", include_hs: bool = False
    ) -> dict[int, list[int]]:
        """Detect enantiomers in ensemble.

        Args:
            thres: RMSD threshold for detecting enantiomers in terms of coordinates
            method: RMSD calculation method: 'obrms-batch', 'obrms-iter', openbabel' or
                'spyrmsd'
            include_hs: Whether to include H atoms when determining enantiomers

        Returns:
            enantiomers: Mapping of enantiomer with conformer id as keys and enantiomer
                ids as values.
        """
        # Add set of inverted conformers
        n_conformers = self.n_conformers
        self.add_inverted()

        # Map conformers to enantiomers
        enantiomer_map: dict[int, set[int]] = {i: set() for i in range(n_conformers)}
        rmsds = self.get_rmsd(method=method, include_hs=include_hs, symmetry=False)
        for i in range(n_conformers):
            # Do test that inverted conformer should have 0 RMSD to be
            # enantiomer
            test = np.where(rmsds[i, n_conformers:] < thres)[0]
            if len(test) > 0:
                if i not in test:
                    enantiomer_map[i].update(test)
                    for j in test:
                        enantiomer_map[j].add(i)
        enantiomers = {key: list(value) for key, value in enantiomer_map.items()}

        # Reset number of conformers
        self.conformers = self.conformers[:n_conformers]

        return enantiomers

    @classmethod
    def from_crest(cls: Type[T], path: str | PathLike) -> T:
        """Generate conformer ensemble from CREST output.

        Args:
            path: Path to CREST folder

        Returns:
            ce: Conformer ensemble object.
        """
        cp = CrestParser(path)
        ce = cls(
            cp.elements,
            cp.conformer_coordinates,
            energies=cp.energies,
            degeneracies=cp.degeneracies,
        )

        return ce

    @classmethod
    def from_rdkit(cls: Type[T], *args: Any, **kwargs: Any) -> T:
        """Generate conformer ensemble from RDKit.

        See the documentation for the function conformers_from_rdkit for
        more information.

        Args:
            *args: Positional arguments
            **kwargs: Keyword arguments

        Returns:
            ce: Conformer ensemble object.
        """
        # Run RDKit conformer search and generate ensemble.
        (
            elements,
            conformer_coordinates,
            energies,
            connectivity_matrix,
            charges,
            mol,
        ) = conformers_from_rdkit(*args, **kwargs)
        ce = cls(
            elements,
            conformer_coordinates,
            energies=energies,
            connectivity_matrix=connectivity_matrix,
            formal_charges=charges,
        )
        ce.mol = mol
        ce.set_multiplicity_from_mol()

        # Set reference CIP label if enantiomerically pure.
        cip_labels = ce.get_cip_labels()
        if len(set(cip_labels)) == 1:
            ce.ref_cip_label = cip_labels[0]

        return ce

    @classmethod
    def from_ob_ga(
        cls: Type[T],
        *args: Any,
        generate_rdkit_mol: bool = False,
        update_charges: bool = True,
        update_connectivity: bool = True,
        update_multiplicity: bool = True,
        **kwargs: Any,
    ) -> T:
        """Generate conformer ensemble from OpenBabel with GA method.

        See the documentation for the function conformers_from_ob_ga
        for more information.

        Args:
            *args: Positional arguments for conformers_from_ob_ga
            generate_rdkit_mol: Generate RDKit mol object for ensemble
            update_charges: Update formal charges from generated RDKit Mol object. Only
                used if generate_rdkit_mol is True.
            update_connectivity: Update connectivity from generated RDKit Mol object.
                Only used if generate_rdkit_mol is True.
            update_multiplicity: Update multiplicity from generated RDKit Mol object.
                Only used if generate_rdkit_mol is True.
            **kwargs: Keyword arguments for conformers_from_ob_ga

        Returns:
            ce: Conformer ensemble.
        """
        # Run Openbabel conformer search and generate ensemble.
        (
            elements,
            conformer_coordinates,
            connectivity_matrix,
            charges,
            ob_mol,
        ) = conformers_from_ob_ga(*args, **kwargs)
        ce = cls(
            elements,
            conformer_coordinates,
            connectivity_matrix=connectivity_matrix,
            formal_charges=charges,
        )

        # Generate RDKit mol object and CIP labels
        if generate_rdkit_mol:
            ce.generate_mol(
                update_charges=update_charges,
                update_connectivity=update_connectivity,
                update_multiplicity=update_multiplicity,
            )

            cip_labels = ce.get_cip_labels()
            if len(set(cip_labels)) == 1:
                ce.ref_cip_label = cip_labels[0]

        return ce

    @classmethod
    def from_ob_ff(
        cls: Type[T],
        *args: Any,
        generate_rdkit_mol: bool = False,
        update_charges: bool = True,
        update_connectivity: bool = True,
        update_multiplicity: bool = True,
        **kwargs: Any,
    ) -> T:
        """Generate conformer ensemble from OpenBabel with FF method.

        See the documentation for the function conformers_from_ob_ff
        for more information.

        Args:
            *args: Positional arguments for conformers_from_ob_ga
            generate_rdkit_mol: Generate RDKit mol object for ensemble
            update_charges: Update formal charges from generated RDKit Mol object. Only
                used if generate_rdkit_mol is True.
            update_connectivity: Update connectivity from generated RDKit Mol object.
                Only used if generate_rdkit_mol is True.
            update_multiplicity: Update multiplicity from generated RDKit Mol object.
                Only used if generate_rdkit_mol is True.
            **kwargs: Keyword arguments for conformers_from_ob_ga

        Returns:
            ce: Conformer ensemble.
        """
        # Run Openbabel conformer search and generate ensemble.
        (
            elements,
            conformer_coordinates,
            connectivity_matrix,
            charges,
            ob_mol,
        ) = conformers_from_ob_ff(*args, **kwargs)
        ce = cls(
            elements,
            conformer_coordinates,
            connectivity_matrix=connectivity_matrix,
            formal_charges=charges,
        )

        # Generate RDKit mol object and CIP labels
        if generate_rdkit_mol:
            ce.generate_mol(
                update_charges=update_charges,
                update_connectivity=update_connectivity,
                update_multiplicity=update_multiplicity,
            )

            cip_labels = ce.get_cip_labels()
            if len(set(cip_labels)) == 1:
                ce.ref_cip_label = cip_labels[0]

        return ce

    @requires_dependency([Import(module="rdkit", item="Chem")], globals())
    def generate_mol(
        self,
        update_charges: bool = True,
        update_connectivity: bool = True,
        update_multiplicity: bool = True,
    ) -> None:
        """Generate RDKit Mol object.

        Args:
            update_charges: Update formal charges from generated RDKit Mol object
            update_connectivity: Update connectivity from generated RDKit Mol object
            update_multiplicity: Update multiplicity from generated RDKit Mol object
        """
        mol = _get_rdkit_mol(
            self.elements,
            self.get_coordinates(),
            self.connectivity_matrix,
            self.formal_charges,
        )
        self.mol = mol

        if update_charges:
            self.formal_charges = np.array(
                [atom.GetFormalCharge() for atom in mol.GetAtoms()]
            )
        if update_connectivity:
            self.connectivity_matrix = np.array(
                Chem.GetAdjacencyMatrix(mol, useBO=True)
            )
        if update_multiplicity:
            self.set_multiplicity_from_mol()

    @requires_dependency([Import(module="rdkit", item="Chem")], globals())
    def get_cip_labels(self) -> list[tuple[str, ...]]:
        """Generate tuples of CIP labels for conformer.

        Returns:
            cip_labels: Tuples of CIP labels for each conformer.
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

    def get_coordinates(self) -> Array3DFloat:
        """Get conformer coordinates.

        Returns:
            conformer_coordinates: Conformer coordinates (Å)
        """
        conformer_coordinates: Array3DFloat = np.array(
            [conformer.coordinates for conformer in self.conformers]
        )

        return conformer_coordinates

    def get_degeneracies(self) -> Array1DInt:
        """Get conformer degeneracies.

        Returns:
            degeneriacies: Degeneracies
        """
        degeneracies: Array1DInt = np.array(
            [conformer.degeneracy for conformer in self.conformers]
        )

        return degeneracies

    def get_energies(self) -> Array1DFloat:
        """Get conformer energies.

        Returns:
            energies: Energy (a.u.)
        """
        energies: Array1DFloat = np.array(
            [conformer.energy for conformer in self.conformers]
        )
        return energies

    def get_properties(self) -> dict[str, Array1DFloat]:
        """Get conformer properties.

        Returns:
            properties: Conformer properties
        """
        properties: dict[str, list[float]] = {}
        for conformer in self.conformers:
            for key, value in conformer.properties.items():
                properties.setdefault(key, []).append(value)
        properties: dict[str, Array1DFloat] = {
            key: np.array(value) for key, value in properties.items()
        }

        return properties

    def get_relative_energies(
        self, unit: str = "kcal/mol", relative: bool = True
    ) -> Array1DFloat:
        """Get conformer energies with choice of units and reference value.

        Args:
            unit: Unit of returned energies: 'hartree', 'kcal/mol' or 'kJ/mol'.
            relative: Return energies relative to lowest-energy conformer.

        Returns:
            energies: Conformer energies.
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

    def get_rmsd(
        self,
        i_s: ArrayLike1D | None = None,
        j_s: ArrayLike1D | None = None,
        include_hs: bool = False,
        symmetry: bool = False,
        method: str = "rdkit",
    ) -> Array2DFloat:
        """Get RSMD between conformers.

        For very small systems 'openbabel' or 'spyrmsd' work well. For larger systems a
        significant speed-up is attained with 'rdkit', 'obrms-batch' or 'obrms-iter'.

        Args:
            i_s: Indices of conformers
            j_s: Indices of conformers
            include_hs: Whether to include H atoms in RMSD calculation. Ignored for
                'obrms-iter' and 'obrms-batch' that only use heavy atoms.
            symmetry: Consider symmetry (requires connectivity matrix). Ignored for
                'obrms-iter' and 'obrms-batch' that always use symmetry.
            method: RMSD calculation method: 'obrms-batch', 'obrms-iter', 'openbabel',
                'rdkit' or 'spyrmsd'.

        Returns:
            rmsds: RSMDs (Å)

        Raises:
            ValueError: If method not supported.
        """
        i_s_: np.ndarray
        if i_s is None:
            i_s_ = np.arange(1, len(self.conformers) + 1)
        else:
            i_s_ = np.array(i_s)
        i_s = i_s_
        j_s_: np.ndarray
        if j_s is None:
            j_s_ = np.arange(1, len(self.conformers) + 1)
        else:
            j_s_ = np.array(j_s)
        j_s = j_s_

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
        else:
            raise ValueError("Method not supported.")
        return rmsds

    @property
    def n_conformers(self) -> int:
        """Number of conformers."""
        n_conformers = len(self.conformers)

        return n_conformers

    def optimize_qc_engine(
        self,
        ids: Sequence[int] | None = None,
        program: str | None = None,
        model: dict[str, Any] | None = None,
        keywords: dict[str, Any] | None = None,
        local_options: dict[str, Any] | None = None,
        procedure: str = "berny",
    ) -> ConformerEnsemble:
        """Optimize conformers with QCEngine interface.

        Args:
            ids: Conformer indices (1-indexed). If None, all are optimized
            program: QCEngine program
            model: QCEngine model
            keywords: QCEngine keywords
            local_options: QCEngine local options
            procedure: QCEngine procedure

        Returns:
            self: Self

        Raises:
            Exception: When RDKit requested without formal charges.
        """
        if program is not None:
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

        return self

    def condense_enantiomeric(self, thres: float | None = None) -> ConformerEnsemble:
        """Condense enantiomers into single enantiomer per pair.

        Args:
            thres: RMSD threshold for assessing enantiomers.

        Returns:
            self: Self

        Raises:
            Exception: If molecule is chiral.
        """
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

        return self

    def prune_enantiomers(
        self, keep: str = "original", ref_label: tuple[str, ...] | None = None
    ) -> ConformerEnsemble:
        """Prune conformers so that only one enantiomer is present in the ensemble.

        Args:
            keep: Which enantiomer to keep: 'original', 'most common' or 'specified'.
                Choice of 'original' requires that the ref_cip_label attribute is set.
                Choice of 'specified' requires ref_label to be given.
            ref_label: Reference CIP labels for all atoms

        Returns:
            self: Self
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
            counter: Counter[tuple[str, ...]] = Counter(cip_labels)
            ref_label = counter.most_common(n=1)[0][0]

        # Prune conformers
        to_keep = []
        for i, cip_label in enumerate(cip_labels):
            if cip_label == ref_label:
                to_keep.append(i)
        self.conformers = [self.conformers[i] for i in to_keep]

        return self

    def prune_energy(
        self, threshold: float = 3.0, unit: str = "kcal/mol"
    ) -> ConformerEnsemble:
        """Prune conformers based on energy compared to minimum energy conformer.

        Args:
            threshold: Energy threshold for pruning
            unit: Unit for energy threshold 'hartree', 'kcal/mol' or 'kJ/mol'

        Returns:
            self: Self
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

        return self

    def prune_rmsd(
        self,
        thres: float = 0.35,
        include_hs: bool = False,
        symmetry: bool = False,
        method: str = "rdkit",
    ) -> ConformerEnsemble:
        """Prune conformers based on RMSD.

        Args:
            thres: Threshold for RSMD pruning (Å)
            include_hs: Whether to include H atoms in RMSD calculation. Ignored for
                'obrms-iter' and 'obrms-batch' which only use heavy atoms.
            symmetry: Consider symmetry (requires connectivity matrix). Ignored for
                'obrms-iter' and 'obrms-batch' which always use symmetry.
            method: RMSD calculation method: 'obrms-batch', 'obrms-iter', 'openbabel',
                'rdkit' or 'spyrmsd'

        Returns:
            self: Self
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
                    candidates + 1,
                    include_hs=include_hs,
                    symmetry=symmetry,
                    method=method,
                )
                candidates = candidates[rmsd[0] > thres]
        elif method == "obrms-batch":
            rmsds = self.get_rmsd(
                include_hs=include_hs, symmetry=symmetry, method=method
            )
            working_array = rmsds
            while len(working_array) > 0:
                keeper = candidates[0]
                keep_list.append(keeper)
                rmsd = working_array[0]
                mask = rmsd > thres
                candidates = candidates[mask]
                working_array = working_array[mask, :][:, mask]

        # Update conformer list
        self.conformers = [
            conformer for i, conformer in enumerate(self.conformers) if i in keep_list
        ]

        return self

    def set_coordinates(self, conformer_coordinates: ArrayLike3D) -> ConformerEnsemble:
        """Set conformer coordinates.

        Args:
            conformer_coordinates: Conformer coordinates (Å)

        Returns:
            self: Self

        Raises:
            ValueError: When number of conformer coordinates is different from number of
                conformers
        """
        conformer_coordinates: Array3DFloat = np.array(conformer_coordinates)
        if len(conformer_coordinates) != self.n_conformers:
            msg = (
                f"Number of coordinates ({len(conformer_coordinates)}) "
                f"!= number of conformers ({self.n_conformers})."
            )
            raise ValueError(msg)
        for conformer, coordinates in zip(self.conformers, conformer_coordinates):
            conformer.coordinates = np.array(coordinates)

        return self

    def set_degeneracies(self, degeneracies: ArrayLike1D) -> ConformerEnsemble:
        """Set degeneracies.

        Args:
            degeneracies: Degeneracies

        Returns:
            self: Self

        Raises:
            ValueError: When number of degeneracies is different from number of
                conformers
        """
        degeneracies: Array1DInt = np.array(degeneracies)
        if len(degeneracies) != self.n_conformers:
            msg = (
                f"Number of degeneracies ({len(degeneracies)}) "
                f"!= number of conformers ({self.n_conformers})."
            )
            raise ValueError(msg)
        for conformer, degeneracy in zip(self.conformers, degeneracies):
            conformer.degeneracy = degeneracy

        return self

    def set_energies(self, energies: ArrayLike1D) -> ConformerEnsemble:
        """Set energies.

        Args:
            energies: Energy (a.u.)

        Returns:
            self: Self

        Raises:
            ValueError: When number of energies is different from number of
                conformers
        """
        energies: Array1DFloat = np.array(energies)
        if len(energies) != self.n_conformers:
            msg = (
                f"Number of energies ({len(energies)}) != number of "
                f"conformers ({self.n_conformers})."
            )
            raise ValueError(msg)
        for conformer, energy in zip(self.conformers, energies):
            conformer.energy = energy

        return self

    @requires_dependency([Import(module="rdkit.Chem", item="Descriptors")], globals())
    def set_multiplicity_from_mol(self) -> ConformerEnsemble:
        """Sets multiplicity based on unpaired electrons in Mol object."""
        num_radical = Descriptors.NumRadicalElectrons(self.mol)

        # Assume maximum spin pairing
        if (num_radical % 2) == 0:
            multiplicity = 1
        else:
            multiplicity = 2

        self.multiplicity = multiplicity

        return self

    def set_properties(self, key: str, values: Iterable[float]) -> ConformerEnsemble:
        """Set conformer properties.

        Args:
            key: Name of property
            values: Property values

        Returns:
            self: Self
        """
        for conformer, value in zip(self.conformers, values):
            conformer.properties[key] = value

        return self

    def sort(self) -> ConformerEnsemble:
        """Sort conformers based on energy."""
        energies = [conformer.energy for conformer in self.conformers]
        if not all([isinstance(energy, float) for energy in energies]):
            raise ValueError("Not all conformers have energies.")
        energies = cast("list[float]", energies)
        indices = np.argsort(energies)
        self.conformers = [self.conformers[i] for i in indices]

        return self

    def sp_qc_engine(
        self,
        ids: Sequence[int] | None = None,
        program: str = "xtb",
        model: dict[str, Any] | None = None,
        keywords: dict[str, Any] | None = None,
        local_options: dict[str, Any] | None = None,
    ) -> ConformerEnsemble:
        """Calculate conformer energies with QCEngine interface.

        Args:
            ids: Conformer indices (1-indexed). If None, all are calculated.
            program: QCEngine program
            model: QCEngine model
            keywords: QCEngine keywords
            local_options: QCEngine local options

        Returns:
            self: Self

        Raises:
            Exception: When trying to use RDKit with formal charges.
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

        return self

    def write_xyz(
        self,
        file: str | PathLike,
        ids: Iterable[int] | None = None,
        unit: str = "kcal/mol",
        relative: bool = True,
        separate: bool = False,
        n_decimals: int = 3,
    ) -> None:
        """Write conformers to xyz file.

        Args:
            file: Filename or path object. Needs filename if `separate=True`
            ids: Conformer indices (1-indexed)
            unit: Output unit for energies in xyz comment field: 'hartree', 'kcal/mol',
                'kJ/mol'
            relative: Whether to give energies relative to lowest energy conformer
            separate: Whether to write conformers to separate xyz files
            n_decimals: Number of decimals for energies

        Raises:
            TypeError: When separate=True and file is not str
        """
        ids_: np.ndarray
        if ids is None:
            ids_ = np.arange(len(self.conformers))
        else:
            ids_ = np.array(ids) - 1
        ids = ids_

        # Retrieve symbols, coordinates and energies
        symbols = convert_elements(self.elements, output="symbols")
        conformer_coordinates = self.get_coordinates()[ids]
        energies = self.get_relative_energies(unit=unit, relative=relative)[ids].round(
            n_decimals
        )

        # Write conformers
        if separate:
            if not isinstance(file, str):
                raise TypeError("file must be str when separate=True")
            for i, coordinates, energy in zip(ids, conformer_coordinates, energies):
                conf_filename = file.split(".")[0] + f"_{i + 1}.xyz"
                write_xyz(conf_filename, symbols, coordinates, comments=[energy])
        else:
            write_xyz(file, symbols, conformer_coordinates, comments=energies)

    def _add_conformers(
        self,
        conformer_coordinates: ArrayLike3D,
        energies: ArrayLike1D | None = None,
        degeneracies: ArrayLike1D | None = None,
        properties: Mapping[str, ArrayLike1D] | None = None,
    ) -> None:
        conformer_coordinates: Array3DFloat = np.array(conformer_coordinates)
        n_conformers = len(conformer_coordinates)

        energies_: np.ndarray
        if energies is None:
            energies_ = np.full(n_conformers, np.nan)
        else:
            energies_ = np.array(energies)
        energies = energies_

        degeneracies_: np.ndarray
        if degeneracies is None:
            degeneracies_ = np.ones(n_conformers)
        else:
            degeneracies_ = np.array(degeneracies)
        degeneracies = degeneracies_

        for coordinates, energy, degeneracy in zip(
            conformer_coordinates, energies, degeneracies
        ):
            conformer = Conformer(self.elements, coordinates, energy, degeneracy)
            self.conformers.append(conformer)
        if properties is not None:
            properties_: dict[str, Array1DFloat] = {
                key: np.array(value) for key, value in properties.items()
            }
            for key, value in properties_.items():
                self.set_properties(key, value)

    @requires_executable(["obrms"])
    def _get_rmsd_obrms_batch(self, i_s: Array1DInt, j_s: Array1DInt) -> Array2DFloat:
        """Calculate RMSD with obrms in batch mode.

        First calculates matrix of all pairwise RMSDs and then takes those of interest.

        Args:
            i_s: Conformer indices (1-indexed)
            j_s: Conformer indices (1-indexed)

        Returns:
            rmsds: RMSDs (Å)
        """
        with tempfile.NamedTemporaryFile(suffix=".xyz") as ref_file:
            p_ref = Path(ref_file.name)
            p_ref.unlink()
            self.write_xyz(p_ref, unit="hartree", relative=False)
            process = subprocess.run(
                f"obrms {p_ref.as_posix()} --cross " "--minimize".split(" "),
                capture_output=True,
            )
        result = np.genfromtxt(process.stdout.splitlines(), delimiter=",")
        # Reshape array to 2D if 1D
        if len(result.shape) == 1:
            result = result.reshape(1, -1)
        rmsds = result[:, 1:]

        rmsds: Array2DFloat = rmsds[i_s - 1, :][:, j_s - 1]

        return rmsds

    @requires_executable(["obrms"])
    def _get_rmsd_obrms_iter(self, i_s: Array1DInt, j_s: Array1DInt) -> Array2DFloat:
        """Calculate RMSD with obrms.

        Does iterative row-wise mode for heavy atoms and without symmetry.

        Args:
            i_s: Conformer indices (1-indexed)
            j_s: Conformer indices (1-indexed)

        Returns:
            rmsds: Conformer RMSDs (Å)
        """
        rmsds = []
        for i in i_s:
            with tempfile.NamedTemporaryFile(
                suffix=".xyz"
            ) as ref_file, tempfile.NamedTemporaryFile(suffix=".xyz") as test_file:
                p_ref = Path(ref_file.name)
                p_ref.unlink()
                p_test = Path(test_file.name)
                p_test.unlink()
                self.write_xyz(p_ref, ids=j_s)
                self.write_xyz(p_test, ids=[i])
                process = subprocess.run(
                    f"obrms {p_ref.as_posix()} "
                    f"{p_test.as_posix()} --minimize".split(" "),
                    capture_output=True,
                )
            row_rmsds = np.genfromtxt(process.stdout.splitlines(), usecols=(-1))  # type: ignore
            rmsds.append(row_rmsds)
        rmsds: Array2DFloat = np.vstack(rmsds)

        return rmsds

    @requires_dependency([Import(module="openbabel.openbabel", alias="ob")], globals())
    def _get_rmsd_openbabel(
        self, i_s: Array1DInt, j_s: Array1DInt, include_hs: bool, symmetry: bool
    ) -> Array2DFloat:
        """Calculate RMSD row-wise with openbabel python interface.

        Args:
            i_s: Conformer indices (1-indexed)
            j_s: Conformer indices (1-indexed)
            include_hs: Whether to include H atoms in calculation
            symmetry: Whether to consider symmetry

        Returns:
            rmsds: Conformer RMSDs (Å)
        """
        rmsds = []
        for i in i_s:
            conformer_1 = self.conformers[i - 1]
            ob_mol_1 = _get_ob_mol(
                self.elements, conformer_1.coordinates, self.connectivity_matrix
            )
            align = ob.OBAlign(include_hs, symmetry)
            align.SetRefMol(ob_mol_1)
            rmsds_row = []
            for j in j_s:
                conformer_2 = self.conformers[j - 1]
                ob_mol_2 = _get_ob_mol(
                    self.elements, conformer_2.coordinates, self.connectivity_matrix
                )
                align.SetTargetMol(ob_mol_2)
                align.Align()
                rmsd = align.GetRMSD()
                rmsds_row.append(rmsd)
            rmsds.append(rmsds_row)
        rmsds: Array2DFloat = np.array(rmsds)

        return rmsds

    @requires_dependency(
        [
            Import(module="rdkit", item="Chem"),
            Import(module="rdkit.Chem", item="AllChem"),
        ],
        globals(),
    )
    def _get_rmsd_rdkit(
        self, i_s: Array1DInt, j_s: Array1DInt, include_hs: bool
    ) -> Array2DFloat:
        """Calculate RMSD row-wise with RDKit.

        Args:
            i_s: Conformer indices (1-indexed)
            j_s: Conformer indices (1-indexed)
            include_hs: Whether to include H atoms in calculation

        Returns:
            rmsds: Conformer RMSDs (Å)
        """
        # Update mol object from conformers and make copy.
        self.update_mol()

        # Construct atom list for rmsd: heavy atoms or all atoms
        if include_hs:
            atom_ids = [atom.GetIdx() for atom in self.mol.GetAtoms()]
        else:
            atom_ids = [
                atom.GetIdx()
                for atom in self.mol.GetAtoms()
                if atom.GetAtomicNum() != 1
            ]

        # Calculated RMSD row-wise with RDKit
        rmsds = []
        conformers = list(self.mol.GetConformers())
        for i in i_s:
            ref_mol = Chem.Mol(self.mol)
            ref_conformer = conformers[i - 1]
            ref_mol.RemoveAllConformers()
            ref_mol.AddConformer(ref_conformer)
            for j in j_s:
                conformer = conformers[j - 1]
                ref_mol.AddConformer(conformer)
            rmsds_row: list[float] = []
            AllChem.AlignMolConformers(ref_mol, atomIds=atom_ids, RMSlist=rmsds_row)
            rmsds.append(rmsds_row)
        rmsds: Array2DFloat = np.array(rmsds)

        return rmsds

    @requires_dependency([Import("spyrmsd"), Import("spyrmsd.rmsd")], globals())
    def _get_rmsd_spyrmsd(
        self, i_s: Array1DInt, j_s: Array1DInt, include_hs: bool, symmetry: bool
    ) -> Array2DFloat:
        """Calculate RMSD row-wise with spyrmsd.

        Args:
            i_s: Conformer indices (1-indexed)
            j_s: Conformer indices (1-indexed)
            include_hs: Whether to include H atoms in calculation
            symmetry: Whether to consider symmetry

        Returns:
            rmsds: Conformer RMSDs (Å)
        """
        # Construct mask for H atoms
        if not include_hs:
            mask = self.elements != 1
        else:
            mask = np.ones(len(self.elements), dtype=bool)

        # Calculate RMSD row-wise with spyrmsd
        rmsds: list[Array1DFloat] = []
        for i in i_s:
            ref_coordinates = self.conformers[i - 1].coordinates[mask]
            if symmetry:
                test_coordinates = [
                    self.conformers[j - 1].coordinates[mask] for j in j_s
                ]
                row_rmsds = spyrmsd.rmsd.symmrmsd(
                    ref_coordinates,
                    test_coordinates,
                    self.elements[mask],
                    self.elements[mask],
                    self.connectivity_matrix[mask, :][:, mask],
                    self.connectivity_matrix[mask, :][:, mask],
                    center=True,
                    minimize=True,
                )
            else:
                row_rmsds = []
                for j in j_s:
                    test_coordinates = self.conformers[j - 1].coordinates[mask]
                    rmsd = spyrmsd.rmsd.rmsd(
                        ref_coordinates,
                        test_coordinates,
                        self.elements[mask],
                        self.elements[mask],
                        center=True,
                        minimize=True,
                    )
                    row_rmsds.append(rmsd)
            rmsds.append(np.array(row_rmsds))
        rmsds: Array2DFloat = np.vstack(rmsds)

        return rmsds

    def update_mol(self) -> ConformerEnsemble:
        """Update Mol object with conformers."""
        self.mol.RemoveAllConformers()
        conformer_coordinates = self.get_coordinates()
        _add_conformers_to_mol(self.mol, conformer_coordinates)

        return self

    def __copy__(self) -> ConformerEnsemble:
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

    def __deepcopy__(self, memo: dict[int, object]) -> ConformerEnsemble:
        # Generate copy where conformers and mol object are new
        cls = type(self)
        ce = cls(
            self.elements,
            conformer_coordinates=self.get_coordinates(),
            energies=self.get_energies(),
            connectivity_matrix=self.connectivity_matrix,
            degeneracies=self.get_degeneracies(),
            properties=self.get_properties(),
            charge=self.charge,
            multiplicity=self.multiplicity,
            formal_charges=self.formal_charges,
            ref_cip_label=self.ref_cip_label,
        )
        ce.mol = deepcopy(ce.mol, memo)
        return ce

    def __delitem__(self, index: int) -> None:
        del self.conformers[index]

    def __getitem__(
        self, index: slice | numbers.Integral
    ) -> ConformerEnsemble | Conformer:
        cls = type(self)
        if isinstance(index, slice):
            # Generate copy of ensemble with selected conformers.
            ce = copy(self)
            ce.conformers = ce.conformers[index]
            return ce
        elif isinstance(index, numbers.Integral):
            conformer = self.conformers[int(index)]
            return conformer
        else:
            raise TypeError(f"{cls.__name__} indices must be integers")

    def __len__(self) -> int:
        return len(self.conformers)

    def __repr__(self) -> str:
        n_conformers = len(self.conformers)
        return f"{self.__class__.__name__}({n_conformers!r} conformers)"


@requires_dependency(
    [
        Import(module="openbabel.openbabel", alias="ob"),
        Import(module="openbabel.pybel", alias="pybel"),
        Import("openbabel"),
    ],
    globals(),
)
def conformers_from_ob_ff(
    mol: str | ob.OBMol,
    n_conformers: int = 30,
    ff: str = "MMFF94",
    method: str = "systematic",
    rings: bool = False,
) -> tuple[Array1DInt, Array3DFloat, Array2DInt, Array1DInt, ob.OBMol]:
    """Generates conformers based on the force field algorithm in OpenBabel.

    Follows the recipe of the command line script obabel --conformer:
    https://github.com/openbabel/openbabel/blob/master/src/ops/conformer.cpp
    If an OBMol object with 3D coordinates is given, the conformer search will
    start from that structure.

    Args:
        mol: Molecule either as SMILES string or OBMol object.
        n_conformers: Maximum number of conformers
        ff: Force field supported by OpenBabel
        method: 'fast', random', 'systematic' or 'weighted'
        rings: Sample ring torsions.

    Returns:
        elements: Elements as atomic numbers
        conformer_coordinates: Conformer coordinates (Å)
        connectivity_matrix: Connectivity matrix
        charges: Formal charges
        mol: OBMol object
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
        ff.RandomRotorSearch(n_conformers, 10, rings)
    elif method == "weighted":
        ff.WeightedRotorSearch(n_conformers, 10, rings)
    ff.GetConformers(ob_mol)

    # Extract information
    (
        elements,
        conformer_coordinates,
        connectivity_matrix,
        charges,
    ) = _extract_from_ob_mol(ob_mol)

    return (elements, conformer_coordinates, connectivity_matrix, charges, ob_mol)


@requires_dependency(
    [
        Import(module="openbabel.openbabel", alias="ob"),
        Import(module="openbabel.pybel", alias="pybel"),
        Import("openbabel"),
    ],
    globals(),
)
def conformers_from_ob_ga(  # noqa: C901
    mol: str | ob.OBMol,
    n_conformers: int | None = None,
    n_children: int | None = None,
    mutability: float | None = None,
    convergence: int | None = None,
    score: str = "rmsd",
    filter_method: str = "steric",
    cutoff: float = 0.8,
    vdw_factor: float = 0.5,
    check_hydrogens: bool = True,
) -> tuple[Array1DInt, Array3DFloat, Array2DInt, Array1DInt, ob.OBMol]:
    """Generates conformers based on the genetic algorithm in OpenBabel.

    Follows the recipe of the command line script obabel --conformer:
    https://github.com/openbabel/openbabel/blob/master/src/ops/conformer.cpp If an OBMol
    object with 3D coordinates is given, the conformer search will start from that
    structure.

    Args:
        mol: Molecule either as SMILES string or OBMol object
        n_conformers: Maximum number of conformers
        n_children: Number of children to generate for each parent
        mutability: Mutation frequency
        convergence: Number of identical generations before convergence is reached
        score: Scoring function: 'rmsd', 'min_rmsd', 'energy', 'min_energy'
        filter_method: Filtering algorithm: 'steric'
        cutoff: Absolute distance in Ånström below which atoms are considered to clash
        vdw_factor: Scale factor applied to van der Waals radii for detecting clashes
        check_hydrogens: Detect clashes with hydrogen atoms

    Returns:
        elements: Elements as atomic numbers
        conformer_coordinates: Conformer coordinates (Å)
        connectivity_matrix: Connectivity matrix
        charges: Formal charges
        mol: OBMol object
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
    if n_conformers is not None:
        conf_search.SetNumConformers(n_conformers)
    if n_children is not None:
        conf_search.SetNumChildren(n_children)
    if convergence is not None:
        conf_search.SetConvergence(convergence)
    if mutability is not None:
        conf_search.SetMutability(mutability)
    if score is not None:
        # Scorers don't work with earlier versions of openbabel
        if not (parse(openbabel.__version__) > parse("3.1.0")):
            warnings.warn(
                "Scorer only works with openbabel version > 3.1.0. "
                "Proceeding without scorer.",
                stacklevel=2,
            )
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
        if not (parse(openbabel.__version__) > parse("3.1.0")):
            warnings.warn(
                "Filter only works with openbabel version > 3.1.0. "
                "Proceeding without filter.",
                stacklevel=2,
            )
        else:
            ob_filter = ob.OBStericConformerFilter(cutoff, vdw_factor, check_hydrogens)
            conf_search.SetFilter(ob_filter)

    # Do conformational search
    conf_search.Search()
    conf_search.GetConformers(ob_mol)

    # Extract information
    (
        elements,
        conformer_coordinates,
        connectivity_matrix,
        charges,
    ) = _extract_from_ob_mol(ob_mol)

    return (elements, conformer_coordinates, connectivity_matrix, charges, ob_mol)


@requires_dependency(
    [
        Import(module="rdkit", item="Chem"),
        Import(module="rdkit.Chem", item="AllChem"),
        Import(module="rdkit.Chem", item="rdDistGeom"),
    ],
    globals(),
)
def conformers_from_rdkit(  # noqa: C901
    mol: str | Chem.Mol,
    n_conformers: int | None = None,
    optimize: str | None = "MMFF94",
    version: int = 2,
    small_rings: bool = True,
    macrocycles: bool = True,
    random_seed: int | None = None,
    rmsd_thres: float | None = 0.35,
    rmsd_symmetry: bool = False,
    n_threads: int = 1,
) -> tuple[
    Array1DInt, Array3DFloat, Array1DFloat | None, Array2DFloat, Array1DInt, Chem.Mol
]:
    """Generates conformers for an RDKit mol object.

    Recipe based on J. Chem. Inf. Modeling 2012, 52, 1146.

    Args:
        mol: Molecule either as SMILES string or RDKit Mol object.
        n_conformers: Number of conformers to generate. The actual number
            of conformers returned may be lower due to RMSD pruning. If None,
            a reasonable number will be set depending on the number of
            rotatable bonds.
        optimize: Force field used for conformer optimization: 'MMFF94', 'MMFF94s' or
            'UFF'. If None, conformers are not optimized.
        version: Version of the experimental torsion-angle preferences
        small_rings: Whether to impose small ring torsion angle preferences
        macrocycles: Whether to impose macrocycle torsion angle preferences
        random_seed: Random seed for conformer generation
        rmsd_thres: Pruning RMSD threshold (Å)
        rmsd_symmetry: Whether to use symmetry for RMSD pruning
        n_threads: Number of threads

    Returns:
        elements: Atomic symbols
        conformer_coordinates: Coordinates for all conformers (Å)
        energies: Conformer energies (a.u.)
        connectivity_matrix: Connectivity matrix with bond orders
        charges: Formal charges
        mol: RDKit Mol object. Only returned when return_mol=True

    Raises:
        Exception: When force field not found
    """
    if optimize not in ("MMFF94", "MMFF94s", "UFF", None):
        raise Exception(
            f"Force field {optimize} not found. Choose one of" "MMFF94, MMFF94s, UFF."
        )
    # Generate mol object
    try:
        mol = Chem.MolFromSmiles(mol)
    except TypeError:
        pass
    mol: Chem.Mol = Chem.AddHs(mol)

    # If n_conformers is not set, set number of conformers based on number of
    # rotatable bonds
    if not n_conformers:
        n_rot_bonds = AllChem.CalcNumRotatableBonds(mol)

        if n_rot_bonds <= 7:
            n_conformers = 50
        elif n_rot_bonds >= 8 and n_rot_bonds <= 12:
            n_conformers = 200
        else:
            n_conformers = 300

    # Generate conformers
    if rmsd_thres is not None:
        rdkit_prune_rmsd = rmsd_thres
    else:
        rdkit_prune_rmsd = -1
    if random_seed is not None:
        rdkit_random_seed = random_seed
    else:
        rdkit_random_seed = -1

    params = AllChem.EmbedParameters()
    params.useSymmetryForPruning = rmsd_symmetry
    params.randomSeed = rdkit_random_seed
    params.useSmallRingTorsions = small_rings
    params.useMacrocycleTorsions = macrocycles
    params.ETversion = version
    params.pruneRmsThresh = rdkit_prune_rmsd
    params.numThreads = n_threads
    rdDistGeom.EmbedMultipleConfs(mol, n_conformers, params)

    # Optimize with force fields
    results = None
    if optimize is not None:
        if optimize == "MMFF94":
            results = AllChem.MMFFOptimizeMoleculeConfs(
                mol, mmffVariant="MMFF94", numThreads=n_threads
            )
        elif optimize == "MMFF94s":
            results = AllChem.MMFFOptimizeMoleculeConfs(
                mol, mmffVariant="MMFF94s", numThreads=n_threads
            )
        if optimize == "UFF":
            results = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=n_threads)

    # Set energies
    if results is not None:
        energies = np.array([i[1] for i in results]) * KCAL_TO_HARTREE
    else:
        energies = None

    # Extract information from mol
    elements, conformer_coordinates, connectivity_matrix, charges = _extract_from_mol(
        mol
    )

    return (
        elements,
        conformer_coordinates,
        energies,
        connectivity_matrix,
        charges,
        mol,
    )


@requires_dependency(
    [
        Import(module="rdkit"),
        Import(module="rdkit", item="Chem"),
    ],
    globals(),
)
def _add_conformers_to_mol(mol: Chem.Mol, conformer_coordinates: ArrayLike3D) -> None:
    """Add conformers to RDKit Mol object.

    Args:
        mol: RDKit mol object
        conformer_coordinates: Conformer coordinates (Å)
    """
    conformer_coordinates: Array3DFloat = np.array(conformer_coordinates)
    if len(conformer_coordinates.shape) == 2:
        conformer_coordinates.reshape(-1, conformer_coordinates.shape[0], 3)

    for coordinates in conformer_coordinates:
        conformer = Chem.Conformer()
        for i, coord in enumerate(coordinates):
            point = rdkit.Geometry.Point3D(*coord)
            conformer.SetAtomPosition(i, point)
        mol.AddConformer(conformer, assignId=True)


@requires_dependency([Import(module="rdkit", item="Chem")], globals())
def _extract_from_mol(
    mol: Chem.Mol,
) -> tuple[list[str], Array3DFloat, Array2DFloat, Array1DInt]:
    """Extract information from RDKit Mol object with conformers."""
    # Take out elements, coordinates and connectivity matrix
    elements: list[str] = [atom.GetSymbol() for atom in mol.GetAtoms()]
    charges: Array1DInt = np.array([atom.GetFormalCharge() for atom in mol.GetAtoms()])

    conformer_coordinates = []
    for conformer in mol.GetConformers():
        coordinates = conformer.GetPositions()
        conformer_coordinates.append(coordinates)
    conformer_coordinates: Array3DFloat = np.array(conformer_coordinates)

    connectivity_matrix: Array2DFloat = Chem.GetAdjacencyMatrix(mol, useBO=True)

    return elements, conformer_coordinates, connectivity_matrix, charges


@requires_dependency(
    [
        Import(module="openbabel.openbabel", alias="ob"),
        Import(module="openbabel.pybel", alias="pybel"),
        Import("openbabel"),
    ],
    globals(),
)
def _extract_from_ob_mol(
    ob_mol: ob.OBMol,
) -> tuple[Array1DInt, Array3DFloat, Array2DInt, Array1DInt]:
    """Extract information from Openbabel OBMol object with conformers."""
    py_mol = pybel.Molecule(ob_mol)
    elements: Array1DInt = np.array([atom.atomicnum for atom in py_mol.atoms])
    charges: Array1DInt = np.array([atom.formalcharge for atom in py_mol.atoms])

    n_atoms = len(py_mol.atoms)
    connectivity_matrix: Array2DInt = np.zeros((n_atoms, n_atoms), dtype=int)
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
        coordinates: Array2DFloat = np.array([atom.coords for atom in py_mol.atoms])
        conformer_coordinates.append(coordinates)
    conformer_coordinates: Array3DFloat = np.array(conformer_coordinates)

    return elements, conformer_coordinates, connectivity_matrix, charges


@requires_dependency([Import(module="openbabel.openbabel", alias="ob")], globals())
def _get_ob_mol(
    elements: Iterable[int] | Iterable[str],
    coordinates: ArrayLike2D,
    connectivity_matrix: ArrayLike2D,
    charges: ArrayLike1D | None = None,
) -> ob.OBMol:
    """Generate OpenBabel OBMol object.

    Args:
        elements: Elements as atomic symbols or numbers.
        coordinates: Coordinates (Å)
        connectivity_matrix: Connectivity matrix with bond orders
        charges: Formal charges

    Returns:
        mol: OpenBabel OBMol object
    """
    elements = convert_elements(elements, output="numbers")

    charges_: np.ndarray
    if charges is None:
        charges_ = np.zeros(len(elements))
    else:
        charges_ = np.array(charges)
    charges = charges_

    connectivity_matrix: Array2DInt = np.array(connectivity_matrix)
    coordinates: Array2DFloat = np.array(coordinates)

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
def _get_rdkit_mol(
    elements: Iterable[int] | Iterable[str],
    conformer_coordinates: ArrayLike3D,
    connectivity_matrix: ArrayLike2D,
    charges: ArrayLike1D | None = None,
) -> Chem.Mol:
    _RDKIT_BOND_TYPES = {
        1.0: Chem.BondType.SINGLE,
        1.5: Chem.BondType.AROMATIC,
        2.0: Chem.BondType.DOUBLE,
        3.0: Chem.BondType.TRIPLE,
        4.0: Chem.BondType.QUADRUPLE,
        5.0: Chem.BondType.QUINTUPLE,
        6.0: Chem.BondType.HEXTUPLE,
    }
    elements = convert_elements(elements, output="symbols")

    charges_: np.ndarray
    if charges is None:
        charges_ = np.zeros(len(elements))
    else:
        charges_ = np.array(charges)
    charges = charges_
    conformer_coordinates: Array3DFloat = np.array(conformer_coordinates)
    connectivity_matrix: Array2DInt = np.array(connectivity_matrix)

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
    _add_conformers_to_mol(mol, conformer_coordinates)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)

    return mol


def cli(smiles: str, generator: str = "rdkit") -> Any:
    """CLI for cone angle.

    Args:
        smiles: SMILES string
        generator: Confomer generator: 'ob-ff', 'ob-ga' or 'rdkit'

    Returns:
        Partially instantiated class
    """
    if generator == "rdkit":
        return functools.partial(ConformerEnsemble.from_rdkit, smiles)
    elif generator == "ob-ff":
        return functools.partial(ConformerEnsemble.from_ob_ff, smiles)
    elif generator == "ob-ga":
        return functools.partial(ConformerEnsemble.from_ob_ga, smiles)
