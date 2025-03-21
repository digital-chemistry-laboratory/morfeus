"""Local force constant code."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
import functools
from os import PathLike
from typing import Any

import numpy as np

from morfeus.data import (
    AFU,
    AMU,
    ANGSTROM,
    ANGSTROM_TO_BOHR,
    atomic_masses,
    BOHR,
    BOHR_TO_ANGSTROM,
    C,
    DYNE,
    HARTREE,
)
from morfeus.geometry import (
    Angle,
    Bond,
    Dihedral,
    InternalCoordinates,
    kabsch_rotation_matrix,
)
from morfeus.io import read_geometry
from morfeus.typing import (
    Array1DFloat,
    Array1DInt,
    Array2DFloat,
    ArrayLike1D,
    ArrayLike2D,
)
from morfeus.utils import convert_elements


class LocalForce:
    """Calculates and stores the results from local force constant calculations.

    The method is described by Cremer in
    10.1002/(SICI)1097-461X(1998)67:1<1::AID-QUA1>3.0.CO;2-Z. Alternatively, the
    compliance matrix method can be used according to 10.1063/1.3413528

    Args:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)

    Attributes:
        internal_coordinates: Internal coordinates
        local_force_constants: Local force constants (mDyne/Å)
        local_frequencies: Local mode frequencies (cm⁻¹)
        n_imag: Number of normal modes with imaginary frequencies
    """

    internal_coordinates: list[Bond | Angle | Dihedral]
    local_force_constants: Array1DFloat
    local_frequencies: Array1DFloat
    n_imag: int
    _B_inv: Array2DFloat
    _B: Array2DFloat
    _coordinates: Array2DFloat
    _D_full: Array2DFloat
    _D: Array2DFloat
    _elements: list[int]
    _fc_matrix: Array2DFloat
    _force_constants: Array1DFloat
    _ifc_matrix: Array2DFloat
    _input_coordinates: Array2DFloat
    _internal_coordinates: InternalCoordinates
    _masses: Array1DFloat
    _normal_modes: Array2DFloat
    _standard_coordinates: Array2DFloat

    def __init__(
        self,
        elements: Iterable[int] | Iterable[str] | None = None,
        coordinates: ArrayLike2D | None = None,
    ) -> None:
        # Set up attributes
        self._internal_coordinates = InternalCoordinates()
        self.internal_coordinates = self._internal_coordinates.internal_coordinates

        if elements is not None:
            elements = convert_elements(elements, output="numbers")
            self._elements = elements
            self._masses: Array1DFloat = np.array([atomic_masses[i] for i in elements])
        else:
            self._elements = []
            self._masses: Array1DFloat = np.empty(0, dtype=float)

        if coordinates is not None:
            self._coordinates: Array2DFloat = np.array(coordinates)
        else:
            self._coordinates: Array2DFloat = np.empty(0, dtype=float)

    def add_internal_coordinate(self, atoms: Sequence[int]) -> "LocalForce":
        """Add internal coordinate.

        Composed of two (bond), three (angle) or four atoms (dihedral).

        Args:
            atoms: Atom indices of internal coordinate

        Returns:
            self: Self
        """
        self._internal_coordinates.add_internal_coordinate(atoms)

        return self

    def compute_compliance(self) -> "LocalForce":
        """Compute local force constants with the compliance matrix method."""
        # Compute B matrix if it does not exists
        if not hasattr(self, "_B"):
            self._B = self._internal_coordinates.get_B_matrix(self._coordinates)

        # Compute compliance matrix and get local force constants
        C = self._B @ np.linalg.pinv(self._fc_matrix) @ self._B.T
        k_s = 1 / np.diag(C) * AFU / BOHR / (DYNE / 1000) * ANGSTROM

        self.local_force_constants = k_s

        return self

    def compute_frequencies(self) -> "LocalForce":
        """Compute local frequencies."""
        # Compute local frequencies
        M: Array1DFloat = np.diag(np.repeat(self._masses, 3))
        G = self._B @ np.linalg.inv(M) @ self._B.T
        frequencies = (
            np.sqrt(
                (
                    np.diag(G)
                    / AMU
                    * (self.local_force_constants / ANGSTROM * DYNE / 1000)
                )
                / (4 * np.pi**2 * C**2)
            )
            / 100
        )

        self.local_frequencies = frequencies

        return self

    def compute_local(  # noqa: C901
        self, project_imag: bool = True, cutoff: float = 1e-3
    ) -> "LocalForce":
        """Compute local force constants with the local modes approach.

        Args:
            project_imag: Whether to project out imaginary frequencies
            cutoff: Cutoff for low force constant (mDyne/Å)

        Returns:
            self: Self
        """
        # Compute D matrix from normal modes and B matrix
        if not hasattr(self, "_D"):
            if not hasattr(self, "_B"):
                self._B = self._internal_coordinates.get_B_matrix(self._coordinates)
            self._D = self._B @ self._normal_modes.T

        # Project out imaginary modes
        if project_imag is True and self.n_imag > 0:
            # Save full D matrix and force constant vector
            self._D_full = self._D

            # Remove force constants with imaginary force constants
            self._force_constants = self._force_constants[self.n_imag :]
            self._D = self._D[:, self.n_imag :]

        # Add cutoff to modes with low force constants
        if cutoff is not None:
            indices_small = np.where(np.array(self._force_constants) < cutoff)[0]
            if len(indices_small) > 0:
                self._force_constants[indices_small] = np.array(cutoff).reshape(-1)

        # Compute local mode force constants
        K: Array1DFloat = np.diag(self._force_constants)

        k_s = []
        for row in self._D:
            if not np.all(np.isclose(row, 0)):
                k = 1 / (row @ np.linalg.inv(K) @ np.conj(row).T)
                k_s.append(k)
            else:
                k_s.append(0.0)
        k_s: Array1DFloat = np.array(k_s)

        # Scale force constants due to projection of imaginary normal modes
        if project_imag and self.n_imag > 0:
            lengths = np.linalg.norm(self._D, axis=1)
            lengths_full = np.linalg.norm(self._D_full, axis=1)
            k_s *= lengths / lengths_full

        self.local_force_constants = k_s

        return self

    def detect_bonds(
        self,
        radii: ArrayLike1D | None = None,
        radii_type: str = "pyykko",
        scale_factor: float = 1.2,
    ) -> "LocalForce":
        """Detect bonds based on scaled sum of covalent radii.

        Args:
            radii: Covalent radii (Å)
            radii_type: Covalent radii type: 'pyykko'
            scale_factor: Scale factor for covalent radii

        Returns:
            self: Self
        """
        self._internal_coordinates.detect_bonds(
            self._coordinates,
            self._elements,
            radii=radii,
            radii_type=radii_type,
            scale_factor=scale_factor,
        )

        return self

    def get_local_force_constant(self, atoms: Sequence[int]) -> float:
        """Return the local force constant between a set of atoms.

        Args:
            atoms: Atom indices of the internal coordinate

        Returns:
            force_constant: Local force constant (mDyne/Å, mDyne Å/rad²)

        Raises:
            ValueError: When no internal coordinate found
        """
        coordinate = _get_internal_coordinate(atoms)
        index = self.internal_coordinates.index(coordinate)
        if index is None:
            raise ValueError(f"No internal coordinate with these atoms: {atoms}")
        force_constant: float = self.local_force_constants[index]

        return force_constant

    def get_local_frequency(self, atoms: Sequence[int]) -> float:
        """Return the local frequency between a set of atoms.

        Args:
            atoms: Atom indices in the internal coordinate

        Returns:
            frequency: Local frequency (cm⁻¹)

        Raises:
            ValueError: When no internal coordinate found
        """
        coordinate = _get_internal_coordinate(atoms)
        index = self.internal_coordinates.index(coordinate)
        if index is None:
            raise ValueError(f"No internal coordinate with these atoms: {atoms}.")
        frequency: float = self.local_frequencies[index]

        return frequency

    def load_file(
        self, file: str | PathLike, program: str, filetype: str
    ) -> "LocalForce":
        """Load data from external file.

        Args:
            file: File
            program: Program used to generate file: 'gaussian', 'unimovib' or 'xtb'
            filetype: Filetype. For 'gaussian': 'fchk' or 'log'. For 'unimovib': 'local'
                or 'log'. For 'xtb': 'hessian' or 'normal_modes'

        Returns:
            self: Self
        """
        choices = {
            "gaussian": {
                "fchk": self._parse_gaussian_fchk,
                "log": self._parse_gaussian_log,
            },
            "xtb": {
                "hessian": self._parse_xtb_hessian,
            },
            "unimovib": {
                "local": self._parse_unimovib_local,
                "log": self._parse_unimovib_log,
                "umv": self._parse_unimovib_umv,
            },
        }
        choices[program][filetype](file)

        return self

    def normal_mode_analysis(
        self,
        hessian: ArrayLike2D | None = None,
        save_hessian: bool = False,
    ) -> "LocalForce":
        """Perform normal mode analysis.

        With projection of translations and vibrations to get normal modes and force
        constants.

        Args:
            hessian: User-supplied Hessian
            save_hessian: Save projected Hessian for use in compliance matrix method.

        Returns:
            self: Self
        """
        # Set up
        coordinates: Array2DFloat = self._coordinates * ANGSTROM_TO_BOHR
        masses = self._masses
        if hessian is None:
            hessian = self._fc_matrix
        else:
            hessian = np.array(hessian)
        n_atoms = len(coordinates)

        # Create mass matrices
        M_minus: Array1DFloat = np.diag(np.repeat(masses, 3) ** (-1 / 2))
        M_plus: Array1DFloat = np.diag(np.repeat(masses, 3) ** (1 / 2))
        m_plus = np.repeat(masses, 3) ** (1 / 2)
        m_minus = np.repeat(masses, 3) ** (-1 / 2)

        # Mass-weight Hessian
        hessian_mw = M_minus @ hessian @ M_minus

        # Find center of mass
        com: Array1DFloat = np.sum(
            masses.reshape(-1, 1) * coordinates, axis=0
        ) / np.sum(masses)

        # Shift origin to center of mass
        coordinates -= com

        # Construct translation vectors
        t_x: Array2DFloat = (np.tile(np.array([1, 0, 0]), n_atoms)).reshape(-1, 3)
        t_y: Array2DFloat = (np.tile(np.array([0, 1, 0]), n_atoms)).reshape(-1, 3)
        t_z: Array2DFloat = (np.tile(np.array([0, 0, 1]), n_atoms)).reshape(-1, 3)

        # Construct mass-weighted rotation vectors
        R_x: Array2DFloat = np.cross(coordinates, t_x).flatten() * m_plus
        R_y: Array2DFloat = np.cross(coordinates, t_y).flatten() * m_plus
        R_z: Array2DFloat = np.cross(coordinates, t_z).flatten() * m_plus

        # Mass-weight translation vectors
        T_x = t_x.flatten() * m_plus
        T_y = t_y.flatten() * m_plus
        T_z = t_z.flatten() * m_plus

        # Remove linear dependencies from translation/rotation space
        TR_vectors = np.vstack([T_x, T_y, T_z, R_x, R_y, R_z])
        Q, R = np.linalg.qr(TR_vectors.T)
        keep_indices = ~np.isclose(np.diag(R), 0, atol=1e-6, rtol=0)
        TR_vectors = Q.T[keep_indices]
        n_tr = len(TR_vectors)

        # Construct P matrix
        P = np.identity(n_atoms * 3)
        for vector in TR_vectors:
            P -= np.outer(vector, vector)

        # Project out translations and rotations
        hessian_proj = P.T @ hessian_mw @ P

        # Diagonalize
        eigenvalues, eigenvectors = np.linalg.eigh(hessian_proj)
        eigenvalues = eigenvalues[n_tr:]
        eigenvectors = eigenvectors[:, n_tr:]

        # Calculate cartesian displacements
        cart = eigenvectors.T * m_minus
        N = 1 / np.linalg.norm(cart, axis=1)
        norm_cart = cart * N.reshape(-1, 1)
        reduced_masses = N**2

        # Calculate frequencies and force constants
        n_imag = int(np.sum(eigenvalues < 0))
        frequencies = (
            np.sqrt(np.abs(eigenvalues) * HARTREE / BOHR**2 / AMU)
            / (2 * np.pi * C)
            / 100
        )
        frequencies[:n_imag] = -frequencies[:n_imag]
        force_constants = (
            4
            * np.pi**2
            * (frequencies * 100) ** 2
            * C**2
            * reduced_masses
            * AMU
            / (DYNE / 1000)
            * ANGSTROM
        )

        # Set up attributes
        self.n_imag = n_imag
        self._force_constants = force_constants
        self._normal_modes = norm_cart
        if save_hessian:
            self._fc_matrix = M_plus @ hessian_proj @ M_plus

        return self

    def print_report(
        self, angles: bool = False, dihedrals: bool = False, angle_units: bool = False
    ) -> None:
        """Print report of results.

        Args:
            angle_units: Wheter to convert angle and dihedral force constants to
                mDyne Å/rad²
            angles: Whether to print angles
            dihedrals: Whether to print dihedrals
        """
        # Print header
        if angle_units:
            unit = "mDyne/Å, mDyne Å/rad²"
        else:
            unit = "mDyne/Å"

        string = f"{'Coordinate':30s}" + f"{'Force constant ' + '(' + unit + ')':>50s}"
        if hasattr(self, "local_frequencies"):
            string += f"{'Frequency (cm⁻¹)':>30s}"
        print(string)

        # Print results for each internal
        sorted_coordinates = sorted(
            self.internal_coordinates, key=lambda x: (len(x.atoms), *x.atoms)
        )
        for coordinate in sorted_coordinates:
            # Check if internal coordinate is angle or dihedral
            if len(coordinate.atoms) == 3 and not angles:
                continue
            if len(coordinate.atoms) == 4 and not dihedrals:
                continue
            index = self.internal_coordinates.index(coordinate)
            force_constant = self.local_force_constants[index]

            # Convert units for angles and dihedrals
            if len(coordinate.atoms) > 2 and angle_units:
                force_constant = force_constant * BOHR_TO_ANGSTROM**2

            # Print out the results
            string = f"{repr(coordinate):30s}" + f"{force_constant:50.3f}"
            if hasattr(self, "local_frequencies"):
                frequency = self.local_frequencies[index]
                string += f"{frequency:30.0f}"
            print(string)

    def reset_internal_coordinates(self) -> "LocalForce":
        """Reset internal coordinate system."""
        self._internal_coordinates = InternalCoordinates()
        self.internal_coordinates = self._internal_coordinates.internal_coordinates

        return self

    def _parse_gaussian_fchk(self, file: str | PathLike) -> None:  # noqa: C901
        # Read fchk file
        with open(file) as f:
            lines = f.readlines()

        # Set up read flags
        read_modes = False
        read_hessian = False
        read_vib_e2 = False
        read_ic = False
        read_atomic_numbers = False
        read_coordinates = False
        read_masses = False

        # Set up containers for reading data
        modes: list[float] = []
        hessian: list[float] = []
        vib_e2: list[float] = []
        internal_coordinates = []
        n_atoms: int
        n_imag: int
        atomic_numbers: list[int] = []
        masses: list[float] = []
        coordinates: list[float] = []

        # Parse fchk file
        for line in lines:
            # Read normal modes
            if read_modes:
                try:
                    split_line = line.strip().split()
                    modes.extend([float(value) for value in split_line])
                except ValueError:
                    read_modes = False
            # Read cartesian force constants
            elif read_hessian:
                try:
                    split_line = line.strip().split()
                    hessian.extend([float(value) for value in split_line])
                except ValueError:
                    read_hessian = False
            # Read normal mode force constants
            elif read_vib_e2:
                try:
                    split_line = line.strip().split()
                    vib_e2.extend([float(value) for value in split_line])
                except ValueError:
                    read_vib_e2 = False
            # Read internal coordinates
            elif read_ic:
                try:
                    split_line = line.strip().split()
                    internal_coordinates.extend([int(value) for value in split_line])
                except ValueError:
                    read_ic = False
            # Read atomic numbers
            elif read_atomic_numbers:
                try:
                    split_line = line.strip().split()
                    atomic_numbers.extend([int(value) for value in split_line])
                except ValueError:
                    read_atomic_numbers = False
            # Read atomic masses
            elif read_masses:
                try:
                    split_line = line.strip().split()
                    masses.extend([float(value) for value in split_line])
                except ValueError:
                    read_masses = False
            # Read coordinates
            elif read_coordinates:
                try:
                    split_line = line.strip().split()
                    coordinates.extend([float(value) for value in split_line])
                except ValueError:
                    read_coordinates = False
            # Read number of atoms
            if "Number of atoms" in line:
                n_atoms = int(line.strip().split()[4])
            # Read number of normal modes
            elif "Number of Normal Modes" in line:
                n_modes = int(line.strip().split()[5])
            # Read number of internal coordinates
            elif "Redundant internal coordinates" in line:
                n_redundant = int(line.strip().split()[5])
            elif "NImag" in line:
                n_imag = int(line.strip().split()[2])
            # Detect when to read data
            elif "Vib-Modes" in line:
                read_modes = True
            elif "Vib-AtMass" in line:
                read_masses = True
            elif "Cartesian Force Constants" in line:
                read_hessian = True
            elif "Vib-E2 " in line:
                read_vib_e2 = True
            elif "Redundant internal coordinate indices" in line:
                read_ic = True
            elif "Atomic numbers" in line:
                read_atomic_numbers = True
            elif "Current cartesian coordinates " in line:
                read_coordinates = True
        # Take out normal mode force constants
        force_constants: Array1DFloat = np.array(vib_e2[n_modes * 2 : n_modes * 3])

        # Construct force constant matrix from lower triangular matrix
        fc_matrix: Array2DFloat = np.zeros((n_atoms * 3, n_atoms * 3), dtype=float)
        fc_matrix[np.tril_indices_from(fc_matrix)] = hessian
        fc_matrix: Array2DFloat = np.triu(fc_matrix.T, 1) + fc_matrix

        # Take out the internal coordinates
        internal_coordinates: Array1DInt = np.array(internal_coordinates)
        coordinate: Array1DInt
        for _, coordinate in enumerate(np.split(internal_coordinates, n_redundant)):
            if all(coordinate >= 0):  # Sort out linear bends
                atoms = [i for i in coordinate if i != 0]
                self.add_internal_coordinate(atoms)

        # Convert coordinates to right Ångström
        coordinates: Array2DFloat = (
            np.array(coordinates).reshape(-1, 3) * BOHR_TO_ANGSTROM
        )

        # Set up attributes
        self._fc_matrix = fc_matrix
        self._normal_modes = np.array(modes).reshape(n_modes, n_atoms * 3)
        self._force_constants = force_constants
        self.n_imag = n_imag
        self._elements = convert_elements(atomic_numbers, output="numbers")
        self._coordinates = coordinates
        self._masses = np.array(masses)

    def _parse_gaussian_log(self, file: str | PathLike) -> None:  # noqa: C901
        # Read the log file
        with open(file) as f:
            lines = f.readlines()

        # Set up read flags
        read_b_atoms = False
        read_b_vectors = False
        read_fc_matrix = False
        read_ifc_matrix = False
        read_input_orientation = False
        read_standard_orientation = False
        read_internal = False
        read_internal_modes = False
        read_hp_modes = False
        read_masses = False

        # Set up containers for reading data
        B_atom_map: dict[int, list[int]] = {}
        B_vectors: dict[int, list[float]] = {}
        normal_modes: list[list[list[float]]] = []
        internal_modes: list[list[float]] = []
        force_constants: list[float] = []
        masses: list[float] = []
        fc_matrix: Array2DFloat = np.empty([])
        ifc_matrix: Array2DFloat = np.empty([])
        input_coordinates: list[float] = []
        standard_coordinates: list[float] = []
        n_imag: int = 0
        n_atoms: int = 0
        internal_indices: dict[frozenset[int], int] = {}
        atomic_numbers: list[int] = []
        coordinates: list[list[float]] = []

        # Parse through log file content
        counter = 0
        internal_names: list[str] = []
        internal_vector: list[float] = []

        values: Any
        value: Any
        for line in lines:
            # Read internal coordinate definitions
            if read_internal:
                if counter > 1:
                    if (
                        " --------------------------------------------------------------------------------"  # noqa: B950
                        in line
                    ):
                        read_internal = False
                        n_internals = len(internal_indices.items())
                    else:
                        split_line = line.strip().split()
                        name = split_line[1]
                        internal_names.append(name)
                        atoms = split_line[2][2:].replace(")", "").split(",")
                        atoms = [int(atom) for atom in atoms]
                        internal_indices[frozenset(atoms)] = counter - 2
                counter += 1
            # Read Cartesian force constant matrix
            elif read_fc_matrix:
                if " FormGI is forming" in line:
                    read_fc_matrix = False
                    fc_matrix = np.triu(fc_matrix.T, 1) + fc_matrix
                elif "          " in line:
                    column_indices = [int(value) for value in line.strip().split()]
                elif "D" in line:
                    split_line = line.strip().split()
                    row_index = int(split_line[0]) - 1
                    values = [
                        float(value.replace("D", "E")) for value in split_line[1:]
                    ]
                    for i, value in enumerate(values):
                        column_index = column_indices[i] - 1
                        fc_matrix[row_index, column_index] = value
            # Read internal force constant matrix
            elif read_ifc_matrix:
                if "Leave Link  716" in line:
                    read_ifc_matrix = False
                    ifc_matrix = np.triu(ifc_matrix.T, 1) + ifc_matrix
                elif "          " in line:
                    column_indices = [int(value) for value in line.strip().split()]
                elif "D" in line:
                    split_line = line.strip().split()
                    row_index = int(split_line[0]) - 1
                    values = [
                        float(value.replace("D", "E")) for value in split_line[1:]
                    ]
                    for i, value in enumerate(values):
                        column_index = column_indices[i] - 1
                        ifc_matrix[row_index, column_index] = value
            # Read atoms for creation of B matrix
            elif read_b_atoms:
                if " B Matrix in FormBX:" in line:
                    read_b_atoms = False
                else:
                    if counter == 0:
                        atoms = [int(value) for value in line.strip().split()]
                        for atom in atoms:
                            B_atom_map[atom] = []
                    if counter > 0 and counter < 5:
                        values = [int(value) for value in line.strip().split()[1:]]
                        for atom, value in zip(atoms, values):
                            B_atom_map[atom].append(value)
                    counter += 1
                    if counter == 5:
                        counter = 0
            # Read values of B matrix
            elif read_b_vectors:
                if (
                    " IB Matrix in Red2BG:" in line
                    or "Iteration" in line
                    or " G Matrix:" in line
                ):
                    read_b_vectors = False
                else:
                    if counter == 0:
                        atoms = [int(value) for value in line.strip().split()]
                        for atom in atoms:
                            B_vectors[atom] = []
                    if counter > 0 and counter < 13:
                        values = [
                            float(value.replace("D", "E"))
                            for value in line.strip().split()[1:]
                        ]
                        for atom, value in zip(atoms, values):
                            B_vectors[atom].append(value)
                    counter += 1
                    if counter == 13:
                        counter = 0
            # Read atomic coordinates in input orientation
            elif read_input_orientation:
                if counter > 3:
                    if (
                        "---------------------------------------------------------------------"  # noqa: B950
                        in line
                    ):
                        read_input_orientation = False
                    else:
                        strip_line = line.strip().split()
                        values = [float(value) for value in strip_line[3:]]
                        input_coordinates.append(values)
                        atomic_numbers.append(int(strip_line[1]))
                counter += 1
            # Read atomic coordinates in standard orientation
            elif read_standard_orientation:
                if counter > 3:
                    if (
                        "---------------------------------------------------------------------"  # noqa: B950
                        in line
                    ):
                        read_standard_orientation = False
                    else:
                        strip_line = line.strip().split()
                        values = [float(value) for value in strip_line[3:]]
                        standard_coordinates.append(values)
                counter += 1
            # Read decomposition of normal modes in internal coordinates
            elif read_internal_modes:
                if counter > 3:
                    if (
                        "--------------------------------------------------------------------------------"  # noqa: B950
                        in line
                    ):
                        read_internal_modes = False
                        internal_modes.append(internal_vector)
                    else:
                        value = float(line.strip().split()[3])
                        internal_vector.append(value)
                counter += 1
            # Read high-precision normal modes
            elif read_hp_modes:
                if counter < n_atoms * 3:
                    strip_line = line.strip().split()
                    values = [float(value) for value in strip_line[3:]]
                    coordinates.append(values)
                if counter == n_atoms * 3:
                    normal_modes.append(coordinates)
                    read_hp_modes = False
                counter += 1
            # Read atomic masses
            elif read_masses:
                if "Molecular mass: " in line:
                    read_masses = False
                elif "and mass" in line:
                    masses.append(float(line.strip().split()[8]))
            # Read number of atoms
            if " NAtoms=" in line and not n_atoms:
                n_atoms = int(line.strip().split()[1])
            # Read normal mode force constants
            elif "Frc consts" in line:
                split_line = line.strip().split()
                values = [float(value) for value in split_line[3:]]
                force_constants.extend(values)
            # Read number of imaginary frequencies
            elif "imaginary frequencies (negative Signs)" in line:
                n_imag = int(line.strip().split()[1])
            # Detect when to read data
            elif (
                "Name  Definition              Value          Derivative Info." in line
            ):
                read_internal = True
                counter = 1
                internal_names = []
            elif "- Thermochemistry -" in line:
                read_masses = True
            elif " IB Matrix in FormBX:" in line:
                read_b_atoms = True
                counter = 0
            elif " B Matrix in FormBX:" in line:
                read_b_vectors = True
                counter = 0
            elif " Force constants in Cartesian coordinates: " in line:
                read_fc_matrix = True
                fc_matrix = np.zeros((3 * n_atoms, 3 * n_atoms))
                counter = 0
            elif " Force constants in internal coordinates: " in line:
                read_ifc_matrix = True
                ifc_matrix = np.zeros((n_internals, n_internals))
                counter = 0
            elif "Input orientation: " in line:
                read_input_orientation = True
                counter = 0
            elif "Standard orientation: " in line:
                read_standard_orientation = True
                counter = 0
            elif "Normal Mode" in line:
                read_internal_modes = True
                counter = 1
                internal_vector = []
            elif " Coord Atom Element:" in line:
                read_hp_modes = True
                coordinates = []
                counter = 0

        # Process internal coordinates
        if len(internal_indices) > 0:
            for name, indices in zip(internal_names, internal_indices):
                if name[0] == "R" and len(indices) == 2:
                    self._internal_coordinates.add_internal_coordinate(list(indices))
                if name[0] == "A" and len(indices) == 3:
                    self._internal_coordinates.add_internal_coordinate(list(indices))
                if name[0] == "D" and len(indices) == 4:
                    self._internal_coordinates.add_internal_coordinate(list(indices))

        # Construct the B matrix from atoms and vectors
        B: np.ndarray
        if B_vectors:
            n_cartesian = n_atoms * 3
            B = np.zeros((n_internals, n_cartesian))
            for i in range(n_internals):
                for j, atom in enumerate(B_atom_map[i + 1]):
                    if atom:
                        B[i][(atom - 1) * 3 : (atom - 1) * 3 + 3] = B_vectors[i + 1][
                            j * 3 : j * 3 + 3
                        ]
            B_inv = np.linalg.pinv(B)
            # Detect whether the internal coordinate system is redundant
            if B.shape[0] == len(force_constants):
                self._redundant = False
            else:
                self._redundant = True
        else:
            B = np.array([])
            B_inv = np.array([])

        # Detect whether the input coordinates have been rotated. If so, rotate
        # B matrix and its inverse.
        input_coordinates: Array2DFloat = np.array(input_coordinates).reshape(-1, 3)
        standard_coordinates: Array2DFloat = np.array(standard_coordinates).reshape(
            -1, 3
        )

        if (
            not np.array_equal(input_coordinates, standard_coordinates)
            and standard_coordinates.size > 0
        ):
            rotation_i_to_s = kabsch_rotation_matrix(
                input_coordinates, standard_coordinates
            )
            if len(B) > 0:
                B = (rotation_i_to_s @ B.reshape(-1, 3).T).T.reshape(n_internals, -1)
                B_inv = np.linalg.pinv(B)

        # Set up attributes
        self._B = B
        self._B_inv = B_inv
        self._fc_matrix = fc_matrix
        if fc_matrix.size == 0 and ifc_matrix.size > 0:
            self._fc_matrix = B.T @ ifc_matrix @ B
        self._ifc_matrix = ifc_matrix
        self._force_constants = np.array(force_constants)
        if len(normal_modes) > 0:
            self._normal_modes = np.hstack(normal_modes).T
        else:
            self._normal_modes = np.array([])
        if len(internal_modes) > 0:
            self._D = np.vstack(internal_modes).T
        self.n_imag = n_imag
        self._masses = np.array(masses)
        self._input_coordinates = input_coordinates
        self._standard_coordinates = standard_coordinates
        self._coordinates = input_coordinates
        self._elements = convert_elements(atomic_numbers, output="numbers")

    def _parse_unimovib_local(self, file: str | PathLike) -> None:  # noqa: C901
        # Read file
        with open(file) as f:
            lines = f.readlines()

        # Set up flags for reading
        read_masses = False
        read_atomic_numbers = False
        read_coordinates = False
        read_hessian = False
        read_normal_modes = False

        # Set up containers for data
        masses = []
        atomic_numbers = []
        coordinates = []
        hessian = []
        normal_modes = []

        # Parse file
        for line in lines:
            # Read atomic masses
            if read_masses:
                if " $ZA  $END" in line:
                    read_masses = False
                else:
                    split_line = line.strip().split()
                    masses.extend(
                        [float(value.replace("D", "E")) for value in split_line]
                    )
            # Read atomic numbers
            if read_atomic_numbers:
                if " $XYZ  $END" in line:
                    read_atomic_numbers = False
                else:
                    split_line = line.strip().split()
                    atomic_numbers.extend(
                        [int(float(value.replace("D", "E"))) for value in split_line]
                    )
            # Read coordinates
            if read_coordinates:
                if " $FFX  $END" in line:
                    read_coordinates = False
                else:
                    split_line = line.strip().split()
                    coordinates.extend(
                        [float(value.replace("D", "E")) for value in split_line]
                    )
            # Read Hessian
            if read_hessian:
                if " $NMMODE  $END" in line:
                    read_hessian = False
                else:
                    split_line = line.strip().split()
                    hessian.extend(
                        [float(value.replace("D", "E")) for value in split_line]
                    )
            # Read normal modes
            if read_normal_modes:
                if " $APT" in line:
                    read_normal_modes = False
                else:
                    split_line = line.strip().split()
                    normal_modes.extend(
                        [float(value.replace("D", "E")) for value in split_line]
                    )
            # Read number of atoms and normal modes
            if " $CONTRL" in line:
                split_line = line.strip().split()
                n_atoms = int(split_line[1].split("=")[1])
                n_modes = int(split_line[2].split("=")[1])
            # Detect when to read data
            if " $AMASS  $END" in line:
                read_masses = True
            if " $ZA  $END" in line:
                read_atomic_numbers = True
            if " $XYZ  $END" in line:
                read_coordinates = True
            if " $FFX  $END" in line:
                read_hessian = True
            if " $NMMODE  $END" in line:
                read_normal_modes = True

        # Convert data to right format
        coordinates = np.array(coordinates).reshape(-1, 3) * BOHR_TO_ANGSTROM
        hessian: Array2DFloat = np.array(hessian).reshape(n_atoms * 3, n_atoms * 3)
        normal_modes = np.array(normal_modes).reshape(-1, n_atoms * 3)[:n_modes]
        elements = convert_elements(atomic_numbers, output="numbers")

        # Set up attributes
        self._normal_modes = normal_modes
        self._masses = np.array(masses)
        self._fc_matrix = hessian
        self._coordinates = coordinates
        self._elements = elements

    def _parse_unimovib_log(self, file: str | PathLike) -> None:  # noqa: C901
        # Read file
        with open(file) as f:
            lines = f.readlines()

        # Set up read flags
        read_coordinates = False
        read_vibrations = False
        read_normal_modes = False

        # Set up data containers
        atomic_numbers = []
        masses = []
        coordinates = []
        force_constants = []
        normal_modes = []
        frequencies = []

        # Parse file
        normal_modes_chunk: list[list[float]] = []
        counter = 0
        n_modes_chunk = 0
        values: Any
        for line in lines:
            # Read coordinates, atomic numbers and atomic masses
            if read_coordinates:
                if counter > 0:
                    if (
                        " ------------------------------------------------------------------------------------------"  # noqa: B950
                        in line
                    ):
                        read_coordinates = False
                    else:
                        split_line = line.strip().split()
                        atomic_numbers.append(int(split_line[2]))
                        coordinates.extend([float(value) for value in split_line[3:6]])
                        masses.append(float(split_line[6]))
                counter += 1
            # Read normal modes and force constants
            if read_vibrations:
                if " Results of translations and rotations:" in line:
                    read_vibrations = False
                if read_normal_modes:
                    split_line = line.strip().split()
                    if len(split_line) == 0:
                        read_normal_modes = False
                        normal_modes.extend(normal_modes_chunk)
                    else:
                        values = [float(value) for value in split_line[2:]]
                        for i in range(n_modes_chunk):
                            normal_modes_chunk[i].append(values[i * 3 : i * 3 + 3])
                elif "Irreps" in line:
                    split_line = line.strip().split()
                    n_modes_chunk = len(split_line) - 1
                    normal_modes_chunk = [[] for i in range(n_modes_chunk)]
                elif "Force constants" in line:
                    split_line = line.strip().split()
                    force_constants.extend([float(value) for value in split_line[2:]])
                elif "Frequencies" in line:
                    split_line = line.strip().split()
                    frequencies.extend([float(value) for value in split_line[2:]])
                elif "        Atom  ZA" in line:
                    read_normal_modes = True
            # Detect when to read data
            if "No.   Atom    ZA" in line:
                read_coordinates = True
                counter = 0
            if "Results of vibrations:" in line:
                read_vibrations = True

        # Convert quantities to right format
        n_atoms = len(masses)
        coordinates: Array2DFloat = np.array(coordinates).reshape(-1, 3)
        normal_modes: Array2DFloat = np.array(normal_modes).reshape(-1, n_atoms * 3)
        elements = convert_elements(atomic_numbers, output="numbers")
        force_constants: Array1DFloat = np.array(force_constants)
        frequencies: Array1DFloat = np.array(frequencies)
        masses: Array1DFloat = np.array(masses)

        # Detect imaginary modes
        self.n_imag = np.sum(frequencies < 0)

        # Set up attributes
        self._elements = elements
        self._coordinates = coordinates
        self._masses = masses
        self._force_constants = force_constants
        self._normal_modes = normal_modes

    def _parse_unimovib_umv(self, file: str | PathLike) -> None:  # noqa: C901
        # Read file
        with open(file) as f:
            lines = f.readlines()

        # Set up flags for reading
        read_n_atoms = False
        read_masses = False
        read_atomic_numbers = False
        read_coordinates = False
        read_hessian = False

        # Set up containers for data
        n_atoms: int
        masses = []
        atomic_numbers = []
        coordinates = []
        hessian = []

        # Parse file
        for line in lines:
            # Read number of atoms
            if read_n_atoms:
                if "AMASS" in line:
                    read_n_atoms = False
                else:
                    n_atoms = int(line.strip())
            # Read atomic masses
            if read_masses:
                if "ZA" in line:
                    read_masses = False
                else:
                    split_line = line.strip().split()
                    masses.extend(
                        [float(value.replace("D", "E")) for value in split_line]
                    )
            # Read atomic numbers
            if read_atomic_numbers:
                if "XYZ" in line:
                    read_atomic_numbers = False
                else:
                    split_line = line.strip().split()
                    atomic_numbers.extend(
                        [int(float(value.replace("D", "E"))) for value in split_line]
                    )
            # Read coordinates
            if read_coordinates:
                if "FFX" in line:
                    read_coordinates = False
                else:
                    split_line = line.strip().split()
                    coordinates.extend(
                        [float(value.replace("D", "E")) for value in split_line]
                    )
            # Read Hessian
            if read_hessian:
                if "APT" in line:
                    read_hessian = False
                else:
                    split_line = line.strip().split()
                    hessian.extend(
                        [float(value.replace("D", "E")) for value in split_line]
                    )
            # Read number of atoms and normal modes
            if "NATM" in line:
                read_n_atoms = True
            # Detect when to read data
            if "AMASS" in line:
                read_masses = True
            if "ZA" in line:
                read_atomic_numbers = True
            if "XYZ" in line:
                read_coordinates = True
            if "FFX" in line:
                read_hessian = True

        # Convert data to right format
        coordinates: Array2DFloat = (
            np.array(coordinates).reshape(-1, 3) * BOHR_TO_ANGSTROM
        )
        hessian: Array2DFloat = np.array(hessian).reshape(n_atoms * 3, n_atoms * 3)
        elements = convert_elements(atomic_numbers, output="numbers")

        # Set up attributes
        self._masses = np.array(masses)
        self._fc_matrix = hessian
        self._coordinates = coordinates
        self._elements = elements

    def _parse_xtb_hessian(self, file: str | PathLike) -> None:
        # Read hessian file
        with open(file) as f:
            lines = f.readlines()

        # Parse file
        hessian = []
        for line in lines:
            try:
                hessian.extend([float(value) for value in line.strip().split()])
            except ValueError:
                pass
        # Set up force constant matrix
        dimension = int(np.sqrt(len(hessian)))
        self._fc_matrix = np.array(hessian).reshape(dimension, dimension)

    def __repr__(self) -> str:
        n_internal = len(self.internal_coordinates)
        return f"{self.__class__.__name__}({n_internal!r} internal coordinates)"


def _get_internal_coordinate(atoms: Sequence[int]) -> Bond | Angle | Dihedral:
    """Returns internal coordinate."""
    # Return bond, angle or dihedral
    if len(atoms) == 2:
        return Bond(*atoms)
    elif len(atoms) == 3:
        return Angle(*atoms)
    elif len(atoms) == 4:
        return Dihedral(*atoms)
    else:
        raise ValueError(f"Sequence of atoms must be 2-4 atoms, not {len(atoms)}.")


def cli(file: str | None = None) -> Any:
    """CLI for local force.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    if file is not None:
        elements, coordinates = read_geometry(file)
        return functools.partial(LocalForce, elements, coordinates)
    else:
        return functools.partial(LocalForce, None, None)
