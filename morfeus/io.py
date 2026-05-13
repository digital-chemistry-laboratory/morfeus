"""Input and output."""

from __future__ import annotations

from collections.abc import Iterable, MutableMapping, Sequence
from os import PathLike
from pathlib import Path
import typing
from typing import Any

import numpy as np

from morfeus.d3_data import r2_r4
from morfeus.data import BOHR_TO_ANGSTROM, KCAL_TO_HARTREE
from morfeus.typing import (
    Array1DAny,
    Array1DFloat,
    Array1DInt,
    Array1DStr,
    Array2DFloat,
    Array2DInt,
    Array3DFloat,
    ArrayLike2D,
    ArrayLike3D,
)
from morfeus.utils import convert_elements, Import, requires_dependency

if typing.TYPE_CHECKING:
    import cclib.io


class CrestParser:
    """Parses folder with CREST output files for conformer information.

    Reads the files: cre_members, crest.energies, crest_conformers.xyz

    Args:
        path: Path to CREST folder

    Attributes:
        conformer_coordinates: Conformer coordinates (Å)
        degeneracies: Degeneracies
        elements: Elements
        energies: Energies (a.u.)
    """

    conformer_coordinates: Array3DFloat
    degeneracies: Array1DInt
    elements: Array1DInt
    energies: Array1DFloat

    def __init__(self, path: str | PathLike) -> None:
        path = Path(path)

        with open(path / "cre_members") as f:
            lines = f.readlines()
        degeneracies = []
        for line in lines:
            strip_line = line.strip().split()
            if len(strip_line) != 1:
                degeneracies.append(int(strip_line[0]))
        degeneracies: Array1DInt = np.array(degeneracies)

        with open(path / "crest.energies") as f:
            lines = f.readlines()

        energies = []
        for line in lines:
            strip_line = line.strip().split()
            energies.append(float(strip_line[1]))
        energies = np.array(energies) * KCAL_TO_HARTREE

        elements, conformer_coordinates = read_xyz(path / "crest_conformers.xyz")
        if conformer_coordinates.ndim == 2:
            conformer_coordinates = conformer_coordinates.reshape(1, -1, 3)

        self.elements = elements
        self.conformer_coordinates = conformer_coordinates
        self.energies = energies
        self.degeneracies = degeneracies


class CubeParser:
    """Parses Gaussian cube file of electron density.

    Args:
        file: Cube file

    Attributes:
        min_x: Minimum x value (Å)
        min_y: Minimum y value (Å)
        min_z: Minimum z value (Å)
        step_x: Step size in x direction (Å)
        step_y: Step size in y direction (Å)
        step_z: Step size in z direction (Å)
        X: 3D array of x values (Å)
        Y: 3D array of y values (Å)
        Z: 3D array of z values (Å)
        S: 3D array of electron density scalars (electorns / Bohr^3)
        n_points_x: Number of points in the x direction
        n_points_y: Number of points in the y direction
        n_points_z: Number of points in the < direction
    """

    min_x: float
    min_y: float
    min_z: float
    step_x: float
    step_y: float
    step_z: float
    X: Array3DFloat
    Y: Array3DFloat
    Z: Array3DFloat
    S: Array3DFloat
    n_points_x: int
    n_points_y: int
    n_points_z: int

    def __init__(self, file: str | PathLike) -> None:
        # Read the lines from the cube file
        with open(file) as file:
            lines = file.readlines()

        # Skip first two lines which are comments
        lines = lines[2:]

        # Get the number of atoms
        n_atoms = int(lines[0].strip().split()[0])

        # Get the minimum values along the axes
        min_x = float(lines[0].strip().split()[1]) * BOHR_TO_ANGSTROM
        min_y = float(lines[0].strip().split()[2]) * BOHR_TO_ANGSTROM
        min_z = float(lines[0].strip().split()[3]) * BOHR_TO_ANGSTROM

        # Get the number of points and step size along each axis
        n_points_x = int(lines[1].strip().split()[0])
        step_x = float(lines[1].strip().split()[1]) * BOHR_TO_ANGSTROM

        n_points_y = int(lines[2].strip().split()[0])
        step_y = float(lines[2].strip().split()[2]) * BOHR_TO_ANGSTROM

        n_points_z = int(lines[3].strip().split()[0])
        step_z = float(lines[3].strip().split()[3]) * BOHR_TO_ANGSTROM

        # Generate grid
        x = min_x + np.arange(0, n_points_x) * step_x
        y = min_y + np.arange(0, n_points_y) * step_y
        z = min_z + np.arange(0, n_points_z) * step_z
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

        # Skips to data lines and read in data
        lines = lines[(n_atoms + 4) :]

        data = []
        for line in lines:
            line_data = [float(datum) for datum in line.strip().split()]
            data.extend(line_data)

        # Create array
        S: Array2DFloat = np.array(data).reshape(X.shape)

        # Set up attributes
        self.X = X
        self.Y = Y
        self.Z = Z
        self.S = S

        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z

        self.step_x = step_x
        self.step_y = step_y
        self.step_z = step_z

        self.n_points_x = n_points_x
        self.n_points_y = n_points_y
        self.n_points_z = n_points_z

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.S.size!r} points)"


class D3Parser:
    """Parses the output of Grimme's D3 program.

    Extracts the C₆ᴬᴬ and C₈ᴬᴬ coefficients

    Args:
        file: File containing output from the D3 program.

    Attributes:
        c6_coefficients: C₆ᴬᴬ coefficients (a.u.)
        c8_coefficients: C₈ᴬᴬ coefficients (a.u.)
    """

    c6_coefficients: Array1DFloat
    c8_coefficients: Array1DFloat

    def __init__(self, file: str | PathLike) -> None:
        # Read the file
        with open(file, encoding="utf-8") as file:
            lines = file.readlines()

        # Parse the file for the coefficients
        c6_coefficients = []
        c8_coefficients = []
        read = False
        for line in lines:
            if read:
                if not line.strip():
                    read = False
                    break
                strip_line = line.strip().split()
                c6 = float(strip_line[7])
                c8 = float(strip_line[8])
                c6_coefficients.append(c6)
                c8_coefficients.append(c8)
            if "C8(AA)" in line:
                read = True

        # Set attributes
        self.c6_coefficients = np.array(c6_coefficients)
        self.c8_coefficients = np.array(c8_coefficients)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}"


class D4Parser:
    """Parses the output of Grimme's D4 program.

    Extracts the C₆ᴬᴬ and C₈ᴬᴬ coefficients

    Args:
        file: File containing output from the D4 program.

    Attributes:
        c6_coefficients: C₆ᴬᴬ coefficients (a.u.)
        c8_coefficients: C₈ᴬᴬ coefficients (a.u.)
    """

    c6_coefficients: Array1DFloat
    c8_coefficients: Array1DFloat

    def __init__(self, file: str | PathLike) -> None:
        # Read the file
        with open(file, encoding="utf-8") as file:
            lines = file.readlines()

        # Parse the file and extract the coefficients
        c6_coefficients = []
        elements = []
        read = False
        for line in lines:
            if read:
                if not line.strip():
                    read = False
                    break
                strip_line = line.strip().split()
                element = int(strip_line[1])
                elements.append(element)
                c6 = float(strip_line[5])
                c6_coefficients.append(c6)
            if "C6AA" in line:
                read = True
        c6_coefficients: Array1DFloat = np.array(c6_coefficients)
        c8_coefficients: Array1DFloat = np.array(
            [
                3 * c6 * r2_r4[element] ** 2
                for c6, element in zip(c6_coefficients, elements)
            ]
        )

        # Set up attributes
        self.c6_coefficients = c6_coefficients
        self.c8_coefficients = c8_coefficients

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}"


class VertexParser:
    """Parses the contents of a Multiwfn vtx.pdb file.

    Extracts the vertices and faces of the surface.

    Args:
        file: Name of file containing the vertices.

    Attributes:
        faces: Faces of surface
        vertices: Vertices of surface
    """

    faces: Array2DInt | None
    vertices: Array1DFloat

    def __init__(self, file: str | PathLike) -> None:  # noqa: C901
        # Parse file to see if it containts connectivity
        with open(file) as file:
            lines = file.readlines()

        # Get the number of vertices
        n_vertices = int(lines[0].strip().split()[5])

        # Parse the vertex positions and their connectivities
        vertices = {}
        vertex_map = {}
        connectivities: dict[int, set[int]] = {
            i: set() for i in range(1, n_vertices + 1)
        }
        vertex_counter = 1
        included_vertex_counter = 1
        for line in lines:
            if "HETATM" in line:
                if line[13] == "C":
                    x = float(line[32:39])
                    y = float(line[40:47])
                    z = float(line[48:55])
                    vertices[vertex_counter] = [x, y, z]
                    vertex_map[vertex_counter] = included_vertex_counter
                    included_vertex_counter += 1
                vertex_counter += 1
            if "CONECT" in line:
                n_entries = int(len(line.strip()) / 6 - 1)
                entries = []
                for i in range(1, n_entries + 1):
                    entry = int(line[i * 6 : i * 6 + 6])
                    entries.append(entry)
                connectivities[entries[0]].update(entries[1:])

        # Establish faces based on connectivity
        # https://stackoverflow.com/questions/1705824/finding-cycle-of-3-nodes-or-triangles-in-a-graph # noqa: B950
        if any(connectivities.values()):
            faces = []
            visited = set()
            for vertex_1 in connectivities:
                temp_visited = set()
                for vertex_2 in connectivities[vertex_1]:
                    if vertex_2 in visited:
                        continue
                    for vertex_3 in connectivities[vertex_2]:
                        if vertex_3 in visited or vertex_3 in temp_visited:
                            continue
                        if vertex_1 in connectivities[vertex_3]:
                            triangle_vertices = [vertex_1, vertex_2, vertex_3]
                            mapped_vertices = [
                                vertex_map[vertex] - 1 for vertex in triangle_vertices
                            ]
                            faces.append(mapped_vertices)
                    temp_visited.add(vertex_2)
                visited.add(vertex_1)
            self.faces = np.array(faces)
        else:
            self.faces = None

        # Set up attributes.
        self.vertices = np.array(list(vertices.values()))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({len(self.vertices)!r} points)"


def read_gjf(file: str | PathLike) -> tuple[Array1DAny, Array1DFloat]:
    """Reads Gaussian gjf/com file and returns elements and coordinates.

    Args:
        file: gjf/com file object.

    Returns:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
    """
    # Read file and split lines
    with open(file) as f:
        text = f.read()
    xyz_lines = text.split("\n\n")[2].splitlines()[1:]
    xyz_lines = [line for line in xyz_lines if line.strip()]

    # Loop over lines and store elements and coordinates
    elements: list[int | str] = []
    coordinates: list[list[float]] = []
    for line in xyz_lines:
        split_line = line.strip().split()
        element = split_line[0]
        if element.isdigit():
            elements.append(int(element))
        else:
            elements.append(element)
        atom_coordinates = [float(i) for i in split_line[1:]]
        coordinates.append(atom_coordinates)
    elements: Array1DInt | Array1DStr = np.array(elements)
    coordinates: Array2DFloat = np.array(coordinates)

    return elements, coordinates


def read_geometry(
    file: str | PathLike,
) -> tuple[Array1DInt | Array1DStr, Array2DFloat | Array3DFloat]:
    """Read geometry file and guess parser based on suffix.

    Args:
        file: Filename or Path object

    Returns:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)

    Raises:
        ValueError: When suffix is not supported
    """
    if isinstance(file, str):
        suffix = file.lower().split(".")[-1]
    elif isinstance(file, PathLike):
        path = Path(file)
        suffix = path.suffix[1:]

    if suffix == "xyz":
        elements, coordinates = read_xyz(file)
    elif suffix in ("gjf", "com"):
        elements, coordinates = read_gjf(file)
    else:
        raise ValueError(f"File suffix {suffix} not supported.")

    return elements, coordinates


def read_xyz(
    file: str | PathLike,
) -> tuple[Array1DInt | Array1DStr, Array2DFloat | Array3DFloat]:
    """Reads xyz file.

    Returns elements as written (atomic numbers or symbols) and coordinates.

    Args:
        file: Filename or Path object

    Returns:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
    """
    # Read file and split lines
    with open(file) as f:
        lines = f.readlines()

    # Loop over lines and store elements and coordinates
    elements: list[int | str] = []
    coordinates: list[list[float]] = []
    n_atoms = int(lines[0].strip())
    line_chunks = zip(*[iter(lines)] * (n_atoms + 2))
    for line_chunk in line_chunks:
        for line in line_chunk[2:]:
            strip_line = line.strip().split()
            atom = strip_line[0]
            if atom.isdigit():
                atom = int(atom)
            elements.append(atom)
            coordinates.append(
                [float(strip_line[1]), float(strip_line[2]), float(strip_line[3])]
            )
    elements = np.array(elements)[:n_atoms]
    coordinates: Array2DFloat | Array3DFloat = np.array(coordinates).reshape(
        -1, n_atoms, 3
    )
    if coordinates.shape[0] == 1:
        coordinates = coordinates[0]

    return elements, coordinates


@requires_dependency([Import("cclib")], globals())
def read_cclib(file: Any) -> tuple[Array1DFloat, Array2DFloat | Array3DFloat]:
    """Reads geometry file with cclib.

    Returns elements as atomic numbers and coordinates.

    Args:
        file: Input file to cclib

    Returns:
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
    """
    data = cclib.io.ccread(file)
    elements = data.atomnos
    coordinates = data.atomcoords

    if coordinates.shape[0] == 1:
        coordinates = coordinates[0]

    return elements, coordinates


def get_xyz_string(
    symbols: Sequence[str], coordinates: Sequence[Sequence[float]], comment: str = ""
) -> str:
    """Return xyz string.

    Args:
        symbols: Atomic symbols
        coordinates: Atomic coordinates (Å)
        comment: Comment

    Returns:
        string: XYZ string
    """
    string = f"{len(symbols)}\n"
    string += f"{comment}\n"
    for s, c in zip(symbols, coordinates):
        string += f"{s:10s}{c[0]:10.5f}{c[1]:10.5f}{c[2]:10.5f}\n"

    return string


def write_xyz(
    file: str | PathLike,
    elements: Iterable[int] | Iterable[str],
    coordinates: ArrayLike2D | ArrayLike3D,
    comments: Sequence[str] | None = None,
) -> None:
    """Writes xyz file from elements and coordinates.

    Args:
        file: xyz file or path object
        elements: Elements as atomic symbols or numbers
        coordinates: Coordinates (Å)
        comments: Comments
    """
    # Convert elements to symbols
    symbols = convert_elements(elements, output="symbols")
    coordinates: Array2DFloat | Array3DFloat = np.array(coordinates).reshape(
        -1, len(symbols), 3
    )
    if comments is None:
        comments = [""] * len(coordinates)

    # Write the xyz file
    with open(file, "w") as f:
        for coord, comment in zip(coordinates, comments):
            xyz_string = get_xyz_string(symbols, coord, comment=comment)
            f.write(xyz_string)


def write_xtb_inp(
    file: str | PathLike,
    instruction_blocks: MutableMapping[str, list[str]],
) -> None:
    """Write input instructions to xtb input file.

    Args:
        file: path to the xtb input file to create
        instruction_blocks: xtb input instructions to write in the file
    """
    string = ""
    for flag, content in instruction_blocks.items():
        string += f"${flag}\n"
        for line in content:
            string += f"   {line}\n"
    string += "$end\n"
    with open(file, "w") as f:
        f.write(string)
