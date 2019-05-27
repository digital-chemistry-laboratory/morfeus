"""Geometry file parsing functions.

Functions:
    read_gjf: Parses Gaussian input files.
    read_xyz: Parses xyz files.
"""
import numpy as np

from steriplus.helpers import convert_elements

class D3Parser:
    """Parses the output of Grimme's D3 program and extracts the C6(AA) and
    C8(AA) coefficients
    
    Args:
        filename (str): File containing output from the D3 program.
    
    Attributes:
        c6_coefficients (list): C6(AA) coefficients (au)
        c8_coefficients (list): C8(AA) coefficients (au)
    """
    def __init__(self, filename):
        # Read the file
        lines = open(filename, encoding="utf-8").readlines()
        
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
        self._filename = filename
        self.c6_coefficients = c6_coefficients
        self.c8_coefficients = c8_coefficients
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self._filename!r})"

class D4Parser:
    """Parses the output of Grimme's D4 program and extracts the C6(AA) and
    C8(AA) coefficients
    
    Args:
        filename (str): Name of file containing output from the D3 program.
    
    Attributes:
        c6_coefficients (list): C6(AA) coefficients (au)
        c8_coefficients (list): C8(AA) coefficients (au)
    """    
    def __init__(self, filename):
        # Read the file
        lines = open(filename, encoding="utf-8").readlines()
        
        # Parse the file and extract the coefficients
        c6_coefficients = []
        c8_coefficients = []
        read = False 
        for line in lines:
            if read:
                if not line.strip():
                    read = False
                    break
                strip_line = line.strip().split()
                c6 = float(strip_line[5])
                c8 = float(strip_line[6])
                c6_coefficients.append(c6)
                c8_coefficients.append(c8)
            if "C6AA" in line:
                read = True
        
        # Set up attributes
        self._filename = filename
        self.c6_coefficients = c6_coefficients
        self.c8_coefficients = c8_coefficients

    def __repr__(self):
        return f"{self.__class__.__name__}({self._filename!r})"

class VertexParser:
    """Parses the contents of a Multiwfn vtx.pdb file and extracts the
    vertices of the surface.

    Args:
        filename (str): Name of file containing the vertices.

    Attributes:
        atom_vertices (dict): Atom indices (starting from 1) as keys and
                              lists of vertices as values
        vertices (list): All vertices
    """
    def __init__(self, filename):
        # Read the datafile with NumPy
        data = np.genfromtxt(filename, 
            delimiter=[6, 5, 3, 8, 4, 12, 8, 8, 6, 6, 12],
            comments="REMARK", skip_footer=1)
        
        # Extract the vertices
        vertices = data[:, [5, 6, 7]]

        # Set up attributes. Dictionary with vertices per atom.
        self.vertices = vertices
        self._filename = filename
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self._filename!r}, " \
        f"{len(self.vertices)!r} points)"

def read_gjf(filename):
    """Reads Gaussian gjf/com file and returns elements as they are written in
    the file (either atomic numbers or symbols) as well as coordinates.

    Args:
        filename (str): Name of gjf/com file

    Returns:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
    """
    # Read file and split lines
    lines = open(filename).readlines()
    lines = [line.strip().split() for line in lines]

    # Loop over lines and store elements and coordinates 
    elements = []
    coordinates = []
    empty_counter = 0
    read_counter = 0
    for line in lines:
        if not line:
            empty_counter += 1
        if empty_counter == 2:
            if read_counter > 1:
                atom = line[0]
                if atom.isdigit():
                    atom = int(atom)
                elements.append(atom)
                coordinates.append([float(line[1]), float(line[2]),
                                    float(line[3])])
            read_counter += 1
    return elements, coordinates

def read_xyz(filename):
    """Reads xyz file and returns elements as they are written in the xyz
    file (either atomic numbers or symbols) as well as coordinates.

    Args:
        filename (str): Name of xyz file

    Returns:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
    """
    # Read file and split lines
    lines = open(filename).readlines()[2:]
    lines = [line.strip().split() for line in lines if line.strip().split()]

    # Loop over lines and store elements and coordinates
    elements = []
    coordinates = []
    for line in lines:
        atom = line[0]
        if atom.isdigit():
            atom = int(atom)
        elements.append(atom)
        coordinates.append([float(line[1]), float(line[2]), float(line[3])])
    
    return elements, coordinates

def write_xyz(xyz_file, elements, coordinates):
    """Writes xyz file from elements and coordinates.
    
    Args:
        elements (list): Elements as atomic symbols or numbers
        coordinates (list): Coordinates (Å)
    """
    # Convert elements to symbols
    elements = convert_elements(elements, output='symbols')

    # Write the xyz file
    with open(xyz_file, 'w') as file:
        file.write(f"{len(elements)}\n")
        file.write("\n")
        for element, coord in zip(elements, coordinates):
            file.write(f"{element:10s}{coord[0]:10.6f}" \
                f"{coord[1]:10.6f}{coord[2]:10.6f}\n")