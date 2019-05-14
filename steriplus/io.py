"""Geometry file parsing functions.

Functions:
    read_gjf: Parses Gaussian input files.
    read_xyz: Parses xyz files.
"""
def read_gjf(filename):
    """Reads Gaussian gjf/com file and returns elements as they are written in the
    file (either atomic numbers or symbols) as well as coordinates.

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