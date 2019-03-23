from steriplus.data import atomic_symbols
try:
    from rdkit import Chem
    rdkit = True
except:
    rdkit = False

def create_rdkit_mol(element_ids, coordinates):
    """Creates a RDKit Mol object from element_ids and coordinates. This object
    has no bonding information so cannot be used for much else than computing
    solvent accesible surface areas.

    Args:
        coordinates (list)  :   List of atomic coordinates (Å)
        element_ids (list)  :   List of atomic numbers or symbols

    Returns:
        mol (object)        :   RDkit Mol object.
    """
    if not rdkit:
        raise Exception("RDKit not available.")

    rdkit_string = f"""\



{len(element_ids):>3d}  0  0  0  0  0  0  0  0  0999 V2000
"""
    for coordinate, element_id in zip(coordinates, element_ids):
        if type(element_id) == int:
            element_id = atomic_symbols[element_id]
        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]
        rdkit_string += f"{x:>10.4f}{y:>10.4f}{z:>10.4f} {element_id}   0  0  0  0  0  0  0  0  0  0  0  0\n"
    rdkit_string += "M  END"
    mol = Chem.MolFromMolBlock(rdkit_string)

    return mol

def read_xyz(filename):
    """Reads xyz file and returns element_ids as they are written in the xyz
    file (either atomic numbers or symbols).

    Args:
        filename (str)          :     Name of xyz file
    Returns:
        element_id_list (list)  :     List of element ids (atomic numbers or
                                      symbols)
        coord_list (list)       :     List of atomic coordinates in Å
    """
    lines = open(filename).readlines()[2:]
    lines = [line.strip().split() for line in lines]
    element_id_list = []
    coord_list = []
    for line in lines:
        atom = line[0]
        if atom.isdigit():
            atom = int(atom)
        element_id_list.append(atom)
        coord_list.append([float(line[1]), float(line[2]), float(line[3])])
    return element_id_list, coord_list

def read_gjf(filename):
    """Reads gjf/com file and returns element_ids as they are written in the
    file (either atomic numbers or symbols).

    Args:
        filename (str)          :     Name of xyz file
    Returns:
        element_id_list (list)  :     List of element ids (atomic numbers or
                                      symbols)
        coord_list (list)       :     List of atomic coordinates in Å
    """
    lines = open(filename).readlines()
    lines = [line.strip().split() for line in lines]
    element_id_list = []
    coord_list = []
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
                element_id_list.append(atom)
                coord_list.append([float(line[1]), float(line[2]), float(line[3])])
            read_counter += 1
    return element_id_list, coord_list
