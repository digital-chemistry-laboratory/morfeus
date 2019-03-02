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
