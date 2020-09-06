"""Command line script to calculate local force constants."""

import argparse

from morfeus import LocalForce, read_xyz

def main():
    # Add arguments
    parser = argparse.ArgumentParser(
        "morfeus script to calculate local force constants")
    parser.add_argument(
        'input_file', type=str, help='Input file')
    parser.add_argument(
        '-c', '--coordinate', type=int, action="append", nargs="+",
        help='Atoms of internal coordinate')
    parser.add_argument(
        '-m', '--method', type=str,
        help='Method: "local" (default) or "compliance"', default='local')
    parser.add_argument(
        '-p', '--program', type=str, help="Quantum-chemical program",
        required=True)
    parser.add_argument(
        '-t', '--type', type=str, help="Filetype", required=True)
    parser.add_argument(
        "-x", '--xyz', type=str, help='xyz file')

    # Parse arguments
    args = parser.parse_args()
    internal_coordinates = args.coordinate
    input_file = args.input_file
    program = args.program
    filetype = args.filetype
    xyz_file = args.xyz
    method = args.method

    # Initialize LocalForce object
    if xyz_file:
        elements, coordinates = read_xyz(xyz_file)
        lf = LocalForce(elements, coordinates)
    else:
        lf = LocalForce()
    
    # Add internal coordinates
    if internal_coordinates:
        for coordinate in internal_coordinates:
            coordinate = [int(i) for i in coordinate]
            lf.add_internal_coordinate(coordinate)
    
    # Load input file and perform calculation
    lf.load_file(input_file, program, filetype)

    if method == "local":
        if len(lf._force_constants) < 1:
            lf.normal_mode_analysis()
        lf.compute_local()
    if method == "compliance":
        lf.compute_compliance()
    lf.compute_frequencies()

    # Print results
    lf.print_report(angles=True, dihedrals=True, angle_units=True)

if __name__ == "__main__":
    main()
