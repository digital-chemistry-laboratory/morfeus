"""Command line script to calculate Sterimol parameters."""
import argparse

from morfeus import read_gjf, read_xyz, Sterimol

def main():
    # Parse the arguments
    parser = argparse.ArgumentParser(
        "morfeus script to calculate Sterimol values")
    parser.add_argument(
        'file', type=str, help='Input file, either .xyz, .gjf or .com')
    parser.add_argument(
        'atom_1', type=int, help='Dummy atom')
    parser.add_argument(
        'atom_2', type=int, help='Atom of substituent connected to dummy atom')
    parser.add_argument(
        '--radii', type=str, help='Radii type: "bondi" or "crc" (default)',
        choices=["bondi", "crc"], default="crc")
    parser.add_argument(
        "--verbose", help='Additional print', action='store_true')

    args = parser.parse_args()
    radii_type = args.radii
    verbose = args.verbose
    atom_1 = args.atom_1
    atom_2 = args.atom_2

    # Parse the geometry file
    file = args.file
    if file[-4:] == ".xyz":
        elements, coordinates = read_xyz(file)
    elif file[-4:] == ".gjf" or file[-4:] == ".com":
        elements, coordinates = read_gjf(file)
    else:
        raise Exception("No valid input file. Use .xyz or .gjf/.com")

    # Run the calculations and print the results
    sterimol = Sterimol(elements, coordinates, atom_1, atom_2,
                        radii_type=radii_type)
    sterimol.print_report(verbose=verbose)

if __name__ == "__main__":
    main()
