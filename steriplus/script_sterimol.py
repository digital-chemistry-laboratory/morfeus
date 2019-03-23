import argparse
from steriplus import Sterimol, read_gjf, read_xyz

def main():
    # Parse the arguments
    parser = argparse.ArgumentParser("Steriplus program to calcaulate Sterimol values")
    parser.add_argument('file', type=str, help='Input file, either .xyz, .gjf or .com')
    parser.add_argument('atom1', type=int, help='Dummy atom',)
    parser.add_argument('atom2', type=int, help='Atom of substituent connected to dummy atom')
    parser.add_argument('--density', type=float, help='Atom of substituent connected to dummy atom', default=0.005)
    parser.add_argument('--radii', type=str, help='Type of radii, either "crc" or "bondi"', choices=["bondi", "crc"], default="crc")
    parser.add_argument('-v', "--verbose", action='store_true')

    args = parser.parse_args()

    file = args.file
    if file[-4:] == ".xyz":
        elements, coordinates = read_xyz(file)
    elif file[-4:] == ".gjf" or file[-4:] == ".com":
        elements, coordinates = read_gjf(file)
    else:
        print("No valid input file. Use .xyz or .gjf/.com")
        return

    radii_type = args.radii
    density = args.density
    verbose = args.verbose

    atom_1 = args.atom1
    atom_2 = args.atom2

    sterimol = Sterimol(elements, coordinates, atom_1, atom_2, radii_type=radii_type, density=density)
    sterimol.print_report(verbose=verbose)

if __name__ == "__main__":
    main()
