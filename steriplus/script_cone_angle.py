import argparse
from steriplus import ConeAngle, read_gjf, read_xyz

def main():
    # Parse the arguments
    parser = argparse.ArgumentParser("Steriplus program to \
                                     calcaulate cone angles")
    parser.add_argument('file', type=str, help='Input file, \
                        either .xyz, .gjf or .com')
    parser.add_argument('atom1', type=int, help='Central atom',)
    parser.add_argument('--radii', type=str, help='Type of radii, either "crc" \
                        (default) or "bondi"', choices=["bondi", "crc"],
                        default="crc")

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
    atom_1 = args.atom1

    cone_angle = ConeAngle(elements, coordinates, atom_1, radii_type=radii_type)
    cone_angle.print_report()

if __name__ == "__main__":
    main()
