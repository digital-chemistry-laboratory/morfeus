"""Command line script to calculate solvent accessible surface area"""

import argparse
from steriplus import SASA, read_gjf, read_xyz

def main():
    # Parse the arguments
    parser = argparse.ArgumentParser(
        "Steriplus script to calcaulate solvent accessible surface areas "
        "and volumes beneath this area")
    parser.add_argument(
        'file', type=str, help='Input file, either .xyz, .gjf or .com')
    parser.add_argument(
        '--radii', type=str, help='Type of radii, either "crc" or "bondi"',
        choices=["bondi", "crc"], default="crc")
    parser.add_argument(
        '-d', '--density', type=float,
        help='Density of points on sphere vdW surface', default=0.01)
    parser.add_argument(
        '-p', '--probe', type=float, help='Probe radius', default=1.4)

    args = parser.parse_args()
    radii_type = args.radii
    density = args.density
    probe_radius = args.probe

    # Parse the geometry file
    file = args.file
    if file[-4:] == ".xyz":
        elements, coordinates = read_xyz(file)
    elif file[-4:] == ".gjf" or file[-4:] == ".com":
        elements, coordinates = read_gjf(file)
    else:
        print("No valid input file. Use .xyz or .gjf/.com")
        return

    # Run calculation and print results
    sasa = SASA(elements, coordinates, radii_type=radii_type, density=density,
                probe_radius=probe_radius)
    sasa.print_report(verbose=True)

if __name__ == "__main__":
    main()
