"""Command line script to calculate local force constants."""
import argparse

from steriplus import LocalForce

def main():
    # Add arguments
    parser = argparse.ArgumentParser(
        "Steriplus script to calculate local force constants")
    parser.add_argument(
        'file', type=str, help='Gaussian log file')
    parser.add_argument(
        '-a', '--atoms', type=int, nargs="+",
        help='Atoms of bond', required=False)

    # Parse arguments
    args = parser.parse_args()
    atoms = args.atoms
    if atoms:
        if len(atoms) != 2:
            raise Exception("Give only two atoms.")
        else:
            atom_1, atom_2 = atoms
    file = args.file

    # Run calculation and print results
    lf = LocalForce(file)
    if atoms:
        force_constant = lf.get_local_force_constant(atom_1, atom_2)
        print(f"{force_constant:.3f}")
    else:
        lf.print_report()

if __name__ == "__main__":
    main()
