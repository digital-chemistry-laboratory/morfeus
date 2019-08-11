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
        '--fchk_file', type=str, help='Gaussian fchk file')
    parser.add_argument(
        '--pes_file', type=str, help='Gaussian PES file')
    parser.add_argument(
        '-a', '--atoms', type=int, nargs="+",
        help='Atoms of bond/internal coordinate')
    parser.add_argument(
        '--no_project_imag', help='Do not project out imaginary modes',
        action='store_false')
    parser.add_argument(
        '--cutoff', type=float, help='Cutoff for low-frequency modes (mDyne/Å)'
        '(default: 0.001)', default=0.001)
    parser.add_argument(
        '--method', type=str, help='Method: "local" (default) or "compliance"',
        default='local')

    # Parse arguments
    args = parser.parse_args()
    atoms = args.atoms
    log_file = args.file
    fchk_file = args.fchk_file
    pes_file = args.pes_file
    project_imag = args.no_project_imag
    cutoff = args.cutoff
    method = args.method

    # Run calculation and print results
    lf = LocalForce(log_file, fchk_file=fchk_file, pes_file=pes_file,
        cutoff=cutoff, project_imag=project_imag, method=method)
    if atoms:
        force_constant = lf.get_local_force_constant(atoms)
        print(f"{force_constant:.3f}")
    else:
        lf.print_report()

if __name__ == "__main__":
    main()
