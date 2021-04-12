"""Command line script to calculate dispersion descriptor."""

import argparse

from morfeus import Dispersion, read_gjf, read_xyz


def main() -> None:  # noqa: C901
    """Calculate dispersion descriptor."""
    # Add arguments
    parser = argparse.ArgumentParser(
        "morfeus script to calcaulate dispersion descriptor."
    )
    parser.add_argument("file", type=str, help="Input file, either .xyz, .gjf or .com")
    parser.add_argument(
        "--density",
        type=float,
        help="Density of points on sphere vdW surface (default: 0.1)",
        default=0.1,
    )
    parser.add_argument("--cube_file", type=str, help="Cube file of electron density.")
    parser.add_argument("--d3_file", type=str, help="Output file of D3 program.")
    parser.add_argument("--d4_file", type=str, help="Output file of D4 program.")
    parser.add_argument("--vertex_file", type=str, help="Vertex file from Multiwfn")
    parser.add_argument(
        "--isodensity",
        type=float,
        help="Isodensity value for  density from cube file.",
        default=0.001,
    )
    parser.add_argument("--verbose", help="Print atom areas", action="store_true")

    # Parse the arguments
    args = parser.parse_args()
    density = args.density
    cube_file = args.cube_file
    d3_file = args.d3_file
    d4_file = args.d4_file
    vertex_file = args.vertex_file
    verbose = args.verbose
    isodensity = args.isodensity

    # Perform checks to ensure no inconsistencies in input
    # TODO fix this
    if cube_file and vertex_file:
        raise Exception("Cannot give both cube and Multiwfn vertex file.")
    if d3_file and d4_file:
        raise Exception("Cannot give both files for D3 and D4.")

    # Parse the geometry file
    file = args.file
    if file[-4:] == ".xyz":
        elements, coordinates = read_xyz(file)
    elif file[-4:] == ".gjf" or file[-4:] == ".com":
        elements, coordinates = read_gjf(file)
    else:
        raise Exception("No valid input file. Use .xyz or .gjf/.com")

    # Check if surface or coefficients files are given
    point_surface = True
    calculate_coefficients = True
    if cube_file or vertex_file:
        point_surface = False
    if d3_file or d4_file:
        calculate_coefficients = False

    # Set up Dispersion object
    dispersion = Dispersion(
        elements,
        coordinates,
        density=density,
        point_surface=point_surface,
        compute_coefficients=calculate_coefficients,
    )

    # Set up surface and coefficients if files are given
    if cube_file:
        dispersion.surface_from_cube(cube_file, isodensity=isodensity)
    elif vertex_file:
        dispersion.surface_from_multiwfn(vertex_file)

    if d3_file:
        dispersion.load_coefficients(d3_file, model="d3")
    elif d4_file:
        dispersion.load_coefficients(d4_file, model="d4")
    else:
        dispersion.compute_coefficients()

    # Perform the calculations and print the results
    if not point_surface or not calculate_coefficients:
        dispersion.compute_p_int()

    # Print report
    dispersion.print_report(verbose=verbose)


if __name__ == "__main__":
    main()
