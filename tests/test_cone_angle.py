"""Test cone angle code."""

import csv
from pathlib import Path

from numpy.testing import assert_almost_equal
import pytest

from morfeus import ConeAngle, read_xyz

DATA_DIR = Path(__file__).parent / "data" / "cone_angle"


def test_PdPMe3():
    """Test PdPMe3."""
    elements, coordinates = read_xyz(DATA_DIR / "pd/PdPMe3.xyz")
    ca = ConeAngle(elements, coordinates, 1, radii_type="bondi")
    assert_almost_equal(ca.cone_angle, 120.4, decimal=1)


def pytest_generate_tests(metafunc):
    """Generate test data from csv file."""
    if "cone_angle_data" in metafunc.fixturenames:
        with open(DATA_DIR / "cone_angles.csv") as f:
            reader = csv.DictReader(f)
            records = [("standard", record) for record in reader]
        with open(DATA_DIR / "cone_angles_max.csv") as f:
            reader = csv.DictReader(f)
            records += [("maximum", record) for record in reader]
        metafunc.parametrize("cone_angle_data", records)


@pytest.mark.benchmark
def test_reference(cone_angle_data):
    """Test against cone angle reference data."""
    for metal in ("pd", "pt", "ni"):
        label, data = cone_angle_data
        cone_angle_ref = float(data[f"{metal}_cone_angle"])
        if label == "standard":
            xyz_path = DATA_DIR / f"{metal}/{data[f'{metal}_xyz']}.xyz"
        elif label == "maximum":
            xyz_path = DATA_DIR / f"{metal}/maximum/{data[f'{metal}_xyz']}.xyz"
        elements, coordinates = read_xyz(xyz_path)
        ca = ConeAngle(elements, coordinates, 1, radii_type="bondi")
        assert_almost_equal(ca.cone_angle, cone_angle_ref, decimal=1)


def test_reference_internal(cone_angle_data):
    """Test the internal algorithm against cone angle reference data.

    The default method's tests only cover the internal algorithm when
    libconeangle is not installed, so it is exercised explicitly here.

    Args:
        cone_angle_data: Reference data record from the csv files
    """
    for metal in ("pd", "pt", "ni"):
        label, data = cone_angle_data
        cone_angle_ref = float(data[f"{metal}_cone_angle"])
        if label == "standard":
            xyz_path = DATA_DIR / f"{metal}/{data[f'{metal}_xyz']}.xyz"
        elif label == "maximum":
            xyz_path = DATA_DIR / f"{metal}/maximum/{data[f'{metal}_xyz']}.xyz"
        elements, coordinates = read_xyz(xyz_path)
        ca = ConeAngle(elements, coordinates, 1, radii_type="bondi", method="internal")
        assert_almost_equal(ca.cone_angle, cone_angle_ref, decimal=1)


def test_degenerate_root_does_not_crash():
    """Test a structure with a numerically degenerate tangency root.

    For this structure the tangency quadratic of one atom triple has a root
    marginally outside [-1, 1], which the looped implementation fed to
    math.acos, raising 'math domain error'. Such a root is not a physical
    cone; it must be discarded and the search must complete.
    """
    elements, coordinates = read_xyz(DATA_DIR / "degenerate.xyz")
    ca = ConeAngle(elements, coordinates, len(elements), method="internal")
    assert_almost_equal(ca.cone_angle, 197.70, decimal=2)
    assert sorted(ca.tangent_atoms) == [1, 52, 75]
