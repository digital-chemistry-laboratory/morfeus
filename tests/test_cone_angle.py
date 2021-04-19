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
