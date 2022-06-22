"""Test solid angle code."""

import csv
from pathlib import Path

from numpy.testing import assert_almost_equal
import pytest

from morfeus import read_xyz, SolidAngle

DATA_DIR = Path(__file__).parent / "data" / "solid_angle"


def test_PdPMe3():
    """Test PdPMe3."""
    elements, coordinates = read_xyz(DATA_DIR / "pd/PdPMe3.xyz")
    sa = SolidAngle(elements, coordinates, 1, radii_type="bondi")
    assert_almost_equal(sa.solid_angle, 2.877080666, decimal=1)
    assert_almost_equal(sa.cone_angle, 114.3468001, decimal=0)


def pytest_generate_tests(metafunc):
    """Generate test data from csv file."""
    if "solid_angle_data" in metafunc.fixturenames:
        with open(DATA_DIR / "solid_angles.csv") as f:
            reader = csv.DictReader(f)
            records = [("standard", record) for record in reader]
        with open(DATA_DIR / "solid_angles_max.csv") as f:
            reader = csv.DictReader(f)
            records += [("maximum", record) for record in reader]
        metafunc.parametrize("solid_angle_data", records)


@pytest.mark.benchmark
def test_reference(solid_angle_data):
    """Test against solid angle reference data."""
    for metal in ("pd", "pt", "ni"):
        label, data = solid_angle_data
        solid_angle_ref = float(data[f"{metal}_solid_angle"])
        cone_angle_ref = float(data[f"{metal}_cone_angle"])
        if label == "standard":
            xyz_path = DATA_DIR / f"{metal}/{data[f'{metal}_xyz']}.xyz"
        elif label == "maximum":
            xyz_path = DATA_DIR / f"{metal}/maximum/{data[f'{metal}_xyz']}.xyz"
        elements, coordinates = read_xyz(xyz_path)
        sa = SolidAngle(elements, coordinates, 1, radii_type="bondi")
        assert_almost_equal(sa.solid_angle, solid_angle_ref, decimal=1)
        assert_almost_equal(sa.cone_angle, cone_angle_ref, decimal=0)
