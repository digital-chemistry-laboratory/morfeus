"""Test Sterimol code."""

import csv
from pathlib import Path

from numpy.testing import assert_almost_equal
import pytest

from morfeus import read_gjf, Sterimol
from morfeus.utils import get_radii

DATA_DIR = Path(__file__).parent / "data" / "sterimol"


def test_H():
    """Test H substituent."""
    elements, coordinates = read_gjf(DATA_DIR / "gjfs" / "H.gjf")
    radii = get_radii(elements, radii_type="bondi")
    radii = [1.09 if radius == 1.20 else radius for radius in radii]
    sterimol = Sterimol(elements, coordinates, 1, 2, radii=radii)
    assert_almost_equal(
        (sterimol.L_value, sterimol.B_1_value, sterimol.B_5_value),
        (2.15, 1.09, 1.09),
        decimal=2,
    )


def pytest_generate_tests(metafunc):
    """Generate test data from csv file."""
    if "sterimol_data" in metafunc.fixturenames:
        with open(DATA_DIR / "reference_data.csv") as f:
            reader = csv.DictReader(f)
            records = list(reader)
        metafunc.parametrize("sterimol_data", records)


@pytest.mark.benchmark
def test_reference(sterimol_data):
    """Test against Sterimol reference data."""
    data = sterimol_data
    name = data["name"]
    L, B_1, B_5 = float(data["L"]), float(data["B_1"]), float(data["B_5"])

    # Use Paton's Bondi radii of 1.09 for H
    elements, coordinates = read_gjf(DATA_DIR / "gjfs" / f"{name}.gjf")
    radii = get_radii(elements, radii_type="bondi")
    radii = [1.09 if radius == 1.20 else radius for radius in radii]
    sterimol = Sterimol(elements, coordinates, 1, 2, radii=radii)

    assert_almost_equal(L, sterimol.L_value, decimal=2)
    assert_almost_equal(B_1, sterimol.B_1_value, decimal=2)
    assert_almost_equal(B_5, sterimol.B_5_value, decimal=2)
