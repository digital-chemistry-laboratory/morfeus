"""Test dispersion code."""

import csv
from pathlib import Path

import numpy as np
from numpy.testing import assert_almost_equal

from morfeus import Dispersion

DATA_DIR = Path(__file__).parent / "data" / "dispersion"


def pytest_generate_tests(metafunc):
    """Generate test data from csv file."""
    if "disp_data" in metafunc.fixturenames:
        with open(DATA_DIR / "periodic_table.csv") as f:
            reader = csv.DictReader(f)
            records = list(reader)
        metafunc.parametrize("disp_data", records)


def test_reference(disp_data):
    """Test against cone angle reference data."""
    data = disp_data
    p_int_ref = float(data["p_int"])
    elements = [data["symbol"]]
    coordinates = [[0.0, 0.0, 0.0]]
    radii = [float(data["radius"])]
    disp = Dispersion(elements, coordinates, radii=radii, compute_coefficients=False)
    disp._c_n_coefficients = {
        6: np.array([float(data["c6"])]),
        8: np.array([float(data["c8"])]),
    }
    disp.compute_p_int()
    assert_almost_equal(disp.p_int, p_int_ref, decimal=1)
