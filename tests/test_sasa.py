"""Test SASA code."""

import csv
from pathlib import Path

import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from morfeus import read_xyz, SASA

DATA_DIR = Path(__file__).parent / "data" / "sasa"
XYZ_DIR = Path(__file__).parent / "data" / "buried_volume"


def pytest_generate_tests(metafunc):
    """Generate test data from csv file."""
    if "sasa_data" in metafunc.fixturenames:
        with open(DATA_DIR / "areas.csv") as f:
            reader = csv.DictReader(f)
            records = list(reader)
            records = [(i + 1, record) for i, record in enumerate(records)]
        metafunc.parametrize("sasa_data", records)


def test_reference(sasa_data):
    """Test against SASA reference data."""
    idx, data = sasa_data
    sasa_ref = float(data["area"])
    atom_areas_ref = np.loadtxt(DATA_DIR / f"atom_areas/{idx}.txt")
    elements, coordinates = read_xyz(XYZ_DIR / f"{idx}.xyz")
    sasa = SASA(elements, coordinates)
    atom_areas = np.array(list(sasa.atom_areas.values()))
    assert_almost_equal(sasa.area, sasa_ref, decimal=0)
    assert_array_almost_equal(atom_areas, atom_areas_ref, decimal=0)
