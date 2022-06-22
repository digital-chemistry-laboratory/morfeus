"""Test SASA code."""

import csv
from pathlib import Path

import numpy as np
from numpy.testing import assert_almost_equal, assert_array_almost_equal
import pytest

from morfeus import read_xyz, SASA
from morfeus.typing import Array1DFloat

DATA_DIR = Path(__file__).parent / "data" / "sasa"


def test_one():
    """Test complex one."""
    atom_areas_ref = np.loadtxt(DATA_DIR / "atom_areas/1.txt")

    elements, coordinates = read_xyz(DATA_DIR / "xyz" / "1.xyz")
    sasa = SASA(elements, coordinates)
    atom_areas: Array1DFloat = np.array(list(sasa.atom_areas.values()))
    assert_almost_equal(sasa.area, 624, decimal=0)
    assert_array_almost_equal(atom_areas, atom_areas_ref, decimal=0)


def pytest_generate_tests(metafunc):
    """Generate test data from csv file."""
    if "sasa_data" in metafunc.fixturenames:
        with open(DATA_DIR / "areas.csv") as f:
            reader = csv.DictReader(f)
            records = list(reader)
            records = [(i + 1, record) for i, record in enumerate(records)]
        metafunc.parametrize("sasa_data", records)


@pytest.mark.benchmark
def test_reference(sasa_data):
    """Test against SASA reference data."""
    idx, data = sasa_data
    sasa_ref = float(data["area"])
    atom_areas_ref = np.loadtxt(DATA_DIR / f"atom_areas/{idx}.txt")
    elements, coordinates = read_xyz(DATA_DIR / "xyz" / f"{idx}.xyz")
    sasa = SASA(elements, coordinates)
    atom_areas: Array1DFloat = np.array(list(sasa.atom_areas.values()))
    assert_almost_equal(sasa.area, sasa_ref, decimal=0)
    assert_array_almost_equal(atom_areas, atom_areas_ref, decimal=0)
