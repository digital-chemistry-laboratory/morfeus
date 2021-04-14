"""Test buried volume code."""

import csv
from pathlib import Path

from numpy.testing import assert_almost_equal

from morfeus import BuriedVolume, read_xyz

DATA_DIR = Path(__file__).parent / "data" / "buried_volume"


def pytest_generate_tests(metafunc):
    """Generate test data from csv file."""
    if "bv_data" in metafunc.fixturenames:
        with open(DATA_DIR / "reference_data.csv") as f:
            reader = csv.DictReader(f)
            records = list(reader)
        metafunc.parametrize("bv_data", records)


def test_reference(bv_data):
    """Test against buried volume reference data."""
    data = bv_data
    idx = int(data["idx"])
    percent_buried_volume_ref = float(data["buried_volume"])
    excluded_atoms = [int(i) for i in data["excluded_atoms"].split()]
    elements, coordinates = read_xyz(DATA_DIR / "xyz" / f"{idx}.xyz")
    bv = BuriedVolume(elements, coordinates, 1, excluded_atoms=excluded_atoms)

    assert_almost_equal(
        bv.percent_buried_volume * 100, percent_buried_volume_ref, decimal=0
    )
