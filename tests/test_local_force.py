"""Test dispersion code."""

import csv
from pathlib import Path

from numpy.testing import assert_almost_equal

from morfeus import LocalForce

DATA_DIR = Path(__file__).parent / "data" / "local_force"


def pytest_generate_tests(metafunc):
    """Generate test data from csv file."""
    if "lf_data" in metafunc.fixturenames:
        with open(DATA_DIR / "reference_data.csv") as f:
            reader = csv.DictReader(f)
            records = list(reader)
        metafunc.parametrize("lf_data", records)


def test_reference(lf_data):
    """Test against local force reference data."""
    # Test criteria are relatively loose as reference geometries don't exist
    data = lf_data

    fc_ref = float(data["force_constant"])
    freq_ref = float(data["frequency"])
    atoms = [int(i) for i in data["atoms"].split()]
    idx = data["label"].split(".")[0]
    for file_type in ["hp", "lm"]:
        lf = LocalForce()
        lf.load_file(DATA_DIR / f"{idx}" / f"freq-{file_type}.log", "gaussian", "log")
        lf.compute_local()
        lf.compute_frequencies()
        fc = lf.get_local_force_constant(atoms)
        freq = lf.get_local_frequency(atoms)

        assert_almost_equal(fc, fc_ref, decimal=1)
        assert_almost_equal(freq, freq_ref, decimal=-2)
