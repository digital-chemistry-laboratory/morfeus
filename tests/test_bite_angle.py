"""Test buried volume code."""

import csv
from pathlib import Path

from numpy.testing import assert_almost_equal
import pytest

from morfeus import BiteAngle, read_xyz

DATA_DIR = Path(__file__).parent / "data" / "bite_angle"


def test_one():
    """Test one complex."""
    _, coordinates = read_xyz(DATA_DIR / "xyz" / "cis-B01_XantPhos.xyz")
    bite_angle = BiteAngle(coordinates, 60, 1, 37)
    assert_almost_equal(bite_angle.angle, 101.91, decimal=2)


def test_inverted_vector():
    """Test inverted bite angle by reference vector method."""
    _, coordinates = read_xyz(DATA_DIR / "xyz" / "inverted.xyz")
    ref_vector = coordinates[5] - coordinates[10]
    bite_angle = BiteAngle(coordinates, 11, 10, 14, ref_vector=ref_vector)
    assert_almost_equal(bite_angle.angle, 187.33, decimal=2)


def test_inverted_atoms():
    """Test inverted bite angle by reference vector method."""
    _, coordinates = read_xyz(DATA_DIR / "xyz" / "inverted.xyz")
    bite_angle = BiteAngle(coordinates, 11, 10, 14, ref_atoms=[6])
    assert_almost_equal(bite_angle.angle, 187.33, decimal=2)


def pytest_generate_tests(metafunc):
    """Generate test data from csv file."""
    if "bite_angle_data" in metafunc.fixturenames:
        with open(DATA_DIR / "reference_data.csv") as f:
            reader = csv.DictReader(f)
            records = list(reader)
        metafunc.parametrize("bite_angle_data", records)


@pytest.mark.benchmark
def test_reference(bite_angle_data):
    """Test against buried volume reference data."""
    data = bite_angle_data
    metal_idx = int(data["idx_metal"])
    ligand_idx_1 = int(data["idx_ligand_1"])
    ligand_idx_2 = int(data["idx_ligand_2"])
    bite_angle_ref = float(data["bite_angle"])
    filename = data["filename"]
    _, coordinates = read_xyz(DATA_DIR / "xyz" / filename)
    bite_angle = BiteAngle(coordinates, metal_idx, ligand_idx_1, ligand_idx_2)

    assert_almost_equal(bite_angle.angle, bite_angle_ref, decimal=2)
