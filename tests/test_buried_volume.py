"""Test buried volume code."""

from pathlib import Path

from numpy.testing import assert_almost_equal
import pytest

from morfeus import BuriedVolume, read_xyz

DATA_DIR = Path(__file__).parent / "data" / "buried_volume"


# Buried volume test data.
@pytest.mark.parametrize(
    "idx,excluded_atoms,expected",
    [
        ((1), (1, 2, 3, 4, 5, 6, 7), (29.6)),
        ((2), (1, 2, 3, 4, 5, 6, 7), (36.1)),
        ((3), (1, 2, 3, 4, 5, 6), (35.9)),
        ((4), (1, 2, 3), (47.0)),
        ((5), (1, 2, 3), (43.8)),
        ((6), (1, 2, 3), (53.9)),
        ((7), (1, 2, 3), (49.1)),
        ((8), (1, 2, 3), (52.8)),
        ((9), (1, 2, 3), (47.7)),
        ((10), (1, 2, 3), (40.2)),
        ((11), (1, 2, 3), (35.9)),
        ((12), (1, 2, 3), (36.0)),
        ((13), (1, 2, 3), (64.2)),
        ((14), (1, 2, 3), (65.5)),
        ((15), (1, 2, 3), (60.4)),
        ((16), (1, 2, 3), (64.9)),
        ((17), (1, 2, 3), (66.8)),
        ((18), (1, 2, 3), (64.3)),
    ],
)
def test_reference(idx, excluded_atoms, expected):
    """Test against buried volume reference data."""
    elements, coordinates = read_xyz(DATA_DIR / f"{idx}.xyz")
    bv = BuriedVolume(elements, coordinates, 1, excluded_atoms=excluded_atoms)

    assert_almost_equal(bv.percent_buried_volume * 100, expected, decimal=0)
