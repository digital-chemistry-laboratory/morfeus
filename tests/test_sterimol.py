from pathlib import Path

from numpy.testing import assert_almost_equal, assert_array_almost_equal
import pytest

from morfeus.helpers import get_radii
from morfeus import Sterimol, read_gjf

DATA_DIR = Path(__file__).parent / "data" / "sterimol"

# Paton's test data.
@pytest.mark.parametrize("test_input,expected", [
    (("H"), (2.15, 1.09, 1.09)),
    (("Me"), (3.2, 1.7, 2.13)),
    (("Et"), (4.2, 1.7, 3.26)),
    (("iPr"), (4.2, 1.99, 3.26)),
    (("nBu"), (6.26, 1.7, 4.63)),
    (("CH2iPr"), (5.25, 1.7, 4.54)),
    (("c-Hexyl"), (6.26, 1.99, 3.58)),
    (("nPr"), (5.25, 1.7, 3.58)),
    (("adamantyl"), (6.26, 3.25, 3.58)),
    (("tBu"), (4.2, 2.85, 3.26)),
    (("CH2tBu"), (5.25, 1.7, 4.54)),
    (("CHEt2"), (5.25, 1.99, 4.54)),
    (("CHiPr2"), (5.25, 2.13, 4.54)),
    (("CHPr2"), (6.26, 1.99, 5.76)),
    (("CEt3"), (5.25, 2.85, 4.54)),
    (("Phenyl"), (6.37, 1.7, 3.2)),
    (("Bn"), (4.62, 1.7, 6.11)),
    (("4-ClC6H5"), (7.69, 1.75, 3.2)),
    (("4MePh"), (7.42, 1.7, 3.2)),
    (("4MeOPh"), (8.29, 1.81, 3.2)),
    (("35diMePh"), (6.37, 1.7, 4.39)),
    (("1Nap"), (6.37, 1.7, 5.59)),
])

def test_reference(test_input, expected):
    coordinates = test_input
    L, B_1, B_5 = expected
    
    # Use Paton's Bondi radii of 1.09 for H
    elements, coordinates = read_gjf(DATA_DIR / f"{test_input}.gjf")
    radii = get_radii(elements, radii_type="bondi")
    radii = [1.09 if radius == 1.20 else radius for radius in radii]    
    sterimol = Sterimol(elements, coordinates, 1, 2, radii=radii)

    assert_almost_equal(L, sterimol.L_value, decimal=2)
    assert_almost_equal(B_1, sterimol.B_1_value, decimal=2)
    assert_almost_equal(B_5, sterimol.B_5_value, decimal=2)