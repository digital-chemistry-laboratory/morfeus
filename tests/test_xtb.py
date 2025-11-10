"""Test XTB code."""

from pathlib import Path

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

from morfeus import read_xyz, XTB

DATA_DIR = Path(__file__).parent / "data" / "xtb"


@pytest.mark.xtb
def test_fukui():
    """Test Fukui coefficients."""
    elements, coordinates = read_xyz(DATA_DIR / "1-penten-3-one.xyz")
    xtb = XTB(elements, coordinates)

    ref_data = np.genfromtxt(DATA_DIR / "fukui.csv", delimiter=",", names=True)
    f_nuc = list(xtb.get_fukui(variety="nucleophilicity").values())
    assert_array_almost_equal(f_nuc, ref_data["f_nuc"], decimal=3)
    f_elec = list(xtb.get_fukui(variety="electrophilicity").values())
    assert_array_almost_equal(f_elec, ref_data["f_elec"], decimal=3)
    f_rad = list(xtb.get_fukui(variety="radical").values())
    assert_array_almost_equal(f_rad, ref_data["f_rad"], decimal=3)
    f_dual = list(xtb.get_fukui(variety="dual").values())
    assert_array_almost_equal(f_dual, ref_data["f_dual"], decimal=3)
    f_loc_nuc = list(xtb.get_fukui(variety="local_nucleophilicity").values())
    assert_array_almost_equal(f_loc_nuc, ref_data["f_loc_nuc"], decimal=3)
    f_loc_elec = list(
        xtb.get_fukui(variety="local_electrophilicity", corrected=False).values()
    )
    assert_array_almost_equal(f_loc_elec, ref_data["f_loc_elec"], decimal=3)
