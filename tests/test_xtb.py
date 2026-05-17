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

    ref_data = np.genfromtxt(DATA_DIR / "xtb_fukui.csv", delimiter=",", names=True)
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


@pytest.mark.xtb
def test_charges():
    """Test charges with different methods and solvents."""
    elements, coordinates = read_xyz(DATA_DIR / "1-penten-3-one.xyz")
    ref_data = np.genfromtxt(DATA_DIR / "xtb_charges.csv", delimiter=",", names=True)

    methods: tuple[int | str, ...] = (2, 1, "ptb")
    for method, name_method in zip(methods, ("gfn2", "gfn1", "ptb")):
        for solvent, name_solvent in zip(
            [None, "water", "hexane"], ["vacuum", "water", "hexane"]
        ):
            if method == "ptb" and solvent is not None:
                continue
            xtb = XTB(elements, coordinates, method=method, solvent=solvent)
            charges = list(xtb.get_charges().values())
            assert_array_almost_equal(
                charges, ref_data[f"{name_method}_{name_solvent}_mulliken"], decimal=3
            )
            if method == 1:
                charges_cm5 = list(xtb.get_charges(model="CM5").values())
                assert_array_almost_equal(
                    charges_cm5,
                    ref_data[f"{name_method}_{name_solvent}_cm5"],
                    decimal=3,
                )
