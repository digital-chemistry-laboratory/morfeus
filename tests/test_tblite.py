"""Test native tblite code."""

import numpy as np
import pytest

from morfeus import TBLite

pytest.importorskip("tblite.interface", exc_type=ImportError)

pytestmark = pytest.mark.tblite

ELEMENTS = ["O", "H", "H"]
COORDINATES = np.array(
    [
        [0.000000, 0.000000, 0.000000],
        [0.758602, 0.000000, 0.504284],
        [-0.758602, 0.000000, 0.504284],
    ]
)


def test_tblite_single_point_properties():
    """Test supported single-point properties."""
    tblite = TBLite(ELEMENTS, COORDINATES)

    energy = tblite.get_energy()
    charges = tblite.get_charges()
    bond_orders = tblite.get_bond_orders()
    dipole = tblite.get_dipole()

    assert isinstance(energy, float)
    assert len(charges) == len(ELEMENTS)
    assert sum(charges.values()) == pytest.approx(0, abs=1e-6)
    assert tblite.get_bond_order(1, 2) == bond_orders[(1, 2)]
    assert tblite.get_bond_order(1, 3) == bond_orders[(1, 3)]
    assert dipole.shape == (3,)
    assert tblite.get_dipole_moment(unit="debye") > 0


def test_tblite_frontier_orbitals():
    """Test frontier orbital access."""
    tblite = TBLite(ELEMENTS, COORDINATES)

    homo = tblite.get_homo()
    lumo = tblite.get_lumo()

    assert homo < lumo
    assert tblite.get_homo_lumo_gap(unit="eV") > 0


def test_tblite_unsupported_descriptors():
    """Test clear errors for xtb CLI-only descriptors."""
    tblite = TBLite(ELEMENTS, COORDINATES)

    with pytest.raises(NotImplementedError, match="xtb CLI backend"):
        tblite.get_ip()
    with pytest.raises(NotImplementedError, match="xtb CLI backend"):
        tblite.get_ea()
    with pytest.raises(NotImplementedError, match="xtb CLI backend"):
        tblite.get_fukui("nucleophilicity")


def test_tblite_validation():
    """Test validation of unsupported options."""
    with pytest.raises(ValueError, match="Method"):
        TBLite(ELEMENTS, COORDINATES, method="ptb")
    with pytest.raises(ValueError, match="Solvation model"):
        TBLite(ELEMENTS, COORDINATES, solvation_model="unknown")
