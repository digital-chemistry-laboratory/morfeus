"""Chemistry-focused integration tests for Multiwfn."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

from morfeus.multiwfn import CommandStep, Multiwfn

MULTIWFN_DATA_DIR = Path(__file__).parent / "data" / "multiwfn"
EXAMPLE_MOLDEN_DIR = MULTIWFN_DATA_DIR / "example_molden"
SINGLET_MOLDEN_FILES = sorted(EXAMPLE_MOLDEN_DIR.glob("*singlet*.molden"))
TRIPLET_MOLDEN_FILES = sorted(EXAMPLE_MOLDEN_DIR.glob("*triplet*.molden"))
ALL_MOLDEN_FILES = SINGLET_MOLDEN_FILES + TRIPLET_MOLDEN_FILES
MOLDEN_FILE_BY_NAME = {file_path.name: file_path for file_path in ALL_MOLDEN_FILES}

BOND_ORDER_PARTITION_MODELS = ("mayer", "wiberg", "fuzzy")
DOUBLE_BOND_RANGE = (1.5, 2.5)
SINGLE_BOND_RANGE = (0.75, 1.3)
DOUBLE_BOND_PAIR = (3, 4)
SINGLE_BOND_PAIRS = ((1, 2), (1, 7), (1, 8), (1, 9))
OXYGEN_INDEX = 4
HYDROGEN_INDICES = (7, 8, 9, 10, 11, 12, 13, 14)
SINGLET_TRIPLET_EQUIVALENT_FILES = (
    ("xtb_singlet.molden", "xtb_triplet.molden"),
    (
        "orca_dft_singlet_molden2aim.molden",
        "orca_dft_triplet_uks_molden2aim.molden",
    ),
    ("pyscf_dft_singlet_rks.molden", "pyscf_dft_triplet_uks.molden"),
    (
        "pyscf_casscf_singlet_rhf_natorb.molden",
        "pyscf_casscf_triplet_rohf_natorb.molden",
    ),
)
MULTIPOLE_KEYS = ("quadrupole_spherical_magnitude", "quadrupole_traceless_magnitude")


def test_parse_custom_commands_accepts_supported_step_forms() -> None:
    """Convert custom command shorthand to CommandStep objects."""
    mwfn = object.__new__(Multiwfn)
    existing_step = CommandStep("4", expect="ready", optional=True)

    steps = mwfn._parse_custom_commands(
        ["1", ("2", "prompt"), ("3", "optional prompt", True), existing_step]
    )

    assert steps == [
        CommandStep("1"),
        CommandStep("2", expect="prompt"),
        CommandStep("3", expect="optional prompt", optional=True),
        existing_step,
    ]


@pytest.mark.parametrize(
    "step",
    [
        (),
        ("1",),
        ("1", "prompt", False, "extra"),
        object(),
    ],
)
def test_parse_custom_commands_rejects_invalid_step_forms(step: object) -> None:
    """Reject custom command shorthands with unsupported shape."""
    mwfn = object.__new__(Multiwfn)

    with pytest.raises(ValueError, match="Invalid custom command"):
        mwfn._parse_custom_commands([step])


@pytest.mark.multiwfn
@pytest.mark.parametrize("file_path", ALL_MOLDEN_FILES, ids=lambda path: path.name)
def test_spin_density(file_path: Path, tmp_path: Path) -> None:
    """Compare integrated spin densities to CSV reference values."""
    _assert_atomic_property_matches_csv(
        file_path,
        tmp_path,
        "spin_density",
        MULTIWFN_DATA_DIR / "multiwfn_spin_density.csv",
    )


@pytest.mark.multiwfn
@pytest.mark.parametrize("file_path", ALL_MOLDEN_FILES, ids=lambda path: path.name)
def test_esp_nuclear(file_path: Path, tmp_path: Path) -> None:
    """Compare integrated nuclear ESP descriptors to CSV reference values."""
    _assert_atomic_property_matches_csv(
        file_path,
        tmp_path,
        "esp_nuclear",
        MULTIWFN_DATA_DIR / "multiwfn_esp_nuclear.csv",
    )


@pytest.mark.multiwfn
@pytest.mark.parametrize("file_path", ALL_MOLDEN_FILES, ids=lambda path: path.name)
def test_charges_hirshfeld(file_path: Path, tmp_path: Path) -> None:
    """Compare Hirshfeld charges to CSV reference values."""
    _assert_atomic_property_matches_csv(
        file_path,
        tmp_path,
        "charges_hirshfeld",
        MULTIWFN_DATA_DIR / "multiwfn_charges_hirshfeld.csv",
    )


@pytest.mark.multiwfn
@pytest.mark.parametrize("file_path", SINGLET_MOLDEN_FILES, ids=lambda path: path.name)
def test_singlet_bond_orders_follow_expected_chemistry(
    file_path: Path, tmp_path: Path
) -> None:
    """Assert expected single- and double-bond patterns for singlet Molden inputs."""
    mwfn = _new_mwfn(file_path, tmp_path, has_spin=False)

    for model in BOND_ORDER_PARTITION_MODELS:
        bond_orders = mwfn.get_bond_order(model=model)

        double_bond_value = _get_bond_order(bond_orders, *DOUBLE_BOND_PAIR)
        assert DOUBLE_BOND_RANGE[0] <= double_bond_value <= DOUBLE_BOND_RANGE[1], (
            f"{file_path.name} ({model}): expected bond {DOUBLE_BOND_PAIR[0]}-"
            f"{DOUBLE_BOND_PAIR[1]} in {DOUBLE_BOND_RANGE}, got {double_bond_value:.3f}"
        )

        for atom_i, atom_j in SINGLE_BOND_PAIRS:
            single_bond_value = _get_bond_order(bond_orders, atom_i, atom_j)
            assert SINGLE_BOND_RANGE[0] <= single_bond_value <= SINGLE_BOND_RANGE[1], (
                f"{file_path.name} ({model}): expected bond {atom_i}-{atom_j} in "
                f"{SINGLE_BOND_RANGE}, got {single_bond_value:.3f}"
            )


@pytest.mark.multiwfn
@pytest.mark.parametrize("file_path", SINGLET_MOLDEN_FILES, ids=lambda path: path.name)
def test_singlet_hirshfeld_charges_pattern(file_path: Path, tmp_path: Path) -> None:
    """Assert expected charge signs for oxygen and hydrogens in singlet files."""
    mwfn = _new_mwfn(file_path, tmp_path, has_spin=False)
    charges = mwfn.get_charges(model="hirshfeld")

    assert min(charges, key=lambda atom_index: charges[atom_index]) == OXYGEN_INDEX
    for hydrogen_index in HYDROGEN_INDICES:
        assert charges[hydrogen_index] > 0.0


@pytest.mark.multiwfn
@pytest.mark.parametrize("file_path", ALL_MOLDEN_FILES, ids=lambda path: path.name)
def test_rho_descriptor_max(file_path: Path, tmp_path: Path) -> None:
    """Assert oxygen has the largest integrated electron density contribution."""
    mwfn = _new_mwfn(file_path, tmp_path)
    rho = mwfn.get_descriptor("rho")

    assert max(rho, key=lambda atom_index: rho[atom_index]) == OXYGEN_INDEX


@pytest.mark.multiwfn
@pytest.mark.parametrize("file_path", ALL_MOLDEN_FILES, ids=lambda path: path.name)
def test_polarity(file_path: Path, tmp_path: Path) -> None:
    """Assert all molecules show nonzero dipole and multipole magnitudes."""
    mwfn = _new_mwfn(file_path, tmp_path)
    moments = mwfn.get_electric_moments()

    assert abs(moments["dipole_magnitude_au"]) > 0.0
    multipoles = [moments[key] for key in MULTIPOLE_KEYS if key in moments]
    assert multipoles, f"{file_path.name}: no multipole value was parsed."
    assert any(abs(value) > 0.0 for value in multipoles)


@pytest.mark.multiwfn
@pytest.mark.parametrize("model", BOND_ORDER_PARTITION_MODELS)
@pytest.mark.parametrize(
    "singlet_name,triplet_name",
    SINGLET_TRIPLET_EQUIVALENT_FILES,
)
def test_bond_order_change(
    singlet_name: str, triplet_name: str, model: str, tmp_path: Path
) -> None:
    """Assert C=O bond order decreases going from singlet to equivalent triplet."""
    singlet_file = MOLDEN_FILE_BY_NAME[singlet_name]
    triplet_file = MOLDEN_FILE_BY_NAME[triplet_name]

    singlet_mwfn = _new_mwfn(singlet_file, tmp_path, has_spin=False)
    triplet_mwfn = _new_mwfn(triplet_file, tmp_path, has_spin=True)

    singlet_bond_orders = singlet_mwfn.get_bond_order(model=model)
    triplet_bond_orders = triplet_mwfn.get_bond_order(model=model)

    singlet_co_bond_order = _get_bond_order(singlet_bond_orders, *DOUBLE_BOND_PAIR)
    triplet_co_bond_order = _get_bond_order(triplet_bond_orders, *DOUBLE_BOND_PAIR)

    assert triplet_co_bond_order < singlet_co_bond_order, (
        f"{model}: expected triplet C-O bond order to be lower than singlet for "
        f"{triplet_name} vs {singlet_name}, but got "
        f"{triplet_co_bond_order:.3f} >= {singlet_co_bond_order:.3f}"
    )


def _assert_atomic_values_match_csv_reference(
    values: dict[int, float],
    *,
    file_path: Path,
    csv_path: Path,
    decimal: int = 3,
) -> None:
    """Assert atomic values match a CSV table column keyed by Molden filename."""
    atom_indices = sorted(values)
    expected_atom_indices, expected_values = _load_atomic_reference_column(
        csv_path, file_path.name
    )

    assert atom_indices == expected_atom_indices, (
        f"{file_path.name}: atom indices {atom_indices} do not match "
        f"{csv_path.name} atom_index column {expected_atom_indices}."
    )

    nan_mask = np.isnan(expected_values)
    if nan_mask.all():
        pytest.skip(
            f"Populate {csv_path} column {file_path.name} with reference values."
        )

    if nan_mask.any():
        missing_atoms = [
            str(atom_idx)
            for atom_idx, is_missing in zip(expected_atom_indices, nan_mask)
            if is_missing
        ]
        pytest.fail(
            f"{csv_path} column {file_path.name} has placeholder values for atom(s): "
            f"{', '.join(missing_atoms)}."
        )

    observed_values = np.array(
        [values[atom_idx] for atom_idx in atom_indices], dtype=float
    )
    assert_array_almost_equal(observed_values, expected_values, decimal=decimal)


def _load_atomic_reference_column(
    csv_path: Path, file_name: str
) -> tuple[list[int], np.ndarray]:
    """Load one Molden-file column from an atomic reference CSV."""
    ref_data = np.genfromtxt(
        csv_path,
        delimiter=",",
        names=True,
        dtype=float,
        filling_values=np.nan,
        autostrip=True,
        encoding="utf-8",
        deletechars="",
    )
    ref_data = np.atleast_1d(ref_data)

    if ref_data.dtype.names is None:
        raise AssertionError(f"{csv_path} is missing a header row.")
    if "atom_index" not in ref_data.dtype.names:
        raise AssertionError(f"{csv_path} must contain an 'atom_index' column.")
    if file_name not in ref_data.dtype.names:
        raise AssertionError(
            f"{csv_path} is missing a column for Molden file {file_name!r}."
        )

    atom_indices = [int(atom_index) for atom_index in ref_data["atom_index"]]
    values = np.asarray(ref_data[file_name], dtype=float)

    expected_atom_indices = list(range(1, len(atom_indices) + 1))
    assert (
        atom_indices == expected_atom_indices
    ), f"{csv_path} atom_index column must list consecutive 1-based atom indices."

    return atom_indices, values


def _assert_atomic_property_matches_csv(
    file_path: Path,
    tmp_path: Path,
    property_name: str,
    csv_path: Path,
) -> None:
    """Compute an atomic property and compare it to a CSV reference column."""
    mwfn = _new_mwfn(file_path, tmp_path)
    values = _get_atomic_property(mwfn, property_name)
    _assert_atomic_values_match_csv_reference(
        values, file_path=file_path, csv_path=csv_path
    )


def _get_atomic_property(mwfn: Multiwfn, property_name: str) -> dict[int, float]:
    """Return scalar atomic values for supported CSV-backed test properties."""
    if property_name == "charges_hirshfeld":
        return mwfn.get_charges(model="hirshfeld")
    return mwfn.get_descriptor(property_name)


def _new_mwfn(
    file_path: Path, tmp_path: Path, *, has_spin: bool | None = None
) -> Multiwfn:
    """Build a Multiwfn instance with the standard per-test run directory layout."""
    if has_spin is None:
        has_spin = _is_triplet_file(file_path)
    return Multiwfn(file_path, run_path=tmp_path / file_path.stem, has_spin=has_spin)


def _is_triplet_file(file_path: Path) -> bool:
    """Infer spin state from the example Molden filename convention."""
    return "triplet" in file_path.stem.lower()


def _get_bond_order(
    bond_orders: dict[tuple[int, int], float], atom_i: int, atom_j: int
) -> float:
    """Return bond order value for a pair, independent of index order."""
    if (atom_i, atom_j) in bond_orders:
        return bond_orders[(atom_i, atom_j)]
    if (atom_j, atom_i) in bond_orders:
        return bond_orders[(atom_j, atom_i)]

    raise AssertionError(f"Bond order for atoms {atom_i}-{atom_j} was not reported.")
