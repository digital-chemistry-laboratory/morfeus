"""Chemistry-focused integration tests for Multiwfn."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

import morfeus.multiwfn as multiwfn
from morfeus.multiwfn import CommandStep, Multiwfn, MultiwfnRunResult

MULTIWFN_DATA_DIR = Path(__file__).parent / "data" / "multiwfn"
EXAMPLE_MOLDEN_DIR = MULTIWFN_DATA_DIR / "example_molden"
SINGLET_MOLDEN_FILES = sorted(EXAMPLE_MOLDEN_DIR.glob("*singlet*.molden"))
TRIPLET_MOLDEN_FILES = sorted(EXAMPLE_MOLDEN_DIR.glob("*triplet*.molden"))
ALL_MOLDEN_FILES = SINGLET_MOLDEN_FILES + TRIPLET_MOLDEN_FILES

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
@pytest.mark.parametrize("persistent", [False, True])
@pytest.mark.parametrize("file_path", ALL_MOLDEN_FILES, ids=lambda path: path.name)
def test_descriptor_workflow(file_path, persistent):
    """Compare a descriptor sequence with references and expected chemistry."""
    triplet = "triplet" in file_path.stem
    with Multiwfn(file_path, persistent=persistent, has_spin=triplet) as mwfn:
        for name in ("spin_density", "esp_nuclear", "charges_hirshfeld"):
            values = (
                mwfn.get_charges("hirshfeld")
                if name == "charges_hirshfeld"
                else mwfn.get_descriptor(name)
            )
            reference = np.genfromtxt(
                MULTIWFN_DATA_DIR / f"multiwfn_{name}.csv",
                delimiter=",",
                names=True,
                deletechars="",
            )
            indices = sorted(values)
            expected = reference[file_path.name]
            assert indices == reference["atom_index"].tolist()
            assert np.isfinite(expected).all(), (file_path.name, name)
            assert_array_almost_equal([values[i] for i in indices], expected, decimal=3)
            if name == "charges_hirshfeld" and not triplet:
                assert min(values, key=values.__getitem__) == OXYGEN_INDEX
                assert all(values[i] > 0 for i in HYDROGEN_INDICES)
        rho = mwfn.get_descriptor("rho")
        assert max(rho, key=rho.__getitem__) == OXYGEN_INDEX
        moments = mwfn.get_electric_moments()
        assert abs(moments["dipole_magnitude_au"]) > 0
        assert any(abs(moments.get(key, 0)) > 0 for key in MULTIPOLE_KEYS)


@pytest.mark.multiwfn
@pytest.mark.parametrize("singlet_name,triplet_name", SINGLET_TRIPLET_EQUIVALENT_FILES)
def test_bond_order_workflow(singlet_name, triplet_name):
    """Check bond patterns and spin-state changes across partitioning models."""
    with (
        Multiwfn(
            EXAMPLE_MOLDEN_DIR / singlet_name, persistent=True, has_spin=False
        ) as singlet,
        Multiwfn(
            EXAMPLE_MOLDEN_DIR / triplet_name, persistent=True, has_spin=True
        ) as triplet,
    ):
        for model in BOND_ORDER_PARTITION_MODELS:
            ground = {
                tuple(sorted(pair)): value
                for pair, value in singlet.get_bond_order(model).items()
            }
            excited = {
                tuple(sorted(pair)): value
                for pair, value in triplet.get_bond_order(model).items()
            }
            assert (
                DOUBLE_BOND_RANGE[0] <= ground[DOUBLE_BOND_PAIR] <= DOUBLE_BOND_RANGE[1]
            ), model
            for pair in SINGLE_BOND_PAIRS:
                assert SINGLE_BOND_RANGE[0] <= ground[pair] <= SINGLE_BOND_RANGE[1], (
                    model,
                    pair,
                )
            assert excited[DOUBLE_BOND_PAIR] < ground[DOUBLE_BOND_PAIR], model


def test_grid_cache_distinguishes_user_functions(tmp_path, monkeypatch):
    """Cache repeated requests while keeping different user functions independent."""
    grid = tmp_path / "density.cub"
    grid.touch()
    mwfn = Multiwfn(grid, run_path=tmp_path)
    calls = []

    def run(commands, subdir=None):
        calls.append(commands)
        value = len(calls)
        return MultiwfnRunResult(
            f"Atomic space Value\n 1(C ) {value}.0\nSumming up", tmp_path
        )

    monkeypatch.setattr(mwfn, "_run_commands", run)
    assert mwfn.grid_to_descriptors(None, userfunction=1) == {1: 1.0}
    assert mwfn.grid_to_descriptors(grid, userfunction=2) == {1: 2.0}
    assert mwfn.grid_to_descriptors(None, userfunction=1) == {1: 1.0}
    assert mwfn.grid_to_descriptors(grid, userfunction=2) == {1: 2.0}
    assert len(calls) == 2


def test_loading_settings_invalidates_cached_results(tmp_path):
    """New settings invalidate computed values without losing collected citations."""
    grid = tmp_path / "density.cub"
    grid.touch()
    settings = tmp_path / "settings.ini"
    settings.write_text("iuserfunc=1\n")
    mwfn = Multiwfn(grid, run_path=tmp_path / "run")
    mwfn._results.charges = {"hirshfeld": {1: 0.1}}
    mwfn._results.grid_descriptors = {(str(grid), 1): {1: 2.0}}
    citations = mwfn.get_citations()
    mwfn.load_settingini(settings)
    expected = type(mwfn._results)(citations=set(citations))
    assert mwfn._results == expected
    assert mwfn.get_citations() == citations


def test_rdg_is_available_for_grids_only(tmp_path, monkeypatch):
    """Export RDG without overwriting files or reusing an earlier run's cube."""
    grid = tmp_path / "input.cub"
    grid.touch()
    sent: list[object] = []
    with Multiwfn(grid) as mwfn:

        def run(commands, subdir=None):
            workdir = mwfn._setup_workdir(subdir)
            if not sent:
                (workdir / "RDG.cub").write_text("generated grid")
            sent.extend(commands)
            return MultiwfnRunResult("", workdir)

        monkeypatch.setattr(mwfn, "_run_commands", run)
        output = tmp_path / "RDG.cub"
        assert mwfn.get_grid("rdg", "low", output) == output
        with pytest.raises(FileExistsError):
            mwfn.get_grid("rdg", "low", output)
        with pytest.raises(FileNotFoundError):
            mwfn.get_grid("rdg", "low", tmp_path / "missing.cub")
        assert CommandStep("13", expect="13 Reduced density gradient") in sent
        with pytest.raises(ValueError, match="requires grid data"):
            mwfn.get_descriptor("rdg")
    assert output.read_text() == "generated grid"
    assert not mwfn._run_path.exists()


@pytest.mark.multiwfn
def test_command_failure_closes_child(tmp_path, monkeypatch):
    """Command and initialization failures release the live process."""
    import morfeus.multiwfn as module

    sessions = []
    original = module._PexpectSession

    def session(*args, **kwargs):
        child = original(*args, **kwargs)
        sessions.append(child)
        return child

    monkeypatch.setattr(module, "_PexpectSession", session)
    mwfn = Multiwfn(
        EXAMPLE_MOLDEN_DIR / "xtb_singlet.molden",
        run_path=tmp_path,
        has_spin=False,
        expect_timeout=0.2,
    )
    with pytest.raises(RuntimeError, match="Expected pattern"):
        mwfn.run_commands([CommandStep("300", expect="THIS_PROMPT_DOES_NOT_EXIST")])
    assert len(sessions) == 1
    assert not sessions[0]._child.isalive()
    assert sessions[0]._child.closed
    xyz = MULTIWFN_DATA_DIR.parent / "xtb" / "1-penten-3-one.xyz"
    with pytest.raises(RuntimeError, match="multiplicity"):
        Multiwfn(xyz, run_path=tmp_path, persistent=True)
    assert sessions[-1]._child.closed


@pytest.mark.multiwfn
@pytest.mark.parametrize("robust", [True, False])
def test_persistent_reuse(tmp_path, monkeypatch, robust):
    """Load once, isolate analysis output, and quit once for a descriptor sequence."""
    path = EXAMPLE_MOLDEN_DIR / "xtb_singlet.molden"
    with Multiwfn(path) as fresh:
        charges = fresh.get_charges("hirshfeld")
        moments = fresh.get_electric_moments()
    sent = []
    send = multiwfn._PexpectSession.send

    def record(session, command):
        sent.append(command)
        send(session, command)

    monkeypatch.setattr(multiwfn._PexpectSession, "send", record)
    with Multiwfn(path, run_path=tmp_path, persistent=True, robust_mode=robust) as m:
        session = m._session
        assert session is not None
        assert m.get_charges("hirshfeld") == charges
        assert m.get_electric_moments() == moments
        assert m._session is session and session._child.isalive()
        assert list((m._run_path / "hirshfeld").glob("run-*/*.chg"))
        with Multiwfn(path, run_path=tmp_path, has_spin=False) as other:
            assert m._run_path != other._run_path
        assert "Atomic space" not in session.stdout
        assert sent.count(str(path.resolve())) == 1
        assert "r" not in sent and "q" not in sent
    assert sent.count("q") == 1
    assert session._child.closed and session._child.exitstatus == 0
    assert m._run_path.is_dir() and other._run_path.is_dir()


@pytest.mark.multiwfn
def test_persistent_input_and_settings(tmp_path):
    """Reload for a cube and back, and restart after settings change."""
    path = EXAMPLE_MOLDEN_DIR / "xtb_singlet.molden"
    cube = MULTIWFN_DATA_DIR / "example_xtb_cub" / "spindensity.cub"
    with Multiwfn(path) as fresh:
        density = fresh.get_descriptor("rho")
        integrated = fresh.grid_to_descriptors(cube)
        moments = fresh.get_electric_moments()
    with Multiwfn(path, run_path=tmp_path / "run", persistent=True) as m:
        session = m._session
        assert session is not None
        m.get_vector("coordinates")
        assert m.get_descriptor("rho") == density
        assert m.grid_to_descriptors(cube) == integrated
        assert m.get_electric_moments() == moments
        assert m._session is session
        settings = tmp_path / "run" / "settings.ini"
        settings.write_text("nthreads=1\n")
        m.load_settingini(settings)
        assert session._child.closed and m._results.electric_moments is None
        assert m.get_electric_moments() == moments
        assert m._session is not session


@pytest.mark.multiwfn
def test_persistent_custom_commands_and_failure(tmp_path):
    """Invalidate custom state; discard a failed session and restart on the next call."""
    path = EXAMPLE_MOLDEN_DIR / "xtb_singlet.molden"
    with Multiwfn(path, run_path=tmp_path, persistent=True) as m:
        moments = m.get_electric_moments()
        result = m.run_commands(["300", "5", "0"])
        assert m._parse_electric_moments(result.stdout) == moments
        assert m._results.electric_moments is None and m._loaded_file is None
        assert m.get_electric_moments() == moments
        session = m._session
        assert session is not None
        session._expect_timeout = 0.2
        with pytest.raises(RuntimeError, match="Expected pattern"):
            m.run_commands([CommandStep("300", expect="THIS_PROMPT_DOES_NOT_EXIST")])
        assert session._child.closed and m._session is None
        assert m.get_electric_moments() == moments


@pytest.mark.multiwfn
@pytest.mark.parametrize(
    "events,expected",
    [
        (
            [(0, "Progress: [#-] 50 %"), (1, "Progress: [##] 100 %")],
            "Progress: [#-] 50 %Progress: [##] 100 %",
        ),
        ([(0, "Progress: [#-] 50 %"), (1, "0 Return")], "Progress: [#-] 50 %0 Return"),
        ([(2, "TIMEOUT")], ""),
        (
            [(0, "Progress: [#-] 50 %"), (1, "TIMEOUT"), (31, "TIMEOUT")],
            "Progress: [#-] 50 %",
        ),
        ([(0, "EOF")], ""),
    ],
)
def test_progress_wait_stops_and_preserves_output(monkeypatch, events, expected):
    """Keep menu output and distinguish initial wait, idle timeout and completion."""
    pexpect = pytest.importorskip("pexpect")
    clock = [0.0]
    pending = iter(events)

    def read(**kwargs):
        clock[0], value = next(pending)
        if value == "TIMEOUT":
            raise pexpect.TIMEOUT("idle")
        if value == "EOF":
            raise pexpect.EOF("exited")
        return value

    monkeypatch.setattr(multiwfn, "time", SimpleNamespace(monotonic=lambda: clock[0]))
    session = object.__new__(multiwfn._PexpectSession)
    session._child = Mock(read_nonblocking=Mock(side_effect=read))
    session._debug = False
    session._transcript = []
    session.wait_for_progress(max_wait=30)
    assert session.stdout == expected
    assert session._child.read_nonblocking.call_count == len(events)


@pytest.mark.multiwfn
@pytest.mark.parametrize(
    "filename,count,reference",
    [
        ("pyscf_dft_singlet_rks.molden", 38, {"restricted": (23, 24, -2.1716, 0.2925)}),
        (
            "pyscf_dft_triplet_uks.molden",
            76,
            {"alpha": (24, 25, -0.2070, 3.6907), "beta": (60, 61, -3.9766, -3.0207)},
        ),
    ],
)
def test_orbitals_and_frontier_gaps(tmp_path, filename, count, reference):
    """Check records and gaps against reference Multiwfn energies in eV."""
    with Multiwfn(EXAMPLE_MOLDEN_DIR / filename, run_path=tmp_path) as mwfn:
        orbitals = mwfn.get_orbitals()
        gaps = mwfn.get_frontier_gaps()
    assert len(orbitals) == count
    assert set(gaps) == set(reference)
    for spin, (homo, lumo, occupied_energy, empty_energy) in reference.items():
        for index, occupation, energy in (
            (homo, 2 if spin == "restricted" else 1, occupied_energy),
            (lumo, 0, empty_energy),
        ):
            assert next(o for o in orbitals if o["index"] == index) == {
                "index": index,
                "spin": spin,
                "occupation": occupation,
                "energy_ev": pytest.approx(energy, abs=5e-5, rel=0),
            }
        assert gaps[spin] == {
            "homo_indices": [homo],
            "lumo_indices": [lumo],
            "homo_energy_ev": pytest.approx(occupied_energy, abs=5e-5, rel=0),
            "lumo_energy_ev": pytest.approx(empty_energy, abs=5e-5, rel=0),
            "gap_ev": pytest.approx(empty_energy - occupied_energy, abs=1e-4, rel=0),
        }


@pytest.mark.multiwfn
def test_fractional_orbitals(tmp_path):
    """Expose fractional occupations without assigning an automatic gap."""
    path = EXAMPLE_MOLDEN_DIR / "pyscf_casscf_singlet_rhf_natorb.molden"
    with Multiwfn(path, run_path=tmp_path) as mwfn:
        assert any(
            abs(o["occupation"] - round(o["occupation"])) > 1e-6
            for o in mwfn.get_orbitals()
        )
        with pytest.raises(ValueError, match="Fractional occupations"):
            mwfn.get_frontier_gaps()


def test_degenerate_frontiers(monkeypatch):
    """Select energy extrema independently of index order and degeneracy grouping."""
    mwfn = object.__new__(Multiwfn)
    records = [
        {
            "index": index,
            "spin": "restricted",
            "occupation": occupation,
            "energy_ev": energy,
        }
        for index, occupation, energy in [
            (4, 2, -5.000001),
            (3, 2, -5),
            (7, 0, 2),
            (6, 0, 2.000001),
        ]
    ]
    monkeypatch.setattr(mwfn, "get_orbitals", lambda: records)
    gap = mwfn.get_frontier_gaps()["restricted"]
    exact = mwfn.get_frontier_gaps(energy_tolerance_ev=0)["restricted"]
    assert gap["homo_indices"] == [3, 4] and gap["lumo_indices"] == [6, 7]
    assert exact["homo_indices"] == [3] and exact["lumo_indices"] == [7]
    assert gap["gap_ev"] == exact["gap_ev"] == 7


def test_signed_frontier_gap(monkeypatch):
    """Retain negative gaps for non-Aufbau occupations."""
    mwfn = object.__new__(Multiwfn)
    monkeypatch.setattr(
        mwfn,
        "get_orbitals",
        lambda: [
            {"index": 1, "spin": "restricted", "occupation": 2, "energy_ev": 1},
            {"index": 2, "spin": "restricted", "occupation": 0, "energy_ev": 0},
        ],
    )
    assert mwfn.get_frontier_gaps()["restricted"]["gap_ev"] == -1


def test_restricted_open_shell_gap(monkeypatch):
    """Require an explicit convention for restricted open-shell orbitals."""
    mwfn = object.__new__(Multiwfn)
    monkeypatch.setattr(
        mwfn,
        "get_orbitals",
        lambda: [
            {"index": i, "spin": "restricted", "occupation": occ, "energy_ev": i}
            for i, occ in enumerate([2, 1, 0], 1)
        ],
    )
    with pytest.raises(ValueError, match="Restricted open-shell"):
        mwfn.get_frontier_gaps()


@pytest.mark.multiwfn
@pytest.mark.parametrize("persistent", [False, True])
def test_grid_export_survives_close(tmp_path, persistent):
    """Export a retained cube from scratch and reuse it as input without collisions."""
    path = EXAMPLE_MOLDEN_DIR / "xtb_singlet.molden"
    output = tmp_path / "export" / "rho.cub"
    with Multiwfn(path, persistent=persistent) as m:
        moments = m.get_electric_moments()
        scratch = m._run_path
        assert m.get_grid("rho", "low", output_path=output) == output
        with pytest.raises(FileExistsError):
            m.get_grid("rho", "low", output_path=output)
        integrated = m.grid_to_descriptors(output)
        m._results.electric_moments = None
        assert m.get_electric_moments() == moments
    assert not scratch.exists() and output.is_file()
    with Multiwfn(output) as fresh:
        assert fresh.grid_to_descriptors(None) == integrated
