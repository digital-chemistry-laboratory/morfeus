"""Chemistry-focused integration tests for Multiwfn."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
from unittest.mock import Mock

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest

from morfeus.io import CubeParser
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


@pytest.mark.multiwfn
@pytest.mark.parametrize("persistent", [False, True])
def test_interaction_workflow(tmp_path, persistent):
    """Coupled fields, contributions and cutoff changes preserve descriptor state."""
    box = {"origin": (-4, -4, -4), "spacing": (1, 1, 1), "shape": (9, 9, 9)}
    fragments = [[7, 8, 9], [4, 12, 13]]
    settings = tmp_path / "settings.ini"
    settings.write_text(" nthreads=1\n RDG_maxrho=0.005\n IGMvdwscl=2\n")
    with Multiwfn(
        EXAMPLE_MOLDEN_DIR / "pyscf_dft_singlet_rks.molden", persistent=persistent
    ) as m:
        m.load_settingini(settings)
        charges = m.get_charges("hirshfeld")
        pid = m._session._child.pid if m._session else None
        fields = {}
        for mode in ("inter", "intra", "both"):
            result = m.get_interaction(
                "igmh",
                tmp_path / mode,
                fragments=fragments,
                mode=mode,
                grid=box,
                contributions=("atom", "pair") if mode == "both" else (),
            )
            fields[mode] = {
                key: CubeParser(path).S for key, path in result["grids"].items()
            }
            assert (
                json.loads((tmp_path / mode / "metadata.json").read_text())
                == result["metadata"]
            )
            assert result["metadata"]["multiwfn_version"].startswith("Version ")
            assert result["metadata"]["grid"] == {
                "origin": [-4, -4, -4],
                "spacing": [1, 1, 1],
                "shape": [9, 9, 9],
            }
            assert all(np.isfinite(values).all() for values in fields[mode].values())
        assert box == {"origin": (-4, -4, -4), "spacing": (1, 1, 1), "shape": (9, 9, 9)}
        both = fields["both"]
        assert set(both) == {
            "sign_lambda2_rho",
            "delta_g_inter",
            "delta_g_intra",
            "delta_g_total",
        }
        for mode in ("inter", "intra"):
            assert set(fields[mode]) == {"sign_lambda2_rho", f"delta_g_{mode}"}
            np.testing.assert_allclose(
                fields[mode][f"delta_g_{mode}"], both[f"delta_g_{mode}"]
            )
        np.testing.assert_allclose(
            both["delta_g_total"],
            both["delta_g_inter"] + both["delta_g_intra"],
            atol=2e-7,
        )
        assert both["delta_g_inter"].max() > 0
        atomic, pairs = result["atomic_contributions"], result["pair_contributions"]
        assert atomic is not None and pairs is not None
        assert set(atomic) == {i for fragment in fragments for i in fragment}
        assert set(pairs) == {(i, j) for i in fragments[0] for j in fragments[1]}
        assert pairs[8, 4] == pytest.approx(0.016178, abs=2e-6)
        assert atomic[4] == pytest.approx(0.035466, abs=2e-6)
        for atom, value in atomic.items():
            assert value == pytest.approx(
                sum(v for ij, v in pairs.items() if atom in ij), abs=3e-6
            )
        for cutoff in (None, 0.02, 0.1):
            nci = m.get_interaction(
                "nci", tmp_path / f"nci-{cutoff}", grid=box, density_cutoff=cutoff
            )
            rdg = CubeParser(nci["grids"]["rdg"]).S
            sl2r = CubeParser(nci["grids"]["sign_lambda2_rho"]).S
            np.testing.assert_allclose(sl2r, both["sign_lambda2_rho"])
            assert nci["metadata"]["density_cutoff"] == cutoff
            if cutoff is None:
                raw = rdg
                assert np.any((abs(sl2r) > 0.1) & (raw < 100))
            else:
                np.testing.assert_allclose(rdg, np.where(abs(sl2r) >= cutoff, 100, raw))
        with (
            Multiwfn(nci["grids"]["rdg"]) as cube,
            pytest.raises(ValueError, match="wavefunction"),
        ):
            cube.get_interaction("nci", tmp_path / "unsupported", grid=box)
        repeat = m.get_interaction(
            "igmh", tmp_path / "repeat", fragments=fragments, grid=box
        )
        np.testing.assert_array_equal(
            CubeParser(repeat["grids"]["delta_g_inter"]).S, both["delta_g_inter"]
        )
        assert repeat["metadata"] == json.loads(
            (tmp_path / "inter/metadata.json").read_text()
        )
        with pytest.raises(FileExistsError):
            m.get_interaction("nci", tmp_path / "repeat", grid=box)
        with pytest.raises(ValueError, match="exceeds"):
            m.get_interaction(
                "igmh", tmp_path / "invalid", fragments=[[1], [15]], grid=box
            )
        assert m.get_charges("hirshfeld") == charges
        parent_grid = m.get_grid(
            "rdg", ["5", "-4,-4,-4", "1,1,1", "9,9,9"], tmp_path / "parent.cub"
        )
        np.testing.assert_allclose(
            CubeParser(parent_grid).S, np.where(abs(sl2r) >= 0.005, 100, raw)
        )
        if persistent:
            assert m._session is not None
            assert m._session._child.pid == pid
    assert all(path.is_file() for path in result["grids"].values())
    assert settings.read_text() == " nthreads=1\n RDG_maxrho=0.005\n IGMvdwscl=2\n"


@pytest.mark.multiwfn
def test_interaction_fragment_modes(tmp_path):
    """One fragment has no interfragment field; multiple fragments decompose total."""
    with Multiwfn(EXAMPLE_MOLDEN_DIR / "pyscf_dft_triplet_uks.molden") as m:
        one = m.get_interaction(
            "igmh", tmp_path / "one", fragments=[[1, 2, 7]], mode="intra", grid="low"
        )
        assert CubeParser(one["grids"]["delta_g_intra"]).S.max() > 0
        three = m.get_interaction(
            "igmh",
            tmp_path / "three",
            fragments=[[7], [4], [12]],
            mode="both",
            grid={"origin": (-4, -4, -4), "spacing": (1, 1, 1), "shape": (9, 9, 9)},
        )
        assert np.allclose(CubeParser(three["grids"]["delta_g_intra"]).S, 0, atol=1e-12)
        np.testing.assert_allclose(
            CubeParser(three["grids"]["delta_g_total"]).S,
            CubeParser(three["grids"]["delta_g_inter"]).S,
        )
        exported = m.run_commands(["100", "2", "4", "wavefunction.wfx", "0"])
        with Multiwfn(exported.workdir / "wavefunction.wfx") as wfx:
            assert wfx._has_spin
            result = wfx.get_interaction(
                "igmh",
                tmp_path / "wfx",
                fragments=three["fragments"],
                mode="both",
                grid=three["metadata"]["grid"],
            )
            for name in three["grids"]:
                np.testing.assert_allclose(
                    CubeParser(result["grids"][name]).S,
                    CubeParser(three["grids"][name]).S,
                    atol=2e-7,
                )


@pytest.mark.parametrize(
    "options",
    [
        {"method": "igm"},
        {"fragments": None},
        {"fragments": []},
        {"fragments": [[1], []]},
        {"fragments": [[1, 1], [2]]},
        {"fragments": [[1], [1, 2]]},
        {"fragments": [[0], [2]]},
        {"fragments": [[1.5], [2]]},
        {"fragments": [[True], [2]]},
        {"fragments": [[1]], "mode": "inter"},
        {"mode": "invalid"},
        {"grid": "invalid"},
        {"grid": {"origin": (0, 0, 0)}},
        {
            "grid": {
                "origin": (0, 0, 0),
                "spacing": (1, 1, 1),
                "shape": (2, 2, 2),
                "units": "bohr",
            }
        },
        {"contributions": ("energy",)},
        {"mode": "intra", "contributions": ("atom",)},
        {"fragments": [[1], [2], [3]], "contributions": ("pair",)},
        {"density_cutoff": 0.05},
        {"method": "nci", "fragments": []},
        {"method": "nci", "fragments": None, "mode": "both"},
        {"method": "nci", "fragments": None, "contributions": ("atom",)},
        *(
            {"method": "nci", "fragments": None, "density_cutoff": x}
            for x in (0, -1, float("nan"), float("inf"))
        ),
    ],
)
def test_interaction_validation(tmp_path, monkeypatch, options):
    """Reject incompatible chemical intent before starting a calculation."""
    path = tmp_path / "input.molden"
    path.touch()
    monkeypatch.setattr(
        Multiwfn, "_run_commands", Mock(side_effect=AssertionError("Must not run"))
    )
    args = {"method": "igmh", "fragments": [[1], [2]], **options}
    with Multiwfn(path, has_spin=False) as m, pytest.raises(ValueError):
        m.get_interaction(output_dir=tmp_path / "out", **args)
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(
    "field,value",
    [
        ("origin", (0, 0)),
        ("origin", (0, 0, float("nan"))),
        ("spacing", (1, 0, 1)),
        ("spacing", (1, 1, float("inf"))),
        ("shape", (2, 1, 2)),
        ("shape", (2, 2, 2.5)),
        ("shape", (True, 2, 2)),
    ],
)
def test_interaction_grid_validation(field, value):
    """Reject grids that cannot represent a finite three-dimensional volume."""
    with pytest.raises(ValueError):
        Multiwfn._validate_interaction_grid(
            {
                "origin": (0, 0, 0),
                "spacing": (1, 1, 1),
                "shape": (2, 2, 2),
                field: value,
            }
        )


@pytest.mark.parametrize("failure", ["missing", "geometry", "copy", "collision"])
def test_interaction_export_failures(tmp_path, monkeypatch, failure):
    """Fail on incomplete cubes and preserve existing destinations during export."""
    source = tmp_path / "generated"
    source.mkdir()
    header = "cube\ncomment\n1 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n"
    (source / "func1.cub").write_text(header)
    if failure != "missing":
        (source / "func2.cub").write_text(
            header.replace("2 1 0 0", "3 1 0 0") if failure == "geometry" else header
        )
    path = tmp_path / "input.molden"
    path.touch()
    destination = tmp_path / "out"

    def run(self, commands, subdir=None):
        if commands and failure == "collision":
            destination.mkdir()
            (destination / "keep").write_text("original")
        return MultiwfnRunResult(
            "Version 3.8, update date: 2026-Jan-7\nAtoms: 14, GTFs: 114", source
        )

    def copy(source, target):
        target.write_text("partial")
        raise OSError("Export failed")

    monkeypatch.setattr(Multiwfn, "_run_commands", run)
    if failure == "copy":
        monkeypatch.setattr(shutil, "copyfile", copy)
    error = {
        "missing": FileNotFoundError,
        "geometry": RuntimeError,
        "copy": OSError,
        "collision": FileExistsError,
    }[failure]
    with Multiwfn(path, has_spin=False) as m, pytest.raises(error):
        m.get_interaction("nci", destination)
    if failure == "collision":
        assert (destination / "keep").read_text() == "original"
    else:
        assert not destination.exists()


def test_interaction_contribution_totals():
    """Allow rounding and omitted zeros, but reject missing or inconsistent data."""
    table = (
        "Atom 7 : 0.100001\nAtom 2 : 0.000000\n"
        "Atom 8 : 0.050001\nAtom 3 : 0.050001\n"
        "Atomic pair delta-g indices\n"
        "7 8 : 0.050001\n7 3 : 0.050001\n"
        "Sum of all atomic pair delta-g indices: 0.100001\n"
    )
    fragments = [[7, 2], [8, 3]]
    parse = Multiwfn._parse_interaction_contributions
    atomic, pairs = parse(table, fragments)
    assert atomic[7] == 0.100001
    assert pairs == {(7, 8): 0.050001, (7, 3): 0.050001, (2, 8): 0, (2, 3): 0}
    for damaged in (
        table.replace("7 8 :", "7,8 :").replace("7 3 :", "7,3 :"),
        table.replace("7 8 : 0.050001\n", ""),
        table.replace("Atom 7 : 0.100001", "Atom 7 : 0.200001"),
        table.replace("indices: 0.100001", "indices: 0.200001"),
        table.split("Sum of all")[0],
    ):
        with pytest.raises(RuntimeError, match="IGMH contribution"):
            parse(damaged, fragments)
    zero = table.replace("7 8 : 0.050001\n7 3 : 0.050001\n", "")
    zero = zero.replace("0.100001", "0.000000").replace("0.050001", "0.000000")
    atomic, pairs = parse(zero, fragments)
    assert set(atomic.values()) == set(pairs.values()) == {0}


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
def test_wfn_charges():
    """Preserve Hirshfeld charges after WFN export and reject Mulliken analysis."""
    with Multiwfn(EXAMPLE_MOLDEN_DIR / "pyscf_dft_singlet_rks.molden") as molden:
        expected = molden.get_charges("hirshfeld")
        result = molden.run_commands(["100", "2", "5", "wavefunction.wfn", "0"])
        with Multiwfn(result.workdir / "wavefunction.wfn") as wfn:
            assert wfn.get_charges("hirshfeld") == pytest.approx(expected, abs=1e-6)
            with pytest.raises(ValueError, match="basis function information"):
                wfn.get_charges("mulliken")


@pytest.mark.multiwfn
@pytest.mark.parametrize("persistent", [False, True])
def test_reactivity_and_surface_workflow(persistent):
    """Check atomic tables and both surface formats against printed reference values."""
    path = EXAMPLE_MOLDEN_DIR / "pyscf_dft_singlet_rks.molden"
    with Multiwfn(path, persistent=persistent) as m:
        fukui = m.get_fukui()
        assert [
            fukui[key][1] for key in ("f_plus", "f_minus", "f_zero", "dd")
        ] == pytest.approx([0.00093, 0.00872, 0.00483, -0.00780], abs=1e-5)
        superdeloc = m.get_superdelocalizabilities()
        assert [
            superdeloc[key][1] for key in ("d_n", "d_e", "d_n_0", "d_e_0")
        ] == pytest.approx([-8.81097, -13.40253, -9.68954, -11.90341], abs=1e-5)
        for model, reference in (
            ("esp", {"area_total": 3.57765, "avg_all": 2.50293, "pi": 1.349486}),
            (
                "alie",
                {"area_total": 3.58434, "avg_all": 10.881895, "var_all": 0.147203},
            ),
        ):
            surface = m.get_surface(model)
            assert len(surface["atomic"]) == 12
            assert {
                key: surface["atomic"][3][key] for key in reference
            } == pytest.approx(reference, abs=1e-3)
            assert surface["global"]["overall_surface_area"] > 0


def test_two_column_superdelocalizabilities():
    """Support older headers and absent optional columns; stop before later tables."""
    stdout = "Atom index DN DE\n 1(Cl ) 1.2e-2 -3.4\nSum of DN\n 2(H ) 9 9"
    m = object.__new__(Multiwfn)
    assert m._parse_superdelocalizabilities(stdout) == {
        "d_n": {1: 0.012},
        "d_e": {1: -3.4},
        "d_n_0": {},
        "d_e_0": {},
    }


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
@pytest.mark.parametrize("ending", ["ready", "timeout", "exit"])
def test_calculation_wait_requires_completion(ending):
    """Send only after the final prompt; a stalled or dead process cannot continue."""
    pexpect = pytest.importorskip("pexpect")
    for robust in (False, True):
        session = object.__new__(multiwfn._PexpectSession)
        session._transcript, session._last_command_pos, session._debug = [], 0, False
        session._expect_timeout = 0.1
        session._child = Mock()
        session._child.read_nonblocking.side_effect = [
            "Progress: [####] 100%\nCalculating next stage...",
            "Progress: [####] 100%\n",
            "3 Output cube",
            {
                "ready": " files\n",
                "timeout": pexpect.TIMEOUT("idle"),
                "exit": pexpect.EOF("exit"),
            }[ending],
            pexpect.TIMEOUT("idle"),
        ]
        step = CommandStep("3", expect="3 Output cube files", max_wait=1)
        if ending == "ready":
            session._process_command(step, robust)
            session._child.sendline.assert_called_once_with("3")
            assert session.stdout.endswith("3 Output cube files\n")
        else:
            with pytest.raises(RuntimeError, match="Expected pattern"):
                session._process_command(step, robust)
            session._child.sendline.assert_not_called()


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


@pytest.mark.multiwfn
@pytest.mark.parametrize(
    "method,args",
    [
        ("get_localization", ()),
        ("get_bond_order", ("fuzzy",)),
        ("get_bond_order", ("laplacian",)),
    ],
)
def test_persistent_partitioning_does_not_leak(method, args):
    """Keep default analyses independent of a preceding Hirshfeld integration."""
    path = EXAMPLE_MOLDEN_DIR / "xtb_singlet.molden"
    with Multiwfn(path) as fresh:
        expected = getattr(fresh, method)(*args)
    with Multiwfn(path, persistent=True) as m:
        m.get_descriptor("rho", partitioning="hirshfeld")
        assert getattr(m, method)(*args) == expected
