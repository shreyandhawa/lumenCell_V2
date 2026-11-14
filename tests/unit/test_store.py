import os
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest

from lumencell.api import binary
from lumencell.clonal.isolation import ISOLATION_VERSION
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from lumencell.engine.store import STORE_VERSION, RunStore, run_id
from tests.conftest import Prepared


@pytest.fixture(scope="module")
def result(bifurcation: Prepared) -> TrajectoryResult:
    return run_method(bifurcation.ctx, "paga_dpt", ClusterRoot("true_label", "trunk_early"))


def test_round_trip_preserves_everything(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    back = store.load("abc")
    back.validate()
    for field in (
        "method",
        "topology_kind",
        "milestones",
        "segments",
        "lineages",
        "divergence",
        "root",
        "params",
    ):
        assert getattr(back, field) == getattr(result, field), field
    np.testing.assert_array_equal(back.pseudotime, result.pseudotime)
    np.testing.assert_array_equal(back.membership, result.membership)
    np.testing.assert_array_equal(back.segment, result.segment)
    np.testing.assert_array_equal(back.progress, result.progress)
    assert (back.milestone_weights != result.milestone_weights).nnz == 0
    assert [asdict(w) for w in back.warnings] == [asdict(w) for w in result.warnings]


def test_cells_file_is_the_wire_format(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    decoded = binary.decode(store.cells_bytes("abc"))
    assert set(decoded.arrays) == {"pseudotime", "membership", "primary_lineage", "segment", "progress"}
    assert decoded.arrays["membership"].shape == (len(result.lineages), result.n_obs)
    np.testing.assert_array_equal(decoded.arrays["primary_lineage"], result.primary_lineage)
    assert decoded.meta == {"n_obs": result.n_obs, "n_lineages": len(result.lineages)}


def test_run_id_is_deterministic_and_input_sensitive() -> None:
    base = run_id("paga_dpt", {"n_dcs": 10}, "cluster:a=b", None, "k1")
    assert base == run_id("paga_dpt", {"n_dcs": 10}, "cluster:a=b", None, "k1")
    assert base == run_id("paga_dpt", {"n_dcs": 10}, "cluster:a=b", [], "k1")
    assert run_id("paga_dpt", {"n_dcs": 10}, "cluster:a=b", ["X", "Y"], "k1") == run_id(
        "paga_dpt", {"n_dcs": 10}, "cluster:a=b", ["Y", "X"], "k1"
    )
    for other in (
        run_id("paga_dpt", {"n_dcs": 8}, "cluster:a=b", None, "k1"),
        run_id("paga_dpt", {"n_dcs": 10}, "cluster:a=c", None, "k1"),
        run_id("paga_dpt", {"n_dcs": 10}, "cluster:a=b", None, "k2"),
        run_id("palantir", {"n_dcs": 10}, "cluster:a=b", None, "k1"),
    ):
        assert other != base


def test_listing_and_existence(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    assert store.list() == [] and not store.exists("abc")
    store.save("abc", result)
    assert store.exists("abc")
    [summary] = store.list()
    assert summary["run_id"] == "abc"
    assert summary["n_lineages"] == len(result.lineages)


def test_the_store_records_who_named_the_terminals(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    store.save("auto", result)
    store.save("named", result, terminals=["a", "b"])
    assert store.graph("auto")["terminals_source"] == "auto"
    assert store.graph("named")["terminals_source"] == "user"


def test_listing_skips_runs_it_cannot_read(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    for rid in ("good", "damaged", "other_version"):
        store.save(rid, result)
    (store.path("damaged") / "result.json").write_text("{not json", encoding="utf-8")
    graph = store.graph("other_version")
    graph["store_version"] = 0
    store.update_graph("other_version", graph)
    assert [r["run_id"] for r in store.list()] == ["good"]


def test_failed_save_leaves_nothing(
    tmp_path: Path, result: TrajectoryResult, monkeypatch: pytest.MonkeyPatch
) -> None:
    store = RunStore(tmp_path)

    def boom(*_a: object, **_k: object) -> None:
        raise OSError("disk full")

    monkeypatch.setattr("lumencell.engine.store.sparse.save_npz", boom)
    with pytest.raises(OSError):
        store.save("abc", result)
    assert not store.exists("abc")
    assert list(store.dir.iterdir()) == []


def test_result_json_records_isolation(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    assert result.obs_keys, "run_method records what steered the run"
    store.save("abc", result, prep_key="prep123")
    graph = store.graph("abc")
    assert graph["store_version"] == STORE_VERSION == 2
    assert graph["prep_key"] == "prep123"
    assert graph["inference_isolation"] == ISOLATION_VERSION == 1
    assert graph["obs_keys"] == result.obs_keys
    assert store.load("abc").obs_keys == result.obs_keys
    store.save("unkeyed", result)
    assert store.graph("unkeyed")["prep_key"] is None


def test_a_run_from_another_store_version_is_absent(
    tmp_path: Path, result: TrajectoryResult, caplog: pytest.LogCaptureFixture
) -> None:
    store = RunStore(tmp_path)
    store.save("old", result)
    store.save("new", result)
    old = store.graph("old")
    old["store_version"] = 1
    del old["obs_keys"]
    store.update_graph("old", old)
    assert store.exists("new") and not store.exists("old"), "so a request recomputes it"
    with pytest.raises(FileNotFoundError, match="no run 'old'"):
        store.load("old")
    with caplog.at_level("WARNING", logger="lumencell.engine.store"):
        assert [r["run_id"] for r in store.list()] == ["new"]
    [record] = caplog.records
    assert "1 saved run(s)" in record.getMessage() and "another store version" in record.getMessage()


TRUST = "clonal-0123456789abcdef"


def test_trust_round_trips_atomically(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    assert store.read_trust("abc", TRUST) is None
    payload = {"score": 0.25, "nested": {"n": [1, 2]}}
    path = store.write_trust("abc", TRUST, payload)
    assert path == store.trust_path("abc", TRUST) == store.path("abc") / "trust" / f"{TRUST}.json"
    assert store.read_trust("abc", TRUST) == payload
    store.write_trust("abc", TRUST, {"score": 0.5})
    assert store.read_trust("abc", TRUST) == {"score": 0.5}
    assert [p.name for p in (store.path("abc") / "trust").iterdir()] == [f"{TRUST}.json"]


@pytest.mark.parametrize(
    "name", ["../x", "clonal-xyz", "Clonal-0123456789abcdef", f"{TRUST}\n", f"{TRUST}/x"]
)
def test_trust_names_are_checked(tmp_path: Path, result: TrajectoryResult, name: str) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    with pytest.raises(ValueError):
        store.trust_path("abc", name)
    with pytest.raises(ValueError):
        store.write_trust("abc", name, {})
    with pytest.raises(ValueError):
        store.read_trust("abc", name)
    assert not (store.path("abc") / "trust").exists()


def test_trust_refuses_nan_and_leaves_nothing(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    store.write_trust("abc", TRUST, {"v": 1})
    with pytest.raises(ValueError):
        store.write_trust("abc", TRUST, {"v": float("nan")})
    assert store.read_trust("abc", TRUST) == {"v": 1}, "a failed write leaves the old file"
    assert [p.name for p in (store.path("abc") / "trust").iterdir()] == [f"{TRUST}.json"]
    other = "clonal-fedcba9876543210"
    with pytest.raises(ValueError):
        store.write_trust("abc", other, {"v": float("inf")})
    assert store.read_trust("abc", other) is None
    assert not any(p.name.endswith(".tmp") for p in (store.path("abc") / "trust").iterdir())


def test_trust_needs_a_complete_run(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    with pytest.raises(FileNotFoundError):
        store.write_trust("nope", TRUST, {})
    assert not store.path("nope").exists(), "nothing is created for a missing run"
    assert store.read_trust("nope", TRUST) is None
    store.save("abc", result)
    (store.path("abc") / ".complete").unlink()
    with pytest.raises(FileNotFoundError):
        store.write_trust("abc", TRUST, {})
    assert not (store.path("abc") / "trust").exists()


def test_saving_a_run_again_keeps_its_trust(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    store.write_trust("abc", TRUST, {"score": 0.25})
    before = store.trust_path("abc", TRUST).read_text(encoding="utf-8")
    old_dir = store.path("abc")
    store.save("abc", result)
    assert store.trust_path("abc", TRUST).read_text(encoding="utf-8") == before
    assert store.read_trust("abc", TRUST) == {"score": 0.25}
    assert store.exists("abc")
    assert [p.name for p in store.dir.iterdir()] == ["abc"], "no temporary run folder is left"
    assert [p.name for p in old_dir.joinpath("trust").iterdir()] == [f"{TRUST}.json"]


def test_a_failed_save_keeps_the_old_runs_trust(
    tmp_path: Path, result: TrajectoryResult, monkeypatch: pytest.MonkeyPatch
) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    store.write_trust("abc", TRUST, {"score": 1})

    real, failed = Path.rename, []

    def boom(self: Path, target: Path) -> Path:
        if self.name.startswith(".abc-") and target.name == "abc" and not failed:
            failed.append(self)  # the swap of the new folder into place fails once
            raise RuntimeError("disk full")
        return real(self, target)

    with monkeypatch.context() as patch:
        patch.setattr(Path, "rename", boom)
        with pytest.raises(RuntimeError):
            store.save("abc", result)
    assert failed, "the swap was attempted and failed"
    assert store.read_trust("abc", TRUST) == {"score": 1}
    assert store.exists("abc")
    assert [p.name for p in store.dir.iterdir()] == ["abc"]


def test_a_damaged_trust_file_reads_as_missing(
    tmp_path: Path, result: TrajectoryResult, capsys: pytest.CaptureFixture[str]
) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    store.write_trust("abc", TRUST, {"v": 1})
    path = store.trust_path("abc", TRUST)
    path.write_text("{not json", encoding="utf-8")
    assert store.read_trust("abc", TRUST) is None
    err = capsys.readouterr().err
    assert err.count("\n") == 1 and path.name in err
    path.write_text("[1, 2]", encoding="utf-8")
    assert store.read_trust("abc", TRUST) is None
    assert path.name in capsys.readouterr().err
    path.write_bytes(b"\xff\xfe")
    assert store.read_trust("abc", TRUST) is None
    assert path.name in capsys.readouterr().err


@pytest.mark.parametrize(
    "name",
    [
        "a-0123456789abcdef",
        "a_b-0123456789abcdef",
        "a" + "b" * 31 + "-0123456789abcdef",
    ],
)
def test_trust_names_accept_the_documented_shape(tmp_path: Path, name: str) -> None:
    assert RunStore(tmp_path).trust_path("abc", name).name == f"{name}.json"


@pytest.mark.parametrize(
    "name",
    [
        "a" + "b" * 32 + "-0123456789abcdef",  # kind of 33 characters
        "_x-0123456789abcdef",  # kind starts with an underscore
        "1x-0123456789abcdef",  # kind starts with a digit
        "a1-0123456789abcdef",  # digit inside the kind
        "-0123456789abcdef",  # empty kind
        "clonal-0123456789abcde",  # 15 hex digits
        "clonal-0123456789abcdef0",  # 17 hex digits
        "clonal-0123456789ABCDEF",  # upper-case hex
        "clonal_0123456789abcdef",  # wrong separator
        "clonal-0123456789abcdef.json",
    ],
)
def test_trust_names_reject_near_misses(tmp_path: Path, name: str) -> None:
    with pytest.raises(ValueError):
        RunStore(tmp_path).trust_path("abc", name)


def test_a_failed_trust_write_leaves_no_temp_file(
    tmp_path: Path, result: TrajectoryResult, monkeypatch: pytest.MonkeyPatch
) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    store.write_trust("abc", TRUST, {"v": 1})

    def killed(*_a: object, **_k: object) -> None:
        raise OSError("killed")

    with monkeypatch.context() as patch:
        patch.setattr(os, "replace", killed)
        with pytest.raises(OSError):
            store.write_trust("abc", TRUST, {"v": 2})
    assert store.read_trust("abc", TRUST) == {"v": 1}
    assert [p.name for p in (store.path("abc") / "trust").iterdir()] == [f"{TRUST}.json"]


def test_a_trust_path_that_cannot_be_read_reads_as_missing(
    tmp_path: Path, result: TrajectoryResult, capsys: pytest.CaptureFixture[str]
) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    path = store.trust_path("abc", TRUST)
    path.parent.mkdir()
    path.mkdir()  # reading a folder raises IsADirectoryError, an OSError that is not "missing"
    assert store.read_trust("abc", TRUST) is None
    assert path.name in capsys.readouterr().err


def test_a_missing_trust_file_is_silent(
    tmp_path: Path, result: TrajectoryResult, capsys: pytest.CaptureFixture[str]
) -> None:
    store = RunStore(tmp_path)
    store.save("abc", result)
    assert store.read_trust("abc", TRUST) is None
    assert store.read_trust("nope", TRUST) is None
    assert capsys.readouterr().err == ""


@pytest.mark.parametrize("rid", ["../x", "a/b", "..", "", "a" * 65, "abc\n"])
def test_trust_paths_refuse_a_run_id_that_is_not_a_plain_name(tmp_path: Path, rid: str) -> None:
    store = RunStore(tmp_path)
    with pytest.raises(ValueError):
        store.trust_path(rid, TRUST)
    with pytest.raises(ValueError):
        store.write_trust(rid, TRUST, {})
    with pytest.raises(ValueError):
        store.read_trust(rid, TRUST)
