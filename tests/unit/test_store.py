from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest

from lumencell.api import binary
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from lumencell.engine.store import RunStore, run_id
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
