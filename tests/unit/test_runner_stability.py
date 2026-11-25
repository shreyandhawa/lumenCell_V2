import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, get_args

import pytest

from lumencell.bench.runner import CaseResult, StabilityRow, stability_row
from lumencell.bench.synth import generate, truth_from_adata
from lumencell.data.loader import load
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.tree_support import Level
from lumencell.prep.pipeline import run_pipeline
from lumencell.prep.stages import PrepParams


def _run(tmp_path: Path) -> Any:
    adata = generate("bifurcation", "realistic", n_cells=600, n_genes=200, seed=0)
    path = tmp_path / "b.h5ad"
    adata.write_h5ad(path)
    ds = load(path)
    pp = PrepParams()
    prep = run_pipeline(ds, pp)
    ctx = InferenceContext.from_prep(ds, prep, pp)
    root = ClusterRoot("true_label", truth_from_adata(adata).meta["root_label"])
    return ds, prep, pp, run_method(ctx, "paga_dpt", root)


def test_rows_always_carry_the_stability_fields() -> None:
    from lumencell.bench.metrics import Scores

    fields = {f: 0.0 for f in Scores.__dataclass_fields__}
    fields.update(n_spurious=0)
    row = CaseResult("c", "m", Scores(**fields), 3, 1.0).row()
    assert row["stability_status"] == "not_measured"
    assert math.isnan(row["stability_lineage"])
    assert {"stability_pt", "tree_likeness", "tree_support"} <= set(row)
    assert row["tree_support"] == "not_measured"


def test_stability_row_on_a_synthetic_run(tmp_path: Path) -> None:
    ds, prep, pp, result = _run(tmp_path)
    row = stability_row(ds, prep, pp, result)
    assert row.status == "ok"
    assert row.verdict in get_args(Level)
    assert 0 <= row.lineage <= 1
    assert 0 <= row.pseudotime <= 1


def test_clonal_columns_reach_every_replicate(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    ds, prep, pp, result = _run(tmp_path)
    seen: list[Any] = []

    def spy(*args: Any, **kwargs: Any) -> Any:
        seen.append(args[4])
        raise RuntimeError("stop here")

    monkeypatch.setattr("lumencell.bench.runner.run_stability", spy)
    stability_row(ds, prep, pp, result, clonal_columns=("clone_id", "donor"))
    assert seen[0].clonal_columns == ("clone_id", "donor")


def test_a_stability_error_is_recorded_as_failed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    ds, prep, pp, result = _run(tmp_path)

    def boom(*args: Any, **kwargs: Any) -> Any:
        raise RuntimeError("no replicates today")

    monkeypatch.setattr("lumencell.bench.runner.run_stability", boom)
    row = stability_row(ds, prep, pp, result)
    assert row.status == "failed"
    assert row.reason == "RuntimeError: no replicates today"
    assert "no replicates today" in capsys.readouterr().err
    assert repr(StabilityRow("failed", reason=row.reason)) == repr(row)


def test_a_stability_failure_reason_is_one_escaped_line_of_bounded_length(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    ds, prep, pp, result = _run(tmp_path)

    def boom(*args: Any, **kwargs: Any) -> Any:
        raise ValueError("bad column 'a\x1b[31m\nb' " + "x" * 1000)

    monkeypatch.setattr("lumencell.bench.runner.run_stability", boom)
    reason = stability_row(ds, prep, pp, result).reason
    assert reason is not None and reason.startswith("ValueError: ")
    assert reason.isprintable() and "\n" not in reason and "\x1b" not in reason
    assert len(reason) <= 300


class _FakeKind:
    def __init__(self, n: int) -> None:
        self.rows = [object()] * n


def test_stability_row_is_wired_to_the_run(monkeypatch: pytest.MonkeyPatch) -> None:
    import numpy as np

    from lumencell.bench import runner

    stab = SimpleNamespace(
        lineage_stability=0.11,
        pseudotime_stability=0.22,
        n_failed=2,
        reseeds=_FakeKind(3),
        subsamples=_FakeKind(4),
    )
    seen: dict[str, Any] = {}

    def fake_run_stability(data: Any, prep: Any, params: Any, result: Any, spec: Any) -> Any:
        seen["spec"] = spec
        return stab

    def fake_tree(conn: Any, tree: Any) -> float:
        seen["tree_args"] = (conn, tree)
        return 0.77

    def fake_cov(pt: Any) -> float:
        seen["cov_arg"] = pt
        return 0.99

    def fake_verdict(cov: Any, tree: Any, evidence: Any, clonal: Any, thresholds: Any) -> Any:
        seen["verdict_args"] = (cov, tree, evidence)
        return SimpleNamespace(level="weak")

    monkeypatch.setattr(runner, "run_stability", fake_run_stability)
    monkeypatch.setattr(runner, "tree_likeness", fake_tree)
    monkeypatch.setattr(runner, "coverage", fake_cov)
    monkeypatch.setattr(runner, "verdict", fake_verdict)
    pt = np.array([0.0, 0.5, 1.0], dtype=np.float32)
    root = SimpleNamespace(spec="cluster:true_label=trunk_early")
    result = SimpleNamespace(method="paga_dpt", params={}, root=root, pseudotime=pt, lineages=[1, 2])
    prep = SimpleNamespace(artifacts={"paga_connectivities": "CONN", "paga_tree": "TREE"})
    row = stability_row(object(), prep, object(), result)  # type: ignore[arg-type]
    cov, tree, ev = seen["verdict_args"]
    assert (ev.n_replicates, ev.n_failed, ev.lineage, ev.pseudotime) == (7, 2, 0.11, 0.22)
    assert (cov, tree) == (0.99, 0.77)
    assert seen["tree_args"] == ("CONN", "TREE")
    assert np.array_equal(seen["cov_arg"], pt) and seen["cov_arg"].dtype == np.float64
    assert seen["spec"].root == "cluster:true_label=trunk_early"
    assert (row.status, row.lineage, row.pseudotime, row.tree_likeness, row.verdict) == (
        "ok",
        0.11,
        0.22,
        0.77,
        "weak",
    )


def test_row_maps_a_measured_stability_row() -> None:
    from lumencell.bench.metrics import Scores

    fields = {f: 0.0 for f in Scores.__dataclass_fields__}
    fields.update(n_spurious=0)
    st = StabilityRow("ok", 0.123456, 0.654321, 0.333333, "weak")  # type: ignore[arg-type]
    row = CaseResult("c", "m", Scores(**fields), 3, 1.0, stability=st).row()
    assert row["stability_status"] == "ok"
    assert row["stability_lineage"] == 0.1235
    assert row["stability_pt"] == 0.6543
    assert row["tree_likeness"] == 0.3333
    assert row["tree_support"] == "weak"
