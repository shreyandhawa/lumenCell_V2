import math
from pathlib import Path
from typing import Any, get_args

import pytest
from typer.testing import CliRunner

from lumencell.bench import report
from lumencell.bench.runner import CaseResult, StabilityRow, stability_row
from lumencell.bench.synth import generate, truth_from_adata
from lumencell.cli import app
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
    row = stability_row(ds, prep, pp, result, None)
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
    stability_row(ds, prep, pp, result, None, clonal_columns=("clone_id", "donor"))
    assert seen[0].clonal_columns == ("clone_id", "donor")


def test_a_stability_error_is_recorded_as_failed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    ds, prep, pp, result = _run(tmp_path)

    def boom(*args: Any, **kwargs: Any) -> Any:
        raise RuntimeError("no replicates today")

    monkeypatch.setattr("lumencell.bench.runner.run_stability", boom)
    row = stability_row(ds, prep, pp, result, None)
    assert row.status == "failed"
    assert row.reason == "RuntimeError: no replicates today"
    assert "no replicates today" in capsys.readouterr().err
    assert repr(StabilityRow("failed", reason=row.reason)) == repr(row)


def test_bench_real_has_a_stability_switch() -> None:
    out = CliRunner().invoke(app, ["bench", "real", "--help"]).output
    assert "--no-stability" in out


def test_the_report_shows_the_verdict(monkeypatch: pytest.MonkeyPatch) -> None:
    _, real = report.latest_real() or (None, {"methods": {}})
    d, base = next(iter(real["methods"]["paga_dpt"].items()))
    row = {**base, "tree_support": "weak"}
    data = {"commit": "abc", "methods": {"paga_dpt": {d: row}}}
    monkeypatch.setattr(
        report, "latest_real", lambda: (report.REPO / "benchmarks" / "results" / "real-x.json", data)
    )
    assert "verdict weak" in report.render()
