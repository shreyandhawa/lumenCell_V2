import math
from pathlib import Path
from types import SimpleNamespace
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


def test_a_stability_failure_reason_is_one_escaped_line_of_bounded_length(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    ds, prep, pp, result = _run(tmp_path)

    def boom(*args: Any, **kwargs: Any) -> Any:
        raise ValueError("bad column 'a\x1b[31m\nb' " + "x" * 1000)

    monkeypatch.setattr("lumencell.bench.runner.run_stability", boom)
    reason = stability_row(ds, prep, pp, result, None).reason
    assert reason is not None and reason.startswith("ValueError: ")
    assert reason.isprintable() and "\n" not in reason and "\x1b" not in reason
    assert len(reason) <= 300


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
    assert "verdict Clones do not back these lineages" in report.render()


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
    result = SimpleNamespace(method="paga_dpt", params={}, pseudotime=pt, lineages=[1, 2])
    prep = SimpleNamespace(artifacts={"paga_connectivities": "CONN", "paga_tree": "TREE"})
    row = stability_row(object(), prep, object(), result, "true_label")  # type: ignore[arg-type]
    cov, tree, ev = seen["verdict_args"]
    assert (ev.n_replicates, ev.n_failed, ev.lineage, ev.pseudotime) == (7, 2, 0.11, 0.22)
    assert (cov, tree) == (0.99, 0.77)
    assert seen["tree_args"] == ("CONN", "TREE")
    assert np.array_equal(seen["cov_arg"], pt) and seen["cov_arg"].dtype == np.float64
    assert seen["spec"].label_key == "true_label"
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


def _fake_real(monkeypatch: pytest.MonkeyPatch, now: list[float]) -> list[Any]:
    import numpy as np

    from lumencell.bench import datasets, runner
    from lumencell.bench.metrics import Scores
    from lumencell.prep import cache

    ds = SimpleNamespace(root_spec="x", label_key="cell_type", clonal=None, report_time_key=None)
    monkeypatch.setattr(datasets, "DATASETS", {"d": ds})
    monkeypatch.setattr(datasets, "is_cached", lambda d: True)
    monkeypatch.setattr(datasets, "local_path", lambda d: "p")
    monkeypatch.setattr(datasets, "truth_for", lambda d, obs: None)
    monkeypatch.setattr("lumencell.engine.roots.parse_root", lambda s: "root")
    monkeypatch.setattr(runner, "load", lambda p: SimpleNamespace(obs=None))
    monkeypatch.setattr(cache, "PrepCache", lambda data: SimpleNamespace(run=lambda pp: SimpleNamespace()))
    ctx = SimpleNamespace(leiden=np.array([0, 1, 2]))
    monkeypatch.setattr(runner.InferenceContext, "from_prep", lambda *a, **k: ctx)
    res = SimpleNamespace(pseudotime=None, membership=None, warnings=[])
    monkeypatch.setattr(runner, "run_method", lambda *a, **k: res)
    fields = {f: 0.0 for f in Scores.__dataclass_fields__}
    fields.update(n_spurious=0)
    monkeypatch.setattr(runner, "score", lambda *a, **k: Scores(**fields))

    def tick() -> float:
        now[0] += 2.0
        return now[0] - 2.0

    monkeypatch.setattr(runner.time, "perf_counter", tick)
    calls: list[Any] = []
    return calls


def test_run_real_case_stability_flag(monkeypatch: pytest.MonkeyPatch) -> None:
    from lumencell.bench import runner

    now = [10.0]
    calls = _fake_real(monkeypatch, now)
    monkeypatch.setattr(
        runner, "stability_row", lambda *a, **k: calls.append(a) or StabilityRow("ok", 0.5, 0.5, 0.5, "weak")
    )
    off = runner.run_real_case("d", stability=False)
    assert calls == [] and off.stability.status == "not_measured"

    now[0] = 10.0
    sentinel = StabilityRow("ok", 0.1, 0.2, 0.3, "weak")  # type: ignore[arg-type]

    def slow_stability(*a: Any, **k: Any) -> StabilityRow:
        calls.append(a)
        now[0] += 100.0  # stability takes a long time
        return sentinel

    monkeypatch.setattr(runner, "stability_row", slow_stability)
    on = runner.run_real_case("d", stability=True)
    assert len(calls) == 1 and calls[0][4] == "cell_type"
    assert on.stability is sentinel
    assert on.runtime_seconds == 2.0  # stops before stability's 100 s


def test_bench_real_forwards_stability_and_prints_verdict(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from lumencell.bench import datasets, runner
    from lumencell.bench.metrics import Scores

    fields = {f: 0.0 for f in Scores.__dataclass_fields__}
    fields.update(n_spurious=0)
    flags: list[bool] = []

    def fake_case(d: str, method: str, params: Any = None, stability: bool = False) -> CaseResult:
        flags.append(stability)
        st = StabilityRow("ok", 0.5, 0.5, 0.5, "weak")  # type: ignore[arg-type]
        return CaseResult(d, method, Scores(**fields), 3, 1.0, stability=st)

    monkeypatch.setattr(runner, "run_real_case", fake_case)
    monkeypatch.setattr(datasets, "is_cached", lambda ds: True)
    out = str(tmp_path / "o.json")
    args = ["bench", "real", "--methods", "paga_dpt", "--datasets", "pancreas", "--out", out]
    default = CliRunner().invoke(app, args)
    CliRunner().invoke(app, [*args, "--no-stability"])
    assert flags == [True, False]
    assert "verdict=weak" in default.output


def _real_data(level: str) -> dict[str, Any]:
    _, real = report.latest_real() or (None, {"methods": {}})
    d, base = next(iter(real["methods"]["paga_dpt"].items()))
    row = {k: v for k, v in base.items() if k != "tree_support"}
    if level:
        row["tree_support"] = level
    return {"commit": "abc", "methods": {"paga_dpt": {d: row}}}


def _render(monkeypatch: pytest.MonkeyPatch, data: dict[str, Any]) -> str:
    monkeypatch.setattr(
        report, "latest_real", lambda: (report.REPO / "benchmarks" / "results" / "real-x.json", data)
    )
    return report.render()


def test_report_hides_not_measured_verdict(monkeypatch: pytest.MonkeyPatch) -> None:
    assert ", verdict" not in _render(monkeypatch, _real_data("not_measured"))


def test_report_adds_footnote_only_with_verdicts(monkeypatch: pytest.MonkeyPatch) -> None:
    assert "The verdict is" not in _render(monkeypatch, _real_data("not_measured"))
    assert "The verdict is" not in _render(monkeypatch, _real_data(""))
    assert "The verdict is" in _render(monkeypatch, _real_data("weak"))


def test_the_report_footnote_says_only_supported_is_capped(monkeypatch: pytest.MonkeyPatch) -> None:
    """Owner ruling H-1: the footnote is the one user-facing statement of the cap."""
    text = _render(monkeypatch, _real_data("provisional"))
    assert "A supported verdict is capped at provisional" in text
    assert "A weak or not_supported verdict is a real failure and is not capped" in text
    assert "supported or weak" not in text
