import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any, get_args

import numpy as np
import pytest

from lumencell.bench.gates import check_gates, load_thresholds
from lumencell.bench.report import render
from lumencell.bench.runner import ClonalRow, ClonalRowStatus, clonal_row, run_synthetic_case
from lumencell.bench.synth import generate, truth_from_adata
from lumencell.bench.synth_clones import DEFAULT_LAYER, ClonalLayer, add_clones, clone_truth
from lumencell.bench.truth import Truth

GATES = load_thresholds()["synthetic"]["gates"]


def _row(tmp_path: Path, topology: str, **kw: Any) -> dict[str, Any]:
    return run_synthetic_case(tmp_path, topology, n_cells=1500, n_genes=1000, **kw).row()


def _clone_truth(layer: ClonalLayer = DEFAULT_LAYER, topology: str = "bifurcation") -> Truth:
    adata = generate(topology, "realistic", n_cells=1500, n_genes=200, seed=0)
    return clone_truth(truth_from_adata(adata), add_clones(adata, layer, seed=0))


def _result(membership: Any, truth: Truth) -> Any:
    return SimpleNamespace(membership=membership, pseudotime=truth.order.astype(np.float32))


def _pca(n: int) -> Any:
    return np.random.default_rng(0).normal(size=(n, 10))


def _gate_failures(row: ClonalRow) -> list[Any]:
    values = {
        "case": "c",
        "clonal_status": row.status,
        "clonal_excess_ratio": row.excess_ratio,
        "clonal_p": row.p,
    }
    return check_gates(values, {k: v for k, v in GATES.items() if k.startswith("clonal")})


def test_branching_cases_carry_a_clonal_score(tmp_path: Path) -> None:
    row = _row(tmp_path, "bifurcation")
    assert row["clonal_status"] in get_args(ClonalRowStatus)
    assert row["clonal_eta_oracle"] >= 0.05
    assert row["clonal_n_clones"] >= 20
    assert row["clonal_n_cells"] >= 2 * row["clonal_n_clones"]
    assert row["clonal_effective_clones"] > 0
    assert math.isfinite(row["clonal_eta_m3"])
    assert 0 < row["clonal_p_m3"] <= 1


def test_linear_cases_are_not_applicable(tmp_path: Path) -> None:
    row = _row(tmp_path, "linear")
    assert row["clonal_status"] == "not_applicable"
    assert math.isnan(row["clonal_excess_ratio"])
    assert row["clonal_n_cells"] == 0
    assert math.isnan(row["clonal_effective_clones"])
    assert math.isnan(row["clonal_eta_m3"])
    assert not [f for f in check_gates(row, GATES) if f.metric.startswith("clonal")]


def test_existing_metrics_do_not_move_with_the_clonal_layer(tmp_path: Path) -> None:
    with_layer = _row(tmp_path, "bifurcation", clonal_layer=DEFAULT_LAYER)
    without = _row(tmp_path, "bifurcation", clonal_layer=None)
    assert without["clonal_status"] == "not_measured"
    skip = {"runtime_seconds"}
    keys = {k for k in with_layer if k not in skip and not k.startswith("clonal_")}
    assert keys
    # repr, because NaN != NaN (terminal_share is NaN on a synthetic case)
    assert repr({k: with_layer[k] for k in sorted(keys)}) == repr({k: without[k] for k in sorted(keys)})


def test_a_single_lineage_method_fails_the_gate() -> None:
    truth = _clone_truth()
    row = clonal_row(truth, _result(truth.membership()[:1], truth), _pca(truth.n_obs), 10)
    assert row.status == "single_lineage"
    assert row.eta_oracle >= 0.05 and row.n_clones >= 20
    assert math.isnan(row.excess_ratio) and math.isnan(row.eta_m3)
    assert {f.metric for f in _gate_failures(row)} == {"clonal_excess_ratio", "clonal_p"}


def test_an_oracle_below_its_floor_is_not_assessable() -> None:
    truth = _clone_truth(ClonalLayer(fate_bias=0.0))
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert row.status == "not_assessable"
    assert math.isfinite(row.eta_method) and math.isfinite(row.eta_m3)
    assert math.isnan(row.excess_ratio)
    assert _gate_failures(row)


def test_a_scored_method_gets_a_ratio() -> None:
    truth = _clone_truth()
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert row.status == "ok"
    assert row.excess_ratio == pytest.approx(row.eta_method / row.eta_oracle)
    assert row.excess_ratio == pytest.approx(1.0)  # the method is the truth
    assert row.p <= 0.01 and not _gate_failures(row)


def test_a_truth_without_clones_is_not_measured_and_fails_the_gates() -> None:
    adata = generate("bifurcation", "realistic", n_cells=600, n_genes=100, seed=0)
    truth = truth_from_adata(adata)
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert row.status == "not_measured"
    assert {f.metric for f in _gate_failures(row)} == {"clonal_excess_ratio", "clonal_p"}


def test_a_truth_with_one_lineage_is_not_applicable() -> None:
    truth = _clone_truth(topology="linear")
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert repr(row) == repr(ClonalRow("not_applicable"))


def test_an_error_is_recorded_as_failed(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def boom(*args: object, **kwargs: object) -> None:
        raise RuntimeError("no scores today")

    monkeypatch.setattr("lumencell.clonal.scores.score_clonal", boom)
    truth = _clone_truth()
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert row.status == "failed"
    assert row.reason == "RuntimeError: no scores today"
    assert "no scores today" in capsys.readouterr().err


def test_the_report_says_the_clonal_gate_is_not_a_branch_gate() -> None:
    assert "not a branch gate" in render()
