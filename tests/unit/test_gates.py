import math

import pytest

from lumencell.bench.gates import check_gates, check_regression, load_thresholds
from lumencell.bench.metrics import Scores
from lumencell.bench.runner import CaseResult, ClonalRow

GATES = load_thresholds()["synthetic"]["gates"]


def _row(**kw: object) -> dict[str, object]:
    base: dict[str, object] = {
        "case": "c",
        "pt_spearman_global": 0.95,
        "pt_spearman_lineage": 0.95,
        "branch_accuracy": 0.95,
        "lineage_recall": 1.0,
        "n_spurious": 0,
        "root_pt_excess": 0.0,
        "clonal_excess_ratio": 0.8,
        "clonal_p": 0.001,
        "clonal_status": "ok",
    }
    base.update(kw)
    return base


def test_thresholds_file_has_the_documented_gates() -> None:
    gates = load_thresholds()["synthetic"]["gates"]
    assert set(gates) == {
        "pt_spearman_global",
        "pt_spearman_lineage",
        "branch_accuracy",
        "lineage_recall",
        "n_spurious",
        "root_pt_excess",
        "clonal_excess_ratio",
        "clonal_p",
    }


def test_passing_row_has_no_failures() -> None:
    assert check_gates(_row(), load_thresholds()["synthetic"]["gates"]) == []


def test_min_and_max_bounds() -> None:
    gates = load_thresholds()["synthetic"]["gates"]
    failed = {f.metric for f in check_gates(_row(branch_accuracy=0.5, n_spurious=3), gates)}
    assert failed == {"branch_accuracy", "n_spurious"}


def test_nan_metrics_are_skipped() -> None:
    gates = load_thresholds()["synthetic"]["gates"]
    assert check_gates(_row(branch_accuracy=math.nan), gates) == []


def test_regression_guard() -> None:
    base = {"pt_spearman_global": 0.95, "branch_accuracy": 0.9}
    assert check_regression(_row(pt_spearman_global=0.93), base, 0.03) == []
    failed = check_regression(_row(pt_spearman_global=0.90), base, 0.03)
    assert [f.metric for f in failed] == ["pt_spearman_global"]


def test_clonal_gate_entries_are_declared() -> None:
    assert GATES["clonal_excess_ratio"] == {"min": 0.5, "nan": "fail", "status": "clonal_status"}
    assert GATES["clonal_p"] == {"max": 0.01, "nan": "fail", "status": "clonal_status"}


def test_clonal_nan_fails_with_its_reason() -> None:
    (f,) = check_gates(_row(clonal_excess_ratio=math.nan), GATES)
    assert (f.metric, f.reason) == ("clonal_excess_ratio", "nan")
    assert str(f).endswith(" [nan]")


@pytest.mark.parametrize(
    "status", ["too_few_clones", "failed", "not_assessable", "single_lineage", "not_measured"]
)
def test_any_status_but_ok_fails_both_clonal_gates(status: str) -> None:
    failed = check_gates(_row(clonal_status=status), GATES)
    assert {(f.metric, f.reason) for f in failed} == {
        ("clonal_excess_ratio", status),
        ("clonal_p", status),
    }


def test_not_applicable_skips_the_clonal_gates() -> None:
    row = _row(clonal_status="not_applicable", clonal_excess_ratio=math.nan, clonal_p=math.nan)
    assert check_gates(row, GATES) == []


def test_a_row_from_before_the_clonal_fields_fails_as_not_measured() -> None:
    old = {k: v for k, v in _row().items() if not k.startswith("clonal_")}
    assert {f.reason for f in check_gates(old, GATES)} == {"not_measured"}


def test_existing_gates_still_skip_nan_and_unknown_policies_raise() -> None:
    assert check_gates(_row(branch_accuracy=math.nan), GATES) == []
    with pytest.raises(ValueError, match="nan"):
        check_gates(_row(), {"branch_accuracy": {"min": 0.9, "nan": "maybe"}})


def test_regression_guard_covers_the_clonal_ratio() -> None:
    failed = check_regression(_row(clonal_excess_ratio=0.6), {"clonal_excess_ratio": 0.7}, 0.03)
    assert [f.metric for f in failed] == ["clonal_excess_ratio"]


def test_case_rows_always_carry_the_clonal_fields() -> None:
    scores = Scores(
        pt_spearman_global=0.0,
        pt_spearman_lineage=0.0,
        pt_spearman_by_lineage={},
        branch_accuracy=0.0,
        branch_accuracy_all=0.0,
        lineage_recall=0.0,
        n_spurious=0,
        n_inferred_lineages=0,
        root_pt_median=0.0,
        root_pt_excess=0.0,
        terminal_share=0.0,
        coverage=0.0,
        matches=[],
    )
    row = CaseResult("c", "m", scores=scores, n_clusters=1, runtime_seconds=0.0).row()
    assert row["clonal_status"] == "not_measured" and math.isnan(row["clonal_excess_ratio"])
    assert {
        "clonal_p",
        "clonal_delta_base",
        "clonal_eta_oracle",
        "clonal_eta_method",
        "clonal_n_clones",
        "clonal_n_cells",
        "clonal_effective_clones",
        "clonal_eta_m3",
        "clonal_p_m3",
    } <= row.keys()
    assert row["clonal_n_cells"] == 0
    assert math.isnan(row["clonal_effective_clones"]) and math.isnan(row["clonal_p_m3"])
    assert "clonal_n_cells" not in load_thresholds()["synthetic"]["gates"]  # reported only


def _scores() -> Scores:
    return Scores(
        pt_spearman_global=0.0,
        pt_spearman_lineage=0.0,
        pt_spearman_by_lineage={},
        branch_accuracy=0.0,
        branch_accuracy_all=0.0,
        lineage_recall=0.0,
        n_spurious=0,
        n_inferred_lineages=0,
        root_pt_median=0.0,
        root_pt_excess=0.0,
        terminal_share=0.0,
        coverage=0.0,
        matches=[],
    )


def test_row_reports_every_clonal_field_with_its_own_value() -> None:
    clonal = ClonalRow(
        "ok",
        excess_ratio=0.812345,
        p=0.0011,
        delta_base=0.123456,
        eta_oracle=0.4,
        eta_method=0.3,
        n_clones=21,
        n_cells=99,
        effective_clones=7.5,
        eta_m3=0.2,
        p_m3=0.02,
        reason="x",
    )
    row = CaseResult("c", "m", _scores(), 1, 0.0, clonal=clonal).row()
    got = {k: v for k, v in row.items() if k.startswith("clonal_")}
    assert got == {
        "clonal_excess_ratio": 0.8123,
        "clonal_p": 0.0011,
        "clonal_delta_base": 0.1235,
        "clonal_status": "ok",
        "clonal_eta_oracle": 0.4,
        "clonal_eta_method": 0.3,
        "clonal_effective_clones": 7.5,
        "clonal_eta_m3": 0.2,
        "clonal_p_m3": 0.02,
        "clonal_n_clones": 21,
        "clonal_n_cells": 99,
        "clonal_reason": "x",
    }
    assert CaseResult("c", "m", _scores(), 1, 0.0).row()["clonal_reason"] is None


def test_status_ok_still_enforces_the_bounds() -> None:
    failed = check_gates(_row(clonal_excess_ratio=0.3, clonal_p=0.02), GATES)
    assert {(f.metric, f.reason) for f in failed} == {("clonal_excess_ratio", None), ("clonal_p", None)}
    assert check_gates(_row(clonal_excess_ratio=0.5, clonal_p=0.01), GATES) == []


def test_nan_fail_without_status_fails_with_reason_nan() -> None:
    gate = {"m": {"min": 0.5, "nan": "fail"}}
    (f,) = check_gates({"case": "c", "m": math.nan}, gate)
    assert (f.metric, f.reason) == ("m", "nan")


def test_status_ok_with_nan_skip_skips_a_nan() -> None:
    gate = {"m": {"min": 0.5, "nan": "skip", "status": "s"}}
    assert check_gates({"case": "c", "m": math.nan, "s": "ok"}, gate) == []


def test_missing_metric_under_status_ok_fails_as_nan() -> None:
    gate = {"m": {"min": 0.5, "nan": "fail", "status": "s"}}
    (f,) = check_gates({"case": "c", "s": "ok"}, gate)
    assert (f.metric, f.reason) == ("m", "nan")


@pytest.mark.parametrize(
    ("gate", "match"),
    [
        ({"min": 0.5, "status": "s"}, "nan"),
        ({"mim": 0.5}, "mim"),
        ({"nan": "fail"}, "min or max"),
    ],
)
def test_malformed_gates_raise(gate: dict[str, object], match: str) -> None:
    with pytest.raises(ValueError, match=match):
        check_gates({"case": "c", "m": 1.0, "s": "ok"}, {"m": gate})
