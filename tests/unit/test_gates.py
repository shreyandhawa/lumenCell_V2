import math

import pytest

from lumencell.bench.gates import check_gates, check_regression, load_thresholds
from lumencell.bench.metrics import Scores
from lumencell.bench.runner import CaseResult

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
