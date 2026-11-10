import math

from lumencell.bench.gates import check_gates, check_regression, load_thresholds


def _row(**kw: float) -> dict[str, object]:
    base: dict[str, object] = {
        "case": "c",
        "pt_spearman_global": 0.95,
        "pt_spearman_lineage": 0.95,
        "branch_accuracy": 0.95,
        "lineage_recall": 1.0,
        "n_spurious": 0,
        "root_pt_excess": 0.0,
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
