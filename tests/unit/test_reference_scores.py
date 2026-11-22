"""The packaged reference rows match the committed results, and say what they say.

The pinned values are the committed real-2026-11-20 results and baseline; when a newer result is
committed, regenerate the file and update these pins with it, never the thresholds.
"""

from lumencell.bench.reference import PACKAGED, render
from lumencell.bench.reference_data import rows_for


def _row(method: str, dataset: str, gate: str) -> dict[str, object]:
    (row,) = [r for r in rows_for(method) if str(r["dataset"]).startswith(dataset) and r["gate"] == gate]
    return row


def test_packaged_file_matches_the_committed_results() -> None:
    assert PACKAGED.read_text(encoding="utf-8") == render(), (
        "reference_scores.json is stale: run `uv run python -m lumencell.bench.reference`"
    )


def test_failures_are_kept_and_passes_are_shown() -> None:
    assert _row("paga_dpt", "Pancreas", "branch_accuracy") == {
        "dataset": "Pancreas (mouse, E15.5)",
        "metric": "Branch accuracy",
        "gate": "branch_accuracy",
        "score": 0.6306,
        "threshold": "≥ 0.8",
        "passed": False,
        "detail": None,
        "source": "benchmarks/results/real-2026-11-20.json @ a7d69c4",
    }
    assert _row("palantir", "Pancreas", "branch_accuracy")["passed"] is True
    assert _row("palantir", "Gastrulation", "terminal_share")["passed"] is False


def test_synthetic_rows_count_the_failing_cases() -> None:
    row = _row("paga_dpt", "Synthetic", "pt_spearman_global")
    assert row["detail"] == "worst case; fails 3 of 12 cases"
    assert row["passed"] is False
    # Single-lineage cases have no branch accuracy: they are not measured, so the denominator is 9.
    assert _row("paga_dpt", "Synthetic", "branch_accuracy")["detail"] == "worst case; fails 3 of 9 cases"


def test_a_dataset_with_no_gated_lineage_is_never_a_pass() -> None:
    for method in ("paga_dpt", "palantir", "scfates_ppt"):
        rows = [r for r in rows_for(method) if str(r["dataset"]).startswith("LCMV")]
        assert rows, method
        for row in rows:
            assert row["passed"] is None and row["score"] is None
            assert str(row["detail"]).startswith("lineage gates: not measured")


def test_every_row_names_its_source_and_threshold() -> None:
    for method in ("paga_dpt", "palantir", "scfates_ppt"):
        for row in rows_for(method):
            assert str(row["source"]).startswith(("benchmarks/baseline.json @ ", "benchmarks/results/real-"))
            assert str(row["threshold"]).startswith(("≥ ", "≤ "))


def test_an_unknown_method_has_no_rows() -> None:
    assert rows_for("no_such_method") == []
