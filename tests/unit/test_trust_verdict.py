"""The served trust mappings in lumencell.trust.verdict, pinned value by value."""

import math

import pytest

from lumencell.api import models
from lumencell.engine.tree_support import Thresholds, Verdict, _check
from lumencell.trust.verdict import (
    clonal_evidence,
    stability_evidence,
    to_api_tree_support,
    tree_support_for_run,
)


@pytest.fixture(autouse=True)
def _verified(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("lumencell.engine.tree_support.VERDICT_VERIFIED", True)


def _summary(kind: str, r: int) -> models.StabilitySummary:
    return models.StabilitySummary(
        kind=kind,  # type: ignore[arg-type]
        r=r,
        fraction=0.8,
        seeds=list(range(r)),
        spearman_median=None,
        spearman_min=None,
        share_matched_median=None,
        share_matched_min=None,
        mean_jaccard_median=None,
        mean_jaccard_min=None,
        n_lineages_median=None,
        n_lineages_min=0,
        replicates=[],
    )


def _stability(
    status: str = "ok",
    lineage: float | None = 0.9,
    pseudotime: float | None = 0.8,
    notes: list[str] | None = None,
) -> models.Stability:
    return models.Stability(
        status=status,  # type: ignore[arg-type]
        match_jaccard=0.5,
        reference_lineages=2,
        n_failed=2,
        reseeds=_summary("reseed", 3),
        subsamples=_summary("subsample", 5),
        lineage_stability=lineage,
        pseudotime_stability=pseudotime,
        notes=notes or [],
    )


def test_an_ok_stability_is_done_with_its_counts_and_values() -> None:
    e = stability_evidence(_stability())
    assert e.state == "done"
    assert (e.n_replicates, e.n_failed) == (8, 2)
    assert (e.lineage, e.pseudotime) == (0.9, 0.8)


def test_an_undefined_stability_value_is_nan() -> None:
    e = stability_evidence(_stability(lineage=None))
    assert math.isnan(e.lineage) and e.pseudotime == 0.8


def test_a_failed_stability_is_unavailable_with_its_first_note() -> None:
    e = stability_evidence(_stability(status="failed", notes=["a", "b"]))
    assert e.state == "unavailable" and e.note == "a"


def test_an_ok_stability_is_judged_by_value_in_the_verdict() -> None:
    graph = {"n_unassigned": 0, "n_obs": 10, "lineages": [1, 2]}
    t = Thresholds(tau_tree=None, tau_lin=0.85, tau_pt=0.75)
    passing = tree_support_for_run(graph, 0.9, _stability(), None, t)
    assert {c.id: c.state for c in passing.checks}["C3"] == "pass"
    assert {c.id: c.state for c in passing.checks}["C4"] == "pass"
    failing = tree_support_for_run(graph, 0.9, _stability(lineage=0.5, pseudotime=0.5), None, t)
    assert {c.id: c.state for c in failing.checks}["C3"] == "fail"
    assert {c.id: c.state for c in failing.checks}["C4"] == "fail"


def _agreement(status: str, notes: list[str]) -> models.ClonalAgreement:
    return models.ClonalAgreement(
        status=status,  # type: ignore[arg-type]
        settings=None,
        primary=None,
        sensitivity_m3=None,
        coupling=[],
        notes=notes,
    )


def test_a_failed_clonal_score_carries_its_note_and_undefined_numbers() -> None:
    e = clonal_evidence(_agreement("failed", ["n"]), 3)
    assert e.state == "scored" and e.status == "failed" and e.note == "n"
    assert math.isnan(e.p) and math.isnan(e.delta_eta)


def test_an_ok_clonal_score_carries_no_note() -> None:
    e = clonal_evidence(_agreement("ok", ["n"]), 3)
    assert e.status == "ok" and e.note == ""


def test_an_undefined_check_value_is_served_as_null() -> None:
    v = Verdict("provisional", (_check("C2", "fail", math.nan, 0.5, "x"),), ("x",))
    assert to_api_tree_support(v).checks[0].value is None
