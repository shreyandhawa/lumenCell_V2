import math

import numpy as np
import pytest
from scipy import sparse

from lumencell.clonal.scores import ClonalResult
from lumencell.engine.tree_support import (
    CHECK_NAMES,
    ClonalEvidence,
    StabilityEvidence,
    Thresholds,
    calibrated_thresholds,
    coverage,
    tree_likeness,
    verdict,
)
from tests.conftest import Prepared

TAU = Thresholds(tau_tree=0.7, tau_lin=0.8, tau_pt=0.85)
GOOD_STAB = StabilityEvidence("done", lineage=0.9, pseudotime=0.95)
OK_CLONAL = ClonalEvidence("scored", 3, status="ok", p=0.001, delta_eta=0.2)
NO_CLONES = ClonalEvidence("not_configured", 0)


def _v(
    cov: float = 1.0,
    tree: float | None = 0.9,
    stab: StabilityEvidence = GOOD_STAB,
    clonal: ClonalEvidence = OK_CLONAL,
    thresholds: Thresholds = TAU,
):  # type: ignore[no-untyped-def]
    return verdict(cov, tree, stab, clonal, thresholds)


def _check(v, cid: str):  # type: ignore[no-untyped-def]
    return next(c for c in v.checks if c.id == cid)


CONN = sparse.csr_matrix(np.array([[0, 0.5, 0.2], [0.5, 0, 0.4], [0.2, 0.4, 0]]))


def test_tree_likeness_by_hand() -> None:
    once = sparse.csr_matrix(np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]], dtype=float))
    assert tree_likeness(CONN, once) == pytest.approx(0.9 / 1.1, abs=1e-12)
    assert tree_likeness(CONN, once + once.T) == pytest.approx(0.9 / 1.1, abs=1e-12)


def test_a_path_graph_scores_one() -> None:
    path = sparse.csr_matrix(np.array([[0, 0.5, 0], [0.5, 0, 0.4], [0, 0.4, 0]]))
    assert tree_likeness(path, path) == pytest.approx(1.0)


def test_a_graph_without_edges_fails_c2_with_a_reason() -> None:
    zero = sparse.csr_matrix((1, 1))
    t = tree_likeness(zero, zero)
    assert math.isnan(t)
    v = _v(tree=t)
    assert v.level == "not_supported"
    c2 = _check(v, "C2")
    assert c2.state == "fail" and c2.reason


def test_tree_likeness_of_the_prep_is_in_the_unit_interval(bifurcation: Prepared) -> None:
    art = bifurcation.prep.artifacts
    t = tree_likeness(art["paga_connectivities"], art["paga_tree"])
    assert 0.0 < t <= 1.0


def test_coverage_is_the_finite_share() -> None:
    assert coverage(np.array([0.1, np.nan, 0.3, np.inf])) == 0.5


PENDING = StabilityEvidence("pending")
UNAVAIL = StabilityEvidence("unavailable", note="no replicate ran")
BAD_C5 = ClonalEvidence("scored", 3, status="ok", p=0.5, delta_eta=0.2)
N_A = ClonalEvidence("scored", 3, status="single_lineage")
LEVEL_CASES = [
    (dict(cov=0.5), "not_supported"),
    (dict(tree=0.1), "not_supported"),
    (dict(tree=None), "provisional"),
    (dict(cov=0.5, tree=None), "not_supported"),
    (dict(stab=StabilityEvidence("done", lineage=0.5, pseudotime=0.95)), "weak"),
    (dict(stab=StabilityEvidence("done", lineage=0.9, pseudotime=0.5)), "weak"),
    (dict(stab=PENDING, clonal=BAD_C5), "weak"),
    (dict(stab=PENDING, clonal=NO_CLONES), "provisional"),
    (dict(stab=UNAVAIL), "provisional"),
    ({}, "supported"),
    (dict(clonal=NO_CLONES), "supported"),
]


@pytest.mark.parametrize(("kw", "level"), LEVEL_CASES)
def test_levels(kw: dict, level: str) -> None:  # type: ignore[type-arg]
    assert _v(**kw).level == level


def test_a_left_out_check_is_excluded_with_a_reason() -> None:
    v = _v(thresholds=Thresholds(0.7, 0.8, None))
    assert v.level == "supported"
    c4 = _check(v, "C4")
    assert c4.state == "excluded" and "left out" in (c4.reason or "")


def test_every_replicate_failing_leaves_stability_unmeasured() -> None:
    allfail = StabilityEvidence("done", lineage=0.0, pseudotime=0.0, n_replicates=5, n_failed=5)
    v = _v(stab=allfail)
    assert v.level == "provisional"
    assert _check(v, "C3").state == "unavailable" and _check(v, "C4").state == "unavailable"
    some = StabilityEvidence("done", lineage=0.0, pseudotime=0.0, n_replicates=5, n_failed=2)
    assert _v(stab=some).level == "weak"


def test_c1_bound_is_inclusive() -> None:
    assert _check(_v(cov=0.95), "C1").state == "pass"
    assert _check(_v(cov=0.9499), "C1").state == "fail"


C5_CASES = [
    (ClonalEvidence("scored", 3, status="ok", p=0.01, delta_eta=0.001), "pass", "supported"),
    (ClonalEvidence("scored", 3, status="ok", p=0.011, delta_eta=0.1), "fail", "weak"),
    (ClonalEvidence("scored", 3, status="ok", p=0.001, delta_eta=0.0), "fail", "weak"),
    (ClonalEvidence("scored", 3, status="ok", p=0.001, delta_eta=math.nan), "fail", "weak"),
    (ClonalEvidence("scored", 3, status="ok", p=math.nan, delta_eta=0.1), "fail", "weak"),
    (ClonalEvidence("scored", 3, status="too_few_clones"), "fail", "weak"),
    (ClonalEvidence("scored", 3, status="refused"), "fail", "weak"),
    (ClonalEvidence("scored", 3, status="failed"), "fail", "weak"),
    (NO_CLONES, "not_applicable", "supported"),
    (ClonalEvidence("pending", 0), "pending", "provisional"),
]


@pytest.mark.parametrize(("clonal", "state", "level"), C5_CASES)
def test_c5_rules(clonal: ClonalEvidence, state: str, level: str) -> None:
    v = _v(clonal=clonal)
    assert _check(v, "C5").state == state
    assert v.level == level
    if clonal.state == "scored" and clonal.status != "ok":
        assert clonal.status in (_check(v, "C5").reason or "")


def test_a_single_lineage_clonal_score_fails_c5() -> None:
    # Spec 5.3 over the plan (ruling F12): L >= 2 with every scored clone in one lineage fails C5.
    v = _v(clonal=N_A)
    c5 = _check(v, "C5")
    assert c5.state == "fail" and "single_lineage" in (c5.reason or "")
    assert v.level == "weak"


def test_c5_reasons_claim_no_direction() -> None:
    for c, _, _ in C5_CASES:
        r = _check(_v(clonal=c), "C5").reason or ""
        for word in ("first", "earlier", "precede", "origin", "direction"):
            assert word not in r.lower()


def test_from_result_reads_the_clonal_result() -> None:
    r = ClonalResult.__new__(ClonalResult)
    for k, val in dict(status="ok", p=0.004, delta_eta=0.3, n_lineages=3).items():
        object.__setattr__(r, k, val)
    e = ClonalEvidence.from_result(r)
    assert (e.state, e.status, e.p, e.delta_eta, e.n_lineages) == ("scored", "ok", 0.004, 0.3, 3)


def test_reasons_one_per_check_that_did_not_pass_or_apply() -> None:
    v = _v(cov=0.5, tree=None, stab=PENDING, clonal=BAD_C5, thresholds=Thresholds(0.7, 0.8, None))
    expect = [c.reason for c in v.checks if c.state not in ("pass", "not_applicable")]
    assert all(expect) and v.reasons == tuple(expect)
    assert [c.id for c in v.checks] == ["C1", "C2", "C3", "C4", "C5"]
    assert all(c.name == CHECK_NAMES[c.id] for c in v.checks)
    clean = _v()
    assert clean.reasons == () and all(c.reason is None for c in clean.checks)


def test_c5_is_not_applicable_below_two_lineages() -> None:
    v = _v(clonal=ClonalEvidence("scored", 1, status="single_lineage"))
    c5 = _check(v, "C5")
    assert c5.state == "not_applicable" and c5.reason
    assert v.level == "supported"


def test_c5_with_two_lineages_in_one_still_fails() -> None:
    c5 = _check(_v(clonal=ClonalEvidence("scored", 2, status="single_lineage")), "C5")
    assert c5.state == "fail"


def test_n_lineages_has_no_default() -> None:
    with pytest.raises(TypeError):
        ClonalEvidence("scored")  # type: ignore[call-arg]


def test_values_exactly_at_the_thresholds_pass() -> None:
    v = _v(tree=0.7, stab=StabilityEvidence("done", lineage=0.8, pseudotime=0.85))
    assert [_check(v, c).state for c in ("C2", "C3", "C4")] == ["pass"] * 3
    assert v.level == "supported"


def test_tree_likeness_symmetrises_by_maximum() -> None:
    upper = sparse.csr_matrix(np.array([[0, 0.5, 0], [0, 0, 0.4], [0, 0, 0]]))
    tree_lower = sparse.csr_matrix(np.array([[0, 0, 0], [0.5, 0, 0], [0, 0.4, 0]]))
    assert tree_likeness(upper, upper) == pytest.approx(1.0)
    assert tree_likeness(upper, tree_lower) == pytest.approx(1.0)
    lower = upper.T.tocsr()
    assert tree_likeness(lower, upper) == pytest.approx(1.0)


def test_excluded_c2_stays_excluded_when_tree_likeness_is_none() -> None:
    thr = Thresholds(None, 0.8, 0.85)
    for tree in (0.9, None):
        c2 = _check(_v(tree=tree, thresholds=thr), "C2")
        assert c2.state == "excluded" and c2.reason


def test_not_supported_wins_over_other_fails() -> None:
    low_stab = StabilityEvidence("done", lineage=0.1, pseudotime=0.95)
    assert _v(cov=0.5, stab=low_stab).level == "not_supported"
    assert _v(tree=0.1, clonal=BAD_C5).level == "not_supported"


def test_calibrated_thresholds_are_the_committed_taus() -> None:
    from lumencell.trust import constants

    t = calibrated_thresholds()
    assert (t.tau_tree, t.tau_lin, t.tau_pt) == (constants.TAU_TREE, constants.TAU_LIN, constants.TAU_PT)
    assert (t.coverage_min, t.clonal_p_max) == (constants.COVERAGE_MIN, constants.CLONAL_P_MAX)
