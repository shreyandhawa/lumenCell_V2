import math

import numpy as np
import pytest

from lumencell.bench.metrics import score, spearman
from lumencell.bench.truth import Truth

N = 2000
FORK = 0.4


def _bifurcation_truth(seed: int = 0) -> tuple[Truth, np.ndarray, np.ndarray]:
    """Hand-built truth: trunk below FORK shared by lineages A and B."""
    rng = np.random.default_rng(seed)
    t = rng.random(N)
    branch = rng.integers(0, 2, N)  # 0 -> A, 1 -> B (only meaningful past the fork)
    trunk = t < FORK
    on_a = trunk | (branch == 0)
    on_b = trunk | (branch == 1)
    truth = Truth(order=t, lineages={"A": on_a, "B": on_b}, root=t < 0.05)
    return truth, t, np.stack([on_a, on_b])


def test_perfect_inference_scores_one() -> None:
    truth, t, membership = _bifurcation_truth()
    s = score(t, membership, truth)
    assert s.pt_spearman_global == pytest.approx(1.0)
    assert s.pt_spearman_lineage == pytest.approx(1.0)
    assert s.branch_accuracy == pytest.approx(1.0)
    assert s.lineage_recall == 1.0
    assert s.n_spurious == 0
    assert s.root_pt_median < 0.05
    assert s.coverage == 1.0


def test_reversed_pseudotime_scores_minus_one() -> None:
    truth, t, membership = _bifurcation_truth()
    s = score(1.0 - t, membership, truth)
    assert s.pt_spearman_global == pytest.approx(-1.0)
    assert s.pt_spearman_lineage == pytest.approx(-1.0)
    assert s.root_pt_median > 0.9


def test_random_pseudotime_scores_near_zero() -> None:
    truth, _t, membership = _bifurcation_truth()
    noise = np.random.default_rng(1).random(N)
    s = score(noise, membership, truth)
    assert abs(s.pt_spearman_global) < 0.1


def test_lineage_ids_do_not_matter() -> None:
    truth, t, membership = _bifurcation_truth()
    s = score(t, membership[::-1], truth)
    assert s.branch_accuracy == pytest.approx(1.0)
    assert {m.true_name for m in s.matches} == {"A", "B"}


def test_random_branch_assignment_is_at_chance() -> None:
    truth, t, _ = _bifurcation_truth()
    rng = np.random.default_rng(2)
    coin = rng.integers(0, 2, N)
    trunk = t < FORK
    random_membership = np.stack([trunk | (coin == 0), trunk | (coin == 1)])
    s = score(t, random_membership, truth)
    assert s.branch_accuracy == pytest.approx(0.5, abs=0.05)


def test_missing_the_branch_scores_zero() -> None:
    truth, t, _ = _bifurcation_truth()
    everyone_on_both = np.ones((2, N), dtype=bool)
    s = score(t, everyone_on_both, truth)
    assert s.branch_accuracy == 0.0
    assert s.lineage_recall == 0.0
    assert s.n_spurious == 2


def test_one_lineage_for_a_bifurcation_has_half_recall() -> None:
    truth, t, membership = _bifurcation_truth()
    only_a = membership[:1]
    s = score(t, only_a, truth)
    assert s.lineage_recall == 0.5
    assert s.n_spurious == 0


def test_one_lineage_merging_both_branches_finds_neither() -> None:
    truth, t, _ = _bifurcation_truth()
    merged = np.ones((1, N), dtype=bool)
    s = score(t, merged, truth)
    assert s.lineage_recall == 0.0
    assert s.n_spurious == 1


def test_extra_lineage_counts_as_spurious() -> None:
    truth, t, membership = _bifurcation_truth()
    extra = np.zeros((1, N), dtype=bool)
    extra[0, :20] = True
    s = score(t, np.vstack([membership, extra]), truth)
    assert s.lineage_recall == 1.0
    assert s.n_spurious == 1
    assert s.n_inferred_lineages == 3


def test_nan_pseudotime_is_excluded_and_reported() -> None:
    truth, t, membership = _bifurcation_truth()
    pt = t.copy()
    pt[:200] = np.nan
    s = score(pt, membership, truth)
    assert s.coverage == pytest.approx(0.9)
    assert s.pt_spearman_global == pytest.approx(1.0)


def test_single_lineage_has_no_branch_accuracy() -> None:
    t = np.linspace(0, 1, 100)
    truth = Truth(order=t, lineages={"only": np.ones(100, dtype=bool)}, root=t < 0.1)
    s = score(t, np.ones((1, 100), dtype=bool), truth)
    assert math.isnan(s.branch_accuracy)
    assert s.lineage_recall == 1.0


def test_shape_validation() -> None:
    truth, t, membership = _bifurcation_truth()
    with pytest.raises(ValueError, match="one value per cell"):
        score(t[:10], membership, truth)
    with pytest.raises(ValueError, match="n_lineages, n_obs"):
        score(t, membership[:, :10], truth)


def test_truth_rejects_empty_root() -> None:
    with pytest.raises(ValueError, match="empty"):
        Truth(order=np.zeros(3), lineages={}, root=np.zeros(3, dtype=bool))


def test_spearman_needs_variation() -> None:
    assert math.isnan(spearman(np.ones(10), np.arange(10.0)))
    assert math.isnan(spearman(np.arange(2.0), np.arange(2.0)))
