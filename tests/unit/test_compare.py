import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata, spearmanr

from lumencell.engine.compare import compare, jaccard_matrix, match_lineages
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from tests.conftest import Prepared


def test_a_shared_mask_ignores_cells_outside_it() -> None:
    ma = np.array([[1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0]], dtype=bool)
    mb = ma.copy()
    mb[0, 4] = mb[1, 5] = True  # extra cells, outside the mask
    shared = np.array([1, 1, 1, 1, 0, 0], dtype=bool)
    with_mask = jaccard_matrix(ma, mb, shared)
    without = jaccard_matrix(ma, mb)
    assert with_mask[0, 0] == 1.0 and with_mask[1, 1] == 1.0
    assert without[0, 0] < 1.0 and without[1, 1] < 1.0
    pairs = match_lineages(ma, mb, shared)
    assert [(p.a, p.b, p.jaccard) for p in pairs] == [(0, 0, 1.0), (1, 1, 1.0)]


def test_a_mask_limits_the_cells_compared(bifurcation: Prepared) -> None:
    r = run_method(bifurcation.ctx, "paga_dpt", ClusterRoot("true_label", "trunk_early"))
    shared = np.zeros(r.n_obs, dtype=bool)
    shared[: r.n_obs // 2] = True
    c = compare(r, r, shared)
    assert c.n_common == int((np.isfinite(r.pseudotime) & shared).sum())
    assert np.all(np.isnan(c.rank_difference[~shared]))


def _today(
    a: TrajectoryResult, b: TrajectoryResult
) -> tuple[int, float | None, list[tuple[int, int, float]], np.ndarray]:
    """compare() as it was before the shared-cell mask, frozen as the reference."""
    both = np.isfinite(a.pseudotime) & np.isfinite(b.pseudotime)
    n = int(both.sum())
    rho = float(spearmanr(a.pseudotime[both], b.pseudotime[both]).statistic) if n >= 3 else None
    diff = np.full(a.n_obs, np.nan, dtype=np.float32)
    if n:
        diff[both] = (np.abs(rankdata(a.pseudotime[both]) - rankdata(b.pseudotime[both])) / n).astype(
            np.float32
        )
    ma, mb = a.membership, b.membership
    jac = np.zeros((ma.shape[0], mb.shape[0]))
    for i in range(ma.shape[0]):
        for j in range(mb.shape[0]):
            union = np.count_nonzero(ma[i] | mb[j])
            jac[i, j] = np.count_nonzero(ma[i] & mb[j]) / union if union else 0.0
    pairs = []
    if jac.size:
        rows, cols = linear_sum_assignment(-jac)
        pairs = [(int(r), int(c), float(jac[r, c])) for r, c in zip(rows, cols, strict=True)]
    return n, rho, pairs, diff


@pytest.mark.parametrize("other_root", ["trunk_early", "A_late"])
def test_without_a_mask_compare_is_unchanged(bifurcation: Prepared, other_root: str) -> None:
    a = run_method(bifurcation.ctx, "paga_dpt", ClusterRoot("true_label", "trunk_early"))
    b = run_method(bifurcation.ctx, "paga_dpt", ClusterRoot("true_label", other_root))
    n, rho, pairs, diff = _today(a, b)
    c = compare(a, b)
    assert c.n_common == n
    assert c.spearman == rho
    assert [(p.a, p.b, p.jaccard) for p in c.pairs] == pairs
    np.testing.assert_array_equal(c.rank_difference, diff)


def test_compare_counts_only_cells_finite_in_both(bifurcation: Prepared) -> None:
    from dataclasses import replace

    a = run_method(bifurcation.ctx, "paga_dpt", ClusterRoot("true_label", "trunk_early"))
    pa, pb = a.pseudotime.copy(), a.pseudotime[::-1].copy()
    pa[:5] = np.nan
    pb[3:8] = np.nan
    ra, rb = replace(a, pseudotime=pa), replace(a, pseudotime=pb)
    both = np.isfinite(pa) & np.isfinite(pb)
    c = compare(ra, rb)
    assert c.n_common == int(both.sum()) < int(np.isfinite(pa).sum())
    assert c.spearman == pytest.approx(float(spearmanr(pa[both], pb[both]).statistic))
