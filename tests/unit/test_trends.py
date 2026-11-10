import numpy as np
import pytest
from scipy import sparse

from lumencell.engine.compare import compare
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.trends import N_BINS, bin_lineage, gene_trends, temporal_heatmap
from tests.conftest import Prepared

N = 4000
rng = np.random.default_rng(0)
PT = rng.random(N)
MEMBERS = np.ones(N, dtype=bool)
MEMBERS[:500] = False
# gene 0 rises with pseudotime, gene 1 is noise, gene 2 peaks mid-way, gene 3 is off
X = sparse.csc_matrix(
    np.column_stack(
        [
            3 * PT + rng.normal(0, 0.1, N),
            rng.random(N),
            2 * np.exp(-0.5 * ((PT - 0.4) / 0.08) ** 2) + rng.normal(0, 0.05, N),
            np.zeros(N),
        ]
    ).clip(0)
)


def test_bins_average_their_cells() -> None:
    b = bin_lineage(PT, MEMBERS)
    assert b.counts.sum() == MEMBERS.sum()
    np.testing.assert_allclose(np.asarray(b.weights.sum(axis=1)).ravel()[b.counts > 0], 1.0)


def test_binned_means_match_brute_force() -> None:
    t = gene_trends(X, [0, 2], PT, MEMBERS)
    cells = np.flatnonzero(MEMBERS)
    bins = np.minimum((PT[cells] * N_BINS).astype(int), N_BINS - 1)
    dense = X.toarray()
    for k, gene in enumerate([0, 2]):
        brute = np.array([dense[cells[bins == i], gene].mean() for i in range(N_BINS)])
        np.testing.assert_allclose(t.means[:, k], brute, atol=1e-5)


def test_r2_separates_trends_from_noise() -> None:
    t = gene_trends(X, [0, 1, 2], PT, MEMBERS)
    assert t.r2[0] > 0.9 and t.r2[2] > 0.8
    assert t.r2[1] < 0.05
    assert np.all(np.diff(t.smooth[5:-5, 0]) > 0), "a rising gene's trend rises"


def test_heatmap_picks_temporal_genes_ordered_by_peak() -> None:
    h = temporal_heatmap(X, PT, MEMBERS, top_n=2)
    assert h.genes == [2, 0], (
        "the mid-way peak comes before the late riser; noise and the silent gene are left out"
    )
    assert h.values.shape == (2, N_BINS)
    assert h.values.min() >= 0 and h.values.max() <= 1


def test_heatmap_for_named_genes() -> None:
    h = temporal_heatmap(X, PT, MEMBERS, genes=[0, 2])
    assert sorted(h.genes) == [0, 2]


def test_comparing_a_run_with_itself(bifurcation: Prepared) -> None:
    r = run_method(bifurcation.ctx, "paga_dpt", ClusterRoot("true_label", "trunk_early"))
    c = compare(r, r)
    assert c.spearman == pytest.approx(1.0)
    assert all(p.jaccard == pytest.approx(1.0) and p.a == p.b for p in c.pairs)
    assert np.nanmax(c.rank_difference) == 0.0


def test_comparison_sees_disagreement(bifurcation: Prepared) -> None:
    r = run_method(bifurcation.ctx, "paga_dpt", ClusterRoot("true_label", "trunk_early"))
    flipped = run_method(bifurcation.ctx, "paga_dpt", ClusterRoot("true_label", "A_late"))
    c = compare(r, flipped)
    assert c.spearman is not None and c.spearman < 0
