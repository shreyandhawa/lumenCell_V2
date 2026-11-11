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


# ---- the same results, block by block ----------------------------------------------------

import tracemalloc  # noqa: E402

from hypothesis import given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

import lumencell.engine.trends as trends_mod  # noqa: E402

# Today's implementation, frozen as the reference the block-wise one must reproduce.


def _ref_r2(x: sparse.csc_matrix, b: trends_mod.Binning, smooth: np.ndarray) -> np.ndarray:
    sub = x[b.cells]
    n = len(b.cells)
    mean = np.asarray(sub.mean(axis=0)).ravel()
    sq = np.asarray(sub.multiply(sub).mean(axis=0)).ravel()
    total = np.maximum(sq - mean**2, 1e-12)
    dev = np.nan_to_num(smooth, nan=0.0) - mean[None, :]
    explained = (b.counts[:, None] * dev**2).sum(axis=0) / max(n, 1)
    return np.clip(explained / total, 0.0, 1.0)


def _ref_gene_trends(expr: sparse.csc_matrix, genes: list[int], pt: np.ndarray, members: np.ndarray) -> tuple:  # type: ignore[type-arg]
    b = trends_mod.bin_lineage(pt, members)
    x = sparse.csc_matrix(expr[:, genes], dtype=np.float64)
    means = np.asarray((b.weights @ x).todense())
    means[b.counts == 0] = np.nan
    smooth = trends_mod._smooth(means, b.counts)
    return means, smooth, _ref_r2(x, b, smooth)


def _ref_heatmap(expr: sparse.csc_matrix, pt: np.ndarray, members: np.ndarray, top_n: int) -> tuple:  # type: ignore[type-arg]
    b = trends_mod.bin_lineage(pt, members)
    x = sparse.csc_matrix(expr, dtype=np.float64)
    means = np.asarray((b.weights @ x).todense())
    means[b.counts == 0] = np.nan
    smooth = trends_mod._smooth(means, b.counts)
    r2 = _ref_r2(x, b, smooth)
    expressed = np.nanmax(np.nan_to_num(smooth, nan=0.0), axis=0) >= trends_mod.MIN_MEAN
    candidates = np.flatnonzero(expressed)
    chosen = candidates[np.argsort(-r2[candidates], kind="stable")[:top_n]]
    block = np.nan_to_num(smooth[:, chosen], nan=0.0)
    lo, hi = block.min(axis=0), block.max(axis=0)
    scaled = (block - lo) / np.maximum(hi - lo, 1e-12)
    order = np.argsort(np.argmax(scaled, axis=0), kind="stable")
    return [int(chosen[i]) for i in order], scaled[:, order].T.astype(np.float32), r2[chosen][order]


def _fixture_inputs(bifurcation: Prepared) -> tuple[sparse.csc_matrix, np.ndarray, np.ndarray]:
    expr = sparse.csc_matrix(bifurcation.prep.artifacts["lognorm"], dtype=np.float32)
    pt = bifurcation.ds.obs["true_pseudotime"].to_numpy(dtype=np.float64)
    pt = (pt - pt.min()) / (pt.max() - pt.min())
    pt[::17] = np.nan  # cells without a pseudotime
    members = bifurcation.ds.obs["true_label"].astype(str).to_numpy() != "branch_b_late"
    return expr, pt, members


@pytest.mark.parametrize("block_nnz", [1, 1000, trends_mod.GENE_BLOCK_NNZ])
@pytest.mark.parametrize("which", ["fixture", "toy"])
def test_trends_and_heatmap_match_the_reference_implementation(
    bifurcation: Prepared, which: str, block_nnz: int
) -> None:
    expr, pt, members = _fixture_inputs(bifurcation) if which == "fixture" else (X, PT, MEMBERS)
    genes = list(range(min(expr.shape[1], 40)))
    means, smooth, r2 = _ref_gene_trends(expr, genes, pt, members)
    t = gene_trends(expr, genes, pt, members, block_nnz=block_nnz)
    np.testing.assert_allclose(t.means, means, rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(t.smooth, smooth, rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(t.r2, r2, rtol=1e-6, atol=1e-9)
    top_n = 50 if which == "fixture" else 2
    want_genes, want_values, want_r2 = _ref_heatmap(expr, pt, members, top_n)
    h = temporal_heatmap(expr, pt, members, top_n=top_n, block_nnz=block_nnz)
    assert h.genes == want_genes
    np.testing.assert_allclose(h.values, want_values, rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(h.r2, want_r2, rtol=1e-6, atol=1e-9)


@st.composite
def _problems(draw: st.DrawFn) -> tuple:  # type: ignore[type-arg]
    n_cells, n_genes = draw(st.integers(1, 60)), draw(st.integers(1, 30))
    rng = np.random.default_rng(draw(st.integers(0, 2**32 - 1)))
    dense = rng.random((n_cells, n_genes)).astype(np.float32)
    dense[rng.random((n_cells, n_genes)) < draw(st.floats(0, 1))] = 0
    pt = rng.random(n_cells)
    pt[rng.random(n_cells) < draw(st.floats(0, 0.5))] = np.nan  # cells without a pseudotime
    members = rng.random(n_cells) < draw(st.floats(0.05, 1))  # lineages with empty bins included
    genes = draw(st.none() | st.lists(st.integers(0, n_genes - 1), min_size=1, max_size=n_genes, unique=True))
    # Fewer bins than the smoothing window breaks _smooth (older than this test); every caller uses 50.
    n_bins = draw(st.integers(trends_mod.SMOOTH_WIDTH, 12))
    block_nnz = draw(st.integers(1, n_cells * n_genes + 1))
    return sparse.csc_matrix(dense), pt, members, genes, n_bins, block_nnz


@settings(max_examples=200, deadline=None)
@given(_problems())
def test_results_do_not_depend_on_the_block_size(problem: tuple) -> None:  # type: ignore[type-arg]
    x, pt, members, genes, n_bins, block_nnz = problem
    whole = 10**9
    h = temporal_heatmap(x, pt, members, top_n=5, genes=genes, n_bins=n_bins, block_nnz=block_nnz)
    ref = temporal_heatmap(x, pt, members, top_n=5, genes=genes, n_bins=n_bins, block_nnz=whole)
    assert h.genes == ref.genes
    np.testing.assert_allclose(h.values, ref.values, rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(h.r2, ref.r2, rtol=1e-9, atol=1e-12)
    cols = genes or list(range(x.shape[1]))
    t = gene_trends(x, cols, pt, members, n_bins=n_bins, block_nnz=block_nnz)
    t_ref = gene_trends(x, cols, pt, members, n_bins=n_bins, block_nnz=whole)
    np.testing.assert_allclose(t.smooth, t_ref.smooth, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(t.r2, t_ref.r2, rtol=1e-9, atol=1e-12)


@pytest.mark.slow
def test_an_uncached_heatmap_allocates_far_less_than_one_matrix_copy() -> None:
    rng = np.random.default_rng(0)
    n_cells, n_genes = 40_000, 2_000
    expr = sparse.random(n_cells, n_genes, density=0.3, format="csc", dtype=np.float32, random_state=rng)
    one_copy = expr.data.nbytes + expr.indices.nbytes
    pt = rng.random(n_cells)
    members = rng.random(n_cells) < 0.6
    tracemalloc.start()
    try:
        temporal_heatmap(expr, pt, members, top_n=50)
        _now, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 0.2 * one_copy, f"peak {peak / 1e6:.0f} MB against {one_copy / 1e6:.0f} MB for one copy"


@given(st.lists(st.integers(0, 6), min_size=1, max_size=40), st.integers(1, 20))
def test_gene_blocks_cover_every_gene_once_within_the_budget(nnz_per_gene: list[int], block_nnz: int) -> None:
    indptr = np.concatenate([[0], np.cumsum(nnz_per_gene)])
    x = sparse.csc_matrix(
        (np.ones(indptr[-1]), np.zeros(indptr[-1], dtype=np.int32), indptr), shape=(1, len(nnz_per_gene))
    )
    blocks = list(trends_mod._column_blocks(x, block_nnz))
    assert [c for b in blocks for c in range(b.start, b.stop)] == list(range(len(nnz_per_gene)))
    for b in blocks:
        assert b.stop - b.start == 1 or indptr[b.stop] - indptr[b.start] <= block_nnz


def test_a_matrix_stored_by_cell_is_refused() -> None:
    with pytest.raises(TypeError, match="by gene"):
        temporal_heatmap(X.tocsr(), PT, MEMBERS)


def test_an_empty_lineage_has_no_trend_rather_than_an_error() -> None:
    # It used to raise ZeroDivisionError (a 500). Now: no trend (null in the API) and r2 0.
    t = gene_trends(X, [0, 2], PT, np.zeros(N, dtype=bool))
    assert (t.counts == 0).all() and np.isnan(t.smooth).all()
    np.testing.assert_array_equal(t.r2, [0.0, 0.0])
