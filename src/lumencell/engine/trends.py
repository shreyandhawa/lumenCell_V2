"""Gene expression along pseudotime, per lineage.

Cells on a lineage are put into equal-width pseudotime bins; a gene's trend is
its mean log-normalised expression per bin, lightly smoothed. One sparse matrix
product (bins x cells) @ (cells x genes) gives every gene's binned means at once,
so ranking all genes for the heatmap costs about as much as one gene.

``r2`` is the share of a gene's variance across the lineage's cells that the
smoothed trend explains: near 0 for genes unrelated to progression, high for
genes that switch on or off along it.

The matrix is worked through a block of genes at a time, each block converted
to float64 on its own. A heatmap over every gene therefore never holds more than
a block's worth of copies, however large the dataset.
"""

from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy import sparse

N_BINS = 50
MIN_MEAN = 0.05  # ignore genes this close to unexpressed when ranking
SMOOTH_WIDTH = 5
GENE_BLOCK_NNZ = 500_000  # nonzeros per block of genes: about 6 MB as float64 with its indices


@dataclass(frozen=True)
class Binning:
    weights: sparse.csr_matrix  # (bins, cells), each row averages the cells in that bin
    centers: npt.NDArray[np.float64]
    counts: npt.NDArray[np.int64]
    cells: npt.NDArray[np.int64]  # the lineage's cells with a pseudotime
    cell_bin: npt.NDArray[np.int64]  # bin of each of those cells


def bin_lineage(
    pseudotime: npt.NDArray[np.floating], members: npt.NDArray[np.bool_], n_bins: int = N_BINS
) -> Binning:
    cells = np.flatnonzero(members & np.isfinite(pseudotime))
    pt = np.clip(np.asarray(pseudotime, dtype=np.float64)[cells], 0.0, 1.0)
    cell_bin = np.minimum((pt * n_bins).astype(np.int64), n_bins - 1)
    counts = np.bincount(cell_bin, minlength=n_bins).astype(np.int64)
    vals = 1.0 / np.maximum(counts[cell_bin], 1)
    weights = sparse.csr_matrix((vals, (cell_bin, cells)), shape=(n_bins, len(members)))
    centers = (np.arange(n_bins) + 0.5) / n_bins
    return Binning(weights, centers, counts, cells, cell_bin)


def _smooth(means: npt.NDArray[np.float64], counts: npt.NDArray[np.int64]) -> npt.NDArray[np.float64]:
    """Count-weighted moving average along bins; empty bins stay NaN. Works on (bins, genes)."""
    w = counts.astype(np.float64)[:, None]
    filled = np.nan_to_num(means) * w
    kernel = np.ones(SMOOTH_WIDTH)
    num = np.apply_along_axis(lambda c: np.convolve(c, kernel, mode="same"), 0, filled)
    den = np.convolve(w[:, 0], kernel, mode="same")[:, None]
    out = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)
    out[counts == 0] = np.nan
    result: npt.NDArray[np.float64] = out
    return result


def _column_blocks(x: sparse.csc_matrix, block_nnz: int) -> Iterator[slice]:
    """Contiguous ranges of columns, each holding at most ``block_nnz`` nonzeros (or one column)."""
    indptr, n = x.indptr, x.shape[1]
    start = 0
    while start < n:
        stop = int(np.searchsorted(indptr, indptr[start] + block_nnz, side="right")) - 1
        stop = min(max(stop, start + 1), n)
        yield slice(start, stop)
        start = stop


@dataclass(frozen=True)
class _Binned:
    means: npt.NDArray[np.float64]  # (bins, genes), NaN for empty bins
    smooth: npt.NDArray[np.float64]  # (bins, genes)
    r2: npt.NDArray[np.float64]  # (genes,)


def _binned(x: sparse.csc_matrix, b: Binning, block_nnz: int) -> _Binned:
    """Binned means, the smoothed trend and r2 for every column of ``x``, a block at a time."""
    if x.format != "csc":  # a block of genes is a range of columns, which only CSC stores together
        raise TypeError(f"expression must be stored by gene (CSC), not as {x.format!r}")
    n_genes = x.shape[1]
    n = len(b.cells)
    means = np.empty((b.weights.shape[0], n_genes))
    mean = np.empty(n_genes)  # over the lineage's cells
    sq = np.empty(n_genes)
    # Weights 1/n on the lineage's cells, 0 elsewhere: a product with it is the lineage mean.
    per_cell = np.zeros(x.shape[0])
    per_cell[b.cells] = 1.0 / max(n, 1)
    for cols in _column_blocks(x, block_nnz):
        xb = sparse.csc_matrix(x[:, cols], dtype=np.float64)
        means[:, cols] = (b.weights @ xb).toarray()
        mean[cols] = xb.T @ per_cell
        sq[cols] = xb.multiply(xb).T @ per_cell
    means[b.counts == 0] = np.nan
    smooth = _smooth(means, b.counts)
    return _Binned(means, smooth, _r2(mean, sq, b, smooth))


def _r2(
    mean: npt.NDArray[np.float64],
    sq: npt.NDArray[np.float64],
    b: Binning,
    smooth: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Variance explained by the smoothed trend, per gene, from its lineage mean and mean square."""
    total = np.maximum(sq - mean**2, 1e-12)
    # The fitted value is constant within a bin, so the explained sum of squares
    # only needs bin counts (no cells x genes matrix, which would not fit at 100k cells).
    dev = np.nan_to_num(smooth, nan=0.0) - mean[None, :]
    explained = (b.counts[:, None] * dev**2).sum(axis=0) / max(len(b.cells), 1)
    out: npt.NDArray[np.float64] = np.clip(explained / total, 0.0, 1.0)
    return out


@dataclass(frozen=True)
class GeneTrends:
    centers: npt.NDArray[np.float64]
    counts: npt.NDArray[np.int64]
    means: npt.NDArray[np.float64]  # (bins, genes), raw binned means, NaN for empty bins
    smooth: npt.NDArray[np.float64]  # (bins, genes)
    r2: npt.NDArray[np.float64]  # (genes,)


def gene_trends(
    expr: sparse.csc_matrix,
    genes: list[int],
    pseudotime: npt.NDArray[np.floating],
    members: npt.NDArray[np.bool_],
    n_bins: int = N_BINS,
    block_nnz: int = GENE_BLOCK_NNZ,
) -> GeneTrends:
    b = bin_lineage(pseudotime, members, n_bins)
    t = _binned(expr[:, genes], b, block_nnz)
    return GeneTrends(b.centers, b.counts, t.means, t.smooth, t.r2)


@dataclass(frozen=True)
class Heatmap:
    genes: list[int]  # gene indices, ordered by when they peak
    centers: npt.NDArray[np.float64]
    values: npt.NDArray[np.float32]  # (genes, bins), each row scaled to [0, 1]
    r2: npt.NDArray[np.float64]


def temporal_heatmap(
    expr: sparse.csc_matrix,
    pseudotime: npt.NDArray[np.floating],
    members: npt.NDArray[np.bool_],
    top_n: int = 50,
    genes: list[int] | None = None,
    n_bins: int = N_BINS,
    block_nnz: int = GENE_BLOCK_NNZ,
) -> Heatmap:
    """The genes that change most along a lineage (or the ones given), ordered by peak time."""
    b = bin_lineage(pseudotime, members, n_bins)
    t = _binned(expr if genes is None else expr[:, genes], b, block_nnz)
    smooth, r2 = t.smooth, t.r2
    if genes is None:
        expressed = np.nanmax(np.nan_to_num(smooth, nan=0.0), axis=0) >= MIN_MEAN
        candidates = np.flatnonzero(expressed)
        chosen = candidates[np.argsort(-r2[candidates], kind="stable")[:top_n]]
        index = [int(i) for i in chosen]
    else:
        chosen = np.arange(len(genes))
        index = list(genes)
    block = np.nan_to_num(smooth[:, chosen], nan=0.0)  # (bins, k)
    lo, hi = block.min(axis=0), block.max(axis=0)
    scaled = (block - lo) / np.maximum(hi - lo, 1e-12)
    order = np.argsort(np.argmax(scaled, axis=0), kind="stable")
    return Heatmap(
        genes=[index[i] for i in order],
        centers=b.centers,
        values=scaled[:, order].T.astype(np.float32),
        r2=r2[chosen][order],
    )
