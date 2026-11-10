"""Gene expression along pseudotime, per lineage.

Cells on a lineage are put into equal-width pseudotime bins; a gene's trend is
its mean log-normalised expression per bin, lightly smoothed. One sparse matrix
product (bins x cells) @ (cells x genes) gives every gene's binned means at once,
so ranking all genes for the heatmap costs about as much as one gene.

``r2`` is the share of a gene's variance across the lineage's cells that the
smoothed trend explains: near 0 for genes unrelated to progression, high for
genes that switch on or off along it.
"""

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy import sparse

N_BINS = 50
MIN_MEAN = 0.05  # ignore genes this close to unexpressed when ranking
SMOOTH_WIDTH = 5


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


def _r2(x: sparse.csc_matrix, b: Binning, smooth: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Variance explained by the smoothed trend, per gene (columns of x)."""
    sub = x[b.cells]
    n = len(b.cells)
    mean = np.asarray(sub.mean(axis=0)).ravel()
    sq = np.asarray(sub.multiply(sub).mean(axis=0)).ravel()
    total = np.maximum(sq - mean**2, 1e-12)
    # The fitted value is constant within a bin, so the explained sum of squares
    # only needs bin counts (no cells x genes matrix, which would not fit at 100k cells).
    dev = np.nan_to_num(smooth, nan=0.0) - mean[None, :]
    explained = (b.counts[:, None] * dev**2).sum(axis=0) / max(n, 1)
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
) -> GeneTrends:
    b = bin_lineage(pseudotime, members, n_bins)
    x = sparse.csc_matrix(expr[:, genes], dtype=np.float64)
    means = np.asarray((b.weights @ x).todense())
    means[b.counts == 0] = np.nan
    smooth = _smooth(means, b.counts)
    return GeneTrends(b.centers, b.counts, means, smooth, _r2(x, b, smooth))


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
) -> Heatmap:
    """The genes that change most along a lineage (or the ones given), ordered by peak time."""
    b = bin_lineage(pseudotime, members, n_bins)
    x = sparse.csc_matrix(expr if genes is None else expr[:, genes], dtype=np.float64)
    means = np.asarray((b.weights @ x).todense())
    means[b.counts == 0] = np.nan
    smooth = _smooth(means, b.counts)
    r2 = _r2(x, b, smooth)
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
