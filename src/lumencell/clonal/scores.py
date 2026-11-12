"""Clonal scores of a run's lineages: entropy and purity, η against three nulls, a lineage-free
k-means baseline, and lineage-pair coupling.

The definitions are the ones declared in ``benchmarks/thresholds.toml`` (``[clonal]``) before any
score was computed; spec §4.3 explains them. In short, over the scored cells (``score_clonal``):

- Each clone k has a distribution p_k over the L lineages, from its cells' (possibly fractional)
  counts, with normalised entropy H_k = -sum_l p_kl log p_kl / log L. The mean entropy is
  H = sum_k n_k H_k / sum_k n_k, and the purity is P = 1 - H, after STARTRAC's transition index
  (Zhang et al., Nature 2018, doi:10.1038/s41586-018-0694-x) and LRT's lineage-bias test
  (Xie et al., PLoS Comput Biol 2023, doi:10.1371/journal.pcbi.1011300).
- η = 1 - H_obs / E_null[H]: 0 is no better than chance, 1 is every clone in one lineage. The
  log L cancels, as in the chance-corrected measures of Vinh, Epps & Bailey (JMLR 11(95):
  2837-2854, 2010).
- Each null shuffles clone identities among the scored cells, B times, within strata: N1 donor,
  N2 donor x pseudotime decile, N3 the declared strata (with donor). p = (1 + #{H_b <= H_obs}) /
  (B + 1), one-sided, never zero (Phipson & Smyth, Stat Appl Genet Mol Biol 2010,
  doi:10.2202/1544-6115.1585).
- The lineage-free baseline (``kmeans_partition``) is scored on the same cells and the same N1
  shuffles: η_base, Δη = η - η_base, and a paired permutation p for Δη.
- Coupling (reported, never gated): the cosine of each lineage pair's raw count columns
  against N1, with z, p and a Benjamini-Hochberg p_adj (``_coupling``, after CoSpar).

Clone sharing says two states are related, never which came first: nothing here claims direction.
"""

from collections.abc import Iterator
from dataclasses import dataclass, replace
from typing import Any, Literal

import numpy as np
import numpy.typing as npt
from scipy.cluster.vq import kmeans2
from scipy.special import entr
from scipy.stats import false_discovery_control, rankdata

from lumencell.clonal.labels import BoolArray, CloneCodes, FloatArray, IntArray
from lumencell.trust.constants import (
    CLONAL_SEED,
    MIN_CLONE_SIZE,
    MIN_CLONES,
    N_PERMUTATIONS,
    PSEUDOTIME_BINS,
)

Mode = Literal["own_cells", "oracle_cells"]
NullName = Literal["N1", "N2", "N3"]
ScoreStatus = Literal["ok", "single_lineage", "too_few_clones"]

F64 = npt.NDArray[np.float64]

_NAN = float("nan")
# Slack on every comparison with an observed statistic, so that rounding never decides a count.
_SLACK = 1e-12
# Shuffles are scored in blocks of about this many array elements.
_BLOCK_ELEMENTS = 1 << 18
_NO_PAIRS = np.empty((0, 2), dtype=np.int64)


@dataclass(frozen=True)
class NullComparison:
    """The observed mean entropy against one null: the null's mean and sd (ddof=1) of H, η, p."""

    null: NullName
    mean_entropy_mean: float
    mean_entropy_sd: float
    eta: float
    p: float


@dataclass(frozen=True)
class Coupling:
    """Lineages a < b: cosine of their raw count columns against N1, with z, p and BH p_adj."""

    a: int
    b: int
    observed: float
    null_mean: float
    z: float
    p: float
    p_adj: float


@dataclass(frozen=True)
class ClonalResult:
    """One clonal score. Undefined values are NaN, never a stand-in number.

    ``nulls`` are N1, N2, then N3 when strata exist. ``eta`` is against N1; ``p`` is the largest
    null p. ``n_donors`` and ``n_strata`` count values among the scored cells. Under
    ``single_lineage`` every float is NaN and ``nulls`` and ``coupling`` are empty; so too when no
    clone is scored. ``too_few_clones`` still reports every value.
    """

    status: ScoreStatus
    mode: Mode
    n_lineages: int
    n_cells: int
    n_clones: int
    effective_clones: float
    n_donors: int
    n_strata: int | None
    min_clone_size: int
    n_permutations: int
    seed: int
    purity: float
    mean_entropy: float
    nulls: tuple[NullComparison, ...]
    eta: float
    p: float
    eta_base: float
    delta_eta: float
    p_delta: float
    clone_weighted_purity: float
    purity_without_largest: float
    coupling: tuple[Coupling, ...]


def baseline_space(pca: npt.NDArray[np.floating[Any]], n_pcs: int) -> F64:
    """The PCs the kNN graph uses (callers pass ``PrepParams.neighbors_pcs``), as float64."""
    return pca[:, : min(n_pcs, pca.shape[1])].astype(np.float64)


def kmeans_partition(space: F64, k: int, seed: int) -> IntArray:
    """A lineage-free partition of every cell: seeded k-means with k = L.

    T-cell clones are transcriptionally coherent (scirpy, Sturm et al., Bioinformatics 2020,
    doi:10.1093/bioinformatics/btaa611; CoNGA, Schattgen et al., Nat Biotechnol 2022,
    doi:10.1038/s41587-021-00989-2), so any expression partition beats the shuffle nulls. Scoring
    this partition on the same cells says how much of η the lineages add beyond that.
    """
    _, labels = kmeans2(space, k, minit="++", seed=seed)
    return np.asarray(labels, dtype=np.int64)


def score_clonal(
    membership: BoolArray,
    codes: CloneCodes,
    pseudotime: FloatArray,
    baseline: IntArray,
    *,
    oracle_cells: BoolArray | None = None,
    min_clone_size: int = MIN_CLONE_SIZE,
    min_clones: int = MIN_CLONES,
    n_permutations: int = N_PERMUTATIONS,
    seed: int = CLONAL_SEED,
    n_bins: int = PSEUDOTIME_BINS,
) -> ClonalResult:
    """Score the ``(L, n)`` lineage masks ``membership`` against the clones in ``codes``.

    ``oracle_cells=None`` is own-cells mode: the candidates are the cells in exactly one lineage.
    Otherwise it is oracle-cells mode: the candidates are ``oracle_cells``, and each counts in the
    method's lineages as spec §3 says (one lineage 1; k lineages 1/k each; none 1/L each). The
    scored cells are the candidates with a clone that has at least ``min_clone_size`` candidates.
    ``pseudotime`` defines N2's deciles (NaN: a stratum of its own); ``baseline`` is the
    lineage-free partition (``kmeans_partition``), scored on the same cells.
    """
    membership, pseudotime, baseline = _checked(
        membership,
        codes,
        pseudotime,
        baseline,
        oracle_cells,
        min_clone_size,
        min_clones,
        n_permutations,
        n_bins,
    )
    n_lineages = membership.shape[0]
    # own-cells mode: the cells in exactly one lineage; oracle-cells mode: the oracle's cells.
    candidates = membership.sum(axis=0) == 1 if oracle_cells is None else np.asarray(oracle_cells, bool)
    cells, clone, sizes = _scored_cells(codes.clone, candidates, min_clone_size)
    n_clones = int(sizes.size)
    status: ScoreStatus = (
        "single_lineage" if n_lineages < 2 else "too_few_clones" if n_clones < min_clones else "ok"
    )
    undefined = ClonalResult(
        status=status,
        mode="own_cells" if oracle_cells is None else "oracle_cells",
        n_lineages=n_lineages,
        n_cells=int(cells.size),
        n_clones=n_clones,
        effective_clones=_NAN,
        n_donors=_n_distinct(codes.donor[cells]),
        n_strata=None if codes.stratum is None else _n_distinct(codes.stratum[cells]),
        min_clone_size=min_clone_size,
        n_permutations=n_permutations,
        seed=seed,
        purity=_NAN,
        mean_entropy=_NAN,
        nulls=(),
        eta=_NAN,
        p=_NAN,
        eta_base=_NAN,
        delta_eta=_NAN,
        p_delta=_NAN,
        clone_weighted_purity=_NAN,
        purity_without_largest=_NAN,
        coupling=(),
    )
    if n_lineages < 2 or n_clones == 0:
        return undefined

    method = _Partition(_lineage_weights(membership, cells), n_lineages)
    base_codes = np.unique(baseline[cells], return_inverse=True)[1].reshape(-1)
    n_base = int(base_codes.max()) + 1
    base = _Partition(np.eye(n_base)[base_codes], max(n_lineages, n_base))

    observed = _counts(clone[None, :], method.weights, n_clones)[0]
    h = method.clone_entropy(observed, sizes)
    h_obs = float(_mean_entropy(h, sizes))
    h_base = float(
        _mean_entropy(base.clone_entropy(_counts(clone[None, :], base.weights, n_clones)[0], sizes), sizes)
    )

    # Coupling is defined for pairs whose count columns are both non-zero; each cell keeps its
    # weights under a shuffle, so a column's total is the same in every shuffle.
    pairs = np.array([(a, b) for a in range(n_lineages) for b in range(a + 1, n_lineages)], dtype=np.int64)
    totals = observed.sum(axis=0)
    defined = (totals[pairs[:, 0]] > 0) & (totals[pairs[:, 1]] > 0)

    donor = codes.donor[cells]
    (h_n1, h_base_n1), cosine_n1 = _draw(
        clone, donor, np.random.default_rng([seed, 1]), n_permutations, sizes, [method, base], pairs[defined]
    )
    stratified: list[tuple[NullName, IntArray]] = [
        ("N2", donor * (n_bins + 1) + _pseudotime_bins(pseudotime[cells], n_bins))
    ]
    if codes.stratum is not None:
        stratified.append(("N3", codes.stratum[cells]))
    nulls = [_compare("N1", h_obs, h_n1)]
    for i, (name, groups) in enumerate(stratified, start=2):
        (h_null,), _ = _draw(
            clone, groups, np.random.default_rng([seed, i]), n_permutations, sizes, [method], _NO_PAIRS
        )
        nulls.append(_compare(name, h_obs, h_null))

    eta = nulls[0].eta
    eta_base = _eta(h_base, h_base_n1)
    delta_eta = eta - eta_base
    return replace(
        undefined,
        effective_clones=1.0 / float(np.sum((sizes / sizes.sum()) ** 2)),
        purity=1.0 - h_obs,
        mean_entropy=h_obs,
        nulls=tuple(nulls),
        eta=eta,
        p=max(n.p for n in nulls),
        eta_base=eta_base,
        delta_eta=delta_eta,
        p_delta=_p_delta(delta_eta, h_n1, h_base_n1),
        clone_weighted_purity=1.0 - float(np.mean(h)),
        purity_without_largest=_purity_without_largest(h, sizes),
        coupling=_coupling(pairs, defined, _cosines(observed, pairs[defined]), cosine_n1),
    )


@dataclass(frozen=True)
class _Partition:
    """A partition of the scored cells to score clones against: lineages, or the baseline.

    ``weights`` is ``(n_cells, C)``, each cell's share of each column. Entropies are divided by
    ``log(n_columns)``; for the baseline any constant serves, because η is a ratio.
    """

    weights: F64
    n_columns: int

    def clone_entropy(self, counts: F64, sizes: IntArray) -> F64:
        """H_k of each clone, ``(..., K)`` from ``(..., K, C)`` counts, clipped to [0, 1] against rounding."""
        plogp: F64 = entr(counts / sizes[:, None])
        h: F64 = np.clip(plogp.sum(axis=-1) / np.log(self.n_columns), 0.0, 1.0)
        return h


def _checked(
    membership: BoolArray,
    codes: CloneCodes,
    pseudotime: FloatArray,
    baseline: IntArray,
    oracle_cells: BoolArray | None,
    min_clone_size: int,
    min_clones: int,
    n_permutations: int,
    n_bins: int,
) -> tuple[BoolArray, F64, IntArray]:
    masks = np.asarray(membership, dtype=bool)
    n = codes.clone.shape[0]
    if masks.ndim != 2 or masks.shape[1] != n:
        raise ValueError(f"membership must be (n_lineages, {n}), one column per cell; got {masks.shape}")
    shaped = {"pseudotime": pseudotime, "baseline": baseline, "oracle_cells": oracle_cells}
    for name, values in shaped.items():
        if values is not None and np.shape(values) != (n,):
            raise ValueError(f"{name} must have shape ({n},); got {np.shape(values)}")
    if min_clone_size < 2:
        raise ValueError("min_clone_size must be at least 2: single-cell clones never contribute")
    for name, value in (("min_clones", min_clones), ("n_permutations", n_permutations), ("n_bins", n_bins)):
        if value < 1:
            raise ValueError(f"{name} must be at least 1; got {value}")
    return masks, np.asarray(pseudotime, dtype=np.float64), np.asarray(baseline, dtype=np.int64)


def _scored_cells(
    clone: IntArray, candidates: BoolArray, min_size: int
) -> tuple[IntArray, IntArray, IntArray]:
    """Scored cells (indices, in their original order), their dense clone ids, and clone sizes."""
    labelled = candidates & (clone >= 0)
    candidate_size = np.bincount(clone[labelled], minlength=int(clone.max(initial=-1)) + 1)
    scored = labelled.copy()
    scored[labelled] = candidate_size[clone[labelled]] >= min_size
    cells = np.flatnonzero(scored).astype(np.int64)
    _, dense, sizes = np.unique(clone[cells], return_inverse=True, return_counts=True)
    return cells, dense.reshape(-1).astype(np.int64), sizes.astype(np.int64)


def _n_distinct(values: IntArray) -> int:
    return int(np.unique(values).size)


def _lineage_weights(membership: BoolArray, cells: IntArray) -> F64:
    """Each scored cell's share of each lineage: 1/k in each of its k lineages, 1/L each in none.

    In own-cells mode every scored cell is in one lineage, so each row is one-hot.
    """
    inside = membership[:, cells].T.astype(np.float64)
    k = inside.sum(axis=1, keepdims=True)
    return np.where(k > 0, inside / np.maximum(k, 1.0), 1.0 / membership.shape[0])


def _pseudotime_bins(pseudotime: F64, n_bins: int) -> IntArray:
    """N2's pseudotime bin of each scored cell (deciles at the default n_bins = 10).

    Over the finite values, by average rank: min(n_bins - 1, floor(n_bins * (rank - 1) / n_finite)).
    Cells without a finite pseudotime form bin ``n_bins``, a bin of their own.
    """
    bins = np.full(pseudotime.shape, n_bins, dtype=np.int64)
    finite = np.isfinite(pseudotime)
    n_finite = int(finite.sum())
    if n_finite:
        rank = rankdata(pseudotime[finite], method="average")
        bins[finite] = np.minimum(n_bins - 1, np.floor(n_bins * (rank - 1) / n_finite)).astype(np.int64)
    return bins


def _counts(ids: IntArray, weights: F64, n_clones: int) -> F64:
    """``(rows, K, C)``: for each row of clone ids, the summed weights of each clone in each column."""
    rows, n = ids.shape
    n_columns = weights.shape[1]
    flat = (
        (np.arange(rows)[:, None] * n_clones + ids)[:, :, None] * n_columns + np.arange(n_columns)
    ).reshape(-1)
    summed = np.bincount(
        flat,
        weights=np.broadcast_to(weights, (rows, n, n_columns)).reshape(-1),
        minlength=rows * n_clones * n_columns,
    )
    return summed.reshape(rows, n_clones, n_columns)


def _mean_entropy(h: F64, sizes: IntArray) -> F64:
    """H = sum_k n_k H_k / sum_k n_k, over the last axis of ``h``."""
    return np.asarray(h @ sizes / sizes.sum(), dtype=np.float64)


def _shuffles(
    clone: IntArray, groups: IntArray, rng: np.random.Generator, n_permutations: int
) -> Iterator[IntArray]:
    """The B shuffles of clone ids within strata, as blocks of rows (one shuffle per row).

    Shuffle b draws ``u = rng.random(n_scored)`` over the scored cells in their original order, and
    permutes clone ids within strata by ``lexsort((u, groups))``: the cells of each stratum, in
    their original order, take the ids of that stratum's cells sorted by u. Drawing a block of
    rows at once takes the same stream as drawing the rows one by one.
    """
    n = clone.size
    grouped = np.argsort(groups, kind="stable")
    block = max(1, _BLOCK_ELEMENTS // n)
    for start in range(0, n_permutations, block):
        u = rng.random((min(block, n_permutations - start), n))
        order = np.lexsort((u, np.broadcast_to(groups, u.shape)), axis=-1)
        shuffled = np.empty(u.shape, dtype=np.int64)
        shuffled[:, grouped] = clone[order]
        yield shuffled


def _draw(
    clone: IntArray,
    groups: IntArray,
    rng: np.random.Generator,
    n_permutations: int,
    sizes: IntArray,
    partitions: list[_Partition],
    pairs: IntArray,
) -> tuple[list[F64], F64]:
    """Score every shuffle of one null: the ``(B,)`` mean entropy of each partition, and the
    ``(B, P)`` cosines of the first partition's columns for each pair. Every partition is scored
    on the same shuffled ids, which pairs the baseline comparison with the lineages'."""
    entropy: list[list[F64]] = [[] for _ in partitions]
    cosine: list[F64] = []
    for ids in _shuffles(clone, groups, rng, n_permutations):
        for j, partition in enumerate(partitions):
            counts = _counts(ids, partition.weights, sizes.size)
            entropy[j].append(_mean_entropy(partition.clone_entropy(counts, sizes), sizes))
            if j == 0 and pairs.size:
                cosine.append(_cosines(counts, pairs))
    cosines = np.concatenate(cosine) if cosine else np.empty((n_permutations, 0))
    return [np.concatenate(e) for e in entropy], cosines


def _eta(observed: float, null: F64) -> float:
    """η = 1 - observed / mean(null), NaN when the null mean is 0 (Vinh et al. 2010)."""
    mean = float(null.mean())
    return 1.0 - observed / mean if mean > 0 else _NAN


def _p(at_least_as_extreme: BoolArray) -> float:
    """(1 + #{as extreme as observed}) / (B + 1): never zero (Phipson & Smyth 2010)."""
    return (1.0 + float(at_least_as_extreme.sum())) / (at_least_as_extreme.size + 1)


def _sd(values: F64) -> float:
    return float(values.std(ddof=1)) if values.size > 1 else _NAN


def _compare(name: NullName, observed: float, null: F64) -> NullComparison:
    return NullComparison(
        null=name,
        mean_entropy_mean=float(null.mean()),
        mean_entropy_sd=_sd(null),
        eta=_eta(observed, null),
        p=_p(null <= observed + _SLACK),
    )


def _p_delta(delta_eta: float, entropy: F64, base_entropy: F64) -> float:
    """Paired permutation p for Δη over the N1 shuffles: (1 + #{η_b - η_base,b >= Δη}) / (B + 1).

    NaN whenever Δη is NaN: an undefined statistic is never given the most significant p.
    """
    if np.isnan(delta_eta):
        return _NAN
    eta_b = 1.0 - entropy / entropy.mean()
    eta_base_b = 1.0 - base_entropy / base_entropy.mean()
    return _p(eta_b - eta_base_b >= delta_eta - _SLACK)


def _purity_without_largest(h: F64, sizes: IntArray) -> float:
    """Purity with the clone of most scored cells left out; on a size tie, the purer one (lower
    H_k) is left out. NaN with one clone."""
    if sizes.size < 2:
        return _NAN
    left_out = np.lexsort((h, -sizes))[0]
    keep = np.arange(sizes.size) != left_out
    return float(1.0 - h[keep] @ sizes[keep] / sizes[keep].sum())


def _cosines(counts: F64, pairs: IntArray) -> F64:
    """``(..., P)`` cosine similarity of count columns a and b for each pair (a, b); both non-zero."""
    gram = np.einsum("...ka,...kb->...ab", counts, counts)
    a, b = pairs[:, 0], pairs[:, 1]
    result: F64 = gram[..., a, b] / np.sqrt(gram[..., a, a] * gram[..., b, b])
    return result


def _coupling(pairs: IntArray, defined: BoolArray, observed: F64, null: F64) -> tuple[Coupling, ...]:
    """Lineage-pair coupling against N1, after CoSpar (Wang et al., Nat Biotechnol 2022,
    doi:10.1038/s41587-022-01209-1) but on the raw scored-cell count columns, not CoSpar's
    normalised default.

    ``observed`` ``(P,)`` and ``null`` ``(B, P)`` hold the cosines of the defined pairs. z uses the
    N1 null sd and is NaN when that sd is 0 (to within the slack). A pair with an all-zero column
    has NaN fields and is left out of the Benjamini-Hochberg adjustment.
    """
    undefined = Coupling(0, 0, _NAN, _NAN, _NAN, _NAN, _NAN)
    if not defined.any():
        return tuple(replace(undefined, a=a, b=b) for a, b in pairs.tolist())
    mean = null.mean(axis=0)
    sd = null.std(axis=0, ddof=1) if null.shape[0] > 1 else np.full(null.shape[1], _NAN)
    spread = sd > _SLACK
    z = np.where(spread, (observed - mean) / np.where(spread, sd, 1.0), _NAN)
    p = (1.0 + (null >= observed - _SLACK).sum(axis=0)) / (null.shape[0] + 1)
    p_adj = false_discovery_control(p, method="bh")
    found = iter(range(observed.size))
    result = []
    for (a, b), ok in zip(pairs.tolist(), defined.tolist(), strict=True):
        if not ok:
            result.append(replace(undefined, a=a, b=b))
            continue
        j = next(found)
        result.append(
            Coupling(a, b, float(observed[j]), float(mean[j]), float(z[j]), float(p[j]), float(p_adj[j]))
        )
    return tuple(result)
