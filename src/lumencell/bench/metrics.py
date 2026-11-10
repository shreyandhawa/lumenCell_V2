"""Accuracy metrics for an inferred trajectory against a :class:`Truth`.

The metrics take plain arrays (per-cell pseudotime and a lineage membership
matrix), so they score any method and never depend on the engine's internals.

* ``pt_spearman_global``: rank correlation of pseudotime with the true order.
* ``pt_spearman_lineage``: the same, within each expected lineage; the gate
  uses the worst lineage, so one reversed branch cannot hide behind a good one.
* ``branch_accuracy``: for cells that truly lie on exactly one lineage and whose
  branch is identifiable (``Truth.resolvable``), the share placed exclusively on
  the matching inferred lineage. Chance level is ``1 / n_lineages``.
  ``branch_accuracy_all`` is the same over every post-fork cell, including
  those just past a fork that no method can place (reported, not gated).
* ``lineage_recall`` / ``n_spurious``: how many true lineages were found, and
  how many inferred lineages correspond to none.

Matching: inferred lineages are paired with true ones by the Hungarian
algorithm on Jaccard overlap, computed only over cells that truly lie past a
fork (and are identifiable). Shared trunk cells say nothing about which branch a lineage follows, so
including them would stop a method that found one clean lineage from ever
matching it. A pair counts as "found" at Jaccard >= 0.6, so an inferred lineage
that merges two equal-sized true branches (Jaccard 0.5 with each) finds neither.
* ``terminal_share`` (real datasets): for each expected lineage, the share of its
  terminal-stage cells held by the single inferred lineage that holds most of
  them; the worst lineage is reported. Low when a terminal population is split
  across several inferred lineages or left unassigned.
* ``root_pt_median``: median pseudotime of the true root cells.
  ``root_pt_excess`` subtracts their true median when the truth order is a real
  pseudotime: the root population spans a stretch of time, so even a perfect
  method has a non-zero median there.
"""

import math
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr

from lumencell.bench.truth import BoolArray, FloatArray, Truth

MIN_MATCH_JACCARD = 0.6


def spearman(a: FloatArray, b: FloatArray) -> float:
    """Spearman correlation over cells where both values are finite; NaN if < 3 cells."""
    ok = np.isfinite(a) & np.isfinite(b)
    if int(ok.sum()) < 3:
        return math.nan
    x, y = a[ok], b[ok]
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return math.nan
    return float(spearmanr(x, y).statistic)


def pt_spearman_global(pseudotime: FloatArray, truth: Truth) -> float:
    return spearman(pseudotime, truth.order)


def pt_spearman_lineages(pseudotime: FloatArray, truth: Truth) -> dict[str, float]:
    return {name: spearman(pseudotime[mask], truth.order[mask]) for name, mask in truth.lineages.items()}


def _exclusive(membership: BoolArray) -> npt.NDArray[np.int64]:
    out = np.full(membership.shape[1], -1, dtype=np.int64)
    if membership.shape[0] == 0:
        return out
    single = membership.sum(axis=0) == 1
    out[single] = np.argmax(membership[:, single], axis=0)
    return out


@dataclass(frozen=True)
class LineageMatch:
    inferred: int
    true_name: str
    jaccard: float


def match_lineages(membership: BoolArray, truth: Truth) -> list[LineageMatch]:
    """Hungarian matching of inferred lineages to true ones (see module docstring).

    Returns one entry per matched pair, including weak matches; callers decide
    what counts as "found" using ``MIN_MATCH_JACCARD``.
    """
    true_ex = truth.exclusive()
    post = truth.branch_cells()
    n_pred, n_true = membership.shape[0], len(truth.names)
    if n_pred == 0 or n_true == 0:
        return []
    jac = np.zeros((n_pred, n_true))
    for i in range(n_pred):
        p = membership[i] & post
        for j in range(n_true):
            t = (true_ex == j) & post
            union = int(np.count_nonzero(p | t))
            jac[i, j] = np.count_nonzero(p & t) / union if union else 0.0
    rows, cols = linear_sum_assignment(-jac)
    return [LineageMatch(int(r), truth.names[c], float(jac[r, c])) for r, c in zip(rows, cols, strict=True)]


def branch_accuracy(
    membership: BoolArray, truth: Truth, matches: list[LineageMatch], resolvable_only: bool = True
) -> float:
    """Share of truly post-fork cells placed exclusively on their matched lineage.

    NaN when the truth has fewer than two lineages (there is no branch to get right).
    """
    if len(truth.names) < 2:
        return math.nan
    true_ex = truth.exclusive()
    post = truth.branch_cells() if resolvable_only else true_ex >= 0
    if not post.any():
        return math.nan
    pred_to_true = np.full(membership.shape[0], -1, dtype=np.int64)
    for m in matches:
        pred_to_true[m.inferred] = truth.names.index(m.true_name)
    pred_ex = _exclusive(membership)
    mapped = np.where(pred_ex >= 0, pred_to_true[np.maximum(pred_ex, 0)], -1)
    return float(np.mean(mapped[post] == true_ex[post]))


@dataclass(frozen=True)
class Scores:
    pt_spearman_global: float
    pt_spearman_lineage: float  # worst lineage
    pt_spearman_by_lineage: dict[str, float]
    branch_accuracy: float
    branch_accuracy_all: float
    lineage_recall: float
    n_spurious: int
    n_inferred_lineages: int
    root_pt_median: float
    root_pt_excess: float
    terminal_share: float
    coverage: float  # share of cells with a finite pseudotime
    matches: list[LineageMatch]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def score(pseudotime: FloatArray, membership: BoolArray, truth: Truth) -> Scores:
    pseudotime = np.asarray(pseudotime, dtype=np.float64)
    membership = np.asarray(membership, dtype=bool)
    if pseudotime.shape != (truth.n_obs,):
        raise ValueError("pseudotime must have one value per cell")
    if membership.ndim != 2 or membership.shape[1] != truth.n_obs:
        raise ValueError("membership must be (n_lineages, n_obs)")

    by_lineage = pt_spearman_lineages(pseudotime, truth)
    finite = [v for v in by_lineage.values() if not math.isnan(v)]
    matches = match_lineages(membership, truth)
    found = [m for m in matches if m.jaccard >= MIN_MATCH_JACCARD]
    n_true = len(truth.names)

    root_pt = pseudotime[truth.root]
    root_pt = root_pt[np.isfinite(root_pt)]
    root_median = float(np.median(root_pt)) if root_pt.size else math.nan
    true_root_median = float(np.median(truth.order[truth.root]))
    excess = root_median - true_root_median if truth.meta.get("order_is_pseudotime") else math.nan
    shares = []
    for mask in (truth.terminals or {}).values():
        if mask.any() and membership.shape[0]:
            shares.append(float((membership[:, mask].sum(axis=1) / mask.sum()).max()))
        elif mask.any():
            shares.append(0.0)

    return Scores(
        pt_spearman_global=pt_spearman_global(pseudotime, truth),
        pt_spearman_lineage=min(finite) if finite else math.nan,
        pt_spearman_by_lineage=by_lineage,
        branch_accuracy=branch_accuracy(membership, truth, matches),
        branch_accuracy_all=branch_accuracy(membership, truth, matches, resolvable_only=False),
        lineage_recall=len(found) / n_true if n_true else math.nan,
        n_spurious=membership.shape[0] - len(found),
        n_inferred_lineages=int(membership.shape[0]),
        root_pt_median=root_median,
        root_pt_excess=excess,
        terminal_share=min(shares) if shares else math.nan,
        coverage=float(np.mean(np.isfinite(pseudotime))),
        matches=matches,
    )
