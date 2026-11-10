"""How much two runs agree: pseudotime rank correlation and matched lineages.

Runs from different methods put lineages in different orders and give them
different names, so lineages are paired by the Hungarian algorithm on the
Jaccard overlap of their cells before they are compared.
"""

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata, spearmanr

from lumencell.engine.schema import TrajectoryResult


@dataclass(frozen=True)
class LineagePair:
    a: int
    b: int
    jaccard: float


@dataclass(frozen=True)
class RunComparison:
    n_common: int
    spearman: float | None
    pairs: list[LineagePair]
    rank_difference: npt.NDArray[np.float32]  # per cell |rank_a - rank_b| / n_common, NaN if not in both


def compare(a: TrajectoryResult, b: TrajectoryResult) -> RunComparison:
    if a.n_obs != b.n_obs:
        raise ValueError("the runs are on different datasets")
    both = np.isfinite(a.pseudotime) & np.isfinite(b.pseudotime)
    n = int(both.sum())
    rho = float(spearmanr(a.pseudotime[both], b.pseudotime[both]).statistic) if n >= 3 else None
    diff = np.full(a.n_obs, np.nan, dtype=np.float32)
    if n:
        ra = rankdata(a.pseudotime[both])
        rb = rankdata(b.pseudotime[both])
        diff[both] = (np.abs(ra - rb) / n).astype(np.float32)

    ma, mb = a.membership, b.membership
    jac = np.zeros((ma.shape[0], mb.shape[0]))
    for i in range(ma.shape[0]):
        for j in range(mb.shape[0]):
            union = np.count_nonzero(ma[i] | mb[j])
            jac[i, j] = np.count_nonzero(ma[i] & mb[j]) / union if union else 0.0
    pairs: list[LineagePair] = []
    if jac.size:
        rows, cols = linear_sum_assignment(-jac)
        pairs = [LineagePair(int(r), int(c), float(jac[r, c])) for r, c in zip(rows, cols, strict=True)]
    return RunComparison(n, rho, pairs, diff)
