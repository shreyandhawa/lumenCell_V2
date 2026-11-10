"""Method-independent lineage post-processing: names and divergence points.

Computing these the same way for every method is what makes results from
different methods comparable side by side.
"""

from itertools import combinations

import numpy as np
import numpy.typing as npt

from lumencell.engine.context import LEIDEN_KEY, InferenceContext
from lumencell.engine.schema import Divergence, TrajectoryResult

TIP_QUANTILE = 0.9
MIN_PURITY = 0.5


def _display(label: str, key: str) -> str:
    return f"Cluster {label}" if key == LEIDEN_KEY else label


def name_lineages(ctx: InferenceContext, result: TrajectoryResult, label_key: str | None = None) -> None:
    """Name each lineage after the majority label of its tip cells (in place).

    Tip cells are a lineage's exclusive members at or above the 90th percentile of
    its pseudotime. Below 50% purity the name says the tip is mixed. Duplicate
    names get a numeric suffix.
    """
    key = label_key or ctx.label_key
    labels = ctx.labels(key)
    pt = result.pseudotime
    exclusive = result.membership.sum(axis=0) == 1
    used: dict[str, int] = {}
    for i, lin in enumerate(result.lineages):
        members = result.membership[i] & np.isfinite(pt)
        pool = members & exclusive if (members & exclusive).any() else members
        if not pool.any():
            lin.terminal_label, lin.terminal_purity = None, None
            base = f"Lineage {i + 1}"
        else:
            cut = np.quantile(pt[pool], TIP_QUANTILE)
            tip = pool & (pt >= cut)
            values, counts = np.unique(labels[tip], return_counts=True)
            order = np.argsort(-counts, kind="stable")
            top, share = str(values[order[0]]), float(counts[order[0]] / counts.sum())
            lin.terminal_label, lin.terminal_purity = top, share
            if share >= MIN_PURITY:
                base = _display(top, key)
            else:
                parts = ", ".join(
                    f"{_display(str(values[j]), key)} {counts[j] / counts.sum():.0%}" for j in order[:2]
                )
                base = f"Lineage {i + 1} (mixed: {parts})"
        used[base] = used.get(base, 0) + 1
        name = base if used[base] == 1 else f"{base} ({used[base]})"
        lin.name = name
        lin.auto_name = name
        lin.n_cells = int(members.sum())


def divergence(membership: npt.NDArray[np.bool_], pseudotime: npt.NDArray[np.floating]) -> list[Divergence]:
    """For each pair of lineages, the pseudotime at which they stop sharing cells.

    Uses the 95th percentile of the pseudotime of cells on both lineages, so a few
    stray shared cells deep in a branch do not move the fork.
    """
    out: list[Divergence] = []
    finite = np.isfinite(pseudotime)
    for a, b in combinations(range(membership.shape[0]), 2):
        shared = membership[a] & membership[b] & finite
        point = float(np.quantile(pseudotime[shared], 0.95)) if shared.any() else 0.0
        out.append(Divergence(a, b, point))
    return out
