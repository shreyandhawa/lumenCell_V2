"""Ground truth for scoring a trajectory, in one shape for synthetic and real data.

Synthetic data knows each cell's exact pseudotime and branch. Real datasets only
know an ordering of annotated stages along each expected lineage. Both reduce to:

* ``order``: a per-cell value that should increase along the trajectory (true
  pseudotime, or a stage ordinal), NaN where unknown;
* ``lineages``: for each expected lineage, which cells lie on it (the shared
  trunk counts towards every lineage that passes through it);
* ``root``: which cells form the true starting population;
* ``resolvable`` (optional): cells whose branch can be told apart at all. Just
  past a fork, branches have not diverged yet and no method can separate them;
  branch metrics are scored only where the answer is identifiable.
"""

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import numpy.typing as npt

BoolArray = npt.NDArray[np.bool_]
FloatArray = npt.NDArray[np.floating[Any]]


@dataclass(frozen=True)
class Truth:
    order: FloatArray
    lineages: dict[str, BoolArray]
    root: BoolArray
    meta: dict[str, Any] = field(default_factory=dict)
    resolvable: BoolArray | None = None
    # Cells of each lineage's terminal stage (real datasets), for terminal_share.
    terminals: dict[str, BoolArray] | None = None
    # Per-cell clone label (NaN where the cell has no receptor) and donor. Used only to score
    # clonal agreement; both are set or neither is.
    clone_ids: npt.NDArray[np.object_] | None = None
    donors: npt.NDArray[np.object_] | None = None

    def __post_init__(self) -> None:
        n = self.order.shape[0]
        for name, mask in self.lineages.items():
            if mask.shape != (n,) or mask.dtype != np.bool_:
                raise ValueError(f"lineage {name!r} mask must be a bool array of length {n}")
        if self.root.shape != (n,) or self.root.dtype != np.bool_:
            raise ValueError("root must be a bool array matching order")
        if not self.root.any():
            raise ValueError("the true root population is empty")
        if self.resolvable is not None and (
            self.resolvable.shape != (n,) or self.resolvable.dtype != np.bool_
        ):
            raise ValueError("resolvable must be a bool array matching order")
        if (self.clone_ids is None) != (self.donors is None):
            raise ValueError("clone_ids and donors are set together or not at all")
        for name, column in (("clone_ids", self.clone_ids), ("donors", self.donors)):
            if column is not None and column.shape != (n,):
                raise ValueError(f"{name} must have shape ({n},)")

    @property
    def n_obs(self) -> int:
        return int(self.order.shape[0])

    @property
    def names(self) -> list[str]:
        return list(self.lineages)

    def membership(self) -> BoolArray:
        """(n_lineages, n_obs) membership matrix, in ``names`` order."""
        if not self.lineages:
            return np.zeros((0, self.n_obs), dtype=bool)
        return np.stack([self.lineages[n] for n in self.names])

    def exclusive(self) -> npt.NDArray[np.int64]:
        """Per cell: index of the only lineage it lies on, or -1 (shared or off-lineage)."""
        m = self.membership()
        out = np.full(self.n_obs, -1, dtype=np.int64)
        if m.shape[0] == 0:
            return out
        single = m.sum(axis=0) == 1
        out[single] = np.argmax(m[:, single], axis=0)
        return out

    def branch_cells(self) -> BoolArray:
        """Cells on exactly one true lineage whose branch is identifiable (see ``resolvable``)."""
        post = self.exclusive() >= 0
        return post & self.resolvable if self.resolvable is not None else post
