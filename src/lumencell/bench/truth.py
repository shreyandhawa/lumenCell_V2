"""Ground truth for scoring a trajectory, in one shape for synthetic and real data.

Synthetic data knows each cell's exact pseudotime and branch. Real datasets only
know an ordering of annotated stages along each expected lineage. Both reduce to:

* ``order``: a per-cell value that should increase along the trajectory (true
  pseudotime, or a stage ordinal), NaN where unknown;
* ``lineages``: for each expected lineage, which cells lie on it (the shared
  trunk counts towards every lineage that passes through it);
* ``root``: which cells form the true starting population.
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

    def __post_init__(self) -> None:
        n = self.order.shape[0]
        for name, mask in self.lineages.items():
            if mask.shape != (n,) or mask.dtype != np.bool_:
                raise ValueError(f"lineage {name!r} mask must be a bool array of length {n}")
        if self.root.shape != (n,) or self.root.dtype != np.bool_:
            raise ValueError("root must be a bool array matching order")
        if not self.root.any():
            raise ValueError("the true root population is empty")

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
