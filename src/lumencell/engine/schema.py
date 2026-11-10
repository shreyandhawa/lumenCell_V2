"""The common trajectory result every method adapter produces.

Modelled on dynverse's milestone network (Saelens et al. 2019), which was built
to represent the output of dozens of trajectory methods in one shape:

* **milestones** are points of the trajectory (root, forks, tips, or clusters);
* **segments** connect a parent milestone to a child milestone;
* **lineages** are root-to-terminal paths through the segments;
* **per-cell arrays** place each cell: pseudotime in [0, 1] (NaN = unassigned),
  lineage membership, the segment it sits on and its progress along it, and
  optionally fate probabilities.

``validate()`` enforces the invariants on every adapter's output, so the API,
the UI and the benchmark can rely on them instead of re-checking.
"""

from dataclasses import dataclass, field
from itertools import pairwise
from typing import Any, Literal

import numpy as np
import numpy.typing as npt
from scipy import sparse

SCHEMA_VERSION = 1

Severity = Literal["info", "warning", "error"]
MilestoneKind = Literal["root", "fork", "tip", "cluster"]
TopologyKind = Literal["tree", "fate_star"]

PRIMARY_SHARED = -1
PRIMARY_UNASSIGNED = -2

_PT_EPS = 1e-5
_FATE_EPS = 1e-3


class InvariantError(ValueError):
    """A method produced a result that breaks the schema."""

    def __init__(self, problems: list[str]) -> None:
        super().__init__("; ".join(problems))
        self.problems = problems


@dataclass(frozen=True)
class Notice:
    code: str
    severity: Severity
    message: str


@dataclass(frozen=True)
class Milestone:
    id: int
    kind: MilestoneKind
    label: str
    pseudotime: float


@dataclass(frozen=True)
class Segment:
    id: int
    source: int
    target: int
    length: float


@dataclass
class Lineage:
    id: int
    name: str
    auto_name: str
    terminal_milestone: int
    terminal_label: str | None
    terminal_purity: float | None
    n_cells: int
    segments: list[int]


@dataclass(frozen=True)
class Divergence:
    a: int
    b: int
    pseudotime: float


@dataclass(frozen=True)
class ResolvedRoot:
    spec: str  # canonical text form, e.g. "cluster:clusters=Ductal"
    cell_index: int
    obs_name: str
    cluster_key: str | None
    cluster_label: str | None
    resolved_by: str
    warnings: tuple[Notice, ...] = ()


@dataclass
class TrajectoryResult:
    method: str
    method_version: str
    topology_kind: TopologyKind
    milestones: list[Milestone]
    segments: list[Segment]
    lineages: list[Lineage]
    divergence: list[Divergence]
    pseudotime: npt.NDArray[np.float32]
    membership: npt.NDArray[np.bool_]  # (n_lineages, n_obs)
    segment: npt.NDArray[np.int32]  # -1 = not on a segment
    progress: npt.NDArray[np.float32]  # NaN where segment == -1
    root: ResolvedRoot
    inference_rep: str
    inference_dims: int
    params: dict[str, Any]
    milestone_weights: sparse.csr_matrix  # (n_obs, n_milestones), for drawing on an embedding
    fate: npt.NDArray[np.float32] | None = None  # (n_lineages, n_obs)
    extras: dict[str, npt.NDArray[np.float32]] = field(default_factory=dict)
    # A fitted curve (e.g. a principal tree): cell weights per curve node, and each
    # lineage's node path from the root. When present, it is drawn instead of
    # binned means, so the line on screen is the one that produced the pseudotime.
    curve_nodes: sparse.csr_matrix | None = None  # (n_obs, n_nodes)
    curve_paths: list[list[int]] | None = None  # one ordered node list per lineage
    warnings: list[Notice] = field(default_factory=list)
    runtime_seconds: float = 0.0

    @property
    def n_obs(self) -> int:
        return int(self.pseudotime.shape[0])

    @property
    def primary_lineage(self) -> npt.NDArray[np.int16]:
        """Lineage index for cells on exactly one lineage; -1 shared (trunk), -2 unassigned."""
        counts = self.membership.sum(axis=0)
        out = np.full(self.n_obs, PRIMARY_UNASSIGNED, dtype=np.int16)
        out[counts > 1] = PRIMARY_SHARED
        single = counts == 1
        if self.membership.shape[0]:
            out[single] = np.argmax(self.membership[:, single], axis=0).astype(np.int16)
        return out

    @property
    def n_unassigned(self) -> int:
        return int(np.count_nonzero(~np.isfinite(self.pseudotime)))

    def validate(self) -> None:
        problems: list[str] = []
        n = self.n_obs
        n_lin = len(self.lineages)

        def check(ok: bool, message: str) -> None:
            if not ok:
                problems.append(message)

        # shapes
        check(self.pseudotime.ndim == 1, "pseudotime must be 1-D")
        check(self.membership.shape == (n_lin, n), f"membership must be ({n_lin}, {n})")
        check(self.segment.shape == (n,), "segment must have one value per cell")
        check(self.progress.shape == (n,), "progress must have one value per cell")
        check(self.milestone_weights.shape == (n, len(self.milestones)), "milestone_weights shape")
        if self.fate is not None:
            check(self.fate.shape == (n_lin, n), "fate must be (n_lineages, n_obs)")
        for name, arr in self.extras.items():
            check(arr.shape == (n,), f"extras[{name!r}] must have one value per cell")
        if (self.curve_nodes is None) != (self.curve_paths is None):
            problems.append("curve_nodes and curve_paths come together")
        if self.curve_nodes is not None and self.curve_paths is not None:
            check(self.curve_nodes.shape[0] == n, "curve_nodes must have one row per cell")
            check(len(self.curve_paths) == n_lin, "curve_paths needs one path per lineage")
            n_nodes = self.curve_nodes.shape[1]
            check(
                all(0 <= k < n_nodes for path in self.curve_paths for k in path),
                "curve path node out of range",
            )
        if problems:
            raise InvariantError(problems)

        # pseudotime
        finite = np.isfinite(self.pseudotime)
        if finite.any():
            lo, hi = float(self.pseudotime[finite].min()), float(self.pseudotime[finite].max())
            check(lo >= -_PT_EPS and hi <= 1 + _PT_EPS, f"pseudotime must be in [0, 1], got [{lo}, {hi}]")
        check(0 <= self.root.cell_index < n, "root cell index out of range")
        if 0 <= self.root.cell_index < n:
            root_pt = float(self.pseudotime[self.root.cell_index])
            check(abs(root_pt) <= _PT_EPS, f"the root cell's pseudotime must be 0, got {root_pt}")

        # membership vs pseudotime
        on_any = self.membership.any(axis=0) if n_lin else np.zeros(n, dtype=bool)
        check(bool(np.all(on_any[finite])), "every cell with a pseudotime must belong to a lineage")
        check(not bool(np.any(on_any[~finite])), "cells without a pseudotime must belong to no lineage")

        # segments and progress
        n_seg = len(self.segments)
        check(bool(np.all((self.segment >= -1) & (self.segment < n_seg))), "segment values out of range")
        on_seg = self.segment >= 0
        prog = self.progress[on_seg]
        check(
            bool(np.all(np.isfinite(prog) & (prog >= -_PT_EPS) & (prog <= 1 + _PT_EPS))), "progress in [0, 1]"
        )

        # fate
        if self.fate is not None and finite.any():
            sums = self.fate[:, finite].sum(axis=0)
            check(bool(np.all(np.abs(sums - 1.0) <= _FATE_EPS)), "fate probabilities must sum to 1")

        # graph
        ids = [m.id for m in self.milestones]
        check(len(set(ids)) == len(ids), "milestone ids must be unique")
        roots = [m for m in self.milestones if m.kind == "root"]
        check(len(roots) == 1, f"exactly one root milestone expected, found {len(roots)}")
        id_set = set(ids)
        seg_by_id = {s.id: s for s in self.segments}
        check(len(seg_by_id) == n_seg, "segment ids must be unique")
        for s in self.segments:
            check(s.source in id_set and s.target in id_set, f"segment {s.id} references a missing milestone")

        # lineages
        names = [lin.name for lin in self.lineages]
        check(len(set(names)) == len(names), "lineage names must be unique")
        check([lin.id for lin in self.lineages] == list(range(n_lin)), "lineage ids must be 0..n-1 in order")
        root_id = roots[0].id if len(roots) == 1 else None
        for lin in self.lineages:
            check(lin.terminal_milestone in id_set, f"lineage {lin.name!r} ends at a missing milestone")
            path = [seg_by_id.get(sid) for sid in lin.segments]
            if any(p is None for p in path):
                problems.append(f"lineage {lin.name!r} uses a missing segment")
                continue
            segs = [p for p in path if p is not None]
            if segs:
                check(segs[0].source == root_id, f"lineage {lin.name!r} must start at the root")
                check(
                    segs[-1].target == lin.terminal_milestone,
                    f"lineage {lin.name!r} must end at its terminal",
                )
                for a, b in pairwise(segs):
                    check(a.target == b.source, f"lineage {lin.name!r} segments are not a path")
            else:
                check(
                    lin.terminal_milestone == root_id, f"empty lineage {lin.name!r} must be the root itself"
                )
        for d in self.divergence:
            check(0 <= d.a < n_lin and 0 <= d.b < n_lin and d.a != d.b, "divergence refers to bad lineages")

        if problems:
            raise InvariantError(problems)
