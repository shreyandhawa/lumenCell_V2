"""scFates SimplePPT (Faure et al. 2023): a principal tree fitted in diffusion space.

A tree of principal points is fitted through the cells; each cell is
soft-assigned to points, projected onto the tree, and its pseudotime is the
distance along the tree from the root tip. Unlike the other methods, the curve
drawn on screen is the fitted tree itself (each point placed at the weighted
mean of its cells), so it is exactly what produced the pseudotime.

Defaults, measured on the synthetic benchmark (2026-11-10):

* **Inference space:** Palantir's multiscale diffusion space from the shared
  prep. On a bifurcation it scored branch accuracy 0.975 against 0.33 when fitted
  on 30 PCs (where the tree grew 14 spurious tips).
* **Pruning:** side branches shorter than 5 points are removed. With scFates'
  default of 3, two-fork cases grew 2 spurious branches on 2 of 3 seeds; 5, 8 and
  12 all removed them with identical results elsewhere.
* **Root:** if the root cell sits on a dead-end arm of the tree, the arm's tip;
  if it sits between forks (stem cells with lineages radiating out), its own
  point. See ``_root_point``.
"""

import contextlib
import io
import time
import warnings
from collections import deque
from collections.abc import Iterator
from itertools import pairwise
from typing import ClassVar, Literal

import numpy as np
import numpy.typing as npt
import pandas as pd
from pydantic import Field
from scipy import sparse

from lumencell.engine.context import LEIDEN_KEY, InferenceContext
from lumencell.engine.methods.base import Availability, MethodParams, ProgressFn
from lumencell.engine.schema import (
    Lineage,
    Milestone,
    MilestoneKind,
    Notice,
    ResolvedRoot,
    Segment,
    TrajectoryResult,
)

INSTALL_HINT = 'Install it with: pip install "lumencell[scfates]"'
R_FLOOR = 1e-3  # soft assignments below this are dropped when storing the curve weights


class ScFatesParams(MethodParams):
    use_rep: Literal["diffmap_ms", "pca"] = Field(
        "diffmap_ms", description="Space the tree is fitted in: multiscale diffusion or 30 PCs."
    )
    nodes: int | None = Field(
        None, ge=20, le=1000, description="Principal points; empty = min(200, cells/10)."
    )
    ppt_lambda: float = Field(1.0, gt=0, description="SimplePPT tree-length penalty.")
    ppt_sigma: float = Field(0.1, gt=0, description="SimplePPT soft-assignment bandwidth.")
    ppt_nsteps: int = Field(50, ge=5, le=500)
    minbranchlength: int = Field(5, ge=1, le=50, description="Prune side branches shorter than this.")
    cluster_key: str = Field(LEIDEN_KEY, description="obs column that terminal names refer to.")
    seed: int = 0


@contextlib.contextmanager
def _quiet() -> Iterator[None]:
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        yield


def _node_path(adj: sparse.csr_matrix, start: int, end: int) -> list[int]:
    parent = {start: -1}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        if node == end:
            break
        for nb in (int(x) for x in adj[node].indices):
            if nb not in parent:
                parent[nb] = node
                queue.append(nb)
    path = [end]
    while parent[path[-1]] != -1:
        path.append(parent[path[-1]])
    return path[::-1]


def _arm_end(adj: sparse.csr_matrix, start: int, first: int) -> tuple[int, int]:
    """Walk from ``start`` through ``first`` along degree-2 points; return (end point, steps)."""
    prev, cur, steps = start, first, 1
    while len(adj[cur].indices) == 2:
        nxt = next(int(n) for n in adj[cur].indices if int(n) != prev)
        prev, cur, steps = cur, nxt, steps + 1
    return cur, steps


def _root_point(adj: sparse.csr_matrix, node: int, tips: set[int]) -> int:
    """Where to root the tree, given the point the root cell is most assigned to.

    * On a dead-end arm (its segment ends at a tip), the trajectory starts at that
      tip: the cells between the tip and the root cell are the earliest ones, and
      rooting mid-arm would give them pseudotime counting backwards.
    * Between forks (stem cells with lineages radiating out), root at the point
      itself: any tip lies beyond other lineages, and their cells would be wrongly
      placed before the root.
    """
    if node in tips:
        return node
    neighbours = [int(n) for n in adj[node].indices]
    if len(neighbours) != 2:
        return node  # the root cell sits on a fork
    ends = [_arm_end(adj, node, n) for n in neighbours]
    tip_ends = [(steps, end) for end, steps in ends if end in tips]
    if not tip_ends:
        return node
    return min(tip_ends)[1]  # on a linear tree, the nearer end


class ScFates:
    name: ClassVar[str] = "scfates_ppt"
    label: ClassVar[str] = "scFates principal tree"
    family: ClassVar[str] = "principal tree"
    description: ClassVar[str] = (
        "A principal tree fitted through the cells in diffusion space; pseudotime is distance along it "
        "(scFates SimplePPT)."
    )
    Params: ClassVar[type[MethodParams]] = ScFatesParams
    experimental: ClassVar[bool] = True
    experimental_note: ClassVar[str | None] = (
        "On the synthetic benchmark it misses one gate: global Spearman on two-fork cases (0.87-0.88), "
        "where branches of different length progress at different speeds. Branch accuracy 0.97+ and "
        "every lineage found on all cases."
    )

    def available(self) -> Availability:
        try:
            from importlib.metadata import version

            import scFates  # noqa: F401

            return Availability(ok=True, version=version("scFates"))
        except ImportError:
            return Availability(ok=False, install_hint=INSTALL_HINT)

    def run(
        self,
        ctx: InferenceContext,
        root: ResolvedRoot,
        terminals: list[str] | None,
        params: MethodParams,
        progress: ProgressFn,
    ) -> TrajectoryResult:
        from importlib.metadata import version

        import anndata as ad
        import scFates as scf

        assert isinstance(params, ScFatesParams)
        start = time.perf_counter()
        notices: list[Notice] = []
        rep: npt.NDArray[np.float64] = (
            np.asarray(ctx.artifacts["diffmap_ms"], dtype=np.float64)
            if params.use_rep == "diffmap_ms"
            else np.asarray(ctx.artifacts["pca"], dtype=np.float64)[:, :30]
        )
        view = ad.AnnData(obs=pd.DataFrame(index=pd.Index([str(n) for n in ctx.obs_names])))
        view.obsm["X_rep"] = rep

        progress("fit tree", 0, 3)
        with _quiet():
            scf.tl.tree(
                view,
                method="ppt",
                Nodes=params.nodes or max(20, min(200, ctx.n_obs // 10)),
                use_rep="X_rep",
                ndims_rep=rep.shape[1],
                ppt_lambda=params.ppt_lambda,
                ppt_sigma=params.ppt_sigma,
                ppt_nsteps=params.ppt_nsteps,
                seed=params.seed,
            )
            scf.tl.cleanup(view, minbranchlength=params.minbranchlength)
        graph = view.uns["graph"]
        R: npt.NDArray[np.float64] = np.asarray(view.obsm["X_R"], dtype=np.float64)  # (cells, nodes)
        adj = sparse.csr_matrix(graph["B"])
        tips = [int(t) for t in graph["tips"]]
        forks = {int(f) for f in graph["forks"]}

        # ---- root and terminals ------------------------------------------------------
        progress("pseudotime", 1, 3)
        root_tip = _root_point(adj, int(np.argmax(R[root.cell_index])), set(tips))
        labels = ctx.labels(params.cluster_key)
        if terminals:
            missing = [t for t in terminals if t not in set(labels.tolist())]
            if missing:
                raise ValueError(f"terminal clusters not found in {params.cluster_key}: {', '.join(missing)}")
            others = [t for t in tips if t != root_tip]
            chosen = [max(others, key=lambda t: float(R[labels == lab, t].mean())) for lab in terminals]
            end_tips = list(dict.fromkeys(chosen))
        else:
            end_tips = [t for t in tips if t != root_tip]
        with _quiet():
            scf.tl.root(view, root_tip)
            scf.tl.pseudotime(view, n_jobs=1, n_map=1, seed=params.seed)
        graph = view.uns["graph"]
        pp_seg: pd.DataFrame = graph["pp_seg"]
        pp_info: pd.DataFrame = graph["pp_info"]
        node_time = {int(n): float(t) for n, t in zip(pp_info.index, pp_info["time"], strict=True)}

        t_raw = np.array(view.obs["t"], dtype=np.float64)
        pt = np.clip(t_raw - t_raw[root.cell_index], 0.0, None)
        scale = max(float(np.nanmax(pt)), 1e-12)
        pt = pt / scale
        cell_seg = np.array(view.obs["seg"].astype(str))

        # ---- graph in schema terms -------------------------------------------------------
        progress("lineages", 2, 3)
        seg_rows = [(str(i), int(r["from"]), int(r["to"])) for i, r in pp_seg.iterrows()]
        ms_nodes = list(dict.fromkeys([root_tip] + [n for _, a, b in seg_rows for n in (a, b)]))
        mid = {n: i for i, n in enumerate(ms_nodes)}
        hard = np.argmax(R, axis=1)

        def kind_of(n: int) -> MilestoneKind:
            if n == root_tip:
                return "root"
            return "fork" if n in forks else "tip"

        def label_of(n: int) -> str:
            near = labels[hard == n]
            if near.size == 0:
                return f"node {n}"
            values, counts = np.unique(near, return_counts=True)
            return str(values[np.argmax(counts)])

        milestones = [
            Milestone(
                mid[n],
                kind_of(n),
                label_of(n),
                float(np.clip((node_time[n] - t_raw[root.cell_index]) / scale, 0, 1)),
            )
            for n in ms_nodes
        ]
        seg_index = {key: i for i, (key, _, _) in enumerate(seg_rows)}
        segments = [
            Segment(i, mid[a], mid[b], max(0.0, (node_time[b] - node_time[a]) / scale))
            for i, (_key, a, b) in enumerate(seg_rows)
        ]
        by_ends = {(a, b): key for key, a, b in seg_rows}

        lineages: list[Lineage] = []
        paths: list[list[int]] = []
        masks = []
        for i, tip in enumerate(end_tips):
            path = _node_path(adj, root_tip, tip)
            on_path = [n for n in path if n in mid]
            keys = [by_ends[(a, b)] for a, b in pairwise(on_path) if (a, b) in by_ends]
            lineages.append(
                Lineage(
                    i,
                    f"Lineage {i + 1}",
                    f"Lineage {i + 1}",
                    mid[tip],
                    None,
                    None,
                    0,
                    [seg_index[k] for k in keys],
                )
            )
            paths.append(path)
            masks.append(np.isin(cell_seg, keys))
        membership = np.array(masks, dtype=bool) & np.isfinite(pt)[None, :]
        off = np.isfinite(pt) & ~membership.any(axis=0)
        if off.any():
            notices.append(
                Notice(
                    "cells_off_lineages",
                    "info",
                    f"{int(off.sum())} cells sit on branches that lead to none of the chosen terminals and "
                    "get no pseudotime.",
                )
            )
        pt[~membership.any(axis=0)] = np.nan
        for lin in lineages:
            lin.n_cells = int(membership[lin.id].sum())

        segment = np.array([seg_index.get(s, -1) for s in cell_seg], dtype=np.int32)
        segment[~np.isfinite(pt)] = -1
        t_from = np.array([node_time[seg_rows[s][1]] if s >= 0 else np.nan for s in segment])
        t_to = np.array([node_time[seg_rows[s][2]] if s >= 0 else np.nan for s in segment])
        with np.errstate(invalid="ignore", divide="ignore"):
            prog = np.clip((t_raw - t_from) / np.where(t_to > t_from, t_to - t_from, 1.0), 0.0, 1.0)
        prog[segment < 0] = np.nan

        R_sparse = sparse.csr_matrix(np.where(R >= R_FLOOR, R, 0.0).astype(np.float32))
        notices.append(
            Notice(
                "tree_fitted",
                "info",
                f"Principal tree with {R.shape[1]} points fitted in {params.use_rep}; "
                f"{len(tips)} tips after pruning branches shorter than {params.minbranchlength} points.",
            )
        )
        return TrajectoryResult(
            method=self.name,
            method_version=version("scFates"),
            topology_kind="tree",
            milestones=milestones,
            segments=segments,
            lineages=lineages,
            divergence=[],
            pseudotime=pt.astype(np.float32),
            membership=membership,
            segment=segment,
            progress=prog.astype(np.float32),
            root=root,
            inference_rep="X_diffmap_ms" if params.use_rep == "diffmap_ms" else "X_pca",
            inference_dims=int(rep.shape[1]),
            params=params.model_dump(),
            milestone_weights=sparse.csr_matrix(R_sparse[:, ms_nodes]),
            warnings=notices,
            runtime_seconds=time.perf_counter() - start,
            curve_nodes=R_sparse,
            curve_paths=paths,
        )
