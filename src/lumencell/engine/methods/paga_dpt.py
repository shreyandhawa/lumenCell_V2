"""PAGA + diffusion pseudotime (scanpy).

* **Pseudotime:** diffusion pseudotime (Haghverdi et al. 2016) from the root
  cell, on the shared diffusion map.
* **Topology:** the PAGA graph (Wolf et al. 2019) over clusters, reduced to its
  maximum-connectivity spanning tree. Weak tree edges (connectivity below
  ``min_connectivity``) are cut, and the tree is oriented from the root's cluster.
* **Lineages:** paths from the root cluster to each terminal cluster (the tree's
  leaves, or the clusters you name).

Clusters are milestones; a cell belongs to every lineage whose path contains its
cluster. That is coarse near forks (a cluster spanning a fork counts as shared),
which the benchmark's branch accuracy measures honestly.
"""

import time
from importlib.metadata import version as package_version
from typing import ClassVar

import numpy as np
import scanpy as sc
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
from lumencell.engine.topology import cluster_tree
from lumencell.prep.stages import neighbors_view

# A child cluster whose median pseudotime is less than this past its parent's is
# more likely a sister branch than a descendant (see possible_sibling_chain).
SIBLING_GAP = 0.05


class PagaDptParams(MethodParams):
    n_dcs: int = Field(10, ge=2, le=14, description="Diffusion components used by DPT.")
    min_connectivity: float = Field(
        0.05, ge=0.0, le=1.0, description="PAGA tree edges weaker than this are cut."
    )
    cluster_key: str = Field(LEIDEN_KEY, description="obs column whose groups become PAGA nodes.")


class PagaDpt:
    name: ClassVar[str] = "paga_dpt"
    label: ClassVar[str] = "PAGA + DPT"
    family: ClassVar[str] = "cluster graph + diffusion pseudotime"
    description: ClassVar[str] = (
        "Diffusion pseudotime from the root cell, with lineages along the PAGA cluster tree (scanpy)."
    )
    Params: ClassVar[type[MethodParams]] = PagaDptParams
    # Misses gates on the synthetic matrix (2026-11-10): branch accuracy on
    # trifurcations (0.81-0.92; clusters that straddle a fork get wired through a
    # sibling branch), global Spearman on binary2 (0.86-0.88; DPT runs at different
    # speeds on branches of different length) and root excess once (0.105).
    # Pseudotime within each lineage is 0.95-0.99 everywhere.
    experimental: ClassVar[bool] = True

    def available(self) -> Availability:
        return Availability(ok=True, version=package_version("scanpy"))

    def run(
        self,
        ctx: InferenceContext,
        root: ResolvedRoot,
        terminals: list[str] | None,
        params: MethodParams,
        progress: ProgressFn,
    ) -> TrajectoryResult:
        assert isinstance(params, PagaDptParams)
        start = time.perf_counter()
        warnings: list[Notice] = []

        # ---- pseudotime ---------------------------------------------------------------
        progress("dpt", 0, 3)
        view = neighbors_view(ctx.n_obs, ctx.artifacts, ctx.prep_params)
        view.obsm["X_diffmap"] = np.asarray(ctx.artifacts["diffmap"])
        view.uns["diffmap_evals"] = np.asarray(ctx.artifacts["diffmap_evals"])
        view.uns["iroot"] = root.cell_index
        sc.tl.dpt(view, n_dcs=params.n_dcs, n_branchings=0)
        # np.array copies: under pandas 3 copy-on-write, to_numpy() can return a read-only view.
        pt = np.array(view.obs["dpt_pseudotime"], dtype=np.float64)
        pt[~np.isfinite(pt)] = np.nan

        # ---- cluster tree -------------------------------------------------------------
        progress("paga", 1, 3)
        tree = cluster_tree(ctx, root.cell_index, params.cluster_key, params.min_connectivity)
        cats, group_of, root_group = tree.cats, tree.group_of, tree.root_group
        parent, order, children = tree.parent, tree.order, tree.children
        unreachable = tree.unreachable
        if unreachable:
            n_cells = int(np.isin(group_of, unreachable).sum())
            warnings.append(
                Notice(
                    "clusters_cut_off",
                    "warning",
                    f"{len(unreachable)} clusters ({n_cells} cells) are only weakly connected to the "
                    f"rest (PAGA connectivity < {params.min_connectivity}) and get no pseudotime.",
                )
            )

        # ---- milestones and segments ----------------------------------------------------
        progress("lineages", 2, 3)
        mid = {g: i for i, g in enumerate(order)}
        med = {
            g: float(np.nanmedian(pt[group_of == g])) if np.isfinite(pt[group_of == g]).any() else 0.0
            for g in order
        }
        med[root_group] = 0.0
        milestones: list[Milestone] = []
        for g in order:
            n_kids = len(children[g])
            kind: MilestoneKind
            if g == root_group:
                kind = "root"
            elif n_kids == 0:
                kind = "tip"
            elif n_kids >= 2:
                kind = "fork"
            else:
                kind = "cluster"
            milestones.append(Milestone(mid[g], kind, cats[g], med[g]))

        segments: list[Segment] = []
        seg_into: dict[int, int] = {}
        for g in order[1:]:
            p = parent[g]
            if med[g] < med[p]:
                warnings.append(
                    Notice(
                        "paga_dpt_disagree",
                        "info",
                        f"PAGA places cluster {cats[g]} after {cats[p]}, but its median DPT is lower "
                        f"({med[g]:.2f} < {med[p]:.2f}).",
                    )
                )
            elif p != root_group and med[g] - med[p] < SIBLING_GAP:
                # Sister branches share many neighbour edges just past a fork, so PAGA's
                # maximum-connectivity tree can chain them (trunk -> A -> B). A real
                # descendant sits clearly later than its parent; a sibling does not.
                warnings.append(
                    Notice(
                        "possible_sibling_chain",
                        "warning",
                        f"PAGA places cluster {cats[g]} after {cats[p]}, but both sit at the same "
                        f"pseudotime ({med[p]:.2f} vs {med[g]:.2f}); they may be sister branches. "
                        "Finer clusters (the Leiden default) or naming the terminals usually helps.",
                    )
                )
            seg_into[g] = len(segments)
            segments.append(Segment(len(segments), mid[p], mid[g], max(0.0, med[g] - med[p])))

        # ---- lineages ----------------------------------------------------------------
        if terminals:
            missing = [t for t in terminals if t not in cats or cats.index(t) not in parent]
            if missing:
                raise ValueError(
                    f"terminal clusters not found among reachable clusters: {', '.join(missing)}"
                )
            terminal_groups = [cats.index(t) for t in terminals]
        else:
            terminal_groups = [g for g in order if not children[g] and g != root_group] or [root_group]

        paths = [tree.path_to(g) for g in terminal_groups]
        finite = np.isfinite(pt)
        membership = np.stack([np.isin(group_of, p) & finite for p in paths])
        on_any = membership.any(axis=0)
        orphan = finite & ~on_any
        if orphan.any():
            warnings.append(
                Notice(
                    "cells_off_lineages",
                    "info",
                    f"{int(orphan.sum())} cells sit on branches that lead to none of the chosen "
                    "terminals and get no pseudotime.",
                )
            )
        pt[~on_any] = np.nan
        lineages = [
            Lineage(
                id=i,
                name=f"Lineage {i + 1}",
                auto_name=f"Lineage {i + 1}",
                terminal_milestone=mid[p[-1]],
                terminal_label=None,
                terminal_purity=None,
                n_cells=int(membership[i].sum()),
                segments=[seg_into[g] for g in p[1:]],
            )
            for i, p in enumerate(paths)
        ]

        # ---- per-cell placement --------------------------------------------------------
        segment = np.full(ctx.n_obs, -1, dtype=np.int32)
        progress_arr = np.full(ctx.n_obs, np.nan, dtype=np.float32)
        for g, s in seg_into.items():
            cells = np.flatnonzero((group_of == g) & np.isfinite(pt))
            if cells.size == 0:
                continue
            segment[cells] = s
            ranks = np.argsort(np.argsort(pt[cells], kind="stable"), kind="stable")
            progress_arr[cells] = ranks / max(cells.size - 1, 1)

        rows = np.flatnonzero(np.isin(group_of, order))
        weights = sparse.csr_matrix(
            (np.ones(rows.size, dtype=np.float32), (rows, [mid[int(g)] for g in group_of[rows]])),
            shape=(ctx.n_obs, len(milestones)),
        )

        return TrajectoryResult(
            method=self.name,
            method_version=package_version("scanpy"),
            topology_kind="tree",
            milestones=milestones,
            segments=segments,
            lineages=lineages,
            divergence=[],
            pseudotime=pt.astype(np.float32),
            membership=membership,
            segment=segment,
            progress=progress_arr,
            root=root,
            inference_rep="X_diffmap",
            inference_dims=params.n_dcs,
            params=params.model_dump(),
            milestone_weights=weights,
            warnings=warnings,
            runtime_seconds=time.perf_counter() - start,
        )
