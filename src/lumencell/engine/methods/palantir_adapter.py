"""Palantir (Setty et al. 2019): pseudotime and fate probabilities.

Palantir runs its own diffusion map on the shared PCA, orders waypoints by
shortest paths from the root, and treats differentiation as a Markov chain:
each cell's fate probabilities are its absorption probabilities into the
terminal states. Lineages here are Palantir's branch masks.

Two choices differ from Palantir's defaults, both measured on the synthetic
benchmark (2026-11-10):

* **The root is used as given** (``use_early_cell_as_start=True``). By default
  Palantir swaps the given cell for the extreme of a diffusion component, which
  on a bifurcation can be a branch tip: the same silent root substitution that
  broke V1.
* **Terminal states come from the user, or else from the PAGA tree's leaves.**
  Palantir's automatic terminal detection picked cells near the fork on most
  synthetic cases at every setting tried (best: 3 of 12 cases failing, some
  with no lineage found). With the cell farthest from the root in each PAGA leaf
  cluster as the terminals, every true lineage was found. The trust panel says
  which source was used.
"""

import contextlib
import io
import time
import warnings
from collections.abc import Iterator
from typing import Any, ClassVar, Literal

import numpy as np
import numpy.typing as npt
import pandas as pd
from pydantic import Field
from scipy import sparse

from lumencell.engine.context import LEIDEN_KEY, InferenceContext
from lumencell.engine.methods.base import Availability, MethodParams, ProgressFn
from lumencell.engine.schema import Lineage, Milestone, Notice, ResolvedRoot, Segment, TrajectoryResult
from lumencell.engine.topology import cluster_tree, leaf_terminal_cells

INSTALL_HINT = 'Install it with: pip install "lumencell[palantir]"'


class PalantirParams(MethodParams):
    n_components: int = Field(10, ge=3, le=30, description="Diffusion components Palantir computes.")
    n_eigs: int | None = Field(
        None,
        ge=2,
        le=30,
        description="Components kept in the multiscale space; empty = Palantir's eigengap rule.",
    )
    knn: int = Field(30, ge=5, le=100, description="Neighbours in Palantir's graph.")
    num_waypoints: int = Field(1200, ge=100, le=5000, description="Waypoints sampled to compute pseudotime.")
    terminal_detection: Literal["paga_leaves", "palantir"] = Field(
        "paga_leaves", description="Where terminal states come from when you name none."
    )
    cluster_key: str = Field(LEIDEN_KEY, description="obs column that terminal names refer to.")
    seed: int = 20


@contextlib.contextmanager
def _quiet() -> Iterator[None]:
    """Palantir prints progress and trips scipy efficiency warnings internally; keep logs clean."""
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore", sparse.SparseEfficiencyWarning)
        yield


def _neighbour_weights(ctx: InferenceContext, cells: list[int]) -> sparse.csr_matrix:
    """(n_obs, len(cells)) weights: each anchor cell plus its kNN neighbours."""
    conn = sparse.csr_matrix(ctx.artifacts["connectivities"])
    cols = []
    for c in cells:
        col = conn[c].toarray().ravel().astype(np.float32)
        col[c] = max(1.0, float(col.max()))
        cols.append(col)
    return sparse.csr_matrix(np.stack(cols, axis=1))


class Palantir:
    name: ClassVar[str] = "palantir"
    label: ClassVar[str] = "Palantir"
    family: ClassVar[str] = "Markov chain fate probabilities"
    description: ClassVar[str] = (
        "Waypoint pseudotime and fate probabilities from a Markov chain on the diffusion map "
        "(Setty et al. 2019)."
    )
    Params: ClassVar[type[MethodParams]] = PalantirParams
    # Misses one gate on the synthetic matrix (2026-11-10): root excess on two
    # bifurcation seeds (0.177, 0.158 > 0.10); Palantir's pseudotime is compressed
    # near the start. Branch accuracy >= 0.90 and recall 1.0 on every case.
    experimental: ClassVar[bool] = True
    experimental_note: ClassVar[str | None] = (
        "On the synthetic benchmark its pseudotime is compressed near the start on 2 of 12 cases (root "
        "excess 0.16-0.18 against 0.10). Branch accuracy 0.90-0.99 and every lineage found on all cases. "
        "On the synthetic clone layer its lineages captured 0.37-0.56 "
        "of the chance-corrected clonal signal of the true lineages (gate: 0.5) on 7 of 9 branching cases."
    )

    def available(self) -> Availability:
        try:
            from importlib.metadata import version

            import palantir  # noqa: F401

            return Availability(ok=True, version=version("palantir"))
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
        import palantir

        assert isinstance(params, PalantirParams)
        start = time.perf_counter()
        warnings: list[Notice] = []
        names = [str(n) for n in ctx.obs_names]

        view = ad.AnnData(obs=pd.DataFrame(index=pd.Index(names)))
        view.obsm["X_pca"] = np.asarray(ctx.artifacts["pca"])

        progress("diffusion map", 0, 4)
        with _quiet():
            palantir.utils.run_diffusion_maps(view, n_components=params.n_components, knn=params.knn, seed=0)
            palantir.utils.determine_multiscale_space(view, n_eigs=params.n_eigs)
        multiscale: pd.DataFrame = view.obsm["DM_EigenVectors_multiscaled"]
        ms_dims = int(multiscale.shape[1])

        # ---- terminal states ---------------------------------------------------------
        progress("terminal states", 1, 4)
        terminal_cells: list[str] | None
        if terminals:
            labels = ctx.labels(params.cluster_key)
            missing = [t for t in terminals if t not in set(labels.tolist())]
            if missing:
                raise ValueError(f"terminal clusters not found in {params.cluster_key}: {', '.join(missing)}")
            view.obs["_terminal_groups"] = pd.Categorical(labels)
            with _quiet():
                found = palantir.utils.find_terminal_states(
                    view, terminals, celltype_column="_terminal_groups"
                )
            terminal_cells = [str(c) for c in found.index]
            source = "the terminals you named"
        elif params.terminal_detection == "paga_leaves":
            tree = cluster_tree(ctx, root.cell_index)
            terminal_cells = [names[c] for c in leaf_terminal_cells(ctx, tree, root.cell_index)] or None
            source = "the leaves of the PAGA cluster tree"
        else:
            terminal_cells = None
            source = "Palantir's automatic detection"

        progress("pseudotime and fates", 2, 4)
        run_kw: dict[str, Any] = {
            "early_cell": root.obs_name,
            "knn": params.knn,
            "num_waypoints": min(params.num_waypoints, ctx.n_obs),
            "seed": params.seed,
            "n_jobs": 1,
            "use_early_cell_as_start": True,
        }
        with _quiet():
            palantir.core.run_palantir(view, terminal_states=terminal_cells, **run_kw)
            fates_found: pd.DataFrame = view.obsm["palantir_fate_probabilities"]
            if fates_found.shape[1] == 0:
                # No terminal state at all (a linear trajectory): end at the latest cell.
                last = names[int(np.nanargmax(view.obs["palantir_pseudotime"].to_numpy()))]
                palantir.core.run_palantir(view, terminal_states=[last], **run_kw)
                source = "the latest cell (no branch was detected)"
            palantir.presults.select_branch_cells(view)

        progress("lineages", 3, 4)
        fate_df: pd.DataFrame = view.obsm["palantir_fate_probabilities"]
        terms = [str(c) for c in fate_df.columns]
        # np.array copies: values from pandas can be read-only views.
        fate: npt.NDArray[np.float64] = np.array(fate_df, dtype=np.float64).T  # (L, n)
        fate = fate / np.maximum(fate.sum(axis=0, keepdims=True), 1e-12)
        masks: npt.NDArray[np.bool_] = np.array(view.obsm["branch_masks"], dtype=bool).T.copy()  # (L, n)

        pt = np.array(view.obs["palantir_pseudotime"], dtype=np.float64)
        pt = np.clip(pt - pt[root.cell_index], 0.0, None)
        finite = np.isfinite(pt)
        pt = pt / max(float(np.nanmax(pt)), 1e-12)
        # Cells in no branch mask (Palantir leaves a few) join their most likely fate.
        orphan = finite & ~masks.any(axis=0)
        if orphan.any():
            best = np.argmax(fate[:, orphan], axis=0)
            masks[best, np.flatnonzero(orphan)] = True
        masks &= finite[None, :]

        idx = {n: i for i, n in enumerate(names)}
        term_idx = [idx[t] for t in terms]
        group_labels = ctx.labels(params.cluster_key)
        root_label = group_labels[root.cell_index]
        milestones = [Milestone(0, "root", str(root_label), 0.0)] + [
            Milestone(i + 1, "tip", str(group_labels[c]), float(pt[c])) for i, c in enumerate(term_idx)
        ]
        segments = [Segment(i, 0, i + 1, float(pt[c])) for i, c in enumerate(term_idx)]
        lineages = [
            Lineage(i, f"Lineage {i + 1}", f"Lineage {i + 1}", i + 1, None, None, int(masks[i].sum()), [i])
            for i in range(len(term_idx))
        ]
        segment = np.where(finite, np.argmax(fate, axis=0), -1).astype(np.int32)
        tip_pt = np.array([max(float(pt[c]), 1e-9) for c in term_idx])
        progress_arr = np.where(
            segment >= 0, np.clip(pt / tip_pt[np.maximum(segment, 0)], 0.0, 1.0), np.nan
        ).astype(np.float32)

        warnings.append(
            Notice(
                "terminals_source",
                "info",
                f"Terminal states came from {source}. Name the terminal clusters to use your own.",
            )
        )
        return TrajectoryResult(
            method=self.name,
            method_version=version("palantir"),
            topology_kind="fate_star",
            milestones=milestones,
            segments=segments,
            lineages=lineages,
            divergence=[],
            pseudotime=pt.astype(np.float32),
            membership=masks,
            segment=segment,
            progress=progress_arr,
            root=root,
            inference_rep="palantir multiscale diffusion map",
            inference_dims=ms_dims,
            params=params.model_dump(),
            milestone_weights=_neighbour_weights(ctx, [root.cell_index, *term_idx]),
            fate=fate.astype(np.float32),
            extras={"entropy": np.asarray(view.obs["palantir_entropy"], dtype=np.float32)},
            warnings=warnings,
            runtime_seconds=time.perf_counter() - start,
        )
