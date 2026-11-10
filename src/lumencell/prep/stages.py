"""The preprocessing stages.

Each stage declares the artifacts it reads, the artifacts it writes and the
parameters it depends on, so a stage's cache key is a hash of exactly those
(see ``prep/cache.py`` in the platform track). Changing ``n_top_genes`` then
re-runs HVG selection and everything downstream, but not normalisation.

Artifacts are plain values: numpy arrays, scipy CSR matrices, or small
JSON-compatible dicts. Stages never mutate their inputs or the dataset.

Inference uses the PCA / kNN / diffusion space built here. UMAP is computed
only for display, and only when the file has no 2-D embedding of its own.
"""

import time
from collections.abc import Callable, Mapping
from dataclasses import asdict, dataclass
from typing import Any

import anndata as ad
import numpy as np
import numpy.typing as npt
import pandas as pd
import scanpy as sc
from scipy import sparse

from lumencell.data.loader import Dataset

PIPELINE_VERSION = "1"

Artifacts = Mapping[str, Any]


@dataclass(frozen=True)
class PrepParams:
    n_top_genes: int = 2000
    n_pcs: int = 50
    n_neighbors: int = 30
    neighbors_pcs: int = 30
    leiden_resolution: float = 1.0
    n_diffmap: int = 15
    seed: int = 0
    recompute_umap: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class Stage:
    name: str
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    params: tuple[str, ...]
    run: Callable[[Dataset, Artifacts, PrepParams], dict[str, Any]]
    # Returns False when the stage should be skipped for this dataset.
    applies: Callable[[Dataset, PrepParams], bool] = lambda _ds, _p: True


# ---- helpers ----------------------------------------------------------------------


def neighbors_view(n_obs: int, art: Artifacts, params: PrepParams) -> ad.AnnData:
    """A minimal AnnData that scanpy's graph-based tools accept.

    It carries only the kNN graph and the PCA, built fresh from artifacts, so
    scanpy calls can never write into shared state.
    """
    view = ad.AnnData(obs=pd.DataFrame(index=pd.RangeIndex(n_obs).astype(str)))
    view.obsm["X_pca"] = np.asarray(art["pca"])
    view.obsp["connectivities"] = art["connectivities"]
    view.obsp["distances"] = art["distances"]
    view.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {
            "n_neighbors": params.n_neighbors,
            "method": "umap",
            "random_state": params.seed,
            "metric": "euclidean",
            "use_rep": "X_pca",
            "n_pcs": min(params.neighbors_pcs, np.asarray(art["pca"]).shape[1]),
        },
    }
    return view


def _counts_view(ds: Dataset) -> ad.AnnData:
    return ad.AnnData(X=ds.matrix.copy(), obs=pd.DataFrame(index=pd.RangeIndex(ds.n_obs).astype(str)))


def multiscale_space(evecs: npt.NDArray[Any], evals: npt.NDArray[Any]) -> npt.NDArray[np.float32]:
    """Palantir's multiscale diffusion space (Setty et al. 2019).

    Drops the trivial first eigenvector, keeps components up to the largest
    eigengap (at least two), and scales each by ``lambda / (1 - lambda)``.

    This follows ``palantir.utils.determine_multiscale_space`` exactly: the gaps
    are taken over *all* eigenvalues (trivial one included), ``n_eigs`` is the
    index after the largest gap (the second largest if that leaves fewer than
    3), and components ``1 .. n_eigs - 1`` are kept.
    """
    vals = np.asarray(evals, dtype=np.float64)
    vecs = np.asarray(evecs, dtype=np.float64)
    order = np.argsort(vals[:-1] - vals[1:], kind="stable")
    n_eigs = int(order[-1]) + 1
    if n_eigs < 3 and len(order) > 1:
        n_eigs = int(order[-2]) + 1
    n_eigs = min(max(3, n_eigs), vals.size)
    use = np.arange(1, n_eigs)
    lam = np.clip(vals[use], 0.0, 1.0 - 1e-9)
    return (vecs[:, use] * (lam / (1.0 - lam))[None, :]).astype(np.float32)


# ---- stages -----------------------------------------------------------------------


def _normalize(ds: Dataset, _art: Artifacts, _p: PrepParams) -> dict[str, Any]:
    if ds.x_kind == "lognorm":
        return {"lognorm": ds.matrix.copy()}
    view = _counts_view(ds)
    sc.pp.normalize_total(view, target_sum=1e4)
    sc.pp.log1p(view)
    return {"lognorm": sparse.csr_matrix(view.X, dtype=np.float32)}


def _hvg(ds: Dataset, art: Artifacts, p: PrepParams) -> dict[str, Any]:
    n_top = min(p.n_top_genes, ds.n_vars)
    if ds.x_kind == "counts":
        view = _counts_view(ds)
        sc.pp.highly_variable_genes(view, flavor="seurat_v3", n_top_genes=n_top)
    else:
        view = ad.AnnData(X=art["lognorm"].copy())
        sc.pp.highly_variable_genes(view, flavor="seurat", n_top_genes=n_top)
    return {"hvg": view.var["highly_variable"].to_numpy(dtype=bool)}


def _pca(ds: Dataset, art: Artifacts, p: PrepParams) -> dict[str, Any]:
    hvg = np.asarray(art["hvg"], dtype=bool)
    n_comps = max(2, min(p.n_pcs, int(hvg.sum()) - 1, ds.n_obs - 1))
    view = ad.AnnData(X=art["lognorm"][:, hvg].copy())
    sc.pp.pca(view, n_comps=n_comps, svd_solver="covariance_eigh", random_state=p.seed)
    pcs: npt.NDArray[np.float32] = np.asarray(view.obsm["X_pca"], dtype=np.float32)
    # Fix each component's sign (largest-magnitude loading positive) so runs agree exactly.
    loadings: npt.NDArray[np.float64] = np.asarray(view.varm["PCs"], dtype=np.float64)
    signs = np.sign(loadings[np.argmax(np.abs(loadings), axis=0), np.arange(loadings.shape[1])])
    signs[signs == 0] = 1.0
    return {
        "pca": pcs * signs[None, :].astype(np.float32),
        "pca_variance_ratio": np.asarray(view.uns["pca"]["variance_ratio"], dtype=np.float64),
    }


def _neighbors(ds: Dataset, art: Artifacts, p: PrepParams) -> dict[str, Any]:
    view = ad.AnnData(obs=pd.DataFrame(index=pd.RangeIndex(ds.n_obs).astype(str)))
    pca: npt.NDArray[np.float32] = np.asarray(art["pca"], dtype=np.float32)
    view.obsm["X_pca"] = pca
    n_pcs = min(p.neighbors_pcs, pca.shape[1])
    sc.pp.neighbors(
        view,
        n_neighbors=min(p.n_neighbors, ds.n_obs - 1),
        n_pcs=n_pcs,
        use_rep="X_pca",
        random_state=p.seed,
    )
    return {
        "connectivities": sparse.csr_matrix(view.obsp["connectivities"], dtype=np.float32),
        "distances": sparse.csr_matrix(view.obsp["distances"], dtype=np.float32),
    }


def _leiden(ds: Dataset, art: Artifacts, p: PrepParams) -> dict[str, Any]:
    view = neighbors_view(ds.n_obs, art, p)
    sc.tl.leiden(
        view,
        resolution=p.leiden_resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=p.seed,
        key_added="leiden",
    )
    labels = view.obs["leiden"].astype(int).to_numpy()
    # Renumber clusters by size (largest = 0) so labels are stable and readable.
    order = np.argsort(-np.bincount(labels), kind="stable")
    remap = np.empty_like(order)
    remap[order] = np.arange(order.size)
    return {"leiden": remap[labels].astype(np.int32)}


def _diffmap(ds: Dataset, art: Artifacts, p: PrepParams) -> dict[str, Any]:
    view = neighbors_view(ds.n_obs, art, p)
    n_comps = min(p.n_diffmap, ds.n_obs - 2)
    sc.tl.diffmap(view, n_comps=n_comps, random_state=p.seed)
    evecs: npt.NDArray[np.float64] = np.asarray(view.obsm["X_diffmap"], dtype=np.float64)
    evals = np.asarray(view.uns["diffmap_evals"], dtype=np.float64)
    # Eigenvector signs are arbitrary; fix them so runs agree exactly.
    signs = np.sign(evecs[np.argmax(np.abs(evecs), axis=0), np.arange(evecs.shape[1])])
    signs[signs == 0] = 1.0
    evecs = evecs * signs[None, :]
    return {
        "diffmap": evecs.astype(np.float32),
        "diffmap_evals": evals,
        "diffmap_ms": multiscale_space(evecs, evals),
    }


def _paga(ds: Dataset, art: Artifacts, p: PrepParams) -> dict[str, Any]:
    view = neighbors_view(ds.n_obs, art, p)
    labels = np.asarray(art["leiden"])
    view.obs["leiden"] = pd.Categorical(
        labels.astype(str), categories=[str(i) for i in range(labels.max() + 1)]
    )
    sc.tl.paga(view, groups="leiden")
    sc.pl.paga(view, plot=False, random_state=p.seed)
    return {
        "paga_connectivities": sparse.csr_matrix(view.uns["paga"]["connectivities"], dtype=np.float32),
        "paga_tree": sparse.csr_matrix(view.uns["paga"]["connectivities_tree"], dtype=np.float32),
        "paga_pos": np.asarray(view.uns["paga"]["pos"], dtype=np.float32),
    }


def _umap(ds: Dataset, art: Artifacts, p: PrepParams) -> dict[str, Any]:
    view = neighbors_view(ds.n_obs, art, p)
    labels = np.asarray(art["leiden"])
    view.obs["leiden"] = pd.Categorical(
        labels.astype(str), categories=[str(i) for i in range(labels.max() + 1)]
    )
    view.uns["paga"] = {
        "connectivities": art["paga_connectivities"],
        "connectivities_tree": art["paga_tree"],
        "groups": "leiden",
        "pos": np.asarray(art["paga_pos"]),
    }
    sc.tl.umap(view, init_pos="paga", random_state=p.seed)
    return {"umap": np.asarray(view.obsm["X_umap"], dtype=np.float32)}


def _needs_umap(ds: Dataset, p: PrepParams) -> bool:
    return p.recompute_umap or not ds.embeddings


STAGES: tuple[Stage, ...] = (
    Stage("normalize", (), ("lognorm",), (), _normalize),
    Stage("hvg", ("lognorm",), ("hvg",), ("n_top_genes",), _hvg),
    Stage("pca", ("lognorm", "hvg"), ("pca", "pca_variance_ratio"), ("n_pcs", "seed"), _pca),
    Stage(
        "neighbors",
        ("pca",),
        ("connectivities", "distances"),
        ("n_neighbors", "neighbors_pcs", "seed"),
        _neighbors,
    ),
    Stage(
        "leiden", ("pca", "connectivities", "distances"), ("leiden",), ("leiden_resolution", "seed"), _leiden
    ),
    Stage(
        "diffmap",
        ("pca", "connectivities", "distances"),
        ("diffmap", "diffmap_evals", "diffmap_ms"),
        ("n_diffmap", "seed"),
        _diffmap,
    ),
    Stage(
        "paga",
        ("pca", "connectivities", "distances", "leiden"),
        ("paga_connectivities", "paga_tree", "paga_pos"),
        ("seed",),
        _paga,
    ),
    Stage(
        "umap",
        ("pca", "connectivities", "distances", "leiden", "paga_connectivities", "paga_tree", "paga_pos"),
        ("umap",),
        ("seed", "recompute_umap"),
        _umap,
        applies=_needs_umap,
    ),
)


@dataclass(frozen=True)
class StageRecord:
    name: str
    state: str  # "done" | "skipped"
    seconds: float


def run_stage(stage: Stage, ds: Dataset, art: Artifacts, params: PrepParams) -> tuple[dict[str, Any], float]:
    missing = [k for k in stage.inputs if k not in art]
    if missing:
        raise RuntimeError(f"stage {stage.name!r} is missing inputs {missing}")
    start = time.perf_counter()
    out = stage.run(ds, {k: art[k] for k in stage.inputs}, params)
    elapsed = time.perf_counter() - start
    if set(out) != set(stage.outputs):
        raise RuntimeError(f"stage {stage.name!r} produced {sorted(out)}, declared {sorted(stage.outputs)}")
    return out, elapsed
