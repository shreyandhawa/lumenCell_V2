"""Write a run's results into a copy of the dataset (never the input itself)."""

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from lumencell.engine.context import LEIDEN_KEY
from lumencell.engine.schema import PRIMARY_SHARED, PRIMARY_UNASSIGNED, TrajectoryResult


def export_run(
    source: Path,
    out: Path,
    result: TrajectoryResult,
    run_id: str,
    leiden: np.ndarray | None = None,
    umap: np.ndarray | None = None,
    lineage_names: dict[int, str] | None = None,
) -> Path:
    if Path(out).resolve() == Path(source).resolve():
        raise ValueError("refusing to overwrite the input file; choose a different output path")
    adata = ad.read_h5ad(source)
    names = {lin.id: (lineage_names or {}).get(lin.id, lin.name) for lin in result.lineages}

    adata.obs["lumencell_pseudotime"] = result.pseudotime.astype(np.float32)
    primary = result.primary_lineage
    labels = np.where(
        primary == PRIMARY_SHARED,
        "shared",
        np.where(
            primary == PRIMARY_UNASSIGNED, "unassigned", np.array([names.get(int(p), "") for p in primary])
        ),
    )
    adata.obs["lumencell_lineage"] = pd.Categorical(
        labels, categories=["shared", *names.values(), "unassigned"]
    )
    for lin in result.lineages:
        adata.obs[f"lumencell_on_{lin.id}"] = result.membership[lin.id]
    if leiden is not None:
        adata.obs[LEIDEN_KEY] = pd.Categorical(leiden.astype(str))
    if umap is not None and "X_umap" not in adata.obsm:
        adata.obsm["X_umap"] = umap
    meta: dict[str, Any] = {
        "run_id": run_id,
        "method": result.method,
        "method_version": result.method_version,
        "root": result.root.spec,
        "root_cell": result.root.obs_name,
        "inference_rep": result.inference_rep,
        "lineages": {str(k): v for k, v in names.items()},
    }
    adata.uns["lumencell_run"] = meta
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    return out
