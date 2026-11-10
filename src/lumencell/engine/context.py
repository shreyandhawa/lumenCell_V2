"""What a method sees: the dataset's annotations plus the shared prep artifacts."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import sparse

from lumencell.data.loader import Dataset
from lumencell.prep.pipeline import PrepResult
from lumencell.prep.stages import PrepParams

LEIDEN_KEY = "lumencell_leiden"

# Annotation columns commonly holding cell types or clusters, in order of preference.
# Used to name lineages and to phrase root suggestions; never for inference.
LABEL_CANDIDATES = (
    "cell_type",
    "celltype",
    "cell_types",
    "cell_ontology_class",
    "annotation",
    "annotations",
    "clusters",
    "cluster",
    "louvain",
    "leiden",
)


def detect_label_key(obs: pd.DataFrame) -> str:
    lowered = {c.lower(): c for c in obs.columns}
    for name in LABEL_CANDIDATES:
        col = lowered.get(name)
        if col is None:
            continue
        n_unique = obs[col].nunique(dropna=True)
        if 2 <= n_unique <= 200:
            return str(col)
    return LEIDEN_KEY


@dataclass(frozen=True)
class InferenceContext:
    n_obs: int
    obs: pd.DataFrame  # the file's obs plus LEIDEN_KEY
    obs_names: pd.Index
    artifacts: Mapping[str, Any]
    prep_params: PrepParams
    counts: sparse.csr_matrix | None  # raw counts when the file has them
    uns_root: int | None
    fingerprint: str
    label_key: str

    @classmethod
    def from_prep(cls, ds: Dataset, prep: PrepResult, params: PrepParams) -> "InferenceContext":
        leiden = np.asarray(prep.artifacts["leiden"])
        obs = ds.obs.copy()
        obs[LEIDEN_KEY] = pd.Categorical(
            leiden.astype(str), categories=[str(i) for i in range(int(leiden.max()) + 1)]
        )
        return cls(
            n_obs=ds.n_obs,
            obs=obs,
            obs_names=ds.obs.index.astype(str),
            artifacts=prep.artifacts,
            prep_params=params,
            counts=ds.matrix if ds.x_kind == "counts" else None,
            uns_root=ds.uns_root,
            fingerprint=ds.fingerprint,
            label_key=detect_label_key(ds.obs),
        )

    def labels(self, key: str) -> npt.NDArray[np.str_]:
        if key not in self.obs.columns:
            raise KeyError(key)
        out: npt.NDArray[np.str_] = self.obs[key].astype(str).to_numpy()
        return out

    @property
    def leiden(self) -> npt.NDArray[np.int32]:
        return np.asarray(self.artifacts["leiden"], dtype=np.int32)

    def diffusion_space(self, n_dcs: int = 10) -> npt.NDArray[np.float64]:
        """Diffusion components 1..n_dcs (the trivial first one dropped), eigenvalue-scaled."""
        dm = np.asarray(self.artifacts["diffmap"], dtype=np.float64)
        evals = np.asarray(self.artifacts["diffmap_evals"], dtype=np.float64)
        k = min(n_dcs, dm.shape[1] - 1)
        out: npt.NDArray[np.float64] = dm[:, 1 : k + 1] * evals[1 : k + 1][None, :]
        return out
