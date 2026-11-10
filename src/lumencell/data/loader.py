"""Load an .h5ad and work out what it contains.

Files in the wild differ in what ``X`` holds. lumenCell needs either raw counts
(preferred: it normalises them itself) or log-normalised values. Scaled data
(z-scores, with negative values) cannot be used for inference or for showing
expression, so it is only accepted when counts exist elsewhere in the file.

Where the matrix comes from, in order of preference:
1. ``X`` if it holds counts;
2. ``layers["counts"]`` if it holds counts;
3. ``raw.X`` if it holds counts (with ``raw.var`` as the gene list);
4. ``X`` if it holds log-normalised values.

The input file is opened read-only and never written.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import anndata as ad
import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import sparse

from lumencell.data.fingerprint import fingerprint

XKind = Literal["counts", "lognorm", "scaled"]
MatrixSource = Literal["X", "layers/counts", "raw"]

# log1p of any realistic normalised count stays far below this.
_LOGNORM_MAX = 30.0
_SAMPLE = 200_000

_ENSEMBL = re.compile(r"^ENS[A-Z]*G\d{6,}")
_SYMBOL_COLUMNS = ("feature_name", "gene_symbols", "gene_symbol", "gene_name", "gene_names", "symbol")


class DatasetError(ValueError):
    """The file cannot be used as it is; the message says what to do."""


def _values_sample(m: Any) -> npt.NDArray[np.float64]:
    """Nonzero values (sparse) or all values (dense), subsampled deterministically."""
    data = m.data if sparse.issparse(m) else np.asarray(m).ravel()
    data = np.asarray(data, dtype=np.float64)
    if data.size > _SAMPLE:
        idx = np.linspace(0, data.size - 1, _SAMPLE).astype(np.int64)
        data = data[idx]
    return data


def classify_matrix(m: Any) -> XKind | None:
    """counts / lognorm / scaled, or None for normalised-but-not-logged data."""
    v = _values_sample(m)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return "counts"
    if v.min() < 0:
        return "scaled"
    if np.all(v == np.round(v)):
        return "counts"
    if v.max() <= _LOGNORM_MAX:
        return "lognorm"
    return None


@dataclass(frozen=True)
class Dataset:
    path: Path
    name: str
    obs: pd.DataFrame
    matrix: sparse.csr_matrix  # cells x genes, float32, in the source's own scale
    x_kind: XKind  # "counts" or "lognorm" (the matrix's scale)
    source: MatrixSource
    var_names: pd.Index
    gene_symbols: npt.NDArray[np.str_]  # what the UI shows; == var_names unless a symbol column exists
    gene_symbol_column: str | None
    embeddings: dict[str, npt.NDArray[np.float32]]  # display embeddings found in the file
    uns_root: int | None  # uns["iroot"], if the file stores one
    fingerprint: str
    original_x_kind: XKind | None = field(default=None)

    @property
    def n_obs(self) -> int:
        return int(self.matrix.shape[0])

    @property
    def n_vars(self) -> int:
        return int(self.matrix.shape[1])


def _gene_symbols(var: pd.DataFrame, var_names: pd.Index) -> tuple[npt.NDArray[np.str_], str | None]:
    names = var_names.astype(str)
    looks_ensembl = np.mean([bool(_ENSEMBL.match(n)) for n in names[:500]]) > 0.5 if len(names) else False
    if looks_ensembl:
        for col in _SYMBOL_COLUMNS:
            if col in var.columns:
                symbols = var[col].astype(str).to_numpy()
                # Fall back to the ID wherever a symbol is missing.
                symbols = np.where((symbols == "") | (symbols == "nan"), names.to_numpy(), symbols)
                return symbols.astype(str), col
    return names.to_numpy().astype(str), None


def _display_embeddings(adata: ad.AnnData) -> dict[str, npt.NDArray[np.float32]]:
    out: dict[str, npt.NDArray[np.float32]] = {}
    for key in adata.obsm:
        arr = np.asarray(adata.obsm[key])
        if arr.ndim != 2 or arr.shape[1] < 2 or arr.shape[0] != adata.n_obs:
            continue
        if not np.issubdtype(arr.dtype, np.number):
            continue
        name = key[2:] if key.startswith("X_") else key
        out[name] = np.ascontiguousarray(arr[:, :2], dtype=np.float32)
    return out


def _uns_root(adata: ad.AnnData) -> int | None:
    value = adata.uns.get("iroot")
    if value is None:
        return None
    try:
        idx = int(np.asarray(value).item())
    except (TypeError, ValueError):
        return None
    return idx if 0 <= idx < adata.n_obs else None


def load(path: Path, x_is: XKind | None = None) -> Dataset:
    """Load ``path``. ``x_is`` overrides the detected kind of ``X`` (CLI ``--x-is``)."""
    path = Path(path)
    if not path.exists():
        raise DatasetError(f"{path} does not exist")
    adata = ad.read_h5ad(path)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise DatasetError(f"{path} has no cells or no genes")

    x_kind = x_is or classify_matrix(adata.X)
    candidates: list[tuple[MatrixSource, Any, pd.DataFrame, pd.Index]] = []
    if x_kind == "counts":
        candidates.append(("X", adata.X, adata.var, adata.var_names))
    if "counts" in adata.layers and classify_matrix(adata.layers["counts"]) == "counts":
        candidates.append(("layers/counts", adata.layers["counts"], adata.var, adata.var_names))
    if adata.raw is not None and classify_matrix(adata.raw.X) == "counts":
        candidates.append(("raw", adata.raw.X, adata.raw.var, adata.raw.var_names))

    if candidates:
        source, m, var, var_names = candidates[0]
        kind: XKind = "counts"
    elif x_kind == "lognorm":
        source, m, var, var_names = "X", adata.X, adata.var, adata.var_names
        kind = "lognorm"
    elif x_kind == "scaled":
        raise DatasetError(
            "X holds scaled values (it has negatives) and the file has no raw counts in "
            "layers['counts'] or .raw. lumenCell needs counts or log-normalised expression."
        )
    else:
        raise DatasetError(
            "X looks normalised but not log-transformed, and no raw counts were found. "
            "Pass --x-is lognorm if X is in fact log-normalised."
        )

    matrix = sparse.csr_matrix(m, dtype=np.float32)
    symbols, symbol_col = _gene_symbols(var, var_names)
    return Dataset(
        path=path,
        name=path.stem,
        obs=adata.obs.copy(),
        matrix=matrix,
        x_kind=kind,
        source=source,
        var_names=pd.Index(var_names.astype(str)),
        gene_symbols=symbols,
        gene_symbol_column=symbol_col,
        embeddings=_display_embeddings(adata),
        uns_root=_uns_root(adata),
        fingerprint=fingerprint(
            path, adata.obs_names.astype(str), var_names.astype(str), adata.n_obs, len(var_names)
        ),
        original_x_kind=x_kind,
    )
