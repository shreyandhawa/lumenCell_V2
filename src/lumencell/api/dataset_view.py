"""What the server serves straight from one dataset: annotations, embeddings and expression.

The dataset itself never changes while the server runs, so anything derived
only from it is computed once. Prep artifacts are read from the on-disk cache
once preprocessing has finished; the server only ever reads them, and a stage
found missing or damaged is sent back to the worker.
"""

import sys
import threading
from collections.abc import Callable
from functools import cached_property
from http import HTTPStatus
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import sparse

from lumencell.api import binary, models
from lumencell.api.errors import ApiError
from lumencell.data.loader import Dataset
from lumencell.engine.context import LEIDEN_KEY, InferenceContext, detect_label_key
from lumencell.prep.cache import CacheMiss, PrepCache, prep_key
from lumencell.prep.pipeline import PrepResult
from lumencell.prep.stages import PrepParams

MAX_CATEGORIES_LISTED = 1000
MAX_GENES_PER_REQUEST = 200


def _is_numeric(s: pd.Series) -> bool:
    """Continuous values: numeric and not boolean or categorical."""
    return bool(pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s))


def _is_boolean(s: pd.Series) -> bool:
    """Plain bool and pandas' nullable boolean alike (a categorical of booleans stays categorical)."""
    return bool(not isinstance(s.dtype, pd.CategoricalDtype) and pd.api.types.is_bool_dtype(s))


class DatasetView:
    def __init__(
        self,
        path: Path,
        ds: Dataset,
        cache: PrepCache,
        params: PrepParams,
        on_cache_miss: Callable[[], object] | None = None,
    ) -> None:
        """``on_cache_miss`` is called when prep turns out not to be in the cache after all
        (the session starts a prepare job)."""
        self.path = path
        self.ds = ds
        self.cache = cache
        self.params = params
        self._on_cache_miss = on_cache_miss
        self._healing: set[str] = set()  # stages already sent back to the worker once
        self.label_key = detect_label_key(ds.obs)
        self._lock = threading.RLock()
        self._ready = False
        self._prep: PrepResult | None = None
        self._ctx: InferenceContext | None = None
        self._expr: sparse.csc_matrix | None = None

    # ---- prep --------------------------------------------------------------------------

    @cached_property
    def prep_key(self) -> str:
        return prep_key(self.ds, self.params)

    def prep_ready(self) -> bool:
        # Asked on nearly every request. Once every stage is in the cache the answer
        # is kept: the artifacts are loaded into memory from then on, and checking
        # again would stat the cache each time.
        if not self._ready:
            self._ready = all(s in ("done", "skipped") for s in self.cache.status(self.params).values())
        return self._ready

    def loaded_prep(self) -> PrepResult | None:
        """The prep artifacts, read from the cache once; None while they are not ready.

        A stage that turns out to be missing or unreadable is not computed here (that
        is the worker's job): the readiness latch is cleared, a prepare job is
        started, and the answer is None until that job has put the stage back. A
        stage that is still unreadable after that is not sent back again, so a disk
        that keeps damaging files cannot keep the worker busy for ever; preparing
        the dataset by hand retries it.
        """
        with self._lock:
            if self._prep is None and self.prep_ready():
                try:
                    self._prep = self.cache.load_all(self.params)
                except CacheMiss as exc:
                    self._ready = False
                    self._cache_missed(exc)
            return self._prep

    def _cache_missed(self, exc: CacheMiss) -> None:
        if exc.stage in self._healing:
            print(
                f"lumencell: {exc}, even after preparing it again; not trying again automatically "
                "(prepare the dataset to retry)",
                file=sys.stderr,
            )
        elif self._on_cache_miss is None:
            print(f"lumencell: {exc}", file=sys.stderr)
        else:
            print(f"lumencell: {exc}; preparing the dataset again in the background", file=sys.stderr)
            self._healing.add(exc.stage)
            self._on_cache_miss()

    def prep(self) -> PrepResult:
        prep = self.loaded_prep()
        if prep is None:
            raise ApiError(
                HTTPStatus.CONFLICT,
                "prep_not_ready",
                "The dataset is still being prepared. Wait for the prepare job to finish.",
            )
        return prep

    def ctx(self) -> InferenceContext:
        with self._lock:
            if self._ctx is None:
                self._ctx = InferenceContext.from_prep(self.ds, self.prep(), self.params)
            return self._ctx

    def obs_columns(self) -> list[models.ObsColumn]:
        out = list(self._file_columns)
        if self.loaded_prep() is not None:
            k = int(self.ctx().leiden.max()) + 1
            out.append(
                models.ObsColumn(name=LEIDEN_KEY, kind="categorical", categories=[str(i) for i in range(k)])
            )
        return out

    @cached_property
    def _file_columns(self) -> list[models.ObsColumn]:
        """The file's own columns. Listing them reads every value, so it is done once."""
        out: list[models.ObsColumn] = []
        for name in self.ds.obs.columns:
            s = self.ds.obs[name]
            colors = self.ds_colors(str(name))
            # The same test obs_payload() uses, so a column is never described one
            # way here and served the other way there.
            if not _is_numeric(s):
                cats = self._categories(s)
                kind = (
                    "boolean"
                    if _is_boolean(s)
                    else ("categorical" if len(cats) <= MAX_CATEGORIES_LISTED else "string")
                )
                out.append(
                    models.ObsColumn(
                        name=str(name),
                        kind=kind,
                        categories=cats if kind != "string" else None,
                        colors=colors if colors and len(colors) == len(cats) else None,
                    )
                )
            else:
                values = s.to_numpy(dtype=np.float64, na_value=np.nan)
                finite = values[np.isfinite(values)]
                out.append(
                    models.ObsColumn(
                        name=str(name),
                        kind="continuous",
                        min=float(finite.min()) if finite.size else None,
                        max=float(finite.max()) if finite.size else None,
                    )
                )
        return out

    @staticmethod
    def _categories(s: pd.Series) -> list[str]:
        if isinstance(s.dtype, pd.CategoricalDtype):
            return [str(c) for c in s.cat.categories]
        if _is_boolean(s):
            return ["False", "True"]
        return sorted({str(v) for v in s.dropna().unique()})

    def ds_colors(self, column: str) -> list[str] | None:
        return self._uns_colors.get(column)

    @cached_property
    def _uns_colors(self) -> dict[str, list[str]]:
        # The loader drops uns, so colours are read once from the file when asked for.
        import h5py

        colors: dict[str, list[str]] = {}
        with h5py.File(self.path, "r") as f:
            if "uns" in f:
                for key in f["uns"]:
                    if key.endswith("_colors"):
                        try:
                            colors[key[: -len("_colors")]] = [
                                v.decode() if isinstance(v, bytes) else str(v) for v in f["uns"][key][()]
                            ]
                        except (TypeError, KeyError):
                            continue
        return colors

    def embeddings(self) -> list[models.Embedding]:
        out = [models.Embedding(name=n, dims=2, source="file") for n in self.ds.embeddings]
        prep = self.loaded_prep()
        if prep is not None:
            art = prep.artifacts
            if "umap" in art and "umap" not in self.ds.embeddings:
                out.append(models.Embedding(name="umap", dims=2, source="computed"))
            out.append(models.Embedding(name="dc12", dims=2, source="diffusion"))
        return out

    def embedding(self, name: str) -> npt.NDArray[np.float32]:
        if name in self.ds.embeddings:
            return self.ds.embeddings[name]
        prep = self.loaded_prep()
        if prep is not None:
            art = prep.artifacts
            if name == "umap" and "umap" in art:
                return np.asarray(art["umap"], dtype=np.float32)
            if name == "dc12":
                return np.ascontiguousarray(np.asarray(art["diffmap"])[:, 1:3], dtype=np.float32)
        raise ApiError(HTTPStatus.NOT_FOUND, "not_found", f"no embedding named {name!r}")

    def default_embedding(self) -> str | None:
        names = [e.name for e in self.embeddings()]
        for preferred in ("umap", "tsne", "draw_graph_fa", "dc12"):
            if preferred in names:
                return preferred
        return names[0] if names else None

    def obs_payload(self, column: str) -> bytes:
        if column == LEIDEN_KEY:
            codes = self.ctx().leiden
            cats = [str(i) for i in range(int(codes.max()) + 1)]
            return binary.encode({"codes": ("i16", codes)}, {"kind": "categorical", "categories": cats})
        if column not in self.ds.obs.columns:
            raise ApiError(HTTPStatus.NOT_FOUND, "not_found", f"no obs column named {column!r}")
        s = self.ds.obs[column]
        if _is_numeric(s):
            values = s.to_numpy(dtype=np.float64, na_value=np.nan)
            return binary.encode({"values": ("f32", values)}, {"kind": "continuous"})
        cats = self._categories(s)
        lookup = {c: i for i, c in enumerate(cats)}
        codes = np.array([lookup.get(str(v), -1) if pd.notna(v) else -1 for v in s], dtype=np.int64)
        dtype = "i16" if len(cats) < 32767 else "i32"
        kind = "boolean" if _is_boolean(s) else "categorical"
        return binary.encode({"codes": (dtype, codes)}, {"kind": kind, "categories": cats})

    def search_genes(self, q: str, limit: int) -> models.GeneSearchResult:
        symbols = self.ds.gene_symbols
        needle = q.strip().lower()
        if not needle:
            matches = sorted(symbols.tolist())
        else:
            low = np.char.lower(symbols.astype(str))
            prefix = np.flatnonzero(np.char.startswith(low, needle))
            contains = np.flatnonzero((np.char.find(low, needle) >= 0) & ~np.char.startswith(low, needle))
            matches = sorted(symbols[prefix].tolist(), key=len) + sorted(symbols[contains].tolist())
        return models.GeneSearchResult(
            genes=matches[:limit], total=len(matches), truncated=len(matches) > limit
        )

    def expression_matrix(self) -> sparse.csc_matrix:
        with self._lock:
            if self._expr is None:
                lognorm = self.prep().artifacts["lognorm"]
                self._expr = sparse.csc_matrix(lognorm, dtype=np.float32)
            return self._expr

    @cached_property
    def _gene_lookup(self) -> dict[str, int]:
        """Lower-cased symbol or id -> column; symbols win, and the first of a duplicate wins."""
        lookup: dict[str, int] = {}
        for i, s in enumerate(self.ds.gene_symbols):
            lookup.setdefault(str(s).lower(), i)
        for i, s in enumerate(self.ds.var_names):
            lookup.setdefault(str(s).lower(), i)
        return lookup

    def gene_indices(self, genes: list[str]) -> list[int]:
        lookup = self._gene_lookup
        missing = [g for g in genes if g.lower() not in lookup]
        if missing:
            raise ApiError(HTTPStatus.NOT_FOUND, "unknown_gene", f"unknown gene(s): {', '.join(missing)}")
        return [lookup[g.lower()] for g in genes]

    def expression_payload(self, genes: list[str]) -> bytes:
        if not genes:
            raise ApiError(HTTPStatus.UNPROCESSABLE_ENTITY, "invalid_request", "ask for at least one gene")
        if len(genes) > MAX_GENES_PER_REQUEST:
            raise ApiError(
                HTTPStatus.UNPROCESSABLE_ENTITY,
                "too_many_genes",
                f"ask for at most {MAX_GENES_PER_REQUEST} genes per request",
            )
        idx = self.gene_indices(genes)
        m = self.expression_matrix()
        out = np.zeros((len(idx), self.ds.n_obs), dtype=np.float32)
        for row, j in enumerate(idx):
            start, end = m.indptr[j], m.indptr[j + 1]
            out[row, m.indices[start:end]] = m.data[start:end]
        symbols = [str(self.ds.gene_symbols[j]) for j in idx]
        return binary.encode({"expression": ("f32", out)}, {"genes": symbols})
