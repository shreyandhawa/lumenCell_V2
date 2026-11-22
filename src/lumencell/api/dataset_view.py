"""What the server serves straight from one dataset: annotations, embeddings and expression.

The dataset itself never changes while the server runs, so anything derived
only from it is computed once. Prep artifacts are read from the on-disk cache
once preprocessing has finished; the server only ever reads them, and a stage
found missing or damaged is sent back to the worker.
"""

import sys
import threading
from collections.abc import Callable, Sequence
from functools import cached_property
from http import HTTPStatus
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import sparse

from lumencell.api import binary, models
from lumencell.api.errors import ApiError
from lumencell.clonal.isolation import column_names
from lumencell.data.loader import Dataset
from lumencell.engine.context import LEIDEN_KEY, InferenceContext, detect_label_key, inference_obs
from lumencell.engine.projection import diffusion_axes
from lumencell.engine.tree_support import tree_likeness
from lumencell.prep.cache import CacheMiss, PrepCache, prep_key
from lumencell.prep.pipeline import PrepResult
from lumencell.prep.stages import Artifacts, PrepParams
from lumencell.trust.cell_cycle import cell_cycle_info

MAX_CATEGORIES_LISTED = 1000
MAX_GENES_PER_REQUEST = 200


def _is_numeric(s: pd.Series) -> bool:
    """Continuous values: numeric and not boolean or categorical."""
    return bool(pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s))


def _is_boolean(s: pd.Series) -> bool:
    """Plain bool and pandas' nullable boolean alike (a categorical of booleans stays categorical)."""
    return bool(not isinstance(s.dtype, pd.CategoricalDtype) and pd.api.types.is_bool_dtype(s))


def column_kind(s: pd.Series, categories: list[str]) -> models.ObsKind:
    """How a column is listed and served, by one rule: ``categories`` are the values that occur."""
    if _is_numeric(s):
        return "continuous"
    if _is_boolean(s):
        return "boolean"
    return "categorical" if len(categories) <= MAX_CATEGORIES_LISTED else "string"


def _components(art: Artifacts) -> int:
    """How many diffusion components there are to draw: DC1..DCn (column 0 is the constant one)."""
    return int(np.asarray(art["diffmap"]).shape[1]) - 1


def paga_layout(
    coords: npt.NDArray[np.float32], labels: npt.NDArray[np.integer], paga: sparse.spmatrix
) -> tuple[list[models.GraphNode], list[models.GraphEdge]]:
    """Each cluster at the mean position of its cells, and PAGA's weighted edges between clusters.

    Cells without a finite position are left out of the means; a cluster with none is not drawn,
    and neither are its edges.
    """
    finite = np.isfinite(coords).all(axis=1)
    k = int(paga.shape[0])
    sizes = np.bincount(labels, minlength=k)
    drawn = np.bincount(labels[finite], minlength=k)
    sums = np.stack([np.bincount(labels[finite], weights=c, minlength=k) for c in coords[finite].T], axis=1)
    nodes = [
        models.GraphNode(
            id=c,
            x=float(sums[c, 0] / drawn[c]),
            y=float(sums[c, 1] / drawn[c]),
            z=float(sums[c, 2] / drawn[c]) if coords.shape[1] > 2 else None,
            n_cells=int(sizes[c]),
        )
        for c in range(k)
        if drawn[c] > 0
    ]
    upper = sparse.triu(paga, k=1).tocoo()
    edges = [
        models.GraphEdge(source=int(a), target=int(b), weight=float(v))
        for a, b, v in zip(upper.row, upper.col, upper.data, strict=True)
        if v > 0 and drawn[a] > 0 and drawn[b] > 0
    ]
    return nodes, edges


def knn_edges(
    connectivities: sparse.spmatrix, max_edges: int, seed: int
) -> tuple[npt.NDArray[np.uint32], npt.NDArray[np.uint32], npt.NDArray[np.float32], int]:
    """Each neighbour edge once (src < dst) with its weight, and how many there are in all.

    Above ``max_edges`` a seeded uniform sample is kept, in the graph's own order, so the same
    request draws the same edges.
    """
    upper = sparse.triu(connectivities.maximum(connectivities.T), k=1).tocoo()
    keep = upper.data > 0
    src, dst, w = upper.row[keep], upper.col[keep], upper.data[keep]
    total = len(src)
    if total > max_edges:
        pick = np.sort(np.random.default_rng(seed).choice(total, size=max_edges, replace=False))
        src, dst, w = src[pick], dst[pick], w[pick]
    return src.astype(np.uint32), dst.astype(np.uint32), w.astype(np.float32), total


class DatasetView:
    def __init__(
        self,
        path: Path,
        ds: Dataset,
        cache: PrepCache,
        params: PrepParams,
        on_cache_miss: Callable[[], object] | None = None,
        clonal_columns: tuple[str, ...] = (),
    ) -> None:
        """``on_cache_miss`` is called when prep turns out not to be in the cache after all
        (the session starts a prepare job)."""
        self.path = path
        self.ds = ds
        self.cache = cache
        self.params = params
        self._on_cache_miss = on_cache_miss
        self._healing: set[str] = set()  # stages already sent back to the worker once
        # Columns that are only for scoring runs besides the receptor-pattern ones: the clone and
        # donor columns the user named in the dataset's clonal settings.
        self.clonal_columns: tuple[str, ...] = column_names(clonal_columns)
        self.label_key = detect_label_key(inference_obs(ds.obs, self.clonal_columns)[0])
        self._lock = threading.RLock()
        self._ready = False
        self._prep: PrepResult | None = None
        # tree_likeness of the prep it was computed from; a different or dropped prep invalidates it.
        self._tree: tuple[PrepResult, float] | None = None
        self._ctx: InferenceContext | None = None

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
                    self._prep = self._by_gene(self.cache.load_all(self.params))
                except CacheMiss as exc:
                    self._ready = False
                    self._cache_missed(exc)
            return self._prep

    def tree_likeness(self) -> float | None:
        """C2 at the default resolution: None while prep is not loaded, NaN when PAGA has no edges.

        NaN is kept so the verdict fails C2 rather than reading it as pending; only the dataset
        route turns it into null.
        """
        prep = self.loaded_prep()
        if prep is None:
            return None
        with self._lock:
            if self._tree is None or self._tree[0] is not prep:
                value = tree_likeness(prep.artifacts["paga_connectivities"], prep.artifacts["paga_tree"])
                self._tree = (prep, value)
            return self._tree[1]

    def cell_cycle_info(self) -> models.CellCycleInfo | None:
        """How the cells were scored for cell cycle; None while prep is not loaded."""
        prep = self.loaded_prep()
        if prep is None:
            return None
        return cell_cycle_info(prep.artifacts["cc_info"], prep.artifacts["cc_phase"])

    @staticmethod
    def _by_gene(prep: PrepResult) -> PrepResult:
        """Keep log-normalised expression once, by gene (CSC), which is how every route reads it.

        The cache stores it by cell; the copy by cell is dropped once this one is made.
        """
        artifacts = dict(prep.artifacts)
        artifacts["lognorm"] = sparse.csc_matrix(artifacts["lognorm"], dtype=np.float32)
        return PrepResult(artifacts=artifacts, stages=prep.stages)

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

    def set_clonal(self, columns: Sequence[str]) -> None:
        """Name the clone and donor columns: the inference context is rebuilt, the label key found again."""
        named = column_names(columns)
        with self._lock:
            self.clonal_columns = named
            self._ctx = None
            self.label_key = detect_label_key(inference_obs(self.ds.obs, named)[0])

    def clonal_snapshot(self) -> tuple[tuple[str, ...], str | None]:
        """The clonal columns and the label key as one pair, never half of an update."""
        with self._lock:
            return self.clonal_columns, self.label_key

    def ctx(self) -> InferenceContext:
        with self._lock:
            if self._ctx is None:
                self._ctx = InferenceContext.from_prep(
                    self.ds, self.prep(), self.params, clonal_columns=self.clonal_columns
                )
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
            # column_kind() is the rule obs_payload() uses too, so a column is never described
            # one way here and served the other way there.
            if not _is_numeric(s):
                cats = self._categories(s)
                kind = column_kind(s, cats)
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
            n = _components(art)
            out.append(models.Embedding(name="dc12", dims=2, source="diffusion", components=n))
            if n >= 3:
                out.append(models.Embedding(name="dc123", dims=3, source="diffusion", components=n))
        return out

    def embedding(self, name: str) -> npt.NDArray[np.float32]:
        if name in self.ds.embeddings:
            return self.ds.embeddings[name]
        prep = self.loaded_prep()
        if prep is not None:
            art = prep.artifacts
            if name == "umap" and "umap" in art:
                return np.asarray(art["umap"], dtype=np.float32)
            axes = diffusion_axes(name)
            if axes is not None:
                n = _components(art)
                if len(set(axes)) < len(axes) or max(axes) > n:
                    raise ApiError(
                        HTTPStatus.NOT_FOUND,
                        "not_found",
                        f"no embedding named {name!r}: name two or three different diffusion components "
                        f"from 1 to {n}",
                    )
                return np.ascontiguousarray(np.asarray(art["diffmap"])[:, list(axes)], dtype=np.float32)
        raise ApiError(HTTPStatus.NOT_FOUND, "not_found", f"no embedding named {name!r}")

    # ---- graphs drawn over an embedding -------------------------------------------------

    def paga_graph(self, embedding: str) -> models.PagaGraph:
        art = self.prep().artifacts  # before the name, so "not prepared yet" is not read as "no such view"
        coords = self.embedding(embedding)
        nodes, edges = paga_layout(coords, np.asarray(art["leiden"]), art["paga_connectivities"])
        notices = []
        if not edges:
            notices.append(
                models.Notice(
                    code="clusters_do_not_touch",
                    severity="info",
                    message="No neighbour edge joins two Leiden clusters, so PAGA has no edges to draw.",
                )
            )
        return models.PagaGraph(
            embedding=embedding, cluster_key=LEIDEN_KEY, nodes=nodes, edges=edges, notices=notices
        )

    def knn_payload(self, max_edges: int) -> bytes:
        src, dst, w, total = knn_edges(self.prep().artifacts["connectivities"], max_edges, self.params.seed)
        return binary.encode(
            {"src": ("u32", src), "dst": ("u32", dst), "w": ("f32", w)},
            {"total_edges": total, "sampled": len(src) < total, "seed": self.params.seed},
        )

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
        return binary.encode({"codes": (dtype, codes)}, {"kind": column_kind(s, cats), "categories": cats})

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
        """Log-normalised expression, cells x genes, stored by gene (the one copy the server holds)."""
        matrix: sparse.csc_matrix = self.prep().artifacts["lognorm"]
        return matrix

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
