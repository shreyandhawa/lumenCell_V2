"""Server-side state for one loaded dataset.

The server process keeps the dataset (for annotations, embeddings and expression)
and reads prep artifacts and runs from the on-disk cache. Heavy work (prep,
inference) goes to the job worker. User choices (the root, lineage renames) live
in a small state file outside the cache, so they survive cache eviction.
"""

import json
import os
import re
import threading
from datetime import datetime
from http import HTTPStatus
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
import platformdirs
from pydantic import ValidationError
from scipy import sparse

from lumencell import __version__
from lumencell.api import binary, models
from lumencell.api.errors import ApiError, root_required
from lumencell.data.loader import Dataset, load
from lumencell.engine.compare import RunComparison, compare
from lumencell.engine.context import LEIDEN_KEY, InferenceContext, detect_label_key
from lumencell.engine.projection import project
from lumencell.engine.registry import METHODS, get_method
from lumencell.engine.roots import (
    CellRoot,
    ClusterRoot,
    FileRoot,
    RootError,
    RootSpec,
    format_root,
    parse_root,
    resolve_root,
    suggest_roots,
)
from lumencell.engine.schema import ResolvedRoot, TrajectoryResult
from lumencell.engine.store import RunStore, run_id
from lumencell.engine.trends import gene_trends, temporal_heatmap
from lumencell.jobs.manager import Job, JobManager
from lumencell.prep.cache import PrepCache, dataset_dir, prep_key
from lumencell.prep.pipeline import PrepResult
from lumencell.prep.stages import STAGES, PrepParams

MAX_CATEGORIES_LISTED = 1000
MAX_GENES_PER_REQUEST = 200
# Run ids name folders in the run store, so anything that is not a plain name is not a run.
_RUN_ID = re.compile(r"[A-Za-z0-9_-]{1,64}")


def _is_numeric(s: pd.Series) -> bool:
    """Continuous values: numeric and not boolean or categorical."""
    return bool(pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s))


def default_state_dir() -> Path:
    override = os.environ.get("LUMENCELL_DATA_DIR")
    return Path(override) if override else Path(platformdirs.user_data_dir("lumencell")) / "state"


# ---- conversions between engine and API types -------------------------------------------


def to_api_spec(spec: RootSpec) -> models.CellRoot | models.ClusterRoot | models.FileRoot:
    if isinstance(spec, CellRoot):
        return models.CellRoot(kind="cell", obs_name=spec.obs_name)
    if isinstance(spec, ClusterRoot):
        return models.ClusterRoot(kind="cluster", key=spec.key, value=spec.value)
    return models.FileRoot(kind="file")


def from_api_spec(spec: models.CellRoot | models.ClusterRoot | models.FileRoot) -> RootSpec:
    if isinstance(spec, models.CellRoot):
        return CellRoot(spec.obs_name)
    if isinstance(spec, models.ClusterRoot):
        return ClusterRoot(spec.key, spec.value)
    return FileRoot()


def to_api_root(root: ResolvedRoot | dict[str, Any]) -> models.ResolvedRoot:
    r = root if isinstance(root, dict) else root.__dict__
    warnings = [w if isinstance(w, dict) else w.__dict__ for w in r.get("warnings", ())]
    return models.ResolvedRoot(
        spec=to_api_spec(parse_root(r["spec"])),
        cell_index=r["cell_index"],
        obs_name=r["obs_name"],
        cluster_key=r.get("cluster_key"),
        cluster_label=r.get("cluster_label"),
        resolved_by=r["resolved_by"],
        warnings=[models.Notice(**w) for w in warnings],
    )


class Session:
    def __init__(
        self,
        path: Path,
        jobs: JobManager,
        x_is: str | None = None,
        prep_params: PrepParams | None = None,
        cache_root: Path | None = None,
        state_dir: Path | None = None,
    ) -> None:
        self.path = Path(path).resolve()
        self.x_is = x_is
        self.jobs = jobs
        self.params = prep_params or PrepParams()
        self.cache_root = cache_root
        self.ds: Dataset = load(self.path, x_is)  # type: ignore[arg-type]
        self.cache = PrepCache(self.ds, cache_root)
        self.store = RunStore(dataset_dir(self.ds.fingerprint, cache_root))
        self.state_path = (state_dir or default_state_dir()) / f"{self.ds.fingerprint}.json"
        self.label_key = detect_label_key(self.ds.obs)
        self._lock = threading.RLock()
        self._prep: PrepResult | None = None
        self._ctx: InferenceContext | None = None
        self._expr: sparse.csc_matrix | None = None
        self._heatmaps: dict[tuple[object, ...], bytes] = {}

    # ---- persisted user state --------------------------------------------------------

    def _state(self) -> dict[str, Any]:
        if self.state_path.exists():
            data: dict[str, Any] = json.loads(self.state_path.read_text(encoding="utf-8"))
            return data
        return {}

    def _save_state(self, state: dict[str, Any]) -> None:
        self.state_path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.state_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(state, indent=2), encoding="utf-8")
        tmp.replace(self.state_path)

    # ---- prep --------------------------------------------------------------------------

    @property
    def prep_key(self) -> str:
        return prep_key(self.ds, self.params)

    def _job_payload(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "x_is": self.x_is,
            "prep_params": self.params,
            "cache_root": str(self.cache_root) if self.cache_root else None,
        }

    def prep_ready(self) -> bool:
        return all(s in ("done", "skipped") for s in self.cache.status(self.params).values())

    def prep_status(self) -> models.PrepStatus:
        status = self.cache.status(self.params)
        job = self.jobs.find(f"prepare:{self.prep_key}")
        running = job is not None and job.state in ("queued", "running")
        manifest = self.cache.manifest() or {}
        seconds = {k: v.get("seconds") for k, v in manifest.get("stages", {}).items()}
        stages = []
        for stage in STAGES:
            state = status[stage.name]
            if running and job is not None and state == "pending" and job.stage == stage.name:
                state = "running"
            if job is not None and job.state == "failed" and state == "pending":
                state = "failed"
            stages.append(models.PrepStage(name=stage.name, state=state, seconds=seconds.get(stage.name)))
        if self.prep_ready():
            overall = "ready"
        elif running:
            overall = "running"
        elif job is not None and job.state == "failed":
            overall = "failed"
        else:
            overall = "missing"
        return models.PrepStatus(
            state=overall,
            stages=stages,
            job_id=job.job_id if job is not None else None,
        )

    def start_prep(self) -> Job:
        return self.jobs.submit("prepare", f"prepare:{self.prep_key}", self._job_payload())

    def prep(self) -> PrepResult:
        with self._lock:
            if self._prep is None:
                if not self.prep_ready():
                    raise ApiError(
                        HTTPStatus.CONFLICT,
                        "prep_not_ready",
                        "The dataset is still being prepared. Wait for the prepare job to finish.",
                    )
                self._prep = self.cache.run(self.params)
            return self._prep

    def ctx(self) -> InferenceContext:
        with self._lock:
            if self._ctx is None:
                self._ctx = InferenceContext.from_prep(self.ds, self.prep(), self.params)
            return self._ctx

    # ---- dataset -------------------------------------------------------------------------

    def _obs_columns(self) -> list[models.ObsColumn]:
        out: list[models.ObsColumn] = []
        for name in self.ds.obs.columns:
            s = self.ds.obs[name]
            colors = self.ds_colors(str(name))
            if (
                isinstance(s.dtype, pd.CategoricalDtype)
                or s.dtype == bool
                or s.dtype == object
                or s.dtype == "string"
            ):
                cats = self._categories(s)
                kind = (
                    "boolean"
                    if s.dtype == bool
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
        if self.prep_ready():
            k = int(self.ctx().leiden.max()) + 1
            out.append(
                models.ObsColumn(name=LEIDEN_KEY, kind="categorical", categories=[str(i) for i in range(k)])
            )
        return out

    @staticmethod
    def _categories(s: pd.Series) -> list[str]:
        if isinstance(s.dtype, pd.CategoricalDtype):
            return [str(c) for c in s.cat.categories]
        if s.dtype == bool:
            return ["False", "True"]
        return sorted({str(v) for v in s.dropna().unique()})

    def ds_colors(self, column: str) -> list[str] | None:
        # The loader drops uns, so colours are read once from the file when asked for.
        if not hasattr(self, "_uns_colors"):
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
            self._uns_colors = colors
        return self._uns_colors.get(column)

    def embeddings(self) -> list[models.Embedding]:
        out = [models.Embedding(name=n, dims=2, source="file") for n in self.ds.embeddings]
        if self.prep_ready():
            art = self.prep().artifacts
            if "umap" in art and "umap" not in self.ds.embeddings:
                out.append(models.Embedding(name="umap", dims=2, source="computed"))
            out.append(models.Embedding(name="dc12", dims=2, source="diffusion"))
        return out

    def embedding(self, name: str) -> npt.NDArray[np.float32]:
        if name in self.ds.embeddings:
            return self.ds.embeddings[name]
        if self.prep_ready():
            art = self.prep().artifacts
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

    def dataset_info(self) -> models.DatasetInfo:
        return models.DatasetInfo(
            schema_version=models.SCHEMA_VERSION,
            name=self.ds.name,
            n_obs=self.ds.n_obs,
            n_vars=self.ds.n_vars,
            x_kind=self.ds.x_kind,
            gene_symbol_column=self.ds.gene_symbol_column,
            obs=self._obs_columns(),
            embeddings=self.embeddings(),
            default_embedding=self.default_embedding(),
            cluster_key=LEIDEN_KEY if self.prep_ready() else None,
            label_key=self.label_key,
            prep=self.prep_status(),
            root=self.root_state().root,
            versions={"lumencell": __version__},
        )

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
        kind = "boolean" if s.dtype == bool else "categorical"
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

    def _expression_matrix(self) -> sparse.csc_matrix:
        with self._lock:
            if self._expr is None:
                lognorm = self.prep().artifacts["lognorm"]
                self._expr = sparse.csc_matrix(lognorm, dtype=np.float32)
            return self._expr

    def gene_indices(self, genes: list[str]) -> list[int]:
        lookup: dict[str, int] = {}
        for i, s in enumerate(self.ds.gene_symbols):
            lookup.setdefault(str(s).lower(), i)
        for i, s in enumerate(self.ds.var_names):
            lookup.setdefault(str(s).lower(), i)
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
        m = self._expression_matrix()
        out = np.zeros((len(idx), self.ds.n_obs), dtype=np.float32)
        for row, j in enumerate(idx):
            start, end = m.indptr[j], m.indptr[j + 1]
            out[row, m.indices[start:end]] = m.data[start:end]
        symbols = [str(self.ds.gene_symbols[j]) for j in idx]
        return binary.encode({"expression": ("f32", out)}, {"genes": symbols})

    # ---- methods -------------------------------------------------------------------------

    def methods(self) -> models.MethodList:
        out = []
        for m in METHODS.values():
            avail = m.available()
            out.append(
                models.MethodInfo(
                    name=m.name,
                    label=m.label,
                    family=m.family,
                    description=m.description,
                    available=avail.ok,
                    install_hint=avail.install_hint,
                    experimental=m.experimental,
                    experimental_note=m.experimental_note,
                    params_schema=m.Params.model_json_schema(),
                    reference_scores={},
                )
            )
        return models.MethodList(methods=out)

    # ---- root ----------------------------------------------------------------------------

    def _resolve(self, spec: RootSpec) -> ResolvedRoot:
        try:
            return resolve_root(self.ctx(), spec)
        except RootError as exc:
            raise ApiError(HTTPStatus.UNPROCESSABLE_ENTITY, "invalid_root", str(exc)) from exc

    def root_state(self) -> models.RootState:
        spec_text = self._state().get("root")
        if not spec_text or not self.prep_ready():
            return models.RootState(root=None)
        try:
            return models.RootState(root=to_api_root(resolve_root(self.ctx(), parse_root(spec_text))))
        except RootError:
            return models.RootState(root=None)

    def resolve_root(
        self, spec: models.CellRoot | models.ClusterRoot | models.FileRoot
    ) -> models.ResolvedRoot:
        return to_api_root(self._resolve(from_api_spec(spec)))

    def set_root(self, spec: models.CellRoot | models.ClusterRoot | models.FileRoot) -> models.ResolvedRoot:
        resolved = self._resolve(from_api_spec(spec))
        state = self._state()
        state["root"] = resolved.spec
        self._save_state(state)
        return to_api_root(resolved)

    def root_suggestions(self) -> models.RootSuggestions:
        out = []
        for s in suggest_roots(self.ctx()):
            spec = to_api_spec(s.spec)
            assert not isinstance(spec, models.CellRoot)
            out.append(
                models.RootSuggestion(
                    spec=spec,
                    source=s.source,
                    score=s.score,
                    evidence=s.evidence,
                    is_extremity=s.is_extremity,
                )
            )
        return models.RootSuggestions(suggestions=out)

    # ---- runs ----------------------------------------------------------------------------

    def create_run(self, req: models.TrajectoryRequest) -> tuple[models.RunCreated, bool]:
        """Returns (response, created). created=False means an identical run already exists."""
        try:
            method = get_method(req.method)
        except KeyError as exc:
            raise ApiError(HTTPStatus.NOT_FOUND, "unknown_method", str(exc.args[0])) from exc
        if not method.available().ok:
            raise ApiError(
                HTTPStatus.CONFLICT,
                "method_unavailable",
                f"{method.label} is not installed. {method.available().install_hint or ''}".strip(),
            )
        if req.root is not None:
            spec_text = format_root(from_api_spec(req.root))
        else:
            spec_text = self._state().get("root") or ""
            if not spec_text:
                raise root_required()
        try:
            params = method.Params.model_validate(req.params).model_dump()
        except ValidationError as exc:
            errors = [
                {"loc": [str(part) for part in e["loc"]], "msg": str(e["msg"])}
                for e in exc.errors(include_url=False, include_context=False, include_input=False)
            ]
            raise ApiError(
                HTTPStatus.UNPROCESSABLE_ENTITY,
                "invalid_params",
                f"The parameters for {method.label} are not valid.",
                {"errors": errors},
            ) from exc
        try:
            # The root travels to the worker as text, so it has to survive the round trip.
            spec = parse_root(spec_text)
        except RootError as exc:
            raise ApiError(HTTPStatus.UNPROCESSABLE_ENTITY, "invalid_root", str(exc)) from exc
        self._resolve(spec)  # fail fast on a bad root
        rid = run_id(req.method, params, spec_text, req.terminals, self.prep_key)
        if self.store.exists(rid):
            return models.RunCreated(run_id=rid, job_id=None, state="done"), False
        payload = {
            **self._job_payload(),
            "method": req.method,
            "root": spec_text,
            "terminals": req.terminals,
            "params": params,
            "run_id": rid,
            "label_key": self.label_key,
        }
        job = self.jobs.submit("infer", f"infer:{rid}", payload)
        state = "running" if job.state == "running" else "queued"
        return models.RunCreated(run_id=rid, job_id=job.job_id, state=state), True

    def _renames(self, rid: str) -> dict[str, str]:
        renames: dict[str, str] = self._state().get("lineage_names", {}).get(rid, {})
        return renames

    def list_runs(self) -> models.RunList:
        runs = [
            models.RunSummary(
                run_id=r["run_id"],
                method=r["method"],
                state="done",
                created_at=datetime.fromisoformat(r["created_at"]),
                n_lineages=r["n_lineages"],
            )
            for r in self.store.list()
        ]
        done = {r.run_id for r in runs}
        for job in self.jobs_for_runs():
            rid = job.payload["run_id"]
            if rid in done:
                continue
            runs.append(
                models.RunSummary(
                    run_id=rid,
                    method=job.payload["method"],
                    state=job.state,
                    created_at=datetime.fromtimestamp(job.created).astimezone(),
                    error=job.error,
                )
            )
        return models.RunList(runs=runs)

    def jobs_for_runs(self) -> list[Job]:
        with self.jobs._lock:
            return [j for j in self.jobs._jobs.values() if j.kind == "infer"]

    def _require_run(self, rid: str) -> dict[str, Any]:
        if not _RUN_ID.fullmatch(rid):
            raise ApiError(HTTPStatus.NOT_FOUND, "not_found", f"no run {rid!r}")
        if not self.store.exists(rid):
            job = self.jobs.find(f"infer:{rid}")
            if job is not None:
                raise ApiError(
                    HTTPStatus.CONFLICT,
                    "run_not_ready",
                    f"run {rid} is {job.state}" + (f": {job.error}" if job.error else ""),
                    {"state": job.state, "job_id": job.job_id},
                )
            raise ApiError(HTTPStatus.NOT_FOUND, "not_found", f"no run {rid!r}")
        return self.store.graph(rid)

    def get_run(self, rid: str) -> models.TrajectoryRun:
        g = self._require_run(rid)
        renames = self._renames(rid)
        lineages = [
            models.Lineage(**{**lin, "name": renames.get(str(lin["id"]), lin["name"])})
            for lin in g["lineages"]
        ]
        root = to_api_root(g["root"])
        terminals = [lin.terminal_label for lin in lineages if lin.terminal_label]
        return models.TrajectoryRun(
            run_id=rid,
            state="done",
            topology_kind=g["topology_kind"],
            milestones=[models.Milestone(**m) for m in g["milestones"]],
            segments=[models.Segment(**s) for s in g["segments"]],
            lineages=lineages,
            divergence=[models.Divergence(**d) for d in g["divergence"]],
            warnings=[models.Notice(**w) for w in g["warnings"]],
            n_unassigned=g["n_unassigned"],
            runtime_seconds=g["runtime_seconds"],
            cell_arrays=["pseudotime", "membership", "primary_lineage", "segment", "progress"]
            + (["fate"] if g["has_fate"] else [])
            + [f"extras.{e}" for e in g["extras"]],
            provenance=models.Provenance(
                method=g["method"],
                method_version=g["method_version"],
                lumencell_version=g["lumencell_version"],
                params=g["params"],
                param_hash=rid,
                dataset_fingerprint=self.ds.fingerprint,
                inference_rep=g["inference_rep"],
                inference_dims=g["inference_dims"],
                root=root,
                terminals=terminals,
                terminals_source="auto",
                created_at=datetime.fromisoformat(g["created_at"]),
            ),
        )

    def run_cells(self, rid: str) -> bytes:
        self._require_run(rid)
        return self.store.cells_bytes(rid)

    def run_curves(self, rid: str, embedding: str) -> models.Curves:
        self._require_run(rid)
        result = self.store.load(rid)
        proj = project(result, self.embedding(embedding), embedding)
        return models.Curves(
            run_id=rid,
            embedding=embedding,
            curve_source=proj.curve_source,
            curves=[
                models.LineageCurve(lineage=c.lineage, points=c.points, fidelity=c.fidelity)
                for c in proj.curves
            ],
            milestones=[models.MilestonePosition(id=k, x=v[0], y=v[1]) for k, v in proj.milestones.items()],
            warnings=[models.Notice(**w.__dict__) for w in proj.warnings],
        )

    def rename_lineage(self, rid: str, lineage_id: int, name: str) -> models.Lineage:
        g = self._require_run(rid)
        lineages = {lin["id"]: lin for lin in g["lineages"]}
        if lineage_id not in lineages:
            raise ApiError(HTTPStatus.NOT_FOUND, "not_found", f"run {rid} has no lineage {lineage_id}")
        current = self._renames(rid)
        taken = {current.get(str(i), lin["name"]) for i, lin in lineages.items() if i != lineage_id}
        clean = name.strip()
        if clean in taken:
            raise ApiError(HTTPStatus.CONFLICT, "name_taken", f"another lineage is already called {clean!r}")
        state = self._state()
        state.setdefault("lineage_names", {}).setdefault(rid, {})[str(lineage_id)] = clean
        self._save_state(state)
        return models.Lineage(**{**lineages[lineage_id], "name": clean})

    # ---- trends, heatmap, comparison --------------------------------------------------------

    def _lineage(self, rid: str, lineage: int) -> tuple[TrajectoryResult, np.ndarray]:
        self._require_run(rid)
        result = self.store.load(rid)
        if not 0 <= lineage < len(result.lineages):
            raise ApiError(HTTPStatus.NOT_FOUND, "not_found", f"run {rid} has no lineage {lineage}")
        return result, result.membership[lineage]

    def trends(self, rid: str, lineage: int, genes: list[str]) -> models.Trends:
        if not genes:
            raise ApiError(HTTPStatus.UNPROCESSABLE_ENTITY, "invalid_request", "ask for at least one gene")
        if len(genes) > 50:
            raise ApiError(HTTPStatus.UNPROCESSABLE_ENTITY, "too_many_genes", "ask for at most 50 genes")
        result, members = self._lineage(rid, lineage)
        idx = self.gene_indices(genes)
        t = gene_trends(self._expression_matrix(), idx, result.pseudotime, members)
        return models.Trends(
            run_id=rid,
            lineage=lineage,
            x=[float(v) for v in t.centers],
            counts=[int(c) for c in t.counts],
            genes=[
                models.GeneTrend(
                    gene=str(self.ds.gene_symbols[j]),
                    y=[None if not np.isfinite(v) else float(v) for v in t.smooth[:, k]],
                    r2=float(t.r2[k]),
                )
                for k, j in enumerate(idx)
            ],
        )

    def heatmap_payload(self, rid: str, req: models.HeatmapRequest) -> bytes:
        result, members = self._lineage(rid, req.lineage)
        key = (rid, req.lineage, req.top_n, tuple(req.genes or ()))
        with self._lock:
            cached = self._heatmaps.get(key)
        if cached is not None:
            return cached
        idx = self.gene_indices(req.genes) if req.genes else None
        h = temporal_heatmap(
            self._expression_matrix(), result.pseudotime, members, top_n=req.top_n, genes=idx
        )
        body = binary.encode(
            {"heatmap": ("f32", h.values)},
            {
                "genes": [str(self.ds.gene_symbols[j]) for j in h.genes],
                "r2": [round(float(v), 4) for v in h.r2],
                "bins": [round(float(c), 4) for c in h.centers],
            },
        )
        with self._lock:
            self._heatmaps[key] = body
        return body

    def _comparison(self, a: str, b: str) -> RunComparison:
        self._require_run(a)
        self._require_run(b)
        return compare(self.store.load(a), self.store.load(b))

    def compare(self, a: str, b: str) -> models.Comparison:
        c = self._comparison(a, b)
        return models.Comparison(
            runs=(a, b),
            n_common=c.n_common,
            spearman=c.spearman,
            lineage_matches=[models.LineageMatch(a=p.a, b=p.b, jaccard=p.jaccard) for p in c.pairs],
        )

    def compare_cells(self, a: str, b: str) -> bytes:
        return binary.encode(
            {"rank_difference": ("f32", self._comparison(a, b).rank_difference)}, {"runs": [a, b]}
        )

    # ---- jobs ----------------------------------------------------------------------------

    @staticmethod
    def job_info(job: Job) -> models.JobInfo:
        kind = job.kind if job.kind in ("prepare", "infer", "heatmap") else "prepare"
        return models.JobInfo(
            job_id=job.job_id,
            kind=kind,
            state=job.state,
            stage=job.stage,
            fraction=job.fraction,
            message=job.message,
            error=job.error,
            result_id=job.result_id or job.payload.get("run_id"),
        )

    def require_job(self, job_id: str) -> Job:
        job = self.jobs.get(job_id)
        if job is None:
            raise ApiError(HTTPStatus.NOT_FOUND, "not_found", f"no job {job_id!r}")
        return job
