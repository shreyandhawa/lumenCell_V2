"""Trajectory runs for one dataset: starting them, describing them, and what is derived from them."""

import re
import threading
from collections import OrderedDict
from collections.abc import Callable
from datetime import datetime
from http import HTTPStatus
from typing import Any

import numpy as np
from pydantic import ValidationError

from lumencell.api import binary, models
from lumencell.api.convert import from_api_spec, to_api_root
from lumencell.api.dataset_view import DatasetView
from lumencell.api.errors import ApiError, root_required
from lumencell.api.state import StateStore
from lumencell.engine.compare import RunComparison, compare
from lumencell.engine.projection import project
from lumencell.engine.registry import get_method
from lumencell.engine.roots import RootError, RootSpec, format_root, parse_root
from lumencell.engine.schema import ResolvedRoot, TrajectoryResult
from lumencell.engine.store import RunStore, run_id
from lumencell.engine.trends import gene_trends, temporal_heatmap
from lumencell.jobs.manager import JobManager

# Run ids name folders in the run store, so anything that is not a plain name is not a run.
_RUN_ID = re.compile(r"[A-Za-z0-9_-]{1,64}")
# Curves, trends, heatmaps and comparisons all start from a whole run read from
# disk. A run never changes once saved, so the last few are kept in memory.
MAX_RUNS_IN_MEMORY = 4
MAX_HEATMAPS_IN_MEMORY = 32


class Runs:
    def __init__(
        self,
        view: DatasetView,
        store: RunStore,
        jobs: JobManager,
        state: StateStore,
        job_payload: Callable[[], dict[str, Any]],
        resolve: Callable[[RootSpec], ResolvedRoot],
    ) -> None:
        self.view = view
        self.store = store
        self.jobs = jobs
        self.state = state
        self._job_payload = job_payload
        self._resolve = resolve
        self._lock = threading.Lock()
        self._loaded: OrderedDict[str, TrajectoryResult] = OrderedDict()
        self._heatmaps: OrderedDict[tuple[object, ...], bytes] = OrderedDict()

    def _load(self, rid: str) -> TrajectoryResult:
        """The whole run. Callers must not change it: it is shared between requests."""
        with self._lock:
            result = self._loaded.get(rid)
            if result is not None:
                self._loaded.move_to_end(rid)
                return result
        result = self.store.load(rid)
        with self._lock:
            self._loaded[rid] = result
            self._loaded.move_to_end(rid)
            while len(self._loaded) > MAX_RUNS_IN_MEMORY:
                self._loaded.popitem(last=False)
        return result

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
            spec_text = self.state.read().get("root") or ""
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
        rid = run_id(req.method, params, spec_text, req.terminals, self.view.prep_key)
        if self.store.exists(rid):
            return models.RunCreated(run_id=rid, job_id=None, state="done"), False
        payload = {
            **self._job_payload(),
            "method": req.method,
            "root": spec_text,
            "terminals": req.terminals,
            "params": params,
            "run_id": rid,
            "label_key": self.view.label_key,
        }
        job = self.jobs.submit("infer", f"infer:{rid}", payload)
        state = "running" if job.state == "running" else "queued"
        return models.RunCreated(run_id=rid, job_id=job.job_id, state=state), True

    def _renames(self, rid: str) -> dict[str, str]:
        renames: dict[str, str] = self.state.read().get("lineage_names", {}).get(rid, {})
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
        for job in self.jobs.of_kind("infer"):
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
                dataset_fingerprint=self.view.ds.fingerprint,
                inference_rep=g["inference_rep"],
                inference_dims=g["inference_dims"],
                root=root,
                terminals=terminals,
                # Runs saved before this was recorded count as "auto".
                terminals_source="user" if g.get("terminals_source") == "user" else "auto",
                created_at=datetime.fromisoformat(g["created_at"]),
            ),
        )

    def run_cells(self, rid: str) -> bytes:
        self._require_run(rid)
        return self.store.cells_bytes(rid)

    def run_curves(self, rid: str, embedding: str) -> models.Curves:
        self._require_run(rid)
        result = self._load(rid)
        proj = project(result, self.view.embedding(embedding), embedding)
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
        if not clean:
            raise ApiError(
                HTTPStatus.UNPROCESSABLE_ENTITY, "invalid_request", "a lineage name cannot be blank"
            )
        if clean in taken:
            raise ApiError(HTTPStatus.CONFLICT, "name_taken", f"another lineage is already called {clean!r}")

        def change(state: dict[str, Any]) -> None:
            state.setdefault("lineage_names", {}).setdefault(rid, {})[str(lineage_id)] = clean

        self.state.update(change)
        return models.Lineage(**{**lineages[lineage_id], "name": clean})

    # ---- trends, heatmap, comparison --------------------------------------------------------

    def _lineage(self, rid: str, lineage: int) -> tuple[TrajectoryResult, np.ndarray]:
        self._require_run(rid)
        result = self._load(rid)
        if not 0 <= lineage < len(result.lineages):
            raise ApiError(HTTPStatus.NOT_FOUND, "not_found", f"run {rid} has no lineage {lineage}")
        return result, result.membership[lineage]

    def trends(self, rid: str, lineage: int, genes: list[str]) -> models.Trends:
        if not genes:
            raise ApiError(HTTPStatus.UNPROCESSABLE_ENTITY, "invalid_request", "ask for at least one gene")
        if len(genes) > 50:
            raise ApiError(HTTPStatus.UNPROCESSABLE_ENTITY, "too_many_genes", "ask for at most 50 genes")
        result, members = self._lineage(rid, lineage)
        idx = self.view.gene_indices(genes)
        t = gene_trends(self.view.expression_matrix(), idx, result.pseudotime, members)
        return models.Trends(
            run_id=rid,
            lineage=lineage,
            x=[float(v) for v in t.centers],
            counts=[int(c) for c in t.counts],
            genes=[
                models.GeneTrend(
                    gene=str(self.view.ds.gene_symbols[j]),
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
                self._heatmaps.move_to_end(key)
        if cached is not None:
            return cached
        idx = self.view.gene_indices(req.genes) if req.genes else None
        h = temporal_heatmap(
            self.view.expression_matrix(), result.pseudotime, members, top_n=req.top_n, genes=idx
        )
        body = binary.encode(
            {"heatmap": ("f32", h.values)},
            {
                "genes": [str(self.view.ds.gene_symbols[j]) for j in h.genes],
                "r2": [round(float(v), 4) for v in h.r2],
                "bins": [round(float(c), 4) for c in h.centers],
            },
        )
        with self._lock:
            self._heatmaps[key] = body
            while len(self._heatmaps) > MAX_HEATMAPS_IN_MEMORY:
                self._heatmaps.popitem(last=False)
        return body

    def _comparison(self, a: str, b: str) -> RunComparison:
        self._require_run(a)
        self._require_run(b)
        return compare(self._load(a), self._load(b))

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
