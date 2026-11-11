"""Server-side state for one loaded dataset.

The server process keeps the dataset (for annotations, embeddings and expression)
and reads prep artifacts and runs from the on-disk cache. Heavy work (prep,
inference) goes to the job worker. The session ties together the parts that do
the work: the dataset view, the runs, and the user's saved choices.
"""

from http import HTTPStatus
from pathlib import Path
from typing import Any

from lumencell import __version__
from lumencell.api import models
from lumencell.api.convert import ApiRootSpec, from_api_spec, to_api_root, to_api_spec
from lumencell.api.dataset_view import DatasetView
from lumencell.api.errors import ApiError
from lumencell.api.runs import Runs
from lumencell.api.state import StateStore, default_state_dir
from lumencell.data.loader import Dataset, load
from lumencell.engine.context import LEIDEN_KEY
from lumencell.engine.registry import METHODS
from lumencell.engine.roots import RootError, RootSpec, parse_root, resolve_root, suggest_roots
from lumencell.engine.schema import ResolvedRoot
from lumencell.engine.store import RunStore
from lumencell.jobs.manager import Job, JobManager
from lumencell.prep.cache import PrepCache, dataset_dir
from lumencell.prep.stages import STAGES, PrepParams


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
        # Without its matrix: expression is served from the prep cache, held once.
        self.ds: Dataset = load(self.path, x_is, keep_matrix=False)  # type: ignore[arg-type]
        self.cache = PrepCache(self.ds, cache_root)
        self.state = StateStore((state_dir or default_state_dir()) / f"{self.ds.fingerprint}.json")
        self.view = DatasetView(self.path, self.ds, self.cache, self.params, on_cache_miss=self.start_prep)
        self.runs = Runs(
            self.view,
            RunStore(dataset_dir(self.ds.fingerprint, cache_root)),
            jobs,
            self.state,
            self._job_payload,
            self._resolve,
        )

    # ---- prep --------------------------------------------------------------------------

    @property
    def prep_key(self) -> str:
        return self.view.prep_key

    def _job_payload(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "x_is": self.x_is,
            "prep_params": self.params,
            "cache_root": str(self.cache_root) if self.cache_root else None,
        }

    def prep_ready(self) -> bool:
        return self.view.prep_ready()

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

    # ---- dataset -------------------------------------------------------------------------

    def dataset_info(self) -> models.DatasetInfo:
        return models.DatasetInfo(
            schema_version=models.SCHEMA_VERSION,
            name=self.ds.name,
            n_obs=self.ds.n_obs,
            n_vars=self.ds.n_vars,
            x_kind=self.ds.x_kind,
            gene_symbol_column=self.ds.gene_symbol_column,
            obs=self.view.obs_columns(),
            embeddings=self.view.embeddings(),
            default_embedding=self.view.default_embedding(),
            cluster_key=LEIDEN_KEY if self.view.loaded_prep() is not None else None,
            label_key=self.view.label_key,
            prep=self.prep_status(),
            root=self.root_state().root,
            versions={"lumencell": __version__},
        )

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
            return resolve_root(self.view.ctx(), spec)
        except RootError as exc:
            raise ApiError(HTTPStatus.UNPROCESSABLE_ENTITY, "invalid_root", str(exc)) from exc

    def root_state(self) -> models.RootState:
        spec_text = self.state.read().get("root")
        if not spec_text or self.view.loaded_prep() is None:
            return models.RootState(root=None)
        try:
            return models.RootState(root=to_api_root(resolve_root(self.view.ctx(), parse_root(spec_text))))
        except RootError:
            return models.RootState(root=None)

    def resolve_root(self, spec: ApiRootSpec) -> models.ResolvedRoot:
        return to_api_root(self._resolve(from_api_spec(spec)))

    def set_root(self, spec: ApiRootSpec) -> models.ResolvedRoot:
        resolved = self._resolve(from_api_spec(spec))

        def change(state: dict[str, Any]) -> None:
            state["root"] = resolved.spec

        self.state.update(change)
        return to_api_root(resolved)

    def root_suggestions(self) -> models.RootSuggestions:
        out = []
        for s in suggest_roots(self.view.ctx()):
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
