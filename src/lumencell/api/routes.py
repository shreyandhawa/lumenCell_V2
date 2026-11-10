"""Route declarations for /api/v1.

Every route of the M1 contract is declared here with its request and response
models, so the OpenAPI document (and the generated TypeScript types) are
complete from the start. Routes that are not built yet raise a structured 501;
they are filled in by the platform track without changing their signatures.
"""

import json
import queue
from collections.abc import Iterator
from http import HTTPStatus
from typing import Annotated, Any

from fastapi import APIRouter, Depends, Query, Request
from fastapi.responses import Response, StreamingResponse

from lumencell import __version__
from lumencell.api.binary import MEDIA_TYPE
from lumencell.api.errors import ApiError
from lumencell.api.models import (
    Comparison,
    Curves,
    DatasetInfo,
    ErrorBody,
    GeneSearchResult,
    Health,
    HeatmapRequest,
    JobInfo,
    Lineage,
    LineageRename,
    MethodList,
    PrepStatus,
    ResolvedRoot,
    RootRequest,
    RootState,
    RootSuggestions,
    RunCreated,
    RunList,
    TrajectoryRequest,
    TrajectoryRun,
    Trends,
)
from lumencell.api.session import Session

router = APIRouter(prefix="/api/v1")


def get_session(request: Request) -> Session:
    session: Session | None = getattr(request.app.state, "session", None)
    if session is None:
        raise ApiError(HTTPStatus.SERVICE_UNAVAILABLE, "no_dataset", "The server has no dataset loaded.")
    return session


S = Annotated[Session, Depends(get_session)]


def lcb1(body: bytes) -> Response:
    return Response(content=body, media_type=MEDIA_TYPE)


def _csv(value: str) -> list[str]:
    return [v.strip() for v in value.split(",") if v.strip()]


_ERRORS: dict[int | str, dict[str, Any]] = {
    404: {"model": ErrorBody, "description": "Unknown dataset column, gene, run or job."},
    422: {"model": ErrorBody, "description": "The request is not valid."},
}


def _lcb1(description: str) -> dict[int | str, dict[str, Any]]:
    """OpenAPI response spec for a route that returns an LCB1 buffer."""
    return {
        200: {
            "description": f"LCB1 buffer. {description}",
            "content": {MEDIA_TYPE: {"schema": {"type": "string", "format": "binary"}}},
        },
        **_ERRORS,
    }


# ---- health -----------------------------------------------------------------------


@router.get("/healthz", tags=["meta"])
def healthz() -> Health:
    return Health(status="ok", version=__version__)


# ---- dataset ----------------------------------------------------------------------


@router.get("/dataset", tags=["dataset"], responses=_ERRORS)
def get_dataset(session: S) -> DatasetInfo:
    return session.dataset_info()


@router.get(
    "/embeddings/{name}",
    tags=["dataset"],
    response_class=Response,
    responses=_lcb1("Array `coords` f32[n_obs, 2]."),
)
def get_embedding(name: str, session: S) -> Response:
    from lumencell.api import binary

    return lcb1(binary.encode({"coords": ("f32", session.embedding(name))}, {"name": name}))


@router.get(
    "/obs/{column}",
    tags=["dataset"],
    response_class=Response,
    responses=_lcb1(
        "Categorical/boolean columns: `codes` i16[n_obs] (i32 above 32,766 categories), -1 = missing, "
        "with `kind` and `categories` in meta; continuous columns: `values` f32[n_obs] (NaN = missing)."
    ),
)
def get_obs(column: str, session: S) -> Response:
    return lcb1(session.obs_payload(column))


@router.get("/genes", tags=["dataset"], responses=_ERRORS)
def search_genes(
    session: S,
    q: Annotated[str, Query(max_length=100)] = "",
    limit: Annotated[int, Query(ge=1, le=500)] = 50,
) -> GeneSearchResult:
    return session.search_genes(q, limit)


@router.get(
    "/expression",
    tags=["dataset"],
    response_class=Response,
    responses=_lcb1(
        "Array `expression` f32[n_genes, n_obs], log-normalised, rows in request order; "
        "meta `genes` holds the resolved symbols."
    ),
)
def get_expression(
    session: S, genes: Annotated[str, Query(description="Comma-separated gene symbols.")]
) -> Response:
    return lcb1(session.expression_payload(_csv(genes)))


# ---- prep -------------------------------------------------------------------------


@router.get("/prepare", tags=["prepare"], responses=_ERRORS)
def get_prepare(session: S) -> PrepStatus:
    return session.prep_status()


@router.post("/prepare", tags=["prepare"], status_code=HTTPStatus.ACCEPTED, responses=_ERRORS)
def start_prepare(session: S) -> JobInfo:
    return session.job_info(session.start_prep())


# ---- methods ----------------------------------------------------------------------


@router.get("/methods", tags=["methods"], responses=_ERRORS)
def list_methods(session: S) -> MethodList:
    return session.methods()


# ---- root -------------------------------------------------------------------------


@router.get("/root", tags=["root"], responses=_ERRORS)
def get_root(session: S) -> RootState:
    return session.root_state()


@router.put("/root", tags=["root"], responses=_ERRORS)
def set_root(body: RootRequest, session: S) -> ResolvedRoot:
    return session.set_root(body.root)


@router.post("/root/resolve", tags=["root"], responses=_ERRORS)
def resolve_root(body: RootRequest, session: S) -> ResolvedRoot:
    return session.resolve_root(body.root)


@router.get("/root/suggestions", tags=["root"], responses=_ERRORS)
def root_suggestions(session: S) -> RootSuggestions:
    return session.root_suggestions()


# ---- trajectories -------------------------------------------------------------------


@router.post(
    "/trajectories",
    tags=["trajectories"],
    status_code=HTTPStatus.ACCEPTED,
    responses={
        200: {"model": RunCreated, "description": "An identical run already exists."},
        409: {"model": ErrorBody, "description": "root_required: no root given and none saved."},
        **_ERRORS,
    },
)
def create_trajectory(body: TrajectoryRequest, session: S, response: Response) -> RunCreated:
    created, is_new = session.create_run(body)
    if not is_new:
        response.status_code = HTTPStatus.OK
    return created


@router.get("/trajectories", tags=["trajectories"], responses=_ERRORS)
def list_trajectories(session: S) -> RunList:
    return session.list_runs()


@router.get("/trajectories/{run_id}", tags=["trajectories"], responses=_ERRORS)
def get_trajectory(run_id: str, session: S) -> TrajectoryRun:
    return session.get_run(run_id)


@router.get(
    "/trajectories/{run_id}/cells",
    tags=["trajectories"],
    response_class=Response,
    responses=_lcb1(
        "Arrays: `pseudotime` f32[n_obs] (NaN = unassigned), `membership` u8[n_lineages, n_obs], "
        "`primary_lineage` i16[n_obs] (-1 shared trunk, -2 unassigned), `segment` i32[n_obs], "
        "`progress` f32[n_obs], optional `fate` f32[n_lineages, n_obs] and `extras.*`."
    ),
)
def get_trajectory_cells(run_id: str, session: S) -> Response:
    return lcb1(session.run_cells(run_id))


@router.get("/trajectories/{run_id}/curves", tags=["trajectories"], responses=_ERRORS)
def get_trajectory_curves(run_id: str, embedding: str, session: S) -> Curves:
    return session.run_curves(run_id, embedding)


@router.patch("/trajectories/{run_id}/lineages/{lineage_id}", tags=["trajectories"], responses=_ERRORS)
def rename_lineage(run_id: str, lineage_id: int, body: LineageRename, session: S) -> Lineage:
    return session.rename_lineage(run_id, lineage_id, body.name)


@router.get("/trajectories/{run_id}/trends", tags=["trajectories"], responses=_ERRORS)
def get_trends(
    run_id: str,
    lineage: int,
    genes: Annotated[str, Query(description="Comma-separated gene symbols.")],
    session: S,
) -> Trends:
    return session.trends(run_id, lineage, _csv(genes))


@router.post(
    "/trajectories/{run_id}/heatmap",
    tags=["trajectories"],
    response_class=Response,
    responses=_lcb1(
        "Array `heatmap` f32[n_genes, n_bins], each row scaled to [0, 1]; meta has `genes` (ordered by "
        "peak), `r2` and `bins` (bin centres)."
    ),
)
def get_heatmap(run_id: str, body: HeatmapRequest, session: S) -> Response:
    return lcb1(session.heatmap_payload(run_id, body))


@router.get("/compare", tags=["trajectories"], responses=_ERRORS)
def compare_runs(a: str, b: str, session: S) -> Comparison:
    return session.compare(a, b)


@router.get(
    "/compare/cells",
    tags=["trajectories"],
    response_class=Response,
    responses=_lcb1(
        "Array `rank_difference` f32[n_obs]: |rank in a - rank in b| / cells in both; NaN if not."
    ),
)
def compare_cells(a: str, b: str, session: S) -> Response:
    return lcb1(session.compare_cells(a, b))


# ---- jobs -------------------------------------------------------------------------


@router.get("/jobs/{job_id}", tags=["jobs"], responses=_ERRORS)
def get_job(job_id: str, session: S) -> JobInfo:
    return session.job_info(session.require_job(job_id))


@router.get(
    "/jobs/{job_id}/events",
    tags=["jobs"],
    response_class=Response,
    responses={
        200: {
            "description": "Server-sent events; each `data:` line is a JobEvent.",
            "content": {"text/event-stream": {"schema": {"type": "string"}}},
        },
        **_ERRORS,
    },
)
def job_events(job_id: str, session: S) -> Response:
    session.require_job(job_id)
    q = session.jobs.subscribe(job_id)

    def stream() -> Iterator[str]:
        try:
            while True:
                try:
                    event = q.get(timeout=15)
                except queue.Empty:
                    yield ": keep-alive\n\n"
                    continue
                yield f"data: {json.dumps(event)}\n\n"
                if event["state"] in ("done", "failed", "cancelled"):
                    return
        finally:
            session.jobs.unsubscribe(job_id, q)

    return StreamingResponse(
        stream(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-store", "X-Accel-Buffering": "no"},
    )


@router.delete("/jobs/{job_id}", tags=["jobs"], responses=_ERRORS)
def cancel_job(job_id: str, session: S) -> JobInfo:
    session.require_job(job_id)
    job = session.jobs.cancel(job_id)
    assert job is not None
    return session.job_info(job)
