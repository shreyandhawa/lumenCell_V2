"""Route declarations for /api/v1.

Every route of the M1 contract is declared here with its request and response
models, so the OpenAPI document (and the generated TypeScript types) are
complete from the start. Routes that are not built yet raise a structured 501;
they are filled in by the platform track without changing their signatures.
"""

from http import HTTPStatus
from typing import Annotated, Any

from fastapi import APIRouter, Query
from fastapi.responses import Response

from lumencell import __version__
from lumencell.api.binary import MEDIA_TYPE
from lumencell.api.errors import not_implemented
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

router = APIRouter(prefix="/api/v1")

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
def get_dataset() -> DatasetInfo:
    raise not_implemented("GET /dataset")


@router.get(
    "/embeddings/{name}",
    tags=["dataset"],
    response_class=Response,
    responses=_lcb1("Array `coords` f32[n_obs, 2]."),
)
def get_embedding(name: str) -> Response:
    raise not_implemented("GET /embeddings/{name}")


@router.get(
    "/obs/{column}",
    tags=["dataset"],
    response_class=Response,
    responses=_lcb1(
        "Categorical/boolean columns: `codes` (u8/u16/i32 by category count, -1 = missing); "
        "continuous columns: `values` f32[n_obs] (NaN = missing)."
    ),
)
def get_obs(column: str) -> Response:
    raise not_implemented("GET /obs/{column}")


@router.get("/genes", tags=["dataset"], responses=_ERRORS)
def search_genes(
    q: Annotated[str, Query(max_length=100)] = "",
    limit: Annotated[int, Query(ge=1, le=500)] = 50,
) -> GeneSearchResult:
    raise not_implemented("GET /genes")


@router.get(
    "/expression",
    tags=["dataset"],
    response_class=Response,
    responses=_lcb1("Array `expression` f32[n_genes, n_obs], log-normalised, rows in request order."),
)
def get_expression(genes: Annotated[str, Query(description="Comma-separated gene symbols.")]) -> Response:
    raise not_implemented("GET /expression")


# ---- prep -------------------------------------------------------------------------


@router.get("/prepare", tags=["prepare"], responses=_ERRORS)
def get_prepare() -> PrepStatus:
    raise not_implemented("GET /prepare")


@router.post("/prepare", tags=["prepare"], status_code=HTTPStatus.ACCEPTED, responses=_ERRORS)
def start_prepare() -> JobInfo:
    raise not_implemented("POST /prepare")


# ---- methods ----------------------------------------------------------------------


@router.get("/methods", tags=["methods"], responses=_ERRORS)
def list_methods() -> MethodList:
    raise not_implemented("GET /methods")


# ---- root -------------------------------------------------------------------------


@router.get("/root", tags=["root"], responses=_ERRORS)
def get_root() -> RootState:
    raise not_implemented("GET /root")


@router.put("/root", tags=["root"], responses=_ERRORS)
def set_root(body: RootRequest) -> ResolvedRoot:
    raise not_implemented("PUT /root")


@router.post("/root/resolve", tags=["root"], responses=_ERRORS)
def resolve_root(body: RootRequest) -> ResolvedRoot:
    raise not_implemented("POST /root/resolve")


@router.get("/root/suggestions", tags=["root"], responses=_ERRORS)
def root_suggestions() -> RootSuggestions:
    raise not_implemented("GET /root/suggestions")


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
def create_trajectory(body: TrajectoryRequest) -> RunCreated:
    raise not_implemented("POST /trajectories")


@router.get("/trajectories", tags=["trajectories"], responses=_ERRORS)
def list_trajectories() -> RunList:
    raise not_implemented("GET /trajectories")


@router.get("/trajectories/{run_id}", tags=["trajectories"], responses=_ERRORS)
def get_trajectory(run_id: str) -> TrajectoryRun:
    raise not_implemented("GET /trajectories/{run_id}")


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
def get_trajectory_cells(run_id: str) -> Response:
    raise not_implemented("GET /trajectories/{run_id}/cells")


@router.get("/trajectories/{run_id}/curves", tags=["trajectories"], responses=_ERRORS)
def get_trajectory_curves(run_id: str, embedding: str) -> Curves:
    raise not_implemented("GET /trajectories/{run_id}/curves")


@router.patch("/trajectories/{run_id}/lineages/{lineage_id}", tags=["trajectories"], responses=_ERRORS)
def rename_lineage(run_id: str, lineage_id: int, body: LineageRename) -> Lineage:
    raise not_implemented("PATCH /trajectories/{run_id}/lineages/{lineage_id}")


@router.get("/trajectories/{run_id}/trends", tags=["trajectories"], responses=_ERRORS)
def get_trends(
    run_id: str,
    lineage: int,
    genes: Annotated[str, Query(description="Comma-separated gene symbols.")],
) -> Trends:
    raise not_implemented("GET /trajectories/{run_id}/trends")


@router.post(
    "/trajectories/{run_id}/heatmap",
    tags=["trajectories"],
    response_class=Response,
    responses=_lcb1("Array `heatmap` f32[n_genes, n_bins]; meta has `genes` (peak order) and `bins`."),
)
def get_heatmap(run_id: str, body: HeatmapRequest) -> Response:
    raise not_implemented("POST /trajectories/{run_id}/heatmap")


@router.get("/compare", tags=["trajectories"], responses=_ERRORS)
def compare_runs(a: str, b: str) -> Comparison:
    raise not_implemented("GET /compare")


# ---- jobs -------------------------------------------------------------------------


@router.get("/jobs/{job_id}", tags=["jobs"], responses=_ERRORS)
def get_job(job_id: str) -> JobInfo:
    raise not_implemented("GET /jobs/{job_id}")


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
def job_events(job_id: str) -> Response:
    raise not_implemented("GET /jobs/{job_id}/events")


@router.delete("/jobs/{job_id}", tags=["jobs"], responses=_ERRORS)
def cancel_job(job_id: str) -> JobInfo:
    raise not_implemented("DELETE /jobs/{job_id}")
