"""The HTTP contract: every JSON body the server sends or accepts.

These models are the single source of truth for the API. They generate
``contracts/openapi.json`` (``make contract``), which in turn generates the
TypeScript types in ``web/src/api/generated``; ``make ci`` fails when the
committed files drift from these models.

Per-cell data never travels as JSON: those routes return LCB1 (see
``api/binary.py``), and the JSON models below only describe the shape.
JSON field names are snake_case on both sides.
"""

from datetime import datetime
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field

SCHEMA_VERSION = 1


class Model(BaseModel):
    # Fields with defaults are always present in responses, so response schemas
    # mark them required; the generated TypeScript types are then not optional.
    model_config = ConfigDict(extra="forbid", json_schema_serialization_defaults_required=True)


# ---- shared ---------------------------------------------------------------------


class ErrorInfo(Model):
    code: str = Field(description="Stable machine-readable code, e.g. root_required.")
    message: str = Field(description="Human-readable explanation, safe to show in the UI.")
    detail: dict[str, Any] | None = None


class ErrorBody(Model):
    error: ErrorInfo


class Notice(Model):
    """A warning or diagnostic attached to a result. Shown in the trust panel."""

    code: str
    severity: Literal["info", "warning", "error"]
    message: str


class Health(Model):
    status: Literal["ok"]
    version: str


# ---- dataset ----------------------------------------------------------------------

ObsKind = Literal["categorical", "continuous", "boolean", "string"]
XKind = Literal["counts", "lognorm", "scaled"]


class ObsColumn(Model):
    name: str
    kind: ObsKind
    categories: list[str] | None = None
    colors: list[str] | None = Field(None, description="One hex colour per category, when the file has them.")
    min: float | None = None
    max: float | None = None


class Embedding(Model):
    name: str
    dims: int
    source: Literal["file", "computed", "diffusion"] = Field(
        description=(
            "file: read from obsm; computed: made by lumenCell's prep; "
            "diffusion: diffusion components (DC1 vs DC2, or DC1-DC3 when dims is 3)."
        )
    )


StageState = Literal["pending", "running", "done", "failed", "skipped"]


class PrepStage(Model):
    name: str
    state: StageState
    seconds: float | None = None
    message: str | None = None


class PrepStatus(Model):
    state: Literal["missing", "running", "ready", "failed"]
    stages: list[PrepStage]
    job_id: str | None = None


# ---- root -------------------------------------------------------------------------


class CellRoot(Model):
    kind: Literal["cell"]
    obs_name: str


class ClusterRoot(Model):
    kind: Literal["cluster"]
    key: str = Field(description="obs column holding the clusters.")
    value: str = Field(description="The cluster label to start from.")


class FileRoot(Model):
    """Use the root stored in the .h5ad (uns['iroot'] or a boolean obs 'root' column)."""

    kind: Literal["file"]


RootSpec = Annotated[CellRoot | ClusterRoot | FileRoot, Field(discriminator="kind")]


class ResolvedRoot(Model):
    spec: RootSpec
    cell_index: int
    obs_name: str
    cluster_key: str | None = None
    cluster_label: str | None = None
    resolved_by: str = Field(description="The rule that picked the cell, in words.")
    warnings: list[Notice] = []


class RootRequest(Model):
    root: RootSpec


class RootState(Model):
    root: ResolvedRoot | None


class RootSuggestion(Model):
    """A candidate start point. Never applied automatically."""

    spec: ClusterRoot | FileRoot
    source: Literal["file", "time_column", "gene_count"]
    score: float
    evidence: str = Field(description="Why this is suggested, e.g. '84% of cells at day 0 (p=1e-40)'.")
    is_extremity: bool = Field(description="Whether the cluster sits at an end of the PAGA graph.")


class RootSuggestions(Model):
    suggestions: list[RootSuggestion]


class DatasetInfo(Model):
    schema_version: int
    name: str
    n_obs: int
    n_vars: int
    x_kind: XKind
    gene_symbol_column: str | None
    obs: list[ObsColumn]
    embeddings: list[Embedding]
    default_embedding: str | None
    cluster_key: str | None
    label_key: str | None
    prep: PrepStatus
    root: ResolvedRoot | None
    versions: dict[str, str]


class GeneSearchResult(Model):
    genes: list[str]
    total: int
    truncated: bool


# ---- methods ----------------------------------------------------------------------


class MethodInfo(Model):
    name: str
    label: str
    family: str
    description: str
    available: bool
    install_hint: str | None = None
    experimental: bool = False
    experimental_note: str | None = Field(
        None, description="Why the method is experimental, from the benchmark."
    )
    params_schema: dict[str, Any] = Field(description="JSON Schema of the method's parameters.")
    reference_scores: dict[str, float] = Field(
        default_factory=dict, description="Benchmark scores of the current defaults, by metric."
    )


class MethodList(Model):
    methods: list[MethodInfo]


# ---- trajectories -------------------------------------------------------------------

RunState = Literal["queued", "running", "done", "failed", "cancelled"]


class TrajectoryRequest(Model):
    method: str
    root: RootSpec | None = Field(None, description="Omit to use the saved root; 409 if none is saved.")
    terminals: list[str] | None = Field(None, description="Terminal cluster labels, when known.")
    params: dict[str, Any] = {}


class RunCreated(Model):
    run_id: str
    job_id: str | None
    state: RunState


class Milestone(Model):
    id: int
    kind: Literal["root", "fork", "tip", "cluster"]
    label: str
    pseudotime: float


class Segment(Model):
    id: int
    source: int
    target: int
    length: float


class Lineage(Model):
    id: int
    name: str
    auto_name: str
    terminal_milestone: int
    terminal_label: str | None
    terminal_purity: float | None
    n_cells: int
    segments: list[int]


class Divergence(Model):
    """The pseudotime at which lineages ``a`` and ``b`` stop sharing cells."""

    a: int
    b: int
    pseudotime: float


class Provenance(Model):
    method: str
    method_version: str
    lumencell_version: str
    params: dict[str, Any]
    param_hash: str
    dataset_fingerprint: str
    inference_rep: str
    inference_dims: int
    root: ResolvedRoot
    terminals: list[str]
    terminals_source: Literal["user", "auto"]
    created_at: datetime


class RunSummary(Model):
    run_id: str
    method: str
    state: RunState
    created_at: datetime
    n_lineages: int | None = None
    error: str | None = None


class RunList(Model):
    runs: list[RunSummary]


class TrajectoryRun(Model):
    run_id: str
    state: RunState
    topology_kind: Literal["tree", "fate_star"]
    milestones: list[Milestone]
    segments: list[Segment]
    lineages: list[Lineage]
    divergence: list[Divergence]
    warnings: list[Notice]
    n_unassigned: int
    runtime_seconds: float
    cell_arrays: list[str] = Field(description="Names of the arrays served by GET .../cells.")
    provenance: Provenance


class LineageCurve(Model):
    lineage: int
    points: list[tuple[float, float] | tuple[float, float, float]] = Field(
        description="Curve points in the embedding's own dimensions (2 or 3 coordinates each)."
    )
    fidelity: float | None = Field(
        None, description="Spearman of pseudotime against position along the drawn curve."
    )


class MilestonePosition(Model):
    id: int
    x: float
    y: float
    z: float | None = Field(None, description="Set on 3-D embeddings only.")


class Curves(Model):
    run_id: str
    embedding: str
    curve_source: Literal["fitted", "binned"]
    curves: list[LineageCurve]
    milestones: list[MilestonePosition]
    warnings: list[Notice]


class LineageRename(Model):
    name: str = Field(min_length=1, max_length=80)


class GeneTrend(Model):
    gene: str
    y: list[float | None] = Field(description="Smoothed mean per bin; null where the bin has no cells.")
    r2: float = Field(description="Share of the gene's variance on this lineage that the trend explains.")


class Trends(Model):
    run_id: str
    lineage: int
    x: list[float] = Field(description="Pseudotime bin centres.")
    counts: list[int] = Field(description="Cells per bin.")
    genes: list[GeneTrend]


class HeatmapRequest(Model):
    lineage: int
    genes: list[str] | None = None
    top_n: int = Field(50, ge=1, le=500)


class LineageMatch(Model):
    a: int
    b: int
    jaccard: float


class Comparison(Model):
    runs: tuple[str, str]
    n_common: int
    spearman: float | None
    lineage_matches: list[LineageMatch]


# ---- jobs -------------------------------------------------------------------------

JobKind = Literal["prepare", "infer", "heatmap"]
JobState = Literal["queued", "running", "done", "failed", "cancelled"]


class JobInfo(Model):
    job_id: str
    kind: JobKind
    state: JobState
    stage: str | None = None
    fraction: float | None = None
    message: str | None = None
    error: str | None = None
    result_id: str | None = None


class JobEvent(Model):
    """One server-sent event on GET /jobs/{id}/events."""

    job_id: str
    state: JobState
    stage: str | None
    i: int
    n: int
    fraction: float
    message: str | None
    elapsed: float
    eta: float | None
