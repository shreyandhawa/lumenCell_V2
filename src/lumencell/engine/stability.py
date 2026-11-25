"""Stability replicates of a finished run, and how far each agrees with it.

A run is re-run R times (``[stability]`` in bench/policy/thresholds.toml):

* **reseeds** rerun Leiden, PAGA and the method on all cells with new seeds, so the PAGA tree
  is never stale; every other prep artifact is reused;
* **subsamples** keep a uniform share of the cells (the root cell always among them), reuse the
  full-data normalisation, HVGs and PCA, and rerun kNN, Leiden, the diffusion maps and PAGA on a
  ``Dataset`` restricted to them.

Every replicate is rooted at ``cell:<obs name>`` of the reference run's root cell, so a
``cluster:`` root is never resolved again on different data. It gets the reference run's
clonal columns and goes through ``run_method``, which keeps scoring-only columns out of
inference. Replicate prep lives in memory only and never enters the stage cache.

Each replicate is compared with the reference on the cells both contain (subsample results
padded to all cells). A replicate that fails is recorded with its error and counts as no
agreement: it is never dropped from a median.
"""

import hashlib
import math
import statistics
import sys
import traceback
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any, Literal

import numpy as np
import numpy.typing as npt
from scipy import sparse

from lumencell.data.loader import Dataset
from lumencell.engine.compare import compare
from lumencell.engine.context import LEIDEN_KEY, InferenceContext
from lumencell.engine.methods.base import ProgressFn, no_progress
from lumencell.engine.registry import get_method, run_method
from lumencell.engine.roots import CellRoot
from lumencell.engine.schema import TrajectoryResult
from lumencell.engine.store import RunSpec
from lumencell.prep.pipeline import PrepResult
from lumencell.prep.stages import STAGE_BY_NAME, PrepParams, StageRecord, run_stage
from lumencell.trust.constants import (
    MATCH_JACCARD,
    N_RESEEDS,
    N_SUBSAMPLES,
    STABILITY_SEED,
    SUBSAMPLE_FRACTION,
)

BoolArray = npt.NDArray[np.bool_]
IntArray = npt.NDArray[np.integer[Any]]
ReplicateKind = Literal["reseed", "subsample"]

RESEED_STAGES = ("leiden", "paga")
SUBSAMPLE_STAGES = ("neighbors", "leiden", "diffmap", "paga")
SUBSAMPLE_REUSED = ("hvg", "pca_variance_ratio")


class StabilityError(ValueError):
    """Stability cannot be computed for this run; the message says why."""


@dataclass(frozen=True)
class StabilityConfig:
    n_reseeds: int = N_RESEEDS
    n_subsamples: int = N_SUBSAMPLES
    fraction: float = SUBSAMPLE_FRACTION
    seed: int = STABILITY_SEED
    match_jaccard: float = MATCH_JACCARD


@dataclass(frozen=True)
class ReplicatePlan:
    kind: ReplicateKind
    index: int
    seed: int  # prep seed (reseeds) or the subsample's draw seed (subsamples)
    method_seed: int | None  # None when the method has no seed parameter


def replicate_plan(
    config: StabilityConfig, prep_seed: int, method_params: Mapping[str, Any]
) -> tuple[ReplicatePlan, ...]:
    """The reseeds, then the subsamples, following the ``[stability]`` seed scheme."""
    raw = method_params.get("seed")
    method_seed = None if raw is None else int(raw)
    reseeds = [
        ReplicatePlan("reseed", i, prep_seed + 1 + i, None if method_seed is None else method_seed + 1 + i)
        for i in range(config.n_reseeds)
    ]
    subsamples = [
        ReplicatePlan("subsample", j, config.seed + 1 + j, method_seed) for j in range(config.n_subsamples)
    ]
    return (*reseeds, *subsamples)


def subsample_cells(n_obs: int, root: int, fraction: float, seed: int) -> npt.NDArray[np.int64]:
    """``round(fraction * n_obs)`` sorted cell indices, always including ``root``."""
    if not 0 <= root < n_obs:
        raise StabilityError(f"root cell {root} is outside the {n_obs} cells")
    k = max(1, min(n_obs, round(fraction * n_obs)))
    others = np.delete(np.arange(n_obs, dtype=np.int64), root)
    drawn = np.random.default_rng(seed).choice(others, size=k - 1, replace=False)
    return np.sort(np.append(drawn, np.int64(root))).astype(np.int64)


def restrict(ds: Dataset, cells: npt.NDArray[np.int64]) -> Dataset:
    """The dataset on ``cells`` only, without its matrix, under a fingerprint of its own."""
    cells = np.asarray(cells, dtype=np.int64)
    position = {int(c): i for i, c in enumerate(cells)}
    uns_root = None if ds.uns_root is None else position.get(int(ds.uns_root))
    digest = hashlib.sha256(f"{ds.fingerprint}:subset:".encode() + cells.tobytes()).hexdigest()
    return replace(
        ds,
        obs=ds.obs.iloc[cells].copy(),
        matrix=None,
        embeddings={k: v[cells] for k, v in ds.embeddings.items()},
        uns_root=uns_root,
        fingerprint=f"subset-{digest}",
        shape=(int(cells.size), ds.n_vars),
        genes_detected=ds.genes_detected[cells],
    )


def _rerun(
    ds: Dataset, artifacts: dict[str, Any], params: PrepParams, names: Sequence[str]
) -> list[StageRecord]:
    records = []
    for name in names:
        out, seconds = run_stage(STAGE_BY_NAME[name], ds, artifacts, params)
        artifacts.update(out)
        records.append(StageRecord(name, "done", seconds))
    return records


def reseed_prep(
    ds: Dataset, prep: PrepResult, params: PrepParams, seed: int
) -> tuple[PrepResult, PrepParams]:
    """Leiden and PAGA rerun with ``seed``; every other artifact reused as it is."""
    reseeded = replace(params, seed=seed)
    artifacts = dict(prep.artifacts)
    records = _rerun(ds, artifacts, reseeded, RESEED_STAGES)
    return PrepResult(artifacts=artifacts, stages=records), reseeded


def subsample_prep(
    ds: Dataset, prep: PrepResult, params: PrepParams, cells: npt.NDArray[np.int64]
) -> tuple[Dataset, PrepResult]:
    """kNN, Leiden, diffusion maps and PAGA on ``cells``, from the full data's HVGs and PCA."""
    sub = restrict(ds, cells)
    artifacts: dict[str, Any] = {k: prep.artifacts[k] for k in SUBSAMPLE_REUSED}
    artifacts["pca"] = np.asarray(prep.artifacts["pca"])[cells]
    records = _rerun(sub, artifacts, params, SUBSAMPLE_STAGES)
    return sub, PrepResult(artifacts=artifacts, stages=records)


def map_leiden_terminals(
    terminals: Sequence[str], reference_leiden: IntArray, replicate_leiden: IntArray, cells: IntArray
) -> list[str]:
    """Each Leiden terminal as the replicate cluster holding most of its reference cells.

    ``replicate_leiden`` is indexed like ``cells`` (the replicate's cells in the full data). A
    tie goes to the lowest label; a terminal with no cells left is dropped; duplicates collapse.
    """
    ref = np.asarray(reference_leiden).astype(str)[np.asarray(cells)]
    rep = np.asarray(replicate_leiden)
    out: list[str] = []
    for t in terminals:
        held = rep[ref == str(t)]
        if held.size == 0:
            continue
        labels, counts = np.unique(held, return_counts=True)  # sorted, so argmax takes the lowest on ties
        label = str(labels[int(np.argmax(counts))])
        if label not in out:
            out.append(label)
    return out


def pad_result(result: TrajectoryResult, cells: IntArray, n_obs: int) -> TrajectoryResult:
    """``result`` on all ``n_obs`` cells: cells outside ``cells`` are unassigned."""
    cells = np.asarray(cells, dtype=np.int64)

    def per_cell(arr: npt.NDArray[Any], fill: Any) -> npt.NDArray[Any]:
        out = np.full((*arr.shape[:-1], n_obs), fill, dtype=arr.dtype)
        out[..., cells] = arr
        return out

    def rows(m: sparse.csr_matrix) -> sparse.csr_matrix:
        coo = sparse.coo_matrix(m)
        return sparse.csr_matrix((coo.data, (cells[coo.row], coo.col)), shape=(n_obs, m.shape[1]))

    padded = replace(
        result,
        pseudotime=per_cell(result.pseudotime, np.nan),
        membership=per_cell(result.membership, False),
        segment=per_cell(result.segment, -1),
        progress=per_cell(result.progress, np.nan),
        milestone_weights=rows(result.milestone_weights),
        fate=None if result.fate is None else per_cell(result.fate, np.nan),
        extras={k: per_cell(v, np.nan) for k, v in result.extras.items()},
        curve_nodes=None if result.curve_nodes is None else rows(result.curve_nodes),
        root=replace(result.root, cell_index=int(cells[result.root.cell_index])),
    )
    padded.validate()
    return padded


@dataclass(frozen=True)
class ReplicateScores:
    n_compared: int
    n_lineages: int
    spearman: float
    share_matched: float
    mean_jaccard: float


def score_replicate(
    reference: TrajectoryResult, replicate: TrajectoryResult, shared: BoolArray, match_jaccard: float
) -> ReplicateScores:
    """Agreement with the reference on ``shared`` cells; a value that cannot be computed is NaN."""
    c = compare(reference, replicate, shared)
    n_ref = len(reference.lineages)
    matched = sum(p.jaccard >= match_jaccard for p in c.pairs)
    return ReplicateScores(
        n_compared=c.n_common,
        n_lineages=len(replicate.lineages),
        spearman=math.nan if c.spearman is None else c.spearman,
        share_matched=matched / n_ref if n_ref else math.nan,
        mean_jaccard=float(np.mean([p.jaccard for p in c.pairs])) if c.pairs else math.nan,
    )


@dataclass(frozen=True)
class ReplicateRow:
    kind: ReplicateKind
    index: int
    seed: int
    method_seed: int | None
    n_cells: int
    n_compared: int
    n_lineages: int
    spearman: float
    share_matched: float
    mean_jaccard: float
    terminals: tuple[str, ...] | None
    error: str | None = None


@dataclass(frozen=True)
class KindSummary:
    kind: ReplicateKind
    r: int
    fraction: float
    seeds: tuple[int, ...]
    spearman_median: float
    spearman_min: float
    share_matched_median: float
    share_matched_min: float
    mean_jaccard_median: float
    mean_jaccard_min: float
    n_lineages_median: float
    n_lineages_min: int
    rows: tuple[ReplicateRow, ...]


def _zeroed(values: Sequence[float]) -> list[float]:
    """NaN (undefined, or a failed replicate) counts as no agreement."""
    return [0.0 if math.isnan(v) else float(v) for v in values]


def _median(values: Sequence[float]) -> float:
    return float(statistics.median(_zeroed(values))) if values else math.nan


def _min(values: Sequence[float]) -> float:
    return min(_zeroed(values)) if values else math.nan


def summarise(kind: ReplicateKind, fraction: float, rows: Sequence[ReplicateRow]) -> KindSummary:
    """Median and minimum per value; NaN counts as 0, a failed row as ``n_lineages = 0``."""
    n_lin = [0 if r.error is not None else r.n_lineages for r in rows]
    return KindSummary(
        kind=kind,
        r=len(rows),
        fraction=fraction,
        seeds=tuple(r.seed for r in rows),
        spearman_median=_median([r.spearman for r in rows]),
        spearman_min=_min([r.spearman for r in rows]),
        share_matched_median=_median([r.share_matched for r in rows]),
        share_matched_min=_min([r.share_matched for r in rows]),
        mean_jaccard_median=_median([r.mean_jaccard for r in rows]),
        mean_jaccard_min=_min([r.mean_jaccard for r in rows]),
        n_lineages_median=float(statistics.median(n_lin)) if n_lin else math.nan,
        n_lineages_min=min(n_lin, default=0),
        rows=tuple(rows),
    )


@dataclass(frozen=True)
class StabilityResult:
    n_obs: int
    reference_lineages: int
    match_jaccard: float
    reseeds: KindSummary
    subsamples: KindSummary
    n_failed: int  # replicates that raised; each counts as 0 in every median and minimum
    lineage_stability: float  # C3: median share_matched over all R replicates, NaN counted as 0
    pseudotime_stability: float  # C4: median spearman over all R replicates, NaN counted as 0


def _root_index(ds: Dataset, name: str) -> int:
    hits = np.flatnonzero(ds.obs.index.astype(str) == name)
    if hits.size > 1:
        raise StabilityError(
            f"the root cell's name {name!r} occurs more than once in the file, "
            "so replicates cannot be rooted at it"
        )
    if hits.size == 0:
        raise StabilityError(
            f"the root cell {name!r} is not in the file, so replicates cannot be rooted at it"
        )
    return int(hits[0])


def _method_params(spec: RunSpec) -> dict[str, Any]:
    # A bad param raises here and aborts the whole job: the reference run would have failed too.
    """The run's params with the method's defaults filled in, so a default seed is reseeded too."""
    resolved = get_method(spec.method).Params.model_validate(spec.params).model_dump()
    return {**resolved, **spec.params}


def _terminal_key_default(method: str) -> str | None:
    field = get_method(method).Params.model_fields.get("cluster_key")
    return None if field is None else field.default


def _replicate(
    ds: Dataset,
    prep: PrepResult,
    params: PrepParams,
    reference: TrajectoryResult,
    spec: RunSpec,
    plan: ReplicatePlan,
    root: int,
    config: StabilityConfig,
) -> ReplicateRow:
    n = ds.n_obs
    cells = np.arange(n, dtype=np.int64)
    if plan.kind == "subsample":
        cells = subsample_cells(n, root, config.fraction, plan.seed)
    terminals: list[str] | None = spec.terminals
    try:
        if plan.kind == "reseed":
            rep_ds, (rep_prep, rep_params) = ds, reseed_prep(ds, prep, params, plan.seed)
        else:
            (rep_ds, rep_prep), rep_params = subsample_prep(ds, prep, params, cells), params
        ctx = InferenceContext.from_prep(rep_ds, rep_prep, rep_params, clonal_columns=spec.clonal_columns)
        if spec.terminals and (
            spec.params.get("cluster_key", _terminal_key_default(spec.method)) == LEIDEN_KEY
        ):
            terminals = map_leiden_terminals(
                spec.terminals, np.asarray(prep.artifacts["leiden"]), ctx.leiden, cells
            )
            if not terminals:
                raise StabilityError("none of the run's Leiden terminals has cells left in this replicate")
        method_params = spec.params if plan.method_seed is None else {**spec.params, "seed": plan.method_seed}
        result = run_method(
            ctx,
            spec.method,
            CellRoot(reference.root.obs_name),
            terminals,
            method_params,
            label_key=reference.label_key,  # the column the run was named in, not one re-detected
        )
        if plan.kind == "subsample":
            result = pad_result(result, cells, n)
        shared = np.zeros(n, dtype=bool)
        shared[cells] = True
        s = score_replicate(reference, result, shared, config.match_jaccard)
    except Exception as exc:
        print(f"stability {plan.kind} {plan.index + 1} failed:", file=sys.stderr)
        traceback.print_exception(exc, file=sys.stderr)
        return ReplicateRow(
            plan.kind,
            plan.index,
            plan.seed,
            plan.method_seed,
            int(cells.size),
            0,
            0,
            math.nan,
            math.nan,
            math.nan,
            None if terminals is None else tuple(terminals),
            error=f"{type(exc).__name__}: {exc}",
        )
    return ReplicateRow(
        plan.kind,
        plan.index,
        plan.seed,
        plan.method_seed,
        int(cells.size),
        s.n_compared,
        s.n_lineages,
        s.spearman,
        s.share_matched,
        s.mean_jaccard,
        None if terminals is None else tuple(terminals),
    )


def run_stability(
    ds: Dataset,
    prep: PrepResult,
    params: PrepParams,
    reference: TrajectoryResult,
    spec: RunSpec,
    config: StabilityConfig = StabilityConfig(),  # noqa: B008 (frozen, so safe to share)
    progress: ProgressFn = no_progress,
) -> StabilityResult:
    """Run every replicate in memory and summarise its agreement with ``reference``."""
    if reference.n_obs != ds.n_obs:
        raise StabilityError("the reference run is on a different dataset")
    root = _root_index(ds, reference.root.obs_name)
    plans = replicate_plan(config, params.seed, _method_params(spec))
    rows: list[ReplicateRow] = []
    for i, plan in enumerate(plans):
        progress(f"{plan.kind} {plan.index + 1}", i, len(plans))
        rows.append(_replicate(ds, prep, params, reference, spec, plan, root, config))
    reseeds = [r for r in rows if r.kind == "reseed"]
    subsamples = [r for r in rows if r.kind == "subsample"]
    return StabilityResult(
        n_obs=ds.n_obs,
        reference_lineages=len(reference.lineages),
        match_jaccard=config.match_jaccard,
        n_failed=sum(r.error is not None for r in rows),
        reseeds=summarise("reseed", 1.0, reseeds),
        subsamples=summarise("subsample", config.fraction, subsamples),
        lineage_stability=_median([r.share_matched for r in rows]),
        pseudotime_stability=_median([r.spearman for r in rows]),
    )
