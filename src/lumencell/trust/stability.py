"""The stability of a run, as the API and the trust file hold it.

``to_api_stability`` maps the engine's result to the API models: every undefined value (a Spearman
on fewer than three shared cells, a failed replicate, an empty kind) becomes ``None``, so the trust
file is strict JSON. The stability worker job stores the result under ``stability_name()``, a hash
of the metric version and the ``[stability]`` parameters from ``trust/constants.py``: a changed
parameter makes a new name, and a stored result for the old one is simply not read.
"""

import hashlib
import json
import math
from typing import Literal

from lumencell.api import models
from lumencell.engine.stability import KindSummary, ReplicateRow, StabilityResult
from lumencell.trust.constants import (
    MATCH_JACCARD,
    N_RESEEDS,
    N_SUBSAMPLES,
    STABILITY_METRIC_VERSION,
    STABILITY_SEED,
    SUBSAMPLE_FRACTION,
)

STABILITY_NOTE = (
    "Reseeds rerun Leiden, PAGA and the method with new seeds on all cells. Subsamples rerun the "
    "neighbour graph, Leiden, diffusion maps, PAGA and the method on 80% of cells, reusing "
    "normalisation, HVG and PCA, which slightly understates subsampling variance."
)


def stability_hash() -> str:
    """16 hex characters over the metric version and the stability parameters."""
    params = {
        "n_reseeds": N_RESEEDS,
        "n_subsamples": N_SUBSAMPLES,
        "subsample_fraction": SUBSAMPLE_FRACTION,
        "seed": STABILITY_SEED,
        "match_jaccard": MATCH_JACCARD,
    }
    text = json.dumps({"metric_version": STABILITY_METRIC_VERSION, "params": params}, sort_keys=True)
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def stability_name() -> str:
    return f"stability-{stability_hash()}"


def _num(x: float) -> float | None:
    """A finite value as a float; NaN and infinity become ``None``."""
    return float(x) if math.isfinite(x) else None


def _row(r: ReplicateRow) -> models.StabilityReplicate:
    return models.StabilityReplicate(
        kind=r.kind,
        index=r.index,
        seed=r.seed,
        method_seed=r.method_seed,
        n_cells=r.n_cells,
        n_compared=r.n_compared,
        n_lineages=r.n_lineages,
        spearman=_num(r.spearman),
        share_matched=_num(r.share_matched),
        mean_jaccard=_num(r.mean_jaccard),
        terminals=None if r.terminals is None else list(r.terminals),
        error=r.error,
    )


def _summary(s: KindSummary) -> models.StabilitySummary:
    return models.StabilitySummary(
        kind=s.kind,
        r=s.r,
        fraction=s.fraction,
        seeds=list(s.seeds),
        spearman_median=_num(s.spearman_median),
        spearman_min=_num(s.spearman_min),
        share_matched_median=_num(s.share_matched_median),
        share_matched_min=_num(s.share_matched_min),
        mean_jaccard_median=_num(s.mean_jaccard_median),
        mean_jaccard_min=_num(s.mean_jaccard_min),
        n_lineages_median=_num(s.n_lineages_median),
        n_lineages_min=s.n_lineages_min,
        replicates=[_row(r) for r in s.rows],
    )


def to_api_stability(r: StabilityResult) -> models.Stability:
    failed = [row for s in (r.reseeds, r.subsamples) for row in s.rows if row.error is not None]
    return models.Stability(
        status="ok",
        match_jaccard=r.match_jaccard,
        reference_lineages=r.reference_lineages,
        n_failed=r.n_failed,
        reseeds=_summary(r.reseeds),
        subsamples=_summary(r.subsamples),
        lineage_stability=_num(r.lineage_stability),
        pseudotime_stability=_num(r.pseudotime_stability),
        notes=[STABILITY_NOTE, *(f"{row.kind} {row.index + 1} failed: {row.error}" for row in failed)],
    )


def not_measured(status: Literal["pending", "failed", "not_available"], note: str) -> models.Stability:
    """Stability that was not computed (or not yet), with the reason."""
    return models.Stability(
        status=status,
        match_jaccard=MATCH_JACCARD,
        reference_lineages=None,
        n_failed=0,
        reseeds=None,
        subsamples=None,
        lineage_stability=None,
        pseudotime_stability=None,
        notes=[note],
    )
