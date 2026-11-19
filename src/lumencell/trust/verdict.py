"""A run's tree-support verdict, from what the API serves about it.

The rules are in ``engine/tree_support.py``; this module only turns the served stability and clonal
agreement into its evidence, and the verdict back into the API model. ``verdict`` is never told
whether the rule is verified: it reads the recorded verification itself, so a failed one caps
a served supported verdict at provisional (owner ruling H-1: weak and not_supported are kept).
"""

import math
from typing import Any

from lumencell.api import models
from lumencell.engine.tree_support import (
    ClonalEvidence,
    StabilityEvidence,
    Thresholds,
    Verdict,
    calibrated_thresholds,
    verdict,
)


def _nan(x: float | None) -> float:
    return math.nan if x is None else x


def _num(x: float) -> float | None:
    return None if math.isnan(x) else x


def stability_evidence(s: models.Stability) -> StabilityEvidence:
    if s.status == "pending":
        return StabilityEvidence("pending")
    if s.status != "ok":
        return StabilityEvidence("unavailable", note=s.notes[0] if s.notes else "")
    n = sum(k.r for k in (s.reseeds, s.subsamples) if k is not None)
    return StabilityEvidence(
        "done",
        lineage=_nan(s.lineage_stability),
        pseudotime=_nan(s.pseudotime_stability),
        n_replicates=n,
        n_failed=s.n_failed,
    )


def clonal_evidence(c: models.ClonalAgreement | None, n_lineages: int) -> ClonalEvidence:
    """``n_lineages`` is the run's own count: a failed score has no primary to read it from."""
    if c is None or c.status == "not_configured":
        return ClonalEvidence("not_configured", n_lineages)
    if c.status == "pending":
        return ClonalEvidence("pending", n_lineages)
    p = c.primary
    return ClonalEvidence(
        "scored",
        n_lineages,
        status=c.status,
        p=_nan(p.p if p else None),
        delta_eta=_nan(p.delta_eta if p else None),
        note=c.notes[0] if c.status != "ok" and c.notes else "",
    )


def to_api_tree_support(v: Verdict) -> models.TreeSupport:
    return models.TreeSupport(
        level=v.level,
        checks=[
            models.SupportCheck(
                id=c.id,
                name=c.name,
                state=c.state,
                value=_num(c.value),
                threshold=c.threshold,
                reason=c.reason,
            )
            for c in v.checks
        ],
        reasons=list(v.reasons),
    )


def tree_support_for_run(
    graph: dict[str, Any],
    tree_likeness: float | None,
    stability: models.Stability,
    clonal: models.ClonalAgreement | None,
    thresholds: Thresholds | None = None,
) -> models.TreeSupport:
    """``tree_likeness`` is None while prep is not loaded (C2 pending) and NaN without edges (C2 fails)."""
    coverage = 1 - graph["n_unassigned"] / graph["n_obs"]
    v = verdict(
        coverage,
        tree_likeness,
        stability_evidence(stability),
        clonal_evidence(clonal, len(graph["lineages"])),
        thresholds if thresholds is not None else calibrated_thresholds(),
    )
    return to_api_tree_support(v)
