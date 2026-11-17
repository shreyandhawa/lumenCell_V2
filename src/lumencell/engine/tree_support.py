"""The tree-support verdict (spec 5.3): five checks, three levels, and provisional.

Pure rules on given values. Nothing here reads a run or a file. ``_level`` is the one place that turns
the checks into a level, so a cap on the level (a failed verification) has a single place to go.
No reason says which state came first; C5 speaks only of clones supporting the lineages.
"""

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt
from scipy import sparse

from lumencell.clonal.scores import ClonalResult
from lumencell.trust.constants import (
    CLONAL_P_MAX,
    COVERAGE_MIN,
    TAU_LIN,
    TAU_PT,
    TAU_TREE,
    VERDICT_VERIFIED,
    VERDICT_VERIFIED_ON,
)

CheckId = Literal["C1", "C2", "C3", "C4", "C5"]
CheckState = Literal["pass", "fail", "not_applicable", "pending", "unavailable", "excluded"]
Level = Literal["supported", "weak", "not_supported", "provisional"]

CHECK_NAMES: dict[CheckId, str] = {
    "C1": "connectivity",
    "C2": "tree-likeness",
    "C3": "lineage stability",
    "C4": "pseudotime stability",
    "C5": "clonal support",
}


def tree_likeness(paga_connectivities: sparse.spmatrix, paga_tree: sparse.spmatrix) -> float:
    """Weight on PAGA's spanning tree over the total edge weight, each edge once. NaN with no edges."""
    conn = sparse.csr_matrix(paga_connectivities)
    sym = conn.maximum(conn.T)
    total = float(sparse.triu(sym, k=1).sum())
    if total == 0:
        return math.nan
    tree = sparse.csr_matrix(paga_tree)
    mask = sparse.triu(tree + tree.T, k=1).tocoo()
    dense = sym.tocsr()
    weight = float(sum(dense[i, j] for i, j in zip(mask.row, mask.col, strict=True)))
    return weight / total


def coverage(pseudotime: npt.NDArray[np.float64]) -> float:
    return float(np.isfinite(pseudotime).mean())


@dataclass(frozen=True)
class Thresholds:
    tau_tree: float | None
    tau_lin: float | None
    tau_pt: float | None
    coverage_min: float = COVERAGE_MIN
    clonal_p_max: float = CLONAL_P_MAX


def calibrated_thresholds() -> Thresholds:
    """The taus from the recorded calibration ([tree_support.calibrated]); None leaves a check out."""
    return Thresholds(tau_tree=TAU_TREE, tau_lin=TAU_LIN, tau_pt=TAU_PT)


@dataclass(frozen=True)
class StabilityEvidence:
    state: Literal["done", "pending", "unavailable"]
    lineage: float = math.nan
    pseudotime: float = math.nan
    note: str = ""
    # Replicates run and replicates that raised. Every one failing means stability was not measured.
    n_replicates: int = 0
    n_failed: int = 0


@dataclass(frozen=True)
class ClonalEvidence:
    state: Literal["not_configured", "pending", "scored"]
    # Required: a default would hide L < 2, where C5 does not apply. Unused until scored.
    n_lineages: int
    status: str = ""
    p: float = math.nan
    delta_eta: float = math.nan
    note: str = ""

    @classmethod
    def from_result(cls, r: ClonalResult) -> "ClonalEvidence":
        return cls("scored", r.n_lineages, status=r.status, p=r.p, delta_eta=r.delta_eta)


@dataclass(frozen=True)
class Check:
    id: CheckId
    name: str
    state: CheckState
    value: float
    threshold: float | None
    reason: str | None


@dataclass(frozen=True)
class Verdict:
    level: Level
    checks: tuple[Check, ...]
    reasons: tuple[str, ...]


def _check(
    cid: CheckId, state: CheckState, value: float, threshold: float | None, reason: str | None = None
) -> Check:
    return Check(cid, CHECK_NAMES[cid], state, value, threshold, reason)


def _c1(cov: float, t: Thresholds) -> Check:
    if cov >= t.coverage_min:
        return _check("C1", "pass", cov, t.coverage_min)
    return _check(
        "C1",
        "fail",
        cov,
        t.coverage_min,
        f"Only {cov:.1%} of cells have a pseudotime; at least {t.coverage_min:.0%} are needed.",
    )


def _c2(tree: float | None, t: Thresholds) -> Check:
    if t.tau_tree is None:
        return _check(
            "C2",
            "excluded",
            math.nan if tree is None else tree,
            None,
            "Tree-likeness was left out: it did not separate trees from the negative controls.",
        )
    if tree is None:
        return _check(
            "C2", "pending", math.nan, t.tau_tree, "Tree-likeness is waiting for the dataset's prep."
        )
    if math.isnan(tree):
        return _check(
            "C2",
            "fail",
            math.nan,
            t.tau_tree,
            "The cluster graph has no edges, so tree-likeness is undefined.",
        )
    if tree >= t.tau_tree:
        return _check("C2", "pass", tree, t.tau_tree)
    return _check(
        "C2",
        "fail",
        tree,
        t.tau_tree,
        f"The cluster graph is not tree-like: {tree:.2f} of its edge weight is on the tree, "
        f"below {t.tau_tree:.2f}.",
    )


def _stability(cid: CheckId, name: str, value: float, tau: float | None, s: StabilityEvidence) -> Check:
    if tau is None:
        return _check(
            cid,
            "excluded",
            math.nan,
            None,
            f"{name.capitalize()} was left out: it did not separate trees from the negative controls.",
        )
    if s.state == "pending":
        return _check(cid, "pending", math.nan, tau, f"{name.capitalize()} is waiting for the stability run.")
    all_failed = s.n_replicates > 0 and s.n_failed >= s.n_replicates
    if s.state == "unavailable" or all_failed:
        why = s.note or (
            "Every stability replicate failed." if all_failed else "Stability could not be measured."
        )
        return _check(cid, "unavailable", math.nan, tau, f"{name.capitalize()} could not be measured: {why}")
    if value >= tau:
        return _check(cid, "pass", value, tau)
    shown = "undefined" if math.isnan(value) else f"{value:.2f}"
    return _check(cid, "fail", value, tau, f"{name.capitalize()} is {shown}, below {tau:.2f}.")


def _c5(c: ClonalEvidence, t: Thresholds) -> Check:
    p_max = t.clonal_p_max
    if c.state == "not_configured":
        return _check("C5", "not_applicable", math.nan, p_max)
    if c.state == "pending":
        return _check("C5", "pending", math.nan, p_max, "Clonal support is still being scored.")
    if c.n_lineages < 2:
        return _check(
            "C5",
            "not_applicable",
            math.nan,
            p_max,
            f"Clonal support needs at least two lineages; there {'is' if c.n_lineages == 1 else 'are'} "
            f"{c.n_lineages}.",
        )
    if c.status == "ok":
        if c.p <= p_max and c.delta_eta > 0:
            return _check("C5", "pass", c.p, p_max)
        return _check(
            "C5",
            "fail",
            c.p,
            p_max,
            "Clones do not support the lineages beyond chance and beyond expression clustering "
            f"(p {c.p:.3g}, needs at most {p_max:g}, and a gain over the clustering baseline).",
        )
    extra = f" {c.note}" if c.note else ""
    return _check(
        "C5",
        "fail",
        c.p,
        p_max,
        f"Clonal support could not be shown: the clonal score has status {c.status}.{extra}",
    )


def _unverified_note() -> str:
    return (
        f"The verdict rule failed its verification on {VERDICT_VERIFIED_ON}: control runs were "
        "verdicted supported. Until it passes, no run is shown above provisional, unless it is not supported."
    )


def _level(checks: tuple[Check, ...], verified: bool) -> tuple[Level, str | None]:
    """The one place that decides the level, and the note when it is capped (a failed verification).

    A failed C1 or C2 stays not_supported: that is a real failure whether or not the rule is verified.
    """
    first_two = checks[:2]  # C1 and C2 are never "unavailable"
    if any(c.state == "fail" for c in first_two):
        return "not_supported", None
    if any(c.state == "pending" for c in first_two):
        return "provisional", None
    rest = checks[2:]
    if any(c.state in ("pending", "unavailable") for c in rest) and not any(c.state == "fail" for c in rest):
        return "provisional", None
    level: Level = "weak" if any(c.state == "fail" for c in rest) else "supported"
    if not verified:
        return "provisional", _unverified_note()
    return level, None


def verdict(
    coverage: float,
    tree_likeness: float | None,
    stability: StabilityEvidence,
    clonal: ClonalEvidence,
    thresholds: Thresholds,
    verified: bool | None = None,
) -> Verdict:
    """``verified`` defaults to the recorded verification; the verification itself passes True."""
    if verified is None:
        verified = VERDICT_VERIFIED
    t = thresholds
    checks = (
        _c1(coverage, t),
        _c2(tree_likeness, t),
        _stability("C3", CHECK_NAMES["C3"], stability.lineage, t.tau_lin, stability),
        _stability("C4", CHECK_NAMES["C4"], stability.pseudotime, t.tau_pt, stability),
        _c5(clonal, t),
    )
    level, note = _level(checks, verified)
    reasons = tuple(c.reason for c in checks if c.reason is not None) + ((note,) if note else ())
    return Verdict(level, checks, reasons)
