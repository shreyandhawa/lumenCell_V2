"""Run a method on a benchmark case and score it against the truth."""

import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

from lumencell.bench.metrics import Scores, score
from lumencell.bench.synth import generate, truth_from_adata
from lumencell.data.loader import load
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.prep.pipeline import run_pipeline
from lumencell.prep.stages import PrepParams

ClonalRowStatus = Literal[
    "ok", "single_lineage", "too_few_clones", "not_assessable", "failed", "not_applicable", "not_measured"
]

_NAN = math.nan


@dataclass(frozen=True)
class ClonalRow:
    """Clonal agreement of one case; see [clonal] in benchmarks/thresholds.toml."""

    status: ClonalRowStatus
    excess_ratio: float = _NAN
    p: float = _NAN
    delta_base: float = _NAN
    eta_oracle: float = _NAN
    eta_method: float = _NAN
    n_clones: int = 0
    n_cells: int = 0
    effective_clones: float = _NAN
    eta_m3: float = _NAN
    p_m3: float = _NAN
    reason: str | None = None


@dataclass(frozen=True)
class CaseResult:
    case: str
    method: str
    scores: Scores
    n_clusters: int
    runtime_seconds: float
    warnings: list[str] = field(default_factory=list)
    clonal: ClonalRow = ClonalRow("not_measured")

    def row(self) -> dict[str, Any]:
        s, c = self.scores, self.clonal
        return {
            "case": self.case,
            "method": self.method,
            "pt_spearman_global": round(s.pt_spearman_global, 4),
            "pt_spearman_lineage": round(s.pt_spearman_lineage, 4),
            "branch_accuracy": round(s.branch_accuracy, 4),
            "branch_accuracy_all": round(s.branch_accuracy_all, 4),
            "lineage_recall": round(s.lineage_recall, 4),
            "n_spurious": s.n_spurious,
            "root_pt_median": round(s.root_pt_median, 4),
            "root_pt_excess": round(s.root_pt_excess, 4),
            "terminal_share": round(s.terminal_share, 4),
            "coverage": round(s.coverage, 4),
            "n_clusters": self.n_clusters,
            "runtime_seconds": round(self.runtime_seconds, 2),
            "clonal_excess_ratio": round(c.excess_ratio, 4),
            "clonal_p": round(c.p, 4),
            "clonal_delta_base": round(c.delta_base, 4),
            "clonal_status": str(c.status),
            "clonal_eta_oracle": round(c.eta_oracle, 4),
            "clonal_eta_method": round(c.eta_method, 4),
            "clonal_effective_clones": round(c.effective_clones, 4),
            "clonal_eta_m3": round(c.eta_m3, 4),
            "clonal_p_m3": round(c.p_m3, 4),
            "clonal_n_clones": int(c.n_clones),
            "clonal_n_cells": int(c.n_cells),
            "clonal_reason": c.reason,
        }


def run_synthetic_case(
    workdir: Path,
    topology: str,
    method: str = "paga_dpt",
    preset: str = "realistic",
    seed: int = 0,
    n_cells: int = 3000,
    n_genes: int = 2000,
    params: dict[str, Any] | None = None,
    prep_params: PrepParams | None = None,
) -> CaseResult:
    """Generate a case, run the real pipeline and method on it, and score the result.

    The method only sees counts and the root. The root is the true starting
    cluster, given the way a user would give it (``cluster:true_label=trunk_early``);
    no other truth reaches inference.
    """
    case = f"{topology}/{preset}/seed{seed}"
    adata = generate(topology, preset, n_cells=n_cells, n_genes=n_genes, seed=seed)
    truth = truth_from_adata(adata)
    path = workdir / f"{topology}-{preset}-{seed}.h5ad"
    adata.write_h5ad(path)

    start = time.perf_counter()
    ds = load(path)
    pp = prep_params or PrepParams()
    prep = run_pipeline(ds, pp)
    ctx = InferenceContext.from_prep(ds, prep, pp)
    result = run_method(ctx, method, ClusterRoot("true_label", truth.meta["root_label"]), params=params)
    elapsed = time.perf_counter() - start

    return CaseResult(
        case=case,
        method=method,
        scores=score(result.pseudotime, result.membership, truth),
        n_clusters=int(ctx.leiden.max()) + 1,
        runtime_seconds=elapsed,
        warnings=[w.code for w in result.warnings],
    )


def run_real_case(
    dataset_id: str, method: str = "paga_dpt", params: dict[str, Any] | None = None
) -> CaseResult:
    """Run a method on a downloaded real dataset (prep is cached) and score it."""
    from lumencell.bench.datasets import DATASETS, is_cached, local_path, truth_for
    from lumencell.engine.roots import parse_root
    from lumencell.prep.cache import PrepCache

    ds = DATASETS[dataset_id]
    if not is_cached(ds):
        raise FileNotFoundError(f"{dataset_id} is not downloaded; run: lumencell bench fetch {dataset_id}")
    start = time.perf_counter()
    data = load(local_path(ds))
    pp = PrepParams()
    prep = PrepCache(data).run(pp)
    ctx = InferenceContext.from_prep(data, prep, pp)
    result = run_method(ctx, method, parse_root(ds.root_spec), params=params, label_key=ds.label_key)
    elapsed = time.perf_counter() - start
    return CaseResult(
        case=dataset_id,
        method=method,
        scores=score(result.pseudotime, result.membership, truth_for(ds, data.obs)),
        n_clusters=int(ctx.leiden.max()) + 1,
        runtime_seconds=elapsed,
        warnings=[w.code for w in result.warnings],
    )
