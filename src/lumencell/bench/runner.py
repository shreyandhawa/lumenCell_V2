"""Run a method on a benchmark case and score it against the truth."""

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from lumencell.bench.metrics import Scores, score
from lumencell.bench.synth import generate, truth_from_adata
from lumencell.data.loader import load
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.prep.pipeline import run_pipeline
from lumencell.prep.stages import PrepParams


@dataclass(frozen=True)
class CaseResult:
    case: str
    method: str
    scores: Scores
    n_clusters: int
    runtime_seconds: float
    warnings: list[str] = field(default_factory=list)

    def row(self) -> dict[str, Any]:
        s = self.scores
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
            "coverage": round(s.coverage, 4),
            "n_clusters": self.n_clusters,
            "runtime_seconds": round(self.runtime_seconds, 2),
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
