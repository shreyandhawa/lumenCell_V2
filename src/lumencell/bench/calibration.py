"""Calibration of the clonal nulls: slow, run once and recorded.

Definitions are in ``[clonal.calibration]`` of ``benchmarks/thresholds.toml``, declared before the
run. Each check draws one p-value per seed and counts rejections (p <= alpha) against the upper
bound of a binomial: ``int(binom.ppf(quantile, n_scored, alpha))``, 18 for 200 seeds. A p-value
that is NaN is kept as ``None``, counted in ``n_nan`` and left out of ``n_scored`` and of the
rejections. A method that raises on a seed is a NaN p and is listed in ``config["errors"]``.

Checks (a rejection is p <= alpha):

- ``null/<method>``: fate_bias 0, own-cells mode on the method's own masks and pseudotime, headline p.
- ``donor_confound/N1``, ``donor_confound/N1_donors_ignored``, ``depth_only/N2``, ``depth_only/N1``:
  oracle scores (true lineages, true order, oracle cells) on a constant baseline.
- Monotonicity: oracle η rises strictly with fate_bias, and the oracle's p at fate_bias 0 is above 0.01.
"""

import dataclasses
import json
import math
import subprocess
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy.stats import binom

from lumencell.bench.synth import generate, truth_from_adata
from lumencell.bench.synth_clones import ClonalLayer, Control, add_clones, clone_truth
from lumencell.bench.truth import Truth
from lumencell.clonal.labels import CloneCodes, encode_clones
from lumencell.clonal.scores import ClonalResult, score_clonal, score_with_baseline
from lumencell.data.loader import load
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.prep.pipeline import run_pipeline
from lumencell.prep.stages import PrepParams

REPO = Path(__file__).resolve().parents[3]
CLONE_COLUMNS = ("clone_id", "donor")
# Monotonicity's second condition: the oracle's p at fate_bias 0 must be above this.
P_AT_ZERO_MIN = 0.01


def rejection_bound(n_scored: int, alpha: float, quantile: float) -> int:
    """The most rejections a calibrated check shows at ``quantile``: ``binom.ppf(quantile, n, alpha)``."""
    return int(binom.ppf(quantile, n_scored, alpha)) if n_scored > 0 else 0


@dataclass(frozen=True)
class RejectionCount:
    """One check: ``n`` seeds, ``n_nan`` of them without a p-value (``None`` in ``p_values``).

    ``bound`` comes from the ``n - n_nan`` seeds that scored. A check that expects calibration
    passes when ``rejections <= bound``; one that expects none (a control with teeth) when
    ``rejections > bound``. A check with no scored seed never passes, whichever it expects.
    """

    name: str
    null: str
    n: int
    n_nan: int
    rejections: int
    bound: int
    expect_calibrated: bool
    p_values: tuple[float | None, ...]

    @property
    def n_scored(self) -> int:
        return self.n - self.n_nan

    @property
    def passed(self) -> bool:
        if self.n_scored == 0:
            return False
        within = self.rejections <= self.bound
        return within if self.expect_calibrated else not within


def count_rejections(
    name: str,
    null: str,
    expect_calibrated: bool,
    p_values: Sequence[float],
    alpha: float,
    quantile: float,
) -> RejectionCount:
    """Count ``p <= alpha``, keeping a NaN p as ``None`` and counting it apart."""
    stored = tuple(None if math.isnan(p) else float(p) for p in p_values)
    scored = [p for p in stored if p is not None]
    return RejectionCount(
        name=name,
        null=null,
        n=len(stored),
        n_nan=len(stored) - len(scored),
        rejections=sum(p <= alpha for p in scored),
        bound=rejection_bound(len(scored), alpha, quantile),
        expect_calibrated=expect_calibrated,
        p_values=stored,
    )


@dataclass(frozen=True)
class MonotonicCase:
    """The oracle's η at each ``fate_bias`` on one topology and seed, and its p at fate_bias 0."""

    topology: str
    seed: int
    fate_bias: tuple[float, ...]
    eta: tuple[float, ...]
    p_at_zero: float

    @property
    def passed(self) -> bool:
        rising = all(a < b for a, b in zip(self.eta, self.eta[1:], strict=False))  # NaN is never rising
        return rising and self.p_at_zero > P_AT_ZERO_MIN


@dataclass(frozen=True)
class CalibrationReport:
    checks: tuple[RejectionCount, ...]
    monotonic: tuple[MonotonicCase, ...]
    config: dict[str, Any]
    commit: str

    @property
    def passed(self) -> bool:
        """Every check and monotonic case passes and no seed raised an error."""
        return (
            not self.config.get("errors")
            and all(c.passed for c in self.checks)
            and all(m.passed for m in self.monotonic)
        )


def run_calibration(cfg: dict[str, Any], methods: Sequence[str], workdir: Path) -> CalibrationReport:
    cal, syn, layer_cfg = cfg["clonal"]["calibration"], cfg["synthetic"], cfg["clonal"]["synthetic"]
    alpha, quantile = cal["alpha"], cal["binomial_quantile"]
    topologies: list[str] = cal["topologies"]
    seeds = range(cal["seed_start"], cal["seed_start"] + cal["n_seeds"])
    errors: list[str] = []

    def gen(topology: str, seed: int) -> ad.AnnData:
        return generate(topology, syn["preset"], n_cells=syn["n_cells"], n_genes=syn["n_genes"], seed=seed)

    def layer(control: Control = "none", fate_bias: float = 0.0) -> ClonalLayer:
        return ClonalLayer(
            fate_bias=fate_bias,
            n_donors=layer_cfg["n_donors"],
            label_rate=layer_cfg["label_rate"],
            control=control,
        )

    ps: dict[str, list[float]] = {f"null/{m}": [] for m in methods}
    for key in ("donor_confound/N1", "donor_confound/N1_donors_ignored", "depth_only/N2", "depth_only/N1"):
        ps[key] = []

    for seed in seeds:
        topology = topologies[(seed - cal["seed_start"]) % len(topologies)]
        base = gen(topology, seed)

        # null/<method>: prep once, shared by the methods.
        for method, p in _null_p_values(base, layer(), seed, methods, workdir, errors).items():
            ps[f"null/{method}"].append(p)

        confound = _oracle(base, layer("donor_confound"), seed)
        ps["donor_confound/N1"].append(confound.nulls[0].p if confound.nulls else math.nan)
        ps["donor_confound/N1_donors_ignored"].append(_donors_ignored_n1(base, layer("donor_confound"), seed))

        depth_base = base if topology == cal["depth_only_topology"] else gen(cal["depth_only_topology"], seed)
        depth = _oracle(depth_base, layer("depth_only"), seed)
        ps["depth_only/N2"].append(_null_p(depth, 1))
        ps["depth_only/N1"].append(_null_p(depth, 0))
        print(f"calibration: seed {seed} done", file=sys.stderr, flush=True)

    # name -> (null the p belongs to, whether calibration is expected)
    controls = {
        "donor_confound/N1": ("N1", True),
        "donor_confound/N1_donors_ignored": ("N1", False),
        "depth_only/N2": ("N2", True),
        "depth_only/N1": ("N1", False),
    }
    controls.update({f"null/{m}": ("headline", True) for m in methods})
    checks = [count_rejections(name, *controls[name], values, alpha, quantile) for name, values in ps.items()]

    monotonic = tuple(
        _monotonic_case(gen(topology, seed), topology, seed, cal["monotonic_fate_bias"], layer)
        for topology in cal["monotonic_topologies"]
        for seed in cal["monotonic_seeds"]
    )
    config = {
        "calibration": cal,
        "synthetic": syn,
        "clonal_synthetic": layer_cfg,
        "methods": list(methods),
        "errors": errors,
    }
    return CalibrationReport(tuple(checks), monotonic, config, _commit())


def write_report(report: CalibrationReport, path: Path) -> None:
    """The report as strict JSON, ``passed`` at the top level; an undefined (NaN) value is ``null``."""
    doc = {
        "passed": report.passed,
        "commit": report.commit,
        "config": report.config,
        "checks": [
            {**dataclasses.asdict(c), "n_scored": c.n_scored, "passed": c.passed} for c in report.checks
        ],
        "monotonic": [{**dataclasses.asdict(m), "passed": m.passed} for m in report.monotonic],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_nan_to_none(doc), indent=2, allow_nan=False) + "\n")


def _nan_to_none(value: Any) -> Any:
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, dict):
        return {k: _nan_to_none(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [_nan_to_none(v) for v in value]
    return value


def _null_p_values(
    base: ad.AnnData,
    layer: ClonalLayer,
    seed: int,
    methods: Sequence[str],
    workdir: Path,
    errors: list[str],
) -> dict[str, float]:
    """Headline own-cells p of each method on one seed (NaN where a method raises or is undefined)."""
    truth = truth_from_adata(base)
    adata = add_clones(base, layer, seed)
    path = workdir / f"null-{seed}.h5ad"
    adata.write_h5ad(path)
    out = {m: math.nan for m in methods}
    try:
        ds = load(path)
        pp = PrepParams()
        prep = run_pipeline(ds, pp)
        ctx = InferenceContext.from_prep(ds, prep, pp, clonal_columns=CLONE_COLUMNS)
        pca = np.asarray(prep.artifacts["pca"])
        codes = encode_clones(
            adata.obs["clone_id"].to_numpy(dtype=object), adata.obs["donor"].to_numpy(dtype=object)
        )
        for method in methods:
            try:
                result = run_method(ctx, method, ClusterRoot("true_label", truth.meta["root_label"]))
                score = score_with_baseline(
                    result.membership,
                    codes,
                    result.pseudotime,
                    pca,
                    pp.neighbors_pcs,
                    int(result.membership.shape[0]),
                )
                out[method] = score.p
            except Exception as exc:
                errors.append(f"seed {seed} {method}: {type(exc).__name__}: {exc}")
                print(errors[-1], file=sys.stderr, flush=True)
    except Exception as exc:
        errors.append(f"seed {seed} prep: {type(exc).__name__}: {exc}")
        print(errors[-1], file=sys.stderr, flush=True)
    finally:
        path.unlink(missing_ok=True)
    return out


def _oracle_inputs(base: ad.AnnData, layer: ClonalLayer, seed: int) -> tuple[Truth, CloneCodes]:
    truth = clone_truth(truth_from_adata(base), add_clones(base, layer, seed))
    assert truth.clone_ids is not None and truth.donors is not None
    return truth, encode_clones(truth.clone_ids, truth.donors)


def _score_oracle(truth: Truth, codes: CloneCodes) -> ClonalResult:
    """The true lineages on their single-lineage cells, true order, constant baseline."""
    return score_clonal(
        truth.membership(),
        codes,
        truth.order,
        np.zeros(truth.n_obs, dtype=np.int64),
        oracle_cells=truth.exclusive() >= 0,
    )


def _oracle(base: ad.AnnData, layer: ClonalLayer, seed: int) -> ClonalResult:
    return _score_oracle(*_oracle_inputs(base, layer, seed))


def _donors_ignored_n1(base: ad.AnnData, layer: ClonalLayer, seed: int) -> float:
    """N1 on the same data with every donor code 0 (clone ids stay keyed by donor)."""
    truth, codes = _oracle_inputs(base, layer, seed)
    result = _score_oracle(truth, dataclasses.replace(codes, donor=np.zeros_like(codes.donor)))
    return _null_p(result, 0)


def _null_p(result: ClonalResult, index: int) -> float:
    return result.nulls[index].p if len(result.nulls) > index else math.nan


def _monotonic_case(
    base: ad.AnnData, topology: str, seed: int, fate_bias: Sequence[float], layer: Callable[..., ClonalLayer]
) -> MonotonicCase:
    results = [_oracle(base, layer(fate_bias=fb), seed) for fb in fate_bias]
    zero = list(fate_bias).index(0.0)
    return MonotonicCase(
        topology=topology,
        seed=seed,
        fate_bias=tuple(float(fb) for fb in fate_bias),
        eta=tuple(r.eta for r in results),
        p_at_zero=results[zero].p,
    )


def _commit() -> str:
    """HEAD, with ``-dirty`` appended when the working tree has uncommitted changes."""
    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True, check=True, timeout=10
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"], cwd=REPO, capture_output=True, text=True, check=True, timeout=10
        ).stdout.strip()
        return f"{head}-dirty" if status else head
    except (OSError, subprocess.SubprocessError):
        return "unknown"
