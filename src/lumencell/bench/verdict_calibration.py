"""Calibration of the tree-verdict taus: slow, run once and recorded.

The rule is ``[tree_support.calibration]`` in ``benchmarks/thresholds.toml``, declared before the run:
``tau = (quantile(matrix, matrix_quantile) + quantile(controls, control_quantile)) / 2`` with
``numpy.quantile(..., method="linear")``. A check whose matrix quantile is not above its control
quantile is non-discriminating and left out of the verdict.

- C3 (lineage stability) and C4 (pseudotime stability) pool every (dataset, method) run.
- C2 (tree-likeness) takes one value per dataset, since methods do not change it; an undefined
  value (a PAGA graph with no edges) counts as 0. It is measured at every Leiden resolution factor
  and kept only if it discriminates at all of them; ``tau_tree`` is the factor-1.0 tau.
- A run whose reference fails is left out of every check and counted in ``n_failed``.

Matrix: the ``[synthetic]`` topologies, preset and size. Controls: every ``[synthetic.controls]``
kind. No clonal layer is added, so C5 does not apply here.

``verify`` (Task 9) applies the calibrated verdict to fresh-seed records and checks that no control
run is verdicted supported.
"""

import dataclasses
import datetime
import json
import math
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Literal

import anndata as ad
import numpy as np

from lumencell.bench.synth import generate
from lumencell.bench.synth_controls import control_root_cell, generate_control
from lumencell.data.loader import Dataset, load
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import get_method, run_method
from lumencell.engine.roots import CellRoot, ClusterRoot, RootSpec
from lumencell.engine.stability import RunSpec, run_stability
from lumencell.engine.tree_support import (
    ClonalEvidence,
    StabilityEvidence,
    Thresholds,
    coverage,
    tree_likeness,
    verdict,
)
from lumencell.prep.pipeline import PrepResult, run_pipeline
from lumencell.prep.stages import STAGE_BY_NAME, PrepParams, run_stage

Group = Literal["matrix", "control"]
GROUPS: tuple[Group, ...] = ("matrix", "control")
C2_STAGES = ("leiden", "paga")


@dataclass(frozen=True)
class RunRecord:
    """One (dataset, method) run. ``dataset`` is ``<topology or kind>/seed<s>``.

    ``tree_likeness`` holds (resolution factor, value) pairs of the dataset's PAGA graph (NaN when
    it has no edges). An errored record has NaN values and its ``error`` set.
    """

    group: Group
    dataset: str
    seed: int
    method: str
    coverage: float
    tree_likeness: tuple[tuple[float, float], ...]
    lineage_stability: float
    pseudotime_stability: float
    n_lineages: int
    error: str | None = None


@dataclass(frozen=True)
class CheckCalibration:
    check: str
    factor: float | None  # the Leiden resolution factor for C2, None for C3 and C4
    matrix_quantile: float
    control_quantile: float
    n_matrix: int
    n_control: int
    discriminating: bool
    tau: float | None


@dataclass(frozen=True)
class Calibration:
    checks: tuple[CheckCalibration, ...]
    tau_tree: float | None
    tau_lin: float | None
    tau_pt: float | None
    non_discriminating: tuple[str, ...]
    n_failed: int


def _tree_likeness_at(
    ds: Dataset, prep: PrepResult, params: PrepParams, factors: Sequence[float]
) -> tuple[tuple[float, float], ...]:
    out: list[tuple[float, float]] = []
    for factor in factors:
        artifacts = prep.artifacts
        if factor != 1.0:
            artifacts = dict(prep.artifacts)
            scaled = replace(params, leiden_resolution=1.0 * factor)
            for name in C2_STAGES:
                produced, _ = run_stage(STAGE_BY_NAME[name], ds, artifacts, scaled)
                artifacts.update(produced)
        value = tree_likeness(artifacts["paga_connectivities"], artifacts["paga_tree"])
        out.append((float(factor), value))
    return tuple(out)


def _datasets(cfg: Mapping[str, Any], seeds: Mapping[Group, Sequence[int]]) -> list[tuple[Group, str, int]]:
    synth = cfg["synthetic"]
    out: list[tuple[Group, str, int]] = []
    for group, names in (("matrix", synth["topologies"]), ("control", synth["controls"]["kinds"])):
        for name in names:
            for seed in seeds.get(group, ()):  # type: ignore[call-overload]
                out.append((group, str(name), int(seed)))  # type: ignore[arg-type]
    return out


def _generate(cfg: Mapping[str, Any], group: Group, name: str, seed: int) -> tuple[ad.AnnData, RootSpec]:
    synth = cfg["synthetic"]
    if group == "matrix":
        adata = generate(name, synth["preset"], n_cells=synth["n_cells"], n_genes=synth["n_genes"], seed=seed)
        return adata, ClusterRoot("true_label", str(adata.uns["lumencell_truth"]["root_label"]))
    controls = synth["controls"]
    adata = generate_control(
        name,  # type: ignore[arg-type]
        n_cells=controls["n_cells"],
        n_genes=controls["n_genes"],
        seed=seed,
        preset=synth["preset"],
    )
    return adata, CellRoot(control_root_cell(adata))


def _failed(
    group: Group, dataset: str, seed: int, method: str, trees: tuple[tuple[float, float], ...], error: str
) -> RunRecord:
    print(f"verdict calibration: {dataset} {method} left out: {error}", file=sys.stderr)
    return RunRecord(group, dataset, seed, method, math.nan, trees, math.nan, math.nan, 0, error)


def collect_records(
    cfg: dict[str, Any], seeds: Mapping[Group, Sequence[int]], methods: Sequence[str], workdir: Path
) -> list[RunRecord]:
    """Every (dataset, installed method) run of the calibration, in a fixed order."""
    factors = [float(f) for f in cfg["tree_support"]["calibration"]["c2_resolution_factors"]]
    installed = []
    for method in methods:
        if get_method(method).available().ok:
            installed.append(method)
        else:
            print(f"verdict calibration: {method} is not installed and is skipped", file=sys.stderr)
    records: list[RunRecord] = []
    for group, name, seed in _datasets(cfg, seeds):
        dataset = f"{name}/seed{seed}"
        adata, root = _generate(cfg, group, name, seed)
        path = workdir / f"{group}-{name}-{seed}.h5ad"
        adata.write_h5ad(path)
        ds = load(path)
        params = PrepParams()
        try:
            prep = run_pipeline(ds, params)
            trees = _tree_likeness_at(ds, prep, params, factors)
        except Exception as exc:
            err = f"prep: {type(exc).__name__}: {exc}"
            records += [_failed(group, dataset, seed, m, (), err) for m in installed]
            continue
        ctx = InferenceContext.from_prep(ds, prep, params)
        for method in installed:
            try:
                result = run_method(ctx, method, root)
            except Exception as exc:
                records.append(_failed(group, dataset, seed, method, trees, f"{type(exc).__name__}: {exc}"))
                continue
            try:
                stab = run_stability(ds, prep, params, result, RunSpec(method, result.params, None, None))
            except Exception as exc:
                err = f"stability: {type(exc).__name__}: {exc}"
                records.append(_failed(group, dataset, seed, method, trees, err))
                continue
            records.append(
                RunRecord(
                    group,
                    dataset,
                    seed,
                    method,
                    coverage(np.asarray(result.pseudotime, dtype=np.float64)),
                    trees,
                    stab.lineage_stability,
                    stab.pseudotime_stability,
                    len(result.lineages),
                )
            )
            print(
                f"verdict calibration: {dataset} {method} C3 {stab.lineage_stability:.3f} "
                f"C4 {stab.pseudotime_stability:.3f} failed replicates {stab.n_failed}",
                file=sys.stderr,
            )
    return records


def _quantile(values: Sequence[float], q: float) -> float:
    return float(np.quantile(np.asarray(values, dtype=float), q, method="linear")) if values else math.nan


def _calibrate(
    check: str,
    factor: float | None,
    matrix: Sequence[float],
    control: Sequence[float],
    cal: Mapping[str, Any],
) -> CheckCalibration:
    mq = _quantile(matrix, float(cal["matrix_quantile"]))
    cq = _quantile(control, float(cal["control_quantile"]))
    discriminating = bool(mq > cq)  # False when either is NaN
    tau = (mq + cq) / 2 if discriminating else None
    return CheckCalibration(check, factor, mq, cq, len(matrix), len(control), discriminating, tau)


def derive(records: Sequence[RunRecord], cal: dict[str, Any]) -> Calibration:
    """The taus by the declared rule, from the records that did not fail."""
    ok = [r for r in records if r.error is None]
    checks: list[CheckCalibration] = []
    c2: list[CheckCalibration] = []
    for factor in [float(f) for f in cal["c2_resolution_factors"]]:
        per_group: dict[Group, list[float]] = {g: [] for g in GROUPS}
        seen: set[tuple[Group, str]] = set()
        for r in ok:
            if (r.group, r.dataset) in seen:
                continue
            seen.add((r.group, r.dataset))
            value = dict(r.tree_likeness)[factor]
            per_group[r.group].append(0.0 if math.isnan(value) else value)
        c2.append(_calibrate("C2", factor, per_group["matrix"], per_group["control"], cal))
    checks += c2
    c3 = _calibrate(
        "C3",
        None,
        [r.lineage_stability for r in ok if r.group == "matrix"],
        [r.lineage_stability for r in ok if r.group == "control"],
        cal,
    )
    c4 = _calibrate(
        "C4",
        None,
        [r.pseudotime_stability for r in ok if r.group == "matrix"],
        [r.pseudotime_stability for r in ok if r.group == "control"],
        cal,
    )
    checks += [c3, c4]
    c2_kept = all(c.discriminating for c in c2)
    tau_tree = next(c.tau for c in c2 if c.factor == 1.0) if c2_kept else None
    non_discriminating = tuple(
        name
        for name, kept in (("C2", c2_kept), ("C3", c3.discriminating), ("C4", c4.discriminating))
        if not kept
    )
    return Calibration(tuple(checks), tau_tree, c3.tau, c4.tau, non_discriminating, len(records) - len(ok))


def expected_runs(cfg: Mapping[str, Any]) -> set[tuple[Group, str, str]]:
    """Every (group, dataset, method) the declared calibration runs."""
    cal = cfg["tree_support"]["calibration"]
    seeds: dict[Group, Sequence[int]] = {"matrix": cal["matrix_seeds"], "control": cal["control_seeds"]}
    return {
        (group, f"{name}/seed{seed}", str(method))
        for group, name, seed in _datasets(cfg, seeds)
        for method in cal["methods"]
    }


def calibration_payload(
    records: Sequence[RunRecord], calibration: Calibration, cfg: dict[str, Any], commit: str
) -> dict[str, Any]:
    present = {(r.group, r.dataset, r.method) for r in records}
    synth = {k: v for k, v in cfg["synthetic"].items() if k in ("preset", "topologies", "n_cells", "n_genes")}
    return {
        "generated": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
        "commit": commit,
        "config": {
            "calibration": cfg["tree_support"]["calibration"],
            "synthetic": synth,
            "controls": cfg["synthetic"]["controls"],
            "prep": dataclasses.asdict(PrepParams()),
        },
        "records": [dataclasses.asdict(r) for r in records],
        "checks": [dataclasses.asdict(c) for c in calibration.checks],
        "taus": {
            "tau_tree": calibration.tau_tree,
            "tau_lin": calibration.tau_lin,
            "tau_pt": calibration.tau_pt,
        },
        "non_discriminating": list(calibration.non_discriminating),
        "n_failed": calibration.n_failed,
        "complete": expected_runs(cfg) <= present,
    }


def _nan_to_none(value: Any) -> Any:
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, dict):
        return {k: _nan_to_none(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [_nan_to_none(v) for v in value]
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Strict JSON: NaN is written as ``null``, and any other non-finite number raises."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_nan_to_none(payload), indent=2, allow_nan=False) + "\n")


def _float(value: float | None) -> float:
    return math.nan if value is None else float(value)


def records_from_payload(doc: Mapping[str, Any]) -> list[RunRecord]:
    """The records of a written payload, ``null`` read back as NaN."""
    return [
        RunRecord(
            group=r["group"],
            dataset=r["dataset"],
            seed=int(r["seed"]),
            method=r["method"],
            coverage=_float(r["coverage"]),
            tree_likeness=tuple((float(f), _float(v)) for f, v in r["tree_likeness"]),
            lineage_stability=_float(r["lineage_stability"]),
            pseudotime_stability=_float(r["pseudotime_stability"]),
            n_lineages=int(r["n_lineages"]),
            error=r["error"],
        )
        for r in doc["records"]
    ]


@dataclass(frozen=True)
class Verification:
    """The verdict applied to fresh-seed records. ``verdicts`` holds (dataset, method, level)."""

    verdicts: tuple[tuple[str, str, str], ...]
    controls_supported: tuple[str, ...]  # "<dataset> <method>" of each control run verdicted supported
    matrix_supported_share: float  # reported, not gated; 0.0 with no matrix verdict
    n_failed: int
    passed: bool


def verify(records: Sequence[RunRecord], thresholds: Thresholds, max_controls_supported: int) -> Verification:
    """Verdict every run that did not fail (factor-1.0 tree-likeness, stability, no clones)."""
    verdicts: list[tuple[str, str, str]] = []
    controls: list[str] = []
    matrix_levels: list[str] = []
    for r in records:
        if r.error is not None:
            continue
        v = verdict(
            r.coverage,
            dict(r.tree_likeness)[1.0],
            StabilityEvidence("done", r.lineage_stability, r.pseudotime_stability),
            ClonalEvidence("not_configured", r.n_lineages),
            thresholds,
        )
        verdicts.append((r.dataset, r.method, v.level))
        if r.group == "control" and v.level == "supported":
            controls.append(f"{r.dataset} {r.method}")
        if r.group == "matrix":
            matrix_levels.append(v.level)
    share = matrix_levels.count("supported") / len(matrix_levels) if matrix_levels else 0.0
    n_failed = sum(r.error is not None for r in records)
    return Verification(
        tuple(verdicts), tuple(controls), share, n_failed, len(controls) <= max_controls_supported
    )


def expected_verification_runs(cfg: Mapping[str, Any]) -> set[tuple[Group, str, str]]:
    """Every (group, dataset, method) the declared verification runs."""
    ver = cfg["tree_support"]["verification"]
    seeds: dict[Group, Sequence[int]] = {"matrix": ver["matrix_seeds"], "control": ver["control_seeds"]}
    methods = cfg["tree_support"]["calibration"]["methods"]
    return {
        (group, f"{name}/seed{seed}", str(method))
        for group, name, seed in _datasets(cfg, seeds)
        for method in methods
    }


def verification_payload(
    records: Sequence[RunRecord],
    verification: Verification,
    thresholds: Thresholds,
    cfg: dict[str, Any],
    commit: str,
) -> dict[str, Any]:
    present = {(r.group, r.dataset, r.method) for r in records}
    group_of = {(r.dataset, r.method): r.group for r in records}
    counts: dict[str, dict[str, int]] = {g: {} for g in GROUPS}
    by_kind: dict[str, int] = {}
    for dataset, method, level in verification.verdicts:
        group = group_of[(dataset, method)]
        counts[group][level] = counts[group].get(level, 0) + 1
        if group == "control" and level == "supported":
            kind = dataset.split("/")[0]
            by_kind[kind] = by_kind.get(kind, 0) + 1
    synth = {k: v for k, v in cfg["synthetic"].items() if k in ("preset", "topologies", "n_cells", "n_genes")}
    ver = cfg["tree_support"]["verification"]
    return {
        "generated": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
        "commit": commit,
        "config": {
            "verification": ver,
            "methods": cfg["tree_support"]["calibration"]["methods"],
            "synthetic": synth,
            "controls": cfg["synthetic"]["controls"],
            "prep": dataclasses.asdict(PrepParams()),
            "thresholds": dataclasses.asdict(thresholds),
        },
        "records": [dataclasses.asdict(r) for r in records],
        "verdicts": [list(v) for v in verification.verdicts],
        "level_counts": counts,
        "controls_supported": list(verification.controls_supported),
        "controls_supported_by_kind": by_kind,
        "matrix_supported_share": verification.matrix_supported_share,
        "n_failed": verification.n_failed,
        "max_controls_supported": ver["max_controls_supported"],
        "passed": verification.passed,
        "complete": expected_verification_runs(cfg) <= present,
    }
