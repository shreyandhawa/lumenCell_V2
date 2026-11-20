"""The case-study numbers: every method on each patient alone, and a pooled run as a diagnostic.

The declaration (``[case_study.<name>]`` in benchmarks/thresholds.toml) is merged before any of this
runs. Each method runs once per unit (the patient), because lumenCell has no batch correction; the
pooled run, when declared, is reported apart and labelled a diagnostic, never the headline.

Every number comes from the phase 1 and phase 2 code, called as the worker and ``GET`` call it:
``agreement_for_run`` (own-cells mode, m = 2 and m = 3, N1-N3, the baseline, coupling),
``run_stability`` with the reference run's params, the verdict with no ``verified`` argument (so a
failed verification caps it), and ``cell_cycle_check``. Nothing is re-implemented here.

Clone sharing says that lineages are related, never which came first; nothing here says otherwise.
Unit values come from the file: they never reach a path (unit files are named by a hash).
"""

import hashlib
import json
import os
import sys
import tempfile
import traceback
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np

from lumencell import __version__
from lumencell.api.models import ClonalSettings
from lumencell.bench.gates import load_thresholds
from lumencell.clonal.labels import is_missing
from lumencell.data.loader import load
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import run_method
from lumencell.engine.roots import RootSpec, format_root, parse_root
from lumencell.engine.stability import RunSpec, StabilityError, run_stability
from lumencell.engine.tree_support import (
    StabilityEvidence,
    calibrated_thresholds,
    coverage,
    tree_likeness,
    verdict,
)
from lumencell.prep.cache import PrepCache
from lumencell.prep.stages import PrepParams
from lumencell.trust import constants
from lumencell.trust.cell_cycle import cell_cycle_check
from lumencell.trust.clonal import agreement_for_run, check_settings, settings_columns
from lumencell.trust.stability import not_measured, to_api_stability
from lumencell.trust.verdict import clonal_evidence, stability_evidence, to_api_tree_support

TREE_LIKENESS_NOTE = "no PAGA edges, so tree-likeness is undefined and C2 fails"
PER_UNIT = "Per patient (primary)"
POOLED = "Pooled diagnostic (not the headline: patient structure is not corrected)"
MISSING = "–"


@dataclass(frozen=True)
class CaseStudyDeclaration:
    settings: ClonalSettings
    root: RootSpec
    methods: list[str]
    unit_key: str
    pooled_diagnostic: bool
    # The TOML section as written, reported with the results; empty for a declaration built in code.
    section: Mapping[str, Any] = field(default_factory=dict, compare=False)


def load_declaration(name: str, cfg: dict[str, Any] | None = None) -> CaseStudyDeclaration:
    """``[case_study.<name>]`` of benchmarks/thresholds.toml; a missing field raises ``KeyError``."""
    s = (cfg if cfg is not None else load_thresholds())["case_study"][name]
    settings = ClonalSettings(
        clone_key=s["clone_key"],
        donor_key=s.get("donor_key"),
        single_donor=s.get("single_donor", False),
        strata_keys=list(s["strata_keys"]),
    )
    return CaseStudyDeclaration(
        settings=settings,
        root=parse_root(s["root"]),
        methods=list(s["methods"]),
        unit_key=s["unit_key"],
        pooled_diagnostic=bool(s["pooled_diagnostic"]),
        section=dict(s),
    )


def _unit_name(value: str) -> str:
    return f"unit-{hashlib.sha256(value.encode()).hexdigest()[:12]}.h5ad"


def _write_h5ad_atomic(adata: ad.AnnData, target: Path) -> None:
    fd, name = tempfile.mkstemp(prefix=".unit-", suffix=".tmp.h5ad", dir=target.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        adata.write_h5ad(tmp)
        os.replace(tmp, target)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def split_units(path: Path, unit_key: str, workdir: Path) -> dict[str, Path]:
    """One ``.h5ad`` per distinct non-missing value of obs ``unit_key``, X, var and obs unchanged.

    Files are named by the value's hash inside ``workdir``; the input file is only read.
    """
    adata = ad.read_h5ad(path)
    values = np.array([str(v) for v in adata.obs[unit_key].to_numpy(dtype=object)], dtype=object)
    present = ~is_missing(adata.obs[unit_key].to_numpy(dtype=object))
    workdir.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    for value in sorted({str(v) for v in values[present]}):
        target = workdir / _unit_name(value)
        _write_h5ad_atomic(adata[present & (values == value)].copy(), target)
        out[value] = target
    return out


def _versions() -> dict[str, str]:
    return {"lumencell": __version__, "scanpy": version("scanpy"), "numpy": np.__version__}


def _method_failed(exc: Exception) -> dict[str, Any]:
    return {
        "error": f"{type(exc).__name__}: {exc}",
        "n_lineages": None,
        "coverage": None,
        "runtime_seconds": None,
        "clonal": None,
        "stability": None,
        "tree_support": None,
        "cell_cycle": None,
    }


def run_case_study(
    path: Path,
    settings: ClonalSettings,
    root: RootSpec,
    methods: Sequence[str],
    include_stability: bool = True,
    *,
    cache_root: Path | None = None,
) -> dict[str, Any]:
    """Every method on one file, with its clonal agreement, stability, verdict and cell-cycle check.

    ``cache_root`` is where the prep cache goes (the user's cache when None). A method that raises is
    recorded with its error and the others still run; the settings not fitting the file raises.
    """
    ds = load(path)
    problems = check_settings(ds.obs, settings)
    if problems:
        raise ValueError(f"the clonal settings do not fit {path.name}: {'; '.join(problems)}")
    params = PrepParams()
    prep = PrepCache(ds, cache_root).run(params)
    columns = settings_columns(settings)
    ctx = InferenceContext.from_prep(ds, prep, params, clonal_columns=columns)
    pca = np.asarray(prep.artifacts["pca"])
    tree = tree_likeness(prep.artifacts["paga_connectivities"], prep.artifacts["paga_tree"])

    out: dict[str, Any] = {}
    for method in methods:
        print(f"case study: {path.name} {method}", file=sys.stderr, flush=True)
        try:
            result = run_method(ctx, method, root)
        except Exception as exc:  # recorded: a method that cannot run on this data is a result
            print(f"case study: {method} failed on {path.name}:", file=sys.stderr)
            traceback.print_exception(exc, file=sys.stderr)
            out[method] = _method_failed(exc)
            continue
        agreement = agreement_for_run(result, ds.obs, settings, pca, params.neighbors_pcs, result.obs_keys)
        if include_stability:
            # The reference run's params, as the worker passes graph["params"], so method seeds reseed.
            spec = RunSpec(method, result.params, None, result.label_key, columns)
            try:
                stability = to_api_stability(run_stability(ds, prep, params, result, spec))
            except StabilityError as exc:
                stability = not_measured("failed", str(exc))
            evidence = stability_evidence(stability)
        else:
            stability, evidence = None, StabilityEvidence("pending")
        support = to_api_tree_support(
            verdict(
                coverage(result.pseudotime.astype(np.float64)),
                tree,
                evidence,
                clonal_evidence(agreement, len(result.lineages)),
                calibrated_thresholds(),
            )
        )
        cycle = cell_cycle_check(
            result.primary_lineage,
            [lin.name for lin in result.lineages],
            prep.artifacts["cc_phase"],
            prep.artifacts["cc_info"],
        )
        out[method] = {
            "n_lineages": len(result.lineages),
            "coverage": coverage(result.pseudotime.astype(np.float64)),
            "runtime_seconds": float(result.runtime_seconds),
            "clonal": agreement.model_dump(mode="json"),
            "stability": None if stability is None else stability.model_dump(mode="json"),
            "tree_support": support.model_dump(mode="json"),
            "cell_cycle": cycle.model_dump(mode="json"),
        }

    finite = bool(np.isfinite(tree))
    return {
        "dataset": {
            "name": path.name,
            "n_cells": int(ds.n_obs),
            "n_genes": int(ds.n_vars),
            "fingerprint": ds.fingerprint,
            "tree_likeness": float(tree) if finite else None,
            "tree_likeness_note": None if finite else TREE_LIKENESS_NOTE,
        },
        "settings": settings.model_dump(mode="json"),
        "root": format_root(root),
        "methods": out,
        "versions": _versions(),
    }


def _declaration_record(decl: CaseStudyDeclaration) -> dict[str, Any]:
    if decl.section:
        return dict(decl.section)
    return {
        **decl.settings.model_dump(mode="json"),
        "root": format_root(decl.root),
        "methods": list(decl.methods),
        "unit_key": decl.unit_key,
        "pooled_diagnostic": decl.pooled_diagnostic,
    }


def run_units(
    path: Path, decl: CaseStudyDeclaration, workdir: Path, include_stability: bool = True
) -> dict[str, Any]:
    """Each unit alone, then (when declared) the pooled run; unit files and prep caches in ``workdir``."""
    files = split_units(path, decl.unit_key, workdir / "units")
    cache = workdir / "cache"

    def run(p: Path) -> dict[str, Any]:
        return run_case_study(p, decl.settings, decl.root, decl.methods, include_stability, cache_root=cache)

    units = {value: run(f) for value, f in files.items()}
    return {
        "declaration": _declaration_record(decl),
        "units": units,
        "pooled_diagnostic": run(path) if decl.pooled_diagnostic else None,
    }


def write_report(report: dict[str, Any], path: Path, commit: str) -> None:
    """Strict JSON (NaN refused before anything is written), replaced atomically."""
    payload = {**report, "generated": datetime.now(UTC).date().isoformat(), "commit": commit}
    text = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.stem}-", suffix=".tmp", dir=path.parent)
    tmp = Path(name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


# ---- markdown -------------------------------------------------------------------------


def _cell(v: Any) -> str:
    """A table cell: a missing value is a dash; text from the file cannot break the table."""
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return MISSING
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return f"{v:.3f}"
    return " ".join(str(v).split()).replace("|", "\\|") or MISSING


def _p(v: Any) -> str:
    if v is None or not np.isfinite(v):
        return MISSING
    return f"{v:.3g}"


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return out


def _runs(report: dict[str, Any], pooled: bool) -> list[tuple[str, str, dict[str, Any]]]:
    if pooled:
        p = report.get("pooled_diagnostic")
        return [] if p is None else [("pooled", m, r) for m, r in p["methods"].items()]
    return [(u, m, r) for u, unit in report["units"].items() for m, r in unit["methods"].items()]


ROW_HEADER = (
    "Patient",
    "Method",
    "L",
    "Verdict",
    "Clonal status",
    "Cells",
    "Clones",
    "Effective clones",
    "Strata",
    "η (N1)",
    "η (N2)",
    "η (N3)",
    "Headline p",
    "η_base",
    "Δη",
    "p_Δη",
    "η (m = 3)",
    "p (m = 3)",
    "Purity",
    "Purity without largest",
)


def _row(unit: str, method: str, m: dict[str, Any], donors: bool) -> list[str]:
    if m.get("error"):
        cells = [_cell(unit), _cell(method), MISSING, "failed", MISSING]
        return cells + [MISSING] * (len(ROW_HEADER) - len(cells) + int(donors))
    c = m["clonal"]
    s = c["primary"] or {}
    m3 = c["sensitivity_m3"] or {}
    eta = {n["null"]: n["eta"] for n in s.get("nulls", [])}
    row = [
        _cell(unit),
        _cell(method),
        _cell(m["n_lineages"]),
        _cell(m["tree_support"]["level"]),
        _cell(c["status"]),
        _cell(s.get("n_cells")),
        _cell(s.get("n_clones")),
        _cell(s.get("effective_clones")),
        _cell(s.get("n_strata")),
        _cell(eta.get("N1")),
        _cell(eta.get("N2")),
        _cell(eta.get("N3")),
        _p(s.get("p")),
        _cell(s.get("eta_base")),
        _cell(s.get("delta_eta")),
        _p(s.get("p_delta")),
        _cell(m3.get("eta")),
        _p(m3.get("p")),
        _cell(s.get("purity")),
        _cell(s.get("purity_without_largest")),
    ]
    if donors:
        row.insert(8, _cell(s.get("n_donors")))
    return row


def _row_table(runs: list[tuple[str, str, dict[str, Any]]], donors: bool) -> list[str]:
    header = list(ROW_HEADER)
    if donors:
        header.insert(8, "Donors scored")
    return _table(header, [_row(u, meth, m, donors) for u, meth, m in runs])


def _coupling_table(runs: list[tuple[str, str, dict[str, Any]]]) -> list[str]:
    rows = [
        [
            _cell(u),
            _cell(meth),
            _cell(k["a"]),
            _cell(k["b"]),
            _cell(k["observed"]),
            _cell(k["z"]),
            _p(k["p_adj"]),
        ]
        for u, meth, m in runs
        if m.get("clonal")
        for k in m["clonal"]["coupling"]
    ]
    return _table(("Patient", "Method", "a", "b", "Observed", "z", "p_adj"), rows)


def _stability_table(runs: list[tuple[str, str, dict[str, Any]]]) -> list[str]:
    rows = []
    for u, meth, m in runs:
        st = m.get("stability")
        if st is None:
            rows.append([_cell(u), _cell(meth), "not measured"] + [MISSING] * 9)
            continue
        for kind in ("reseeds", "subsamples"):
            k = st[kind]
            if k is None:
                rows.append([_cell(u), _cell(meth), _cell(f"{kind}: {st['status']}")] + [MISSING] * 9)
                continue
            failed = sum(r["error"] is not None for r in k["replicates"])
            rows.append(
                [
                    _cell(u),
                    _cell(meth),
                    kind,
                    _cell(k["r"]),
                    _cell(failed),
                    _cell(k["share_matched_median"]),
                    _cell(k["share_matched_min"]),
                    _cell(k["spearman_median"]),
                    _cell(k["spearman_min"]),
                    _cell(k["mean_jaccard_median"]),
                    _cell(k["mean_jaccard_min"]),
                    _cell(k["n_lineages_min"]),
                ]
            )
    header = (
        "Patient",
        "Method",
        "Kind",
        "R",
        "Failed",
        "Share matched median (C3)",
        "Share matched min",
        "Spearman median (C4)",
        "Spearman min",
        "Mean Jaccard median",
        "Mean Jaccard min",
        "Lineages min",
    )
    return _table(header, rows)


def _cell_cycle_table(runs: list[tuple[str, str, dict[str, Any]]]) -> list[str]:
    rows = []
    for u, meth, m in runs:
        cc = m.get("cell_cycle")
        if cc is None:
            rows.append([_cell(u), _cell(meth), MISSING, MISSING, MISSING])
            continue
        warnings = "; ".join(w["message"] for w in cc["warnings"])
        rows.append([_cell(u), _cell(meth), _cell(cc["status"]), _cell(cc["dataset_share"]), _cell(warnings)])
    return _table(("Patient", "Method", "Cell-cycle status", "Dataset share cycling", "Warnings"), rows)


def failed_methods(report: dict[str, Any]) -> list[tuple[str, str, str]]:
    """(patient, method, exception class name) for every unit or pooled method that raised.

    The class name only: the raw message (in the JSON, for diagnosis) may carry paths.
    """
    return [
        (u, meth, str(m["error"]).split(":", 1)[0])
        for pooled in (False, True)
        for u, meth, m in _runs(report, pooled)
        if m.get("error")
    ]


def _failed_table(runs: list[tuple[str, str, dict[str, Any]]]) -> list[str]:
    rows = [
        [_cell(u), _cell(meth), _cell(str(m["error"]).split(":", 1)[0])]
        for u, meth, m in runs
        if m.get("error")
    ]
    if not rows:
        return []
    return ["", "#### Failed methods", "", *_table(("Patient", "Method", "Exception"), rows)]


def _footnote() -> list[str]:
    if constants.VERDICT_VERIFIED:
        return []
    return [
        "",
        f"Verdict: the rule failed its verification on {constants.VERDICT_VERIFIED_ON}, so supported is "
        "never shown; a passed check is shown as provisional. Tree-likeness and both stability checks are "
        "left out (no calibrated threshold), so the verdict reads coverage and clonal support only. A "
        "failed check is shown as weak or not_supported.",
    ]


def _section(heading: str, runs: list[tuple[str, str, dict[str, Any]]], pooled: bool) -> list[str]:
    return [
        f"### {heading}",
        "",
        *_row_table(runs, donors=pooled),
        *([] if pooled else _footnote()),
        "",
        "#### Lineage-pair coupling (Benjamini–Hochberg)",
        "",
        *_coupling_table(runs),
        "",
        "#### Stability",
        "",
        *_stability_table(runs),
        "",
        "#### Cell cycle",
        "",
        *_cell_cycle_table(runs),
        *_failed_table(runs),
    ]


def render_markdown(report: dict[str, Any]) -> str:
    """Tables only, from a ``run_units`` report: per patient first, the pooled diagnostic after."""
    out = _section(PER_UNIT, _runs(report, pooled=False), pooled=False)
    if report.get("pooled_diagnostic") is not None:
        out += ["", *_section(POOLED, _runs(report, pooled=True), pooled=True)]
    return "\n".join(out) + "\n"
