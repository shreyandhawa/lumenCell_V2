"""Benchmark tables for the README, generated from the committed result files."""

import json
import math
from pathlib import Path
from typing import Any

from lumencell.bench.datasets import DATASETS, lineage_gates_note
from lumencell.bench.gates import BASELINE, REPO, check_gates, load_thresholds
from lumencell.trust.constants import VERDICT_VERIFIED

README = REPO / "docs" / "benchmarks.md"  # the generated tables; the name is kept so the CLI is untouched
START = "<!-- bench:start -->"
END = "<!-- bench:end -->"
LEVEL_LABELS = {
    "not_supported": "Not supported: cells left without a pseudotime",
    "weak": "Clones do not back these lineages",
    "provisional": "Not shown to be a tree",
    "supported": "Supported" if VERDICT_VERIFIED else "Not shown to be a tree",
}
LABELS = {"paga_dpt": "PAGA + DPT", "palantir": "Palantir", "scfates_ppt": "scFates"}


def _fmt(v: Any) -> str:
    return "–" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.2f}"


def _floor3(v: float) -> str:
    # Rounded down, so a value under a ">=" bound can never print as the bound (0.4996 is 0.499).
    return f"{math.floor(round(v * 1000, 6)) / 1000:.3f}"


def _range(rows: list[dict[str, Any]], key: str, floor: bool = False) -> str:
    # A row without the key (a baseline from before the metric existed) shows the missing mark.
    cells = [r.get(key, math.nan) for r in rows]
    vals = [v for v in cells if not (isinstance(v, float) and math.isnan(v))]
    if not vals:
        return "–"
    if floor:
        return f"{_floor3(min(vals))}–{_floor3(max(vals))}"
    return f"{min(vals):.2f}–{max(vals):.2f}"


def _shown(v: Any) -> str:
    if isinstance(v, str):
        return v
    return str(v) if isinstance(v, int) and not isinstance(v, bool) else _fmt(v)


def _clonal_entry(r: dict[str, Any]) -> str:
    """The clonal part of a real-table cell; nothing for a dataset without clones."""
    status = r.get("clonal_status", "not_applicable")
    if status == "not_applicable":
        return ""
    if status == "ok":
        return f", clonal: {_fmt(r.get('clonal_excess_ratio'))}"
    if status == "not_assessable":
        return f", clonal: not assessable ({r.get('clonal_reason') or 'no reason recorded'})"
    return f", clonal: {status}"


_CLONAL_REPORTED = (
    ("clonal_n_cells", "cells"),
    ("clonal_effective_clones", "effective clones"),
    ("clonal_eta_m3", "η at m = 3"),
    ("clonal_p_m3", "p at m = 3"),
)


def _reported_tables(data: dict[str, Any], datasets: list[str]) -> list[str]:
    out: list[str] = []
    for d in datasets:
        rows = {m: cases[d] for m, cases in data["methods"].items() if d in cases}
        fields = sorted({k for r in rows.values() for k in r if k.startswith("reported_")})
        if not fields:
            continue
        extra = (
            [c for c in _CLONAL_REPORTED]
            if any(r.get("clonal_status", "not_applicable") != "not_applicable" for r in rows.values())
            else []
        )
        columns = [(k, k.removeprefix("reported_").replace("_", " ")) for k in fields] + extra
        out.append(f"\nReported, not gated ({d}):\n")
        out.append("| Method | " + " | ".join(name for _, name in columns) + " |")
        out.append("|---|" + "---|" * len(columns))
        for method, r in rows.items():
            cells = " | ".join(_shown(r.get(k)) for k, _ in columns)
            out.append(f"| {LABELS.get(method, method)} | {cells} |")
    return out


def latest_real() -> tuple[Path, dict[str, Any]] | None:
    files = sorted((REPO / "benchmarks" / "results").glob("real-*.json"))
    if not files:
        return None
    return files[-1], json.loads(files[-1].read_text(encoding="utf-8"))


def render() -> str:
    cfg = load_thresholds()
    out: list[str] = []
    base = json.loads(BASELINE.read_text(encoding="utf-8"))
    out.append(
        f"**Synthetic** ({len(cfg['synthetic']['topologies'])} topologies × {len(cfg['synthetic']['seeds'])} "
        f"seeds, `realistic` preset; baseline from commit `{base['commit']}`). Ranges over all cases:\n"
    )
    out.append(
        "| Method | Pseudotime (global) | Pseudotime (per lineage) | Branch accuracy | "
        "Clonal excess ratio | Gates failed |"
    )
    out.append("|---|---|---|---|---|---|")
    by_gate: list[str] = []
    for method, cases in base["methods"].items():
        rows = list(cases.values())
        failures = [check_gates(r, cfg["synthetic"]["gates"]) for r in rows]
        fails = sum(bool(f) for f in failures)
        counts: dict[str, int] = {}
        for case_failures in failures:
            for f in case_failures:
                counts[f.metric] = counts.get(f.metric, 0) + 1
        by_gate.append(
            f"{LABELS.get(method, method)}: "
            + (", ".join(f"`{m}` {n}" for m, n in sorted(counts.items())) or "none")
        )
        out.append(
            f"| {LABELS.get(method, method)} | {_range(rows, 'pt_spearman_global')} | "
            f"{_range(rows, 'pt_spearman_lineage')} | {_range(rows, 'branch_accuracy')} | "
            f"{_range(rows, 'clonal_excess_ratio', floor=True)} | "
            f"{fails} of {len(rows)} cases |"
        )
    out.append("\nCases failing each gate (a case can fail several): " + "; ".join(by_gate) + ".")
    out.append(
        "\nThe clonal gate is not a branch gate: a method that merges two true lineages can score a "
        "clonal ratio above 1, and lineage recall catches that."
    )
    real = latest_real()
    if real:
        path, data = real
        out.append(
            f"\n**Real datasets** (`{path.relative_to(REPO)}`, commit `{data['commit']}`). "
            "Spearman is against stage order (embryo time for gastrulation):\n"
        )
        datasets = sorted({d for cases in data["methods"].values() for d in cases})
        out.append("| Method | " + " | ".join(datasets) + " |")
        out.append("|---|" + "---|" * len(datasets))
        for method, cases in data["methods"].items():
            cells = []
            for d in datasets:
                r = cases.get(d)
                if r is None:
                    cells.append("–")
                    continue
                failed = check_gates(r, cfg["real"][d]["gates"])
                note = lineage_gates_note(DATASETS[d]) if d in DATASETS else None
                mark = "fail: " + ", ".join(f.metric for f in failed) if failed else (note or "pass")
                level = r.get("tree_support", "not_measured")
                shown = "" if level == "not_measured" else f", verdict {LEVEL_LABELS.get(level, level)}"
                cells.append(
                    f"Spearman {_fmt(r['pt_spearman_lineage'])}, branch {_fmt(r['branch_accuracy'])}, "
                    f"terminal {_fmt(r['terminal_share'])} ({mark}){shown}{_clonal_entry(r)}"
                )
            out.append(f"| {LABELS.get(method, method)} | " + " | ".join(cells) + " |")
        out.extend(_reported_tables(data, datasets))
        if any(
            r.get("tree_support", "not_measured") != "not_measured"
            for cases in data["methods"].values()
            for r in cases.values()
        ):
            out.append(
                "\nA supported verdict is capped at provisional, because the verdict's verification "
                "failed on fresh seeds (2026-11-17): at x1.0 the tree-likeness check (C2) separates blob "
                "and sheet from the matrix in the recorded calibration but not loop or two_groups, and the "
                "stability checks do not separate trees from the controls. A weak or not_supported "
                "verdict is a real failure and is not capped. The verdict is reported, not gated, "
                "and is not evidence of a tree."
            )
    return "\n".join(out) + "\n"


CASE_STUDY_START = "<!-- case-study:start -->"
CASE_STUDY_END = "<!-- case-study:end -->"


def replace_section(path: Path, text: str, start: str, end: str) -> None:
    """Replace what lies between the ``start`` and ``end`` markers of ``path`` with ``text``."""
    content = path.read_text(encoding="utf-8")
    if start not in content or end not in content:
        raise ValueError(f"{path} has no {start} ... {end} section")
    head, rest = content.split(start, 1)
    _old, tail = rest.split(end, 1)
    path.write_text(f"{head}{start}\n{text}{end}{tail}", encoding="utf-8")


def update_readme(text: str, path: Path = README) -> None:
    replace_section(path, text, START, END)
