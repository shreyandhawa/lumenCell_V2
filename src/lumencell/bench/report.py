"""Benchmark tables for the README, generated from the committed result files."""

import json
import math
from pathlib import Path
from typing import Any

from lumencell.bench.gates import BASELINE, REPO, check_gates, load_thresholds

README = REPO / "README.md"
START = "<!-- bench:start -->"
END = "<!-- bench:end -->"
LABELS = {"paga_dpt": "PAGA + DPT", "palantir": "Palantir", "scfates_ppt": "scFates"}


def _fmt(v: Any) -> str:
    return "–" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.2f}"


def _range(rows: list[dict[str, Any]], key: str) -> str:
    # A row without the key (a baseline from before the metric existed) shows the missing mark.
    cells = [r.get(key, math.nan) for r in rows]
    vals = [v for v in cells if not (isinstance(v, float) and math.isnan(v))]
    return "–" if not vals else f"{min(vals):.2f}–{max(vals):.2f}"


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
    for method, cases in base["methods"].items():
        rows = list(cases.values())
        fails = sum(bool(check_gates(r, cfg["synthetic"]["gates"])) for r in rows)
        out.append(
            f"| {LABELS.get(method, method)} | {_range(rows, 'pt_spearman_global')} | "
            f"{_range(rows, 'pt_spearman_lineage')} | {_range(rows, 'branch_accuracy')} | "
            f"{_range(rows, 'clonal_excess_ratio')} | "
            f"{fails} of {len(rows)} cases |"
        )
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
                mark = "pass" if not failed else "fail: " + ", ".join(f.metric for f in failed)
                cells.append(
                    f"Spearman {_fmt(r['pt_spearman_lineage'])}, branch {_fmt(r['branch_accuracy'])}, "
                    f"terminal {_fmt(r['terminal_share'])} ({mark})"
                )
            out.append(f"| {LABELS.get(method, method)} | " + " | ".join(cells) + " |")
    return "\n".join(out) + "\n"


def update_readme(text: str, path: Path = README) -> None:
    content = path.read_text(encoding="utf-8")
    if START not in content or END not in content:
        raise ValueError(f"{path} has no {START} ... {END} section")
    head, rest = content.split(START, 1)
    _old, tail = rest.split(END, 1)
    path.write_text(f"{head}{START}\n{text}{END}{tail}", encoding="utf-8")
