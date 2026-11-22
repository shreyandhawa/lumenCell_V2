"""Writes the reference rows the trust panel shows, from the committed benchmark results.

Run ``uv run python -m lumencell.bench.reference`` after a new result or a threshold change.
``tests/unit/test_reference_scores.py`` fails while the packaged file differs from what this
module builds.

Each row is one gate on one dataset: its worst or measured score, its threshold from
``thresholds.toml`` and whether ``check_gates`` passed it. Synthetic rows cover all cases of the
baseline (the worst value, and how many cases fail). A gate that was not measured is shown as
not measured, never as a pass.
"""

import json
import math
from pathlib import Path
from typing import Any

from lumencell.bench.datasets import DATASETS, lineage_gates_note
from lumencell.bench.gates import BASELINE, check_gates, load_thresholds
from lumencell.bench.report import latest_real

PACKAGED = Path(__file__).with_name("reference_scores.json")

METRIC_LABELS = {
    "pt_spearman_global": "Pseudotime, global",
    "pt_spearman_lineage": "Pseudotime, per lineage",
    "branch_accuracy": "Branch accuracy",
    "lineage_recall": "Lineage recall",
    "terminal_share": "Terminal share",
    "n_spurious": "Spurious branches",
    "root_pt_excess": "Root pseudotime excess",
    "clonal_excess_ratio": "Clonal excess ratio",
    "clonal_p": "Clonal p-value",
}
REAL_LABELS = {
    "pancreas": "Pancreas (mouse, E15.5)",
    "bonemarrow": "Bone marrow (human CD34+)",
    "gastrulation_erythroid": "Gastrulation (mouse, erythroid)",
    "daniel2022": "LCMV CD8+ T cells (mouse)",
}


def _threshold(bound: dict[str, Any]) -> str:
    parts = ([f"≥ {bound['min']}"] if "min" in bound else []) + (
        [f"≤ {bound['max']}"] if "max" in bound else []
    )
    return " and ".join(parts)


def _skipped(case: dict[str, Any], gate: str, bound: dict[str, Any]) -> bool:
    """True when ``check_gates`` leaves this case out: not applicable, or NaN under the default policy."""
    status = bound.get("status")
    if status is not None:
        state = case.get(status, "not_measured")
        if state == "not_applicable":
            return True
        if state != "ok":
            return False
    return math.isnan(float(case[gate])) and bound.get("nan", "skip") == "skip"


def _row(dataset: str, gate: str, bound: dict[str, Any], source: str, **result: Any) -> dict[str, Any]:
    return {
        "dataset": dataset,
        "metric": METRIC_LABELS[gate],
        "gate": gate,
        "score": result["score"],
        "threshold": _threshold(bound),
        "passed": result["passed"],
        "detail": result["detail"],
        "source": source,
    }


def _synthetic_rows(
    cases: dict[str, dict[str, Any]], gates: dict[str, Any], source: str
) -> list[dict[str, Any]]:
    n = len(cases)
    rows = []
    for gate, bound in gates.items():
        measured = [c for c in cases.values() if not _skipped(c, gate, bound)]
        failed = [c for c in measured if check_gates(c, {gate: bound})]
        values = [float(c[gate]) for c in measured if not math.isnan(float(c[gate]))]
        if not measured:
            result: dict[str, Any] = {"score": None, "passed": None, "detail": "not measured on any case"}
        else:
            score = (min(values) if "min" in bound else max(values)) if values else None
            if failed:
                detail = f"worst case; fails {len(failed)} of {len(measured)} cases"
            else:
                detail = f"worst case; passes all {len(measured)} cases"
            result = {"score": score, "passed": not failed, "detail": detail}
        rows.append(_row(f"Synthetic, {n} cases", gate, bound, source, **result))
    return rows


def _real_rows(
    dataset: str, case: dict[str, Any], gates: dict[str, Any], source: str
) -> list[dict[str, Any]]:
    note = lineage_gates_note(DATASETS[dataset])
    rows = []
    for gate, bound in gates.items():
        value = float(case[gate])
        if math.isnan(value):
            result: dict[str, Any] = {"score": None, "passed": None, "detail": note or "not measured"}
        else:
            result = {"score": value, "passed": not check_gates(case, {gate: bound}), "detail": None}
        rows.append(_row(REAL_LABELS[dataset], gate, bound, source, **result))
    return rows


def build() -> dict[str, Any]:
    cfg = load_thresholds()
    base = json.loads(BASELINE.read_text(encoding="utf-8"))
    found = latest_real()
    if found is None:
        raise RuntimeError("no benchmarks/results/real-*.json to build the reference rows from")
    real_path, real = found
    synthetic_source = f"benchmarks/baseline.json @ {base['commit']}"
    real_source = f"benchmarks/results/{real_path.name} @ {real['commit']}"
    methods: dict[str, list[dict[str, Any]]] = {}
    for method in dict.fromkeys([*base["methods"], *real["methods"]]):
        rows: list[dict[str, Any]] = []
        if method in base["methods"]:
            rows += _synthetic_rows(base["methods"][method], cfg["synthetic"]["gates"], synthetic_source)
        for dataset, case in real["methods"].get(method, {}).items():
            rows += _real_rows(dataset, case, cfg["real"][dataset]["gates"], real_source)
        methods[method] = rows
    return {"methods": methods}


def render() -> str:
    return json.dumps(build(), indent=2, ensure_ascii=False) + "\n"


def main() -> None:
    PACKAGED.write_text(render(), encoding="utf-8")
    print(f"wrote {PACKAGED}")


if __name__ == "__main__":
    main()
