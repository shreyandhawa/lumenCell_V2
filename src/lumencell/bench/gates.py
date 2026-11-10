"""Evaluate benchmark results against benchmarks/thresholds.toml and the baseline."""

import json
import math
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[3]
THRESHOLDS = REPO / "benchmarks" / "thresholds.toml"
BASELINE = REPO / "benchmarks" / "baseline.json"

# Metrics where higher is better; the regression guard checks drops on these.
HIGHER_IS_BETTER = ("pt_spearman_global", "pt_spearman_lineage", "branch_accuracy", "lineage_recall")


@dataclass(frozen=True)
class GateFailure:
    case: str
    metric: str
    value: float
    bound: str

    def __str__(self) -> str:
        return f"{self.case}: {self.metric} = {self.value:.3f} (needs {self.bound})"


def load_thresholds(path: Path = THRESHOLDS) -> dict[str, Any]:
    with path.open("rb") as fh:
        return tomllib.load(fh)


def check_gates(row: dict[str, Any], gates: dict[str, dict[str, float]]) -> list[GateFailure]:
    failures: list[GateFailure] = []
    for metric, bound in gates.items():
        value = float(row[metric])
        if math.isnan(value):
            continue
        if "min" in bound and value < bound["min"]:
            failures.append(GateFailure(row["case"], metric, value, f">= {bound['min']}"))
        if "max" in bound and value > bound["max"]:
            failures.append(GateFailure(row["case"], metric, value, f"<= {bound['max']}"))
    return failures


def load_baseline(path: Path = BASELINE) -> dict[str, dict[str, dict[str, float]]]:
    """{method: {case: {metric: value}}}"""
    if not path.exists():
        return {}
    data: dict[str, dict[str, dict[str, float]]] = json.loads(path.read_text(encoding="utf-8"))["methods"]
    return data


def check_regression(row: dict[str, Any], baseline: dict[str, float], max_drop: float) -> list[GateFailure]:
    failures: list[GateFailure] = []
    for metric in HIGHER_IS_BETTER:
        base, value = baseline.get(metric), float(row[metric])
        if base is None or math.isnan(float(base)) or math.isnan(value):
            continue
        if value < float(base) - max_drop:
            failures.append(GateFailure(row["case"], metric, value, f">= baseline {base:.3f} - {max_drop}"))
    return failures
