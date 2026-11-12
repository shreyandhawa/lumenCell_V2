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
HIGHER_IS_BETTER = (
    "pt_spearman_global",
    "pt_spearman_lineage",
    "branch_accuracy",
    "lineage_recall",
    "clonal_excess_ratio",
)


@dataclass(frozen=True)
class GateFailure:
    case: str
    metric: str
    value: float
    bound: str
    reason: str | None = None

    def __str__(self) -> str:
        text = f"{self.case}: {self.metric} = {self.value:.3f} (needs {self.bound})"
        return f"{text} [{self.reason}]" if self.reason else text


def load_thresholds(path: Path = THRESHOLDS) -> dict[str, Any]:
    with path.open("rb") as fh:
        return tomllib.load(fh)


_GATE_KEYS = frozenset({"min", "max", "nan", "status"})


def _bound_text(bound: dict[str, Any]) -> str:
    parts = [f">= {bound['min']}"] if "min" in bound else []
    return " and ".join(parts + ([f"<= {bound['max']}"] if "max" in bound else []))


def check_gates(row: dict[str, Any], gates: dict[str, dict[str, Any]]) -> list[GateFailure]:
    """Gates that fail on this row.

    A gate has ``min`` and/or ``max`` and may set ``nan`` ("skip", the default, or "fail")
    and ``status``, the name of a row field. Any other key, no bound, or an unknown ``nan``
    raises ValueError, for every gate. A gate with ``status`` must say its ``nan`` policy
    explicitly: a status gate that silently skipped NaN would pass on a missing metric.
    With ``status``: "not_applicable" skips the gate, anything but "ok" fails it with the
    status as the reason (a missing field counts as "not_measured", a missing metric as
    NaN), and a NaN under ``nan = "fail"`` fails with reason "nan". Without ``status`` a
    missing metric raises KeyError and NaN follows ``nan``.
    """
    failures: list[GateFailure] = []
    for metric, bound in gates.items():
        for key in bound.keys() - _GATE_KEYS:
            raise ValueError(f"gate {metric}: unknown key {key!r}")
        if "min" not in bound and "max" not in bound:
            raise ValueError(f"gate {metric}: needs min or max")
        status_key = bound.get("status")
        if status_key is not None and "nan" not in bound:
            raise ValueError(f"gate {metric}: a gate with status must set nan explicitly")
        policy = bound.get("nan", "skip")
        if policy not in ("skip", "fail"):
            raise ValueError(f"gate {metric}: nan must be 'skip' or 'fail', not {policy!r}")
        text = _bound_text(bound)
        if status_key is None:
            value = float(row[metric])
        else:
            status = row.get(status_key, "not_measured")
            if status == "not_applicable":
                continue
            value = float(row.get(metric, math.nan))
            if status != "ok":
                failures.append(GateFailure(row["case"], metric, value, text, reason=str(status)))
                continue
        if math.isnan(value):
            if policy == "fail":
                failures.append(GateFailure(row["case"], metric, value, text, reason="nan"))
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
