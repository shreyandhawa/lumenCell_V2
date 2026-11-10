"""Run every preprocessing stage in memory (no cache).

The platform track wraps the same stages with an on-disk cache and a background
job; this function is the reference behaviour they must reproduce.
"""

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

from lumencell.data.loader import Dataset
from lumencell.prep.stages import STAGES, PrepParams, StageRecord, run_stage

ProgressFn = Callable[[str, int, int], None]


@dataclass(frozen=True)
class PrepResult:
    artifacts: dict[str, Any]
    stages: list[StageRecord] = field(default_factory=list)

    @property
    def display_embeddings(self) -> list[str]:
        return ["umap"] if "umap" in self.artifacts else []


def run_pipeline(
    ds: Dataset, params: PrepParams | None = None, progress: ProgressFn | None = None
) -> PrepResult:
    params = params or PrepParams()
    artifacts: dict[str, Any] = {}
    records: list[StageRecord] = []
    for i, stage in enumerate(STAGES):
        if progress:
            progress(stage.name, i, len(STAGES))
        if not stage.applies(ds, params):
            records.append(StageRecord(stage.name, "skipped", 0.0))
            continue
        out, seconds = run_stage(stage, ds, artifacts, params)
        artifacts.update(out)
        records.append(StageRecord(stage.name, "done", seconds))
    if progress:
        progress("done", len(STAGES), len(STAGES))
    return PrepResult(artifacts=artifacts, stages=records)
