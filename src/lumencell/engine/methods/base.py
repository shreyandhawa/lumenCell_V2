"""The interface every method adapter implements."""

from collections.abc import Callable
from dataclasses import dataclass
from typing import ClassVar, Protocol

from pydantic import BaseModel, ConfigDict

from lumencell.engine.context import InferenceContext
from lumencell.engine.schema import ResolvedRoot, TrajectoryResult

ProgressFn = Callable[[str, int, int], None]


def no_progress(_stage: str, _i: int, _n: int) -> None:
    return None


class MethodParams(BaseModel):
    model_config = ConfigDict(extra="forbid")


@dataclass(frozen=True)
class Availability:
    ok: bool
    install_hint: str | None = None
    version: str | None = None


class TrajectoryMethod(Protocol):
    name: ClassVar[str]
    label: ClassVar[str]
    family: ClassVar[str]
    description: ClassVar[str]
    Params: ClassVar[type[MethodParams]]
    # True while the method misses any gate in benchmarks/thresholds.toml.
    # tests/golden checks that this matches the measured results.
    experimental: ClassVar[bool]
    # One sentence on which gates it misses, shown in the UI when experimental.
    experimental_note: ClassVar[str | None]

    def available(self) -> Availability: ...

    def run(
        self,
        ctx: InferenceContext,
        root: ResolvedRoot,
        terminals: list[str] | None,
        params: MethodParams,
        progress: ProgressFn,
    ) -> TrajectoryResult: ...
