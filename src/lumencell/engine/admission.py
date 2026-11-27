"""Run admission: a request for a run, turned into a canonical run spec and its run id, or a refusal.

The API and the CLI admit every run they start, so a request gets the same id and the same
refusal whichever way it arrives. The infer job runs the spec the API admitted.
"""

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from pydantic import ValidationError

from lumencell.clonal.isolation import scoring_only_message
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import UnknownMethodError, get_method
from lumencell.engine.roots import ClusterRoot, RootError, RootSpec, format_root, parse_root, resolve_root
from lumencell.engine.store import RunSpec, run_id

ROOT_REQUIRED = (
    "Pick where the trajectory starts before running a method. "
    "lumenCell never guesses the root, because a wrong guess silently reverses pseudotime."
)


class RunRefused(Exception):
    """A request that cannot become a run. ``code`` is stable; ``details`` is what the code carries."""

    def __init__(self, code: str, message: str, details: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details


@dataclass(frozen=True)
class AdmittedRun:
    spec: RunSpec
    run_id: str


def admit_run(
    method: str,
    params: dict[str, Any] | None,
    root_text: str | None,
    terminals: list[str] | None,
    context: Callable[[], InferenceContext],
    prep_key: str,
) -> AdmittedRun:
    """The run a request asks for, or :class:`RunRefused` for its first problem.

    Checked in this order: ``unknown_method``, ``method_unavailable``, ``root_required``,
    ``invalid_params``, ``invalid_root`` (the text does not parse), ``scoring_only_column``
    (naming every such column), then ``invalid_root`` (the root does not resolve). ``context``
    is called only once the request itself is sound, so a dataset still being prepared is
    reported after the request's own problems.
    The spec holds the root in canonical text, the params with defaults filled in and the
    context's clonal columns, so the id is the one its run record gives back.
    """
    try:
        impl = get_method(method)
    except UnknownMethodError as exc:
        raise RunRefused("unknown_method", str(exc.args[0])) from exc
    availability = impl.available()
    if not availability.ok:
        hint = availability.install_hint or ""
        raise RunRefused("method_unavailable", f"{impl.label} is not installed. {hint}".strip())
    if not root_text:
        raise RunRefused("root_required", ROOT_REQUIRED)
    try:
        filled = impl.Params.model_validate(params or {}).model_dump()
    except ValidationError as exc:
        errors = [
            {"loc": [str(part) for part in e["loc"]], "msg": str(e["msg"])}
            for e in exc.errors(include_url=False, include_context=False, include_input=False)
        ]
        raise RunRefused(
            "invalid_params", f"The parameters for {impl.label} are not valid.", {"errors": errors}
        ) from exc
    try:
        root = parse_root(root_text)
    except RootError as exc:
        raise RunRefused("invalid_root", str(exc)) from exc
    ctx = context()
    _refuse_scoring_only(ctx, root, filled)
    try:
        resolve_root(ctx, root)
    except RootError as exc:
        raise RunRefused("invalid_root", str(exc)) from exc
    spec = RunSpec(method, filled, format_root(root), terminals, ctx.clonal_columns)
    return AdmittedRun(spec, run_id(spec, prep_key))


def _refuse_scoring_only(ctx: InferenceContext, root: RootSpec, params: dict[str, Any]) -> None:
    """Checked before the root is resolved, so one refusal names every such column."""
    asked: list[tuple[str | None, str]] = [
        (root.key if isinstance(root, ClusterRoot) else None, "root column"),
        (params.get("cluster_key"), "cluster column"),
    ]
    refused = [(column, role) for column, role in asked if column is not None and column in ctx.scoring_only]
    if refused:
        column, role = refused[0]
        raise RunRefused(
            "scoring_only_column",
            scoring_only_message(column, role),
            {"columns": list(dict.fromkeys(c for c, _ in refused))},
        )
