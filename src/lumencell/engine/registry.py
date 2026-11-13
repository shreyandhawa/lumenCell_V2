"""The method registry and the single entry point for running inference."""

import time
from typing import Any

from lumencell.engine.context import InferenceContext, check_inference_column
from lumencell.engine.lineages import divergence, name_lineages
from lumencell.engine.methods.base import ProgressFn, TrajectoryMethod, no_progress
from lumencell.engine.methods.paga_dpt import PagaDpt
from lumencell.engine.methods.palantir_adapter import Palantir
from lumencell.engine.methods.scfates_adapter import ScFates
from lumencell.engine.roots import RootSpec, interior_root_warning, resolve_root
from lumencell.engine.schema import TrajectoryResult

METHODS: dict[str, TrajectoryMethod] = {m.name: m for m in (PagaDpt(), Palantir(), ScFates())}


class UnknownMethodError(KeyError):
    pass


def get_method(name: str) -> TrajectoryMethod:
    try:
        return METHODS[name]
    except KeyError:
        raise UnknownMethodError(f"unknown method {name!r}; available: {', '.join(METHODS)}") from None


def run_method(
    ctx: InferenceContext,
    method: str,
    root: RootSpec,
    terminals: list[str] | None = None,
    params: dict[str, Any] | None = None,
    progress: ProgressFn = no_progress,
    label_key: str | None = None,
) -> TrajectoryResult:
    """Resolve the root, run a method, add method-independent post-processing, validate.

    There is no default root: ``root`` is required, so a pseudotime can never be
    computed from a guessed start. A scoring-only column (``ctx.scoring_only``) is
    refused as the cluster column, the label column or the root before any work.
    """
    impl = get_method(method)
    availability = impl.available()
    if not availability.ok:
        raise RuntimeError(f"{impl.label} is not installed. {availability.install_hint or ''}".strip())
    parsed = impl.Params.model_validate(params or {})
    cluster_key: str | None = getattr(parsed, "cluster_key", None)
    # When you name terminals in a column, lineages are named in that column too;
    # otherwise in the detected label column.
    used_label_key = label_key or (cluster_key if terminals else None) or ctx.label_key
    if cluster_key:
        check_inference_column(ctx, cluster_key, "cluster column")
    check_inference_column(ctx, used_label_key, "label column")
    resolved = resolve_root(ctx, root)
    start = time.perf_counter()
    result = impl.run(ctx, resolved, terminals, parsed, progress)

    name_lineages(ctx, result, used_label_key)
    result.divergence = divergence(result.membership, result.pseudotime)
    result.warnings = [*resolved.warnings, *result.warnings]
    interior = interior_root_warning(ctx, resolved, result.pseudotime)
    if interior:
        result.warnings.append(interior)
    result.runtime_seconds = time.perf_counter() - start
    result.obs_keys = sorted({k for k in (resolved.cluster_key, cluster_key, used_label_key) if k})
    result.validate()
    return result
