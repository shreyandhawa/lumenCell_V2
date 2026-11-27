"""The method registry and the single entry point for running inference."""

import time
from typing import TYPE_CHECKING, Any

from lumencell.data.loader import Dataset
from lumencell.engine.context import InferenceContext, check_inference_column
from lumencell.engine.lineages import divergence, name_lineages
from lumencell.engine.methods.base import ProgressFn, TrajectoryMethod, no_progress
from lumencell.engine.methods.paga_dpt import PagaDpt
from lumencell.engine.methods.palantir_adapter import Palantir
from lumencell.engine.methods.scfates_adapter import ScFates
from lumencell.engine.roots import RootSpec, interior_root_warning, parse_root, resolve_root
from lumencell.engine.schema import TrajectoryResult
from lumencell.engine.store import RunStore, run_id
from lumencell.prep.cache import prep_key
from lumencell.prep.pipeline import PrepResult
from lumencell.prep.stages import PrepParams

if TYPE_CHECKING:
    from lumencell.engine.admission import AdmittedRun  # admission imports this module

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

    Lineages are named in the detected label column (``ctx.label_key``), also when terminals
    are named in another column. ``label_key`` overrides it only to replay a stored run.
    """
    impl = get_method(method)
    availability = impl.available()
    if not availability.ok:
        raise RuntimeError(f"{impl.label} is not installed. {availability.install_hint or ''}".strip())
    parsed = impl.Params.model_validate(params or {})
    cluster_key: str | None = getattr(parsed, "cluster_key", None)
    used_label_key = label_key or ctx.label_key
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
    result.label_key = used_label_key
    result.obs_keys = sorted({k for k in (resolved.cluster_key, cluster_key, used_label_key) if k})
    result.validate()
    return result


def infer_and_save(
    ds: Dataset,
    prep: PrepResult,
    prep_params: PrepParams,
    admitted: "AdmittedRun",
    store: RunStore,
    progress: ProgressFn = no_progress,
) -> TrajectoryResult:
    """Run an admitted run on a prepared dataset and save it under its run id.

    The id is first recomputed from the spec and this prep: a spec admitted by another lumencell
    version, or on another prep, fails before any work. The inference context and the run record
    come from the same spec.
    """
    spec = admitted.spec
    key = prep_key(ds, prep_params)
    rid = run_id(spec, key)
    if rid != admitted.run_id:
        raise ValueError(f"the run spec gives run id {rid}, not {admitted.run_id}; nothing was run")
    ctx = InferenceContext.from_prep(ds, prep, prep_params, clonal_columns=spec.clonal_columns)
    result = run_method(ctx, spec.method, parse_root(spec.root), spec.terminals, spec.params, progress)
    store.save(rid, result, spec, key)
    return result
