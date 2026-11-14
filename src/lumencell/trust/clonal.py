"""Clonal settings for a dataset, and the clonal agreement of a run.

The settings only say which obs columns carry clone labels, donors and strata. Candidates are
suggested by name and never applied; the user names the columns.

The agreement (``agreement_for_run``) scores a finished run's lineages against the clones in
own-cells mode and maps the result to the API models, NaN to ``None``. The score_clonal worker job
stores it as a per-run trust file named by ``trust_name``, a hash: column names never reach a path.
Clone sharing says two lineages are related, never which came first.
"""

import hashlib
import json
import math
import re
import sys
from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from pydantic import ValidationError

from lumencell.api.models import (
    ClonalAgreement,
    ClonalScore,
    ClonalSettings,
    ClonalStatus,
    ClonalSummary,
    LineageCoupling,
    NullComparison,
)
from lumencell.clonal.isolation import scoring_only_columns
from lumencell.clonal.labels import encode_clones, is_missing
from lumencell.clonal.scores import ClonalResult, Coupling, score_with_baseline
from lumencell.engine.schema import TrajectoryResult
from lumencell.trust.constants import (
    CLONAL_METRIC_VERSION,
    CLONAL_SEED,
    MAX_DONORS,
    MIN_CLONE_SIZE,
    MIN_CLONES,
    N_PERMUTATIONS,
    PSEUDOTIME_BINS,
    SENSITIVITY_CLONE_SIZE,
)

CANDIDATE = re.compile(r"^(?:airr:)?(?:clone_id|raw_clonotype_id|clonotype)$", re.IGNORECASE)


def clonal_candidates(columns: Iterable[object]) -> list[str]:
    """Columns whose name looks like a clone column, in obs order. Suggestions only."""
    return [str(c) for c in columns if CANDIDATE.match(str(c))]


def settings_columns(settings: ClonalSettings) -> tuple[str, ...]:
    """The clone and donor columns, which become scoring-only. Never the strata keys."""
    return (settings.clone_key,) if settings.donor_key is None else (settings.clone_key, settings.donor_key)


def _is_label_column(s: pd.Series) -> bool:
    if isinstance(s.dtype, pd.CategoricalDtype | pd.StringDtype):
        return True
    return s.dtype == object and pd.api.types.infer_dtype(s, skipna=True) in {"string", "empty"}


def _present(s: pd.Series) -> npt.NDArray[np.object_]:
    """Values as an object array: a categorical gives its values, never its unused categories."""
    out: npt.NDArray[np.object_] = s.to_numpy(dtype=object)
    return out


def check_settings(obs: pd.DataFrame, settings: ClonalSettings) -> list[str]:
    """Problems with the settings against ``obs``; empty when they are valid."""
    problems: list[str] = []
    named = [settings.clone_key, *([settings.donor_key] if settings.donor_key else []), *settings.strata_keys]
    missing = [c for c in dict.fromkeys(named) if c not in obs.columns]
    problems += [f"there is no obs column named {c!r}" for c in missing]
    if settings.donor_key is not None and settings.donor_key == settings.clone_key:
        problems.append("the clone column and the donor column must be different columns")
    if len(set(settings.strata_keys)) != len(settings.strata_keys):
        problems.append("the strata columns must be distinct")
    for key in dict.fromkeys(settings.strata_keys):
        if key == settings.clone_key or key == settings.donor_key:
            problems.append(f"strata column {key!r} is also the clone or donor column")

    clone_ok = False
    if settings.clone_key in obs.columns:
        clone = obs[settings.clone_key]
        if not _is_label_column(clone):
            problems.append(
                f"the clone column {settings.clone_key!r} must hold text labels "
                f"(categorical or string), not {clone.dtype}"
            )
        elif bool(is_missing(_present(clone)).all()):
            problems.append(f"the clone column {settings.clone_key!r} has no clone label in any cell")
        else:
            clone_ok = True
    donors_ok = settings.donor_key is None
    if settings.donor_key is not None and settings.donor_key in obs.columns:
        donors = _present(obs[settings.donor_key])
        n = len({str(v) for v in donors[~is_missing(donors)]})
        if not 1 <= n <= MAX_DONORS:
            problems.append(
                f"the donor column {settings.donor_key!r} has {n} distinct values; "
                f"it needs between 1 and {MAX_DONORS}"
            )
        else:
            donors_ok = True
    if clone_ok and donors_ok and settings.donor_key is not None:
        labelled = ~is_missing(_present(obs[settings.clone_key]))
        if not bool((labelled & ~is_missing(_present(obs[settings.donor_key]))).any()):
            problems.append("no cell has both a clone label and a donor")
    return problems


def settings_summary(obs: pd.DataFrame, settings: ClonalSettings) -> ClonalSummary:
    """Clone sizes over every labelled cell with a donor, singletons included."""
    labels = _present(obs[settings.clone_key])
    donors = None if settings.donor_key is None else _present(obs[settings.donor_key])
    codes = encode_clones(labels, donors)
    sizes = np.bincount(codes.clone[codes.clone >= 0])
    present_donors = np.unique(codes.donor[codes.donor >= 0])

    if not sizes.size:
        # check_settings refuses settings with no labelled cell that has a donor.
        raise ValueError(
            "no cell has both a clone label and a donor; check_settings should have refused this"
        )

    return ClonalSummary(
        settings=settings,
        n_cells_labelled=int(sizes.sum()),
        n_clones=int(sizes.size),
        n_clones_2plus=int((sizes >= 2).sum()),
        n_donors=int(present_donors.size),
        size_median=float(np.median(sizes)),
        size_p90=float(np.percentile(sizes, 90)),
        size_p99=float(np.percentile(sizes, 99)),
        size_max=int(sizes.max()),
    )


def settings_from_state(state: dict[str, Any]) -> ClonalSettings | None:
    """The saved settings; ``None`` when there are none or the saved value is invalid."""
    if "clonal" not in state or state["clonal"] is None:
        return None
    try:
        return ClonalSettings.model_validate(state["clonal"])
    except ValidationError as exc:
        print(
            f"lumencell: ignoring the saved clonal settings, which are not valid: {exc.errors()[0]['msg']}",
            file=sys.stderr,
        )
        return None


# ---- the agreement of a run -----------------------------------------------------------

COUPLING_NOTE = (
    "Coupling is the cosine similarity of raw scored-cell counts per clone, not CoSpar's normalised "
    "default. A high z says two lineages share more clones than chance; it does not say which came first."
)


def trust_hash(settings: ClonalSettings) -> str:
    """16 hex characters over the settings, the metric version and the scoring parameters."""
    dumped = settings.model_dump(mode="json")
    dumped["strata_keys"] = sorted(dumped["strata_keys"])
    params = {
        "min_clone_size": MIN_CLONE_SIZE,
        "sensitivity_clone_size": SENSITIVITY_CLONE_SIZE,
        "min_clones": MIN_CLONES,
        "n_permutations": N_PERMUTATIONS,
        "seed": CLONAL_SEED,
        "pseudotime_bins": PSEUDOTIME_BINS,
    }
    text = json.dumps(
        {"settings": dumped, "metric_version": CLONAL_METRIC_VERSION, "params": params}, sort_keys=True
    )
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def trust_name(settings: ClonalSettings) -> str:
    """The trust file name for these settings; it never holds a column name."""
    return f"clonal-{trust_hash(settings)}"


def _num(x: float) -> float | None:
    """A finite value as a float; NaN and infinity become ``None``, so the file stays strict JSON."""
    return float(x) if math.isfinite(x) else None


def to_api_score(r: ClonalResult) -> ClonalScore:
    delta_eta = _num(r.delta_eta)
    return ClonalScore(
        status=r.status,
        mode=r.mode,
        n_lineages=r.n_lineages,
        min_clone_size=r.min_clone_size,
        n_cells=r.n_cells,
        n_clones=r.n_clones,
        effective_clones=_num(r.effective_clones),
        n_donors=r.n_donors,
        n_strata=r.n_strata,
        n_permutations=r.n_permutations,
        seed=r.seed,
        purity=_num(r.purity),
        mean_entropy=_num(r.mean_entropy),
        nulls=[
            NullComparison(
                null=n.null,
                mean_entropy_mean=_num(n.mean_entropy_mean),
                mean_entropy_sd=_num(n.mean_entropy_sd),
                eta=_num(n.eta),
                p=_num(n.p),
            )
            for n in r.nulls
        ],
        eta=_num(r.eta),
        p=_num(r.p),
        eta_base=_num(r.eta_base),
        delta_eta=delta_eta,
        # A p for an undefined difference would be a stand-in number.
        p_delta=None if delta_eta is None else _num(r.p_delta),
        clone_weighted_purity=_num(r.clone_weighted_purity),
        purity_without_largest=_num(r.purity_without_largest),
    )


def to_api_coupling(c: Coupling) -> LineageCoupling:
    return LineageCoupling(
        a=c.a,
        b=c.b,
        observed=_num(c.observed),
        null_mean=_num(c.null_mean),
        z=_num(c.z),
        p=_num(c.p),
        p_adj=_num(c.p_adj),
    )


def not_scored(status: ClonalStatus, settings: ClonalSettings | None, note: str) -> ClonalAgreement:
    """An agreement with no score: the status and the one note that says why."""
    return ClonalAgreement(
        status=status, settings=settings, primary=None, sensitivity_m3=None, coupling=[], notes=[note]
    )


def agreement_for_run(
    result: TrajectoryResult,
    obs: pd.DataFrame,
    settings: ClonalSettings,
    pca: npt.NDArray[np.floating[Any]],
    neighbors_pcs: int,
    obs_keys: Sequence[str],
) -> ClonalAgreement:
    """Own-cells clonal agreement of ``result`` with the clones named by ``settings``.

    ``obs_keys`` are the obs columns that steered the run; a run steered by a scoring-only column
    is refused, and settings that do not fit ``obs`` (``check_settings``) fail with the problems in
    the note. Scores at m = 2 (primary) and m = 3 (sensitivity), each against a k-means baseline
    with k = the number of lineages, on ``baseline_space(pca, neighbors_pcs)``.
    """
    problems = check_settings(obs, settings)
    if problems:
        return not_scored(
            "failed", settings, f"The clonal settings do not fit this file: {'; '.join(problems)}."
        )
    scoring_only = scoring_only_columns(obs.columns, settings_columns(settings))
    steered = [k for k in dict.fromkeys(obs_keys) if k in scoring_only]
    if steered:
        named = ", ".join(repr(k) for k in steered)
        return not_scored(
            "refused",
            settings,
            f"The run used {named}, which lumenCell keeps only for scoring; run it again without it.",
        )
    n_lineages = int(result.membership.shape[0])
    if n_lineages < 2:
        return not_scored(
            "single_lineage", settings, "The run has one lineage; clones are scored over two or more."
        )

    codes = encode_clones(
        _present(obs[settings.clone_key]),
        None if settings.donor_key is None else _present(obs[settings.donor_key]),
        [_present(obs[k]) for k in settings.strata_keys],
    )

    def score(min_clone_size: int) -> ClonalResult:
        return score_with_baseline(
            result.membership,
            codes,
            result.pseudotime,
            pca,
            neighbors_pcs,
            n_lineages,
            min_clone_size=min_clone_size,
            seed=CLONAL_SEED,
        )

    primary, m3 = score(MIN_CLONE_SIZE), score(SENSITIVITY_CLONE_SIZE)
    coupling = [to_api_coupling(c) for c in primary.coupling]
    # The coupling note only goes with a coupling.
    notes = [COUPLING_NOTE] if coupling else []
    if primary.status == "single_lineage":
        notes.append("The scored clones fall in only one lineage; clones are scored over two or more.")
    if codes.n_missing_donor > 0:
        notes.append(f"{codes.n_missing_donor} cells with a receptor have no donor and were not scored.")
    if primary.status == "too_few_clones":
        notes.append(
            f"Only {primary.n_clones} clones were scored; "
            f"at least {MIN_CLONES} are needed before the score counts."
        )
    return ClonalAgreement(
        status=primary.status,
        settings=settings,
        primary=to_api_score(primary),
        sensitivity_m3=to_api_score(m3),
        coupling=coupling,
        notes=notes,
    )
