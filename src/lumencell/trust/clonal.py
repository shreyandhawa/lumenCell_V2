"""Dataset-level clonal settings: what a user may name, whether it is usable, and what it holds.

The settings only say which obs columns carry clone labels, donors and strata. Candidates are
suggested by name and never applied; the user names the columns.
"""

import re
import sys
from collections.abc import Iterable
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from pydantic import ValidationError

from lumencell.api.models import ClonalSettings, ClonalSummary
from lumencell.clonal.labels import encode_clones, is_missing
from lumencell.trust.constants import MAX_DONORS

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

    if settings.clone_key in obs.columns:
        clone = obs[settings.clone_key]
        if not _is_label_column(clone):
            problems.append(
                f"the clone column {settings.clone_key!r} must hold text labels "
                f"(categorical or string), not {clone.dtype}"
            )
        elif bool(is_missing(_present(clone)).all()):
            problems.append(f"the clone column {settings.clone_key!r} has no clone label in any cell")
    if settings.donor_key is not None and settings.donor_key in obs.columns:
        donors = _present(obs[settings.donor_key])
        n = len({str(v) for v in donors[~is_missing(donors)]})
        if not 1 <= n <= MAX_DONORS:
            problems.append(
                f"the donor column {settings.donor_key!r} has {n} distinct values; "
                f"it needs between 1 and {MAX_DONORS}"
            )
    return problems


def settings_summary(obs: pd.DataFrame, settings: ClonalSettings) -> ClonalSummary:
    """Clone sizes over every labelled cell with a donor, singletons included."""
    labels = _present(obs[settings.clone_key])
    donors = None if settings.donor_key is None else _present(obs[settings.donor_key])
    codes = encode_clones(labels, donors)
    sizes = np.bincount(codes.clone[codes.clone >= 0])
    present_donors = np.unique(codes.donor[codes.donor >= 0])

    def pct(q: float) -> float:
        return float(np.percentile(sizes, q)) if sizes.size else 0.0

    return ClonalSummary(
        settings=settings,
        n_cells_labelled=int(sizes.sum()),
        n_clones=int(sizes.size),
        n_clones_2plus=int((sizes >= 2).sum()),
        n_donors=int(present_donors.size),
        size_median=float(np.median(sizes)) if sizes.size else 0.0,
        size_p90=pct(90),
        size_p99=pct(99),
        size_max=int(sizes.max()) if sizes.size else 0,
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
