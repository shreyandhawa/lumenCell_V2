"""Clonal agreement and reported (never gated) fields of a real dataset's run.

The clonal result is scored in oracle-cells mode, the same as the synthetic benchmark, through
``clonal_row``. A declared ``not_assessable`` result is shown with its reason and scores nothing.
The reported fields describe the run next to the gates: convergence of receptors across captures,
agreement with 10x clonotypes, the TRB-only clone definition, and pseudotime against real time.
Clone sharing relates states; no field here says which state came first.
"""

import math
import re
from dataclasses import replace
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import spearmanr

from lumencell.bench.datasets import RealDataset
from lumencell.bench.runner import ClonalRow, clonal_row
from lumencell.bench.truth import Truth
from lumencell.clonal.labels import is_missing
from lumencell.engine.schema import TrajectoryResult

_NUMBER = re.compile(r"[-+]?\d*\.?\d+")


def adjusted_rand_index(a: npt.NDArray[np.int64], b: npt.NDArray[np.int64]) -> float:
    """Hubert-Arabie adjusted Rand index of two partitions given as integer codes.

    NaN with fewer than two cells. The denominator is zero (expected index equals its maximum)
    only when both partitions are one cluster or both are all singletons, and then they are
    identical up to relabelling: 1.0.
    """
    n = int(a.shape[0])
    if n < 2:
        return math.nan
    _, ai = np.unique(a, return_inverse=True)
    _, bi = np.unique(b, return_inverse=True)
    table = np.zeros((int(ai.max()) + 1, int(bi.max()) + 1), dtype=np.int64)
    np.add.at(table, (ai, bi), 1)

    def pairs(x: npt.NDArray[np.int64]) -> float:
        return float((x * (x - 1) // 2).sum())

    index = pairs(table.ravel())
    rows, cols = pairs(table.sum(axis=1)), pairs(table.sum(axis=0))
    total = n * (n - 1) / 2
    expected = rows * cols / total
    maximum = (rows + cols) / 2
    if maximum == expected:
        return 1.0
    return (index - expected) / (maximum - expected)


def clone_recurrence(clones: npt.NDArray[np.object_], captures: npt.NDArray[np.object_]) -> tuple[float, int]:
    """Share of distinct clone strings found in two or more captures, and the number of strings.

    Over cells with a non-missing clone string and a non-missing capture; ``(nan, 0)`` when there
    are none. A cell without a capture is left out, never counted as a capture called "nan". The
    namespace is ignored on purpose: it estimates convergent receptors if the captures hold
    different mice.
    """
    keep = ~is_missing(clones) & ~is_missing(captures)
    if not keep.any():
        return math.nan, 0
    frame = pd.DataFrame({"clone": clones[keep].astype(str), "capture": captures[keep].astype(str)})
    n_captures = frame.groupby("clone")["capture"].nunique()
    return float((n_captures >= 2).mean()), len(n_captures)


def real_clonal_row(
    ds: RealDataset, truth: Truth, result: TrajectoryResult, pca: Any, neighbors_pcs: int
) -> ClonalRow:
    if ds.clonal is None:
        return ClonalRow("not_applicable")
    if not ds.clonal.assessable:
        return ClonalRow("not_assessable", reason=f"declared {ds.clonal.declared}, {ds.clonal.reason}")
    return clonal_row(truth, result, pca, neighbors_pcs)


def _numbers(values: pd.Series) -> npt.NDArray[np.float64]:
    """Day values as numbers: a numeric column is used as is (NaN stays NaN); only a string or
    object column is read by its first number, so ``D21`` and ``21`` both give 21."""
    if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
        return np.asarray(pd.to_numeric(values), dtype=np.float64)
    return np.array(
        [float(m.group()) if (m := _NUMBER.search(str(v))) else math.nan for v in values], dtype=np.float64
    )


def _codes(values: npt.NDArray[np.object_]) -> npt.NDArray[np.int64]:
    codes: npt.NDArray[np.int64] = pd.factorize(pd.Series(values.astype(str)))[0].astype(np.int64)
    return codes


def reported_metrics(
    ds: RealDataset,
    obs: pd.DataFrame,
    result: TrajectoryResult,
    truth: Truth,
    pca: Any,
    neighbors_pcs: int,
) -> dict[str, float | str]:
    """The ``reported_*`` fields of one run. None of them is gated."""
    out: dict[str, float | str] = {}
    if ds.report_time_key is not None:
        pt = np.asarray(result.pseudotime, dtype=np.float64)
        day = _numbers(obs[ds.report_time_key])
        ok = np.isfinite(pt) & np.isfinite(day)
        rho = math.nan
        if len(np.unique(day[ok])) >= 2 and len(np.unique(pt[ok])) >= 2:
            rho = float(spearmanr(pt[ok], day[ok]).statistic)
        out["reported_day_spearman"] = rho
    c = ds.clonal
    if c is None:
        return out
    clone = obs[c.clone_key].to_numpy(dtype=object)
    capture = obs[c.capture_key].to_numpy(dtype=object)
    share, n_strings = clone_recurrence(clone, capture)
    out["reported_clone_recurrence"] = share
    out["reported_n_clone_strings"] = n_strings

    raw = obs[c.raw_clonotype_key].to_numpy(dtype=object)
    has_clone, has_raw = ~is_missing(clone), ~is_missing(raw)
    both = has_clone & has_raw
    keyed = np.array([f"{p}\x1f{q}" for p, q in zip(capture[both], clone[both], strict=True)], dtype=object)
    out["reported_raw_clonotype_ari"] = adjusted_rand_index(_codes(keyed), _codes(raw[both]))
    # Declared text, literally: cells with a clone_id over cells with a raw_clonotype_id. It can
    # exceed 1 when cells carry a clone_id but no raw id.
    out["reported_labelled_share_of_raw"] = (
        float(has_clone.sum() / has_raw.sum()) if has_raw.any() else math.nan
    )

    trb_truth = replace(
        truth,
        clone_ids=obs[c.trb_only_key].to_numpy(dtype=object),
        donors=obs[c.donor_key].to_numpy(dtype=object),
    )
    trb = real_clonal_row(ds, trb_truth, result, pca, neighbors_pcs)
    out["reported_trb_status"] = trb.status
    out["reported_trb_excess_ratio"] = trb.excess_ratio
    out["reported_trb_p"] = trb.p
    return out
