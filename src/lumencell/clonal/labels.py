"""Raw clone labels, donors and strata to integer codes.

Clone identity is the pair (donor, label): the same label string in two donors is two clones.
A missing label means the cell has no receptor. A labelled cell whose donor is missing is left
out of scoring and counted, never pooled into a "nan" donor. ``is_missing`` decides what missing
means for labels, donors and strata alike.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

IntArray = npt.NDArray[np.int64]
BoolArray = npt.NDArray[np.bool_]
FloatArray = npt.NDArray[np.floating[Any]]
ObjectArray = npt.NDArray[np.object_]

_MISSING_STRINGS = frozenset({"", "nan", "None"})


@dataclass(frozen=True)
class CloneCodes:
    """Per-cell codes, each ``(n,)`` int64.

    ``clone`` is -1 for a cell without a label or without a donor, otherwise a dense id over the
    distinct (donor, label) pairs. ``donor`` is -1 where the donor is missing. ``stratum`` codes
    (donor, each strata value) when strata are given, and is None otherwise. ``n_missing_donor``
    counts labelled cells whose donor is missing.
    """

    clone: IntArray
    donor: IntArray
    stratum: IntArray | None
    n_missing_donor: int


def is_missing(values: ObjectArray) -> BoolArray:
    """True where a value is None, a float NaN, ``pd.NA``, or ``""``, ``"nan"`` or ``"None"``."""
    return np.fromiter((_missing(v) for v in values), dtype=np.bool_, count=len(values))


def _missing(value: object) -> bool:
    if value is None or value is pd.NA:
        return True
    if isinstance(value, float | np.floating) and bool(np.isnan(value)):
        return True
    return str(value) in _MISSING_STRINGS


def encode_clones(
    labels: ObjectArray, donors: ObjectArray | None, strata: Sequence[ObjectArray] = ()
) -> CloneCodes:
    """Codes for one clone column, an optional donor column and any strata columns.

    ``donors=None`` means the file holds one donor (every code 0). Callers pass
    ``obs[col].to_numpy(dtype=object)``, which keeps only the categories that occur.
    """
    labels = _column(labels, "labels", None)
    n = labels.shape[0]
    labelled = ~is_missing(labels)
    donor = np.zeros(n, dtype=np.int64) if donors is None else _codes(_column(donors, "donors", n))
    keep = labelled & (donor >= 0)
    clone = np.full(n, -1, dtype=np.int64)
    if keep.any():
        clone[keep] = _group(np.stack([donor[keep], _codes(labels[keep])], axis=1))
    stratum = None
    if strata:
        keys = [donor]
        # A missing strata value codes as -1: a level of its own, apart from every string.
        keys += [_codes(_column(column, f"strata[{i}]", n)) for i, column in enumerate(strata)]
        stratum = _group(np.stack(keys, axis=1))
    return CloneCodes(clone, donor, stratum, int((labelled & (donor < 0)).sum()))


def _column(values: ObjectArray, name: str, n: int | None) -> ObjectArray:
    column = np.asarray(values, dtype=object)
    if column.ndim != 1 or (n is not None and column.shape[0] != n):
        raise ValueError(f"{name} must be one value per cell, shape ({n if n is not None else 'n'},)")
    return column


def _codes(values: ObjectArray) -> IntArray:
    """Dense codes over the distinct strings of the values; -1 where a value is missing."""
    missing = is_missing(values)
    codes = np.full(values.shape[0], -1, dtype=np.int64)
    # Python strings in an object array: a fixed-width numpy string would drop trailing "\0"s.
    present = np.array([str(v) for v in values[~missing]], dtype=object)
    if present.size:
        codes[~missing] = np.unique(present, return_inverse=True)[1].reshape(-1)
    return codes


def _group(keys: IntArray) -> IntArray:
    """Dense codes over the distinct rows of an ``(n, k)`` key array."""
    return np.unique(keys, axis=0, return_inverse=True)[1].reshape(-1).astype(np.int64)
