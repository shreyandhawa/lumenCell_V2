"""Keeping clone identity out of inference.

Clone labels are used only for scoring. Cells of one clone share a receptor, so receptor segment
genes (TRAV, TRBV, TRAC, ...) would let a method recover clones from expression. They are kept
out of HVG selection, and with them out of the PCA and everything built on it.

Clone and receptor annotations are kept out of what methods see in the same spirit: they are
scoring-only, so no root, cluster column or label column may be one of them.
"""

import re
from collections.abc import Iterable

import numpy as np
import numpy.typing as npt

# Recorded with every run; raise it whenever what counts as scoring-only changes.
ISOLATION_VERSION = 1

# Annotation columns holding clone or receptor data (clone_id, IR_VJ_1_junction_aa,
# raw_clonotype_id, chain_pairing, ...), matched by name, ignoring case. A donor column is
# scoring-only only when the caller names it.
RECEPTOR_COLUMN = re.compile(
    r"^(?:clone_id|clonotype|cdr3|ir_|airr:|receptor_)|^(?:raw_clonotype_id|chain_pairing)$", re.IGNORECASE
)

# Spec amendment 1: the segment families V, D, J (followed by a digit) and the constant regions
# TRAC, TRBC1, TRBC2, TRDC, TRGC1, ... Case-insensitive so mouse symbols (Trbc1) match. TRADD,
# TRAT1, TRAPPC1, TRDN and TRDMT1 are not receptor genes and do not match.
RECEPTOR_GENE = re.compile(r"^TR[ABDG]([VDJ]\d|C\d*$)", re.IGNORECASE)


def receptor_gene_mask(symbols: npt.NDArray[np.str_]) -> npt.NDArray[np.bool_]:
    """True where the gene symbol is a T-cell receptor segment gene."""
    return np.fromiter(
        (RECEPTOR_GENE.match(str(s)) is not None for s in symbols), dtype=np.bool_, count=len(symbols)
    )


def column_names(named: Iterable[str]) -> tuple[str, ...]:
    """``named`` as a tuple of names. A bare string is refused: it would be read one character at a time."""
    if isinstance(named, str):
        raise TypeError(f"clonal columns must be a collection of column names, not the string {named!r}")
    return tuple(named)


def scoring_only_columns(columns: Iterable[object], named: Iterable[str] = ()) -> frozenset[str]:
    """Every column whose name matches RECEPTOR_COLUMN, plus every name in ``named``."""
    return frozenset(str(c) for c in columns if RECEPTOR_COLUMN.match(str(c))) | frozenset(
        column_names(named)
    )


def scoring_only_message(column: str, role: str) -> str:
    """Why a scoring-only column cannot be the ``role`` (e.g. "root column"), in its own words."""
    if RECEPTOR_COLUMN.match(column):
        held = "holds clone or receptor data"
    else:
        held = "is set as a clone or donor column"
    return f"{column!r} {held}, which lumenCell uses only to score runs; it cannot be the {role}"
