"""Keeping clone identity out of inference.

Clone labels are used only for scoring. Cells of one clone share a receptor, so receptor segment
genes (TRAV, TRBV, TRAC, ...) would let a method recover clones from expression. They are kept
out of HVG selection, and with them out of the PCA and everything built on it.
"""

import re

import numpy as np
import numpy.typing as npt

# Spec amendment 1: the segment families V, D, J (followed by a digit) and the constant regions
# TRAC, TRBC1, TRBC2, TRDC, TRGC1, ... Case-insensitive so mouse symbols (Trbc1) match. TRADD,
# TRAT1, TRAPPC1, TRDN and TRDMT1 are not receptor genes and do not match.
RECEPTOR_GENE = re.compile(r"^TR[ABDG]([VDJ]\d|C\d*$)", re.IGNORECASE)


def receptor_gene_mask(symbols: npt.NDArray[np.str_]) -> npt.NDArray[np.bool_]:
    """True where the gene symbol is a T-cell receptor segment gene."""
    return np.fromiter(
        (RECEPTOR_GENE.match(str(s)) is not None for s in symbols), dtype=np.bool_, count=len(symbols)
    )
