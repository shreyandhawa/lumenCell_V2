"""The per-run cell-cycle check: lineages made mostly of cycling cells (spec section on cell cycle).

Report only. A lineage is flagged when its share of S or G2M cells, over the cells that sit on that
lineage alone, is at least ``CELL_CYCLE_RATIO_MIN`` times the dataset's share and at least
``CELL_CYCLE_SHARE_MIN``. The notice text carries the lineage name the user gave; it is served
and never used in a path or a log line.
"""

from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt

from lumencell.api import models
from lumencell.trust.constants import CELL_CYCLE_RATIO_MIN, CELL_CYCLE_SHARE_MIN

CYCLING = frozenset({"S", "G2M"})


def _cycling(cc_phase: npt.NDArray[np.str_]) -> npt.NDArray[np.bool_]:
    return np.isin(cc_phase, sorted(CYCLING))


def cell_cycle_info(cc_info: dict[str, Any], cc_phase: npt.NDArray[np.str_]) -> models.CellCycleInfo:
    scored = cc_info["status"] == "scored"
    return models.CellCycleInfo(
        status=cc_info["status"],
        reason=cc_info["reason"],
        requested_species=cc_info["requested_species"],
        species=cc_info["species"],
        gene_set=cc_info["gene_set"],
        genes_found=sum(cc_info["genes_found"].values()),
        genes_listed=sum(cc_info["genes_listed"].values()),
        source=cc_info["source"] or "",
        share_cycling=float(_cycling(cc_phase).mean()) if scored and len(cc_phase) else None,
    )


def cell_cycle_check(
    primary_lineage: npt.NDArray[np.int16],
    names: Sequence[str],
    cc_phase: npt.NDArray[np.str_],
    cc_info: dict[str, Any],
) -> models.CellCycleCheck:
    if cc_info["status"] != "scored":
        return models.CellCycleCheck(
            status="not_scored",
            dataset_share=None,
            ratio_min=CELL_CYCLE_RATIO_MIN,
            share_min=CELL_CYCLE_SHARE_MIN,
            lineages=[],
            warnings=[],
            note=cc_info["reason"],
        )
    cycling = _cycling(cc_phase)
    dataset_share = float(cycling.mean()) if len(cycling) else 0.0
    lineages: list[models.LineageCycling] = []
    warnings: list[models.Notice] = []
    for lineage, name in enumerate(names):
        mine = primary_lineage == lineage
        n = int(mine.sum())
        share = float(cycling[mine].mean()) if n else None
        flagged = (
            share is not None
            and share >= CELL_CYCLE_RATIO_MIN * dataset_share
            and share >= CELL_CYCLE_SHARE_MIN
        )
        lineages.append(
            models.LineageCycling(lineage=lineage, name=name, n_cells=n, share_cycling=share, flagged=flagged)
        )
        if flagged:
            label = name if name.startswith("Lineage ") else f"Lineage {name}"
            warnings.append(
                models.Notice(
                    code="lineage_mostly_cycling",
                    severity="warning",
                    message=(
                        f"{label} is mostly cycling cells; it may be a cell-cycle state rather than a fate."
                    ),
                )
            )
    return models.CellCycleCheck(
        status="scored",
        dataset_share=dataset_share,
        ratio_min=CELL_CYCLE_RATIO_MIN,
        share_min=CELL_CYCLE_SHARE_MIN,
        lineages=lineages,
        warnings=warnings,
        note=None,
    )
