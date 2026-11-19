"""The per-run cell-cycle check: which lineages are made mostly of cycling cells."""

from __future__ import annotations

from typing import Any

import numpy as np

from lumencell.trust.cell_cycle import CYCLING, cell_cycle_check, cell_cycle_info

SCORED: dict[str, Any] = {
    "status": "scored",
    "reason": None,
    "requested_species": "auto",
    "species": "human",
    "gene_set": "tirosh2016",
    "genes_found": {"s": 40, "g2m": 50},
    "genes_listed": {"s": 43, "g2m": 54},
    "source": "Tirosh et al. 2016",
}


def _phase(n: int, cycling: list[int]) -> np.ndarray:
    out = np.full(n, "G1", dtype="<U8")
    out[cycling] = "S"
    return out


def _lineage(*values: int) -> np.ndarray:
    return np.array(values, dtype=np.int16)


def test_a_mostly_cycling_lineage_is_flagged() -> None:
    primary = _lineage(0, 0, 0, 0, 1, 1, 1, 1, -1, -2)
    check = cell_cycle_check(primary, ["Alpha", "Beta"], _phase(10, [0, 1, 2]), SCORED)
    assert check.status == "scored"
    assert check.dataset_share == 0.3
    a, b = check.lineages
    assert (a.n_cells, a.share_cycling, a.flagged) == (4, 0.75, True)
    assert (b.n_cells, b.share_cycling, b.flagged) == (4, 0.0, False)
    assert [w.message for w in check.warnings] == [
        "Lineage Alpha is mostly cycling cells; it may be a cell-cycle state rather than a fate."
    ]
    assert check.warnings[0].code == "lineage_mostly_cycling"
    assert check.warnings[0].severity == "warning"


def test_the_bounds_are_inclusive() -> None:
    primary = _lineage(0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    phase = _phase(16, [0, 1, 4, 5])  # lineage 0: 2/4 = 0.5; dataset 4/16 = 0.25
    check = cell_cycle_check(primary, ["A", "B"], phase, SCORED)
    assert check.dataset_share == 0.25
    assert check.lineages[0].share_cycling == 0.5
    assert check.lineages[0].flagged is True


def test_a_high_share_below_twice_the_dataset_is_not_flagged() -> None:
    primary = _lineage(*([0] * 5 + [1] * 5))
    phase = _phase(10, [0, 1, 2, 5, 6])  # lineage 0: 0.6; dataset 0.5 (needs 1.0)
    check = cell_cycle_check(primary, ["A", "B"], phase, SCORED)
    assert check.dataset_share == 0.5
    assert not any(lin.flagged for lin in check.lineages)
    assert check.warnings == []


def test_shared_and_unassigned_cells_count_only_for_the_dataset() -> None:
    primary = _lineage(0, 0, -1, -1, -2, -2)
    phase = _phase(6, [2, 3, 4, 5])
    check = cell_cycle_check(primary, ["A"], phase, SCORED)
    assert check.dataset_share == 4 / 6
    assert check.lineages[0].n_cells == 2
    assert check.lineages[0].share_cycling == 0.0


def test_a_lineage_without_single_lineage_cells_has_no_share() -> None:
    primary = _lineage(0, 0, -1)
    check = cell_cycle_check(primary, ["A", "B"], _phase(3, [0]), SCORED)
    b = check.lineages[1]
    assert (b.n_cells, b.share_cycling, b.flagged) == (0, None, False)


def test_names_starting_with_lineage_are_not_doubled() -> None:
    primary = _lineage(0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
    check = cell_cycle_check(primary, ["Lineage 2", "B"], _phase(10, [0, 1, 2, 3]), SCORED)
    assert check.warnings[0].message.startswith("Lineage 2 is mostly")


def test_not_scored_passes_the_reason_on() -> None:
    info = {**SCORED, "status": "not_scored", "reason": "only 0 of 97 Tirosh 2016 genes were found"}
    phase = np.full(4, "unscored", dtype="<U8")
    check = cell_cycle_check(_lineage(0, 0, 1, 1), ["A", "B"], phase, info)
    assert check.status == "not_scored"
    assert check.note == info["reason"]
    assert check.lineages == [] and check.warnings == []
    assert check.dataset_share is None


def test_dataset_info_sums_both_lists() -> None:
    info = cell_cycle_info(SCORED, _phase(10, [0, 1]))
    assert (info.status, info.genes_found, info.genes_listed) == ("scored", 90, 97)
    assert info.share_cycling == 0.2
    assert info.species == "human" and info.gene_set == "tirosh2016"
    off = cell_cycle_info(
        {**SCORED, "status": "not_scored", "reason": "r", "species": None},
        np.full(10, "unscored", dtype="<U8"),
    )
    assert off.share_cycling is None and off.reason == "r" and off.species is None
    assert frozenset({"S", "G2M"}) == CYCLING
