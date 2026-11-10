from collections.abc import Callable
from dataclasses import replace

import numpy as np
import pytest
from scipy import sparse

from lumencell.engine.schema import (
    PRIMARY_SHARED,
    PRIMARY_UNASSIGNED,
    InvariantError,
    Lineage,
    Milestone,
    ResolvedRoot,
    Segment,
    TrajectoryResult,
)

N = 12


def _valid() -> TrajectoryResult:
    """A small Y: root(0) -> fork(1) -> tips 2 and 3.

    Cells 0-3 are the trunk, 4-7 lead to tip 2, 8-10 to tip 3, and cell 11 is unassigned.
    """
    pt = np.array([0.0, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 1.0, 0.5, 0.7, 0.9, np.nan], dtype=np.float32)
    trunk = np.zeros(N, dtype=bool)
    trunk[:4] = True
    a = trunk.copy()
    a[4:8] = True
    b = trunk.copy()
    b[8:11] = True
    seg = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, -1], dtype=np.int32)
    prog = np.array([0, 0.33, 0.66, 1, 0, 0.33, 0.66, 1, 0, 0.5, 1, np.nan], dtype=np.float32)
    weights = sparse.csr_matrix(
        (np.ones(11), (np.arange(11), [0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3])), shape=(N, 4)
    )
    return TrajectoryResult(
        method="toy",
        method_version="0",
        topology_kind="tree",
        milestones=[
            Milestone(0, "root", "r", 0.0),
            Milestone(1, "fork", "f", 0.3),
            Milestone(2, "tip", "A", 1.0),
            Milestone(3, "tip", "B", 0.9),
        ],
        segments=[Segment(0, 0, 1, 0.3), Segment(1, 1, 2, 0.7), Segment(2, 1, 3, 0.6)],
        lineages=[
            Lineage(0, "A", "A", 2, "A", 1.0, 8, [0, 1]),
            Lineage(1, "B", "B", 3, "B", 1.0, 7, [0, 2]),
        ],
        divergence=[],
        pseudotime=pt,
        membership=np.stack([a, b]),
        segment=seg,
        progress=prog,
        root=ResolvedRoot("cell:c0", 0, "c0", None, None, "given cell"),
        inference_rep="X_diffmap",
        inference_dims=10,
        params={},
        milestone_weights=weights,
    )


def test_valid_result_passes() -> None:
    _valid().validate()


def test_primary_lineage() -> None:
    primary = _valid().primary_lineage
    assert primary.dtype == np.int16
    assert list(primary[:4]) == [PRIMARY_SHARED] * 4
    assert list(primary[4:8]) == [0] * 4
    assert list(primary[8:11]) == [1] * 3
    assert primary[11] == PRIMARY_UNASSIGNED
    assert _valid().n_unassigned == 1


def _break_root_pt(r: TrajectoryResult) -> None:
    r.pseudotime[0] = 0.2


def _break_range(r: TrajectoryResult) -> None:
    r.pseudotime[5] = 1.5


def _break_orphan(r: TrajectoryResult) -> None:
    r.membership[:, 5] = False


def _break_ghost(r: TrajectoryResult) -> None:
    r.membership[0, 11] = True


def _break_progress(r: TrajectoryResult) -> None:
    r.progress[4] = np.nan


def _break_segment_index(r: TrajectoryResult) -> None:
    r.segment[4] = 9


def _break_path(r: TrajectoryResult) -> None:
    r.lineages[0].segments = [1]


def _break_names(r: TrajectoryResult) -> None:
    r.lineages[1].name = "A"


def _break_two_roots(r: TrajectoryResult) -> None:
    r.milestones[1] = Milestone(1, "root", "f", 0.3)


def _break_fate(r: TrajectoryResult) -> None:
    r.fate = np.full((2, N), 0.3, dtype=np.float32)


def _break_shape(r: TrajectoryResult) -> None:
    r.segment = r.segment[:5]


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (_break_root_pt, "root cell's pseudotime must be 0"),
        (_break_range, r"pseudotime must be in \[0, 1\]"),
        (_break_orphan, "must belong to a lineage"),
        (_break_ghost, "must belong to no lineage"),
        (_break_progress, "progress"),
        (_break_segment_index, "segment values out of range"),
        (_break_path, "must start at the root"),
        (_break_names, "names must be unique"),
        (_break_two_roots, "exactly one root"),
        (_break_fate, "sum to 1"),
        (_break_shape, "one value per cell"),
    ],
)
def test_each_invariant_is_enforced(mutate: Callable[[TrajectoryResult], None], message: str) -> None:
    result = _valid()
    mutate(result)
    with pytest.raises(InvariantError, match=message):
        result.validate()


def test_valid_fate_passes() -> None:
    r = _valid()
    fate = np.zeros((2, N), dtype=np.float32)
    fate[0, :] = 0.5
    fate[1, :] = 0.5
    r.fate = fate
    r.validate()


def test_all_problems_are_reported_together() -> None:
    r = _valid()
    r.pseudotime[0] = 0.2
    r.lineages[1] = replace(r.lineages[1], name="A")
    with pytest.raises(InvariantError) as info:
        r.validate()
    assert len(info.value.problems) >= 2
