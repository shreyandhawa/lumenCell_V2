"""N2's lineage-onset split (the N2a amendment, declared under [clonal] in benchmarks/thresholds.toml)."""

from collections.abc import Sequence

import numpy as np
import pytest

from lumencell.clonal.labels import CloneCodes, encode_clones
from lumencell.clonal.scores import _lineage_weights, _n2_groups, _onset_segments, score_clonal

NAN = float("nan")


def codes(labels: Sequence[object], donors: Sequence[object] | None = None) -> CloneCodes:
    return encode_clones(
        np.array(labels, dtype=object), None if donors is None else np.array(donors, dtype=object)
    )


def test_onsets_are_the_earliest_whole_weight_cell_of_each_lineage() -> None:
    # Lineage 0 alone from 0.1; lineage 1 alone from 0.6. Cell 2 is split, so it sets no onset.
    weights = np.array([[1, 0], [1, 0], [0.5, 0.5], [0, 1], [0, 1]], dtype=np.float64)
    pt = np.array([0.1, 0.5, 0.05, 0.6, 0.9])
    assert _onset_segments(weights, pt).tolist() == [1, 1, 0, 2, 2]


def test_a_lineage_without_a_whole_weight_cell_adds_no_cut() -> None:
    # Lineage 2 holds only shared cells (1/3 each): two cuts, at 0.2 and 0.7, never three.
    third = 1 / 3
    weights = np.array([[1, 0, 0], [third, third, third], [0, 1, 0], [1, 0, 0]], dtype=np.float64)
    pt = np.array([0.2, 0.0, 0.7, 0.9])
    assert _onset_segments(weights, pt).tolist() == [1, 0, 2, 2]


def test_nan_pseudotime_never_sets_an_onset_and_keeps_its_own_stratum() -> None:
    weights = np.array([[1, 0], [0, 1], [0, 1], [1, 0]], dtype=np.float64)
    pt = np.array([0.3, NAN, 0.8, NAN])
    assert _onset_segments(weights, pt).tolist() == [1, 0, 2, 0]
    groups = _n2_groups(np.zeros(4, dtype=np.int64), pt, weights, 10, "lineage_onset")
    # The two NaN cells share a stratum that no finite cell is in.
    assert groups[1] == groups[3] and groups[1] not in {groups[0], groups[2]}


def test_near_equal_onsets_are_not_merged() -> None:
    weights = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1]], dtype=np.float64)
    pt = np.array([0.3, 0.6513, 0.6527, 0.9])
    assert _onset_segments(weights, pt).tolist() == [1, 2, 3, 3]


def test_none_keeps_the_phase_1_strata() -> None:
    weights = np.array([[1, 0], [0, 1], [1, 0]], dtype=np.float64)
    donor = np.array([0, 1, 1])
    pt = np.array([0.1, 0.5, NAN])
    groups = _n2_groups(donor, pt, weights, 10, "none")
    assert groups.tolist() == [0 * 11 + 0, 1 * 11 + 5, 1 * 11 + 10]


def step(n_per_side: int = 120, seed: int = 0) -> tuple[np.ndarray, CloneCodes, np.ndarray]:
    """A binary2-like step inside one decile, no clonal signal beyond depth.

    Before t = 0.65 every cell is lineage 0 (B); after it, lineages 0, 1 and 2 in equal shares.
    All pseudotime lies in [0.6, 0.7), and every clone sits on one side of the step, so the
    within-decile shuffle mixes pure-B clones with mixed-side cells.
    """
    rng = np.random.default_rng(seed)
    pt = np.concatenate([rng.uniform(0.60, 0.649, n_per_side), rng.uniform(0.65, 0.70, n_per_side)])
    lineage = np.concatenate([np.zeros(n_per_side, int), rng.integers(3, size=n_per_side)])
    # Clones of 3 cells, each within one side of the step, labels shuffled within a side.
    before = rng.permutation(n_per_side) // 3
    after = rng.permutation(n_per_side) // 3 + n_per_side // 3
    clone = np.concatenate([before, after])
    membership = np.eye(3, dtype=bool)[lineage].T
    return membership, codes([f"c{k}" for k in clone]), pt


def test_a_step_inside_one_decile_is_split_and_the_null_holds() -> None:
    membership, c, pt = step()
    baseline = np.zeros(pt.size, dtype=np.int64)
    weights = _lineage_weights(membership, np.arange(pt.size))
    # Lineage 0's onset is the first cell; 1 and 2 start after the step, so no stratum spans it.
    segments = _onset_segments(weights, pt)
    later = min(pt[membership[1]].min(), pt[membership[2]].min())
    assert later >= 0.65
    assert set(segments[pt < later].tolist()) == {1}
    assert set(segments[pt >= later].tolist()) <= {2, 3}
    deciles = score_clonal(membership, c, pt, baseline, n_bins=1, n_permutations=199, n2_split="none")
    split = score_clonal(membership, c, pt, baseline, n_bins=1, n_permutations=199, n2_split="lineage_onset")
    assert deciles.nulls[1].null == split.nulls[1].null == "N2"
    # One bin: the decile N2 is N1 and rejects; the split removes the step and does not.
    assert deciles.nulls[1].p <= 0.01
    assert split.nulls[1].p > 0.05
    # N1 is untouched by the split.
    assert split.nulls[0] == deciles.nulls[0]


def test_the_strata_never_read_clone_labels() -> None:
    membership, c, pt = step()
    baseline = np.zeros(pt.size, dtype=np.int64)
    weights = _lineage_weights(membership, np.arange(pt.size))
    donor = c.donor
    groups = _n2_groups(donor, pt, weights, 10, "lineage_onset")
    relabelled = CloneCodes(c.clone[::-1].copy(), c.donor, c.stratum, c.n_missing_donor)
    assert np.array_equal(groups, _n2_groups(relabelled.donor, pt, weights, 10, "lineage_onset"))
    # And through score_clonal: the same N2 null distribution whatever the labels' names.
    renamed = codes([f"x{k}" for k in c.clone])
    a = score_clonal(membership, c, pt, baseline, n_permutations=99, n2_split="lineage_onset")
    b = score_clonal(membership, renamed, pt, baseline, n_permutations=99, n2_split="lineage_onset")
    assert a.nulls[1].p == b.nulls[1].p
    assert a.nulls[1].mean_entropy_mean == pytest.approx(b.nulls[1].mean_entropy_mean, abs=1e-12)


def test_a_split_that_leaves_the_strata_unchanged_gives_the_phase_1_n2() -> None:
    # Both lineages start at the earliest cell: every finite cell is in segment 2.
    membership, c, pt = step()
    membership = membership[:2].copy()
    membership[1, :] = ~membership[0, :]
    first = int(np.argmin(pt))
    second = int(np.argsort(pt)[1])
    membership[:, first] = [True, False]
    membership[:, second] = [False, True]
    pt = pt.copy()
    pt[second] = pt[first]
    baseline = np.zeros(pt.size, dtype=np.int64)
    a = score_clonal(membership, c, pt, baseline, n_permutations=99, n2_split="none")
    b = score_clonal(membership, c, pt, baseline, n_permutations=99, n2_split="lineage_onset")
    assert a.nulls == b.nulls and a.p == b.p


def test_own_and_oracle_modes_use_the_same_rule() -> None:
    membership, c, pt = step()
    baseline = np.zeros(pt.size, dtype=np.int64)
    # Oracle mode with every cell an oracle cell and one-hot masks is own-cells mode in all but name.
    own = score_clonal(membership, c, pt, baseline, n_permutations=99, n2_split="lineage_onset")
    oracle = score_clonal(
        membership,
        c,
        pt,
        baseline,
        oracle_cells=np.ones(pt.size, dtype=bool),
        n_permutations=99,
        n2_split="lineage_onset",
    )
    assert own.mode == "own_cells" and oracle.mode == "oracle_cells"
    assert own.nulls == oracle.nulls


def test_in_oracle_mode_only_whole_weight_cells_set_onsets() -> None:
    # The method puts the early cells in both lineages (1/2 each): they set no onset, so the
    # earliest cut is lineage 0's first exclusive cell.
    membership = np.array([[1, 1, 1, 0, 1, 0], [1, 1, 0, 1, 0, 1]], dtype=bool)
    pt = np.array([0.0, 0.1, 0.4, 0.5, 0.8, 0.9])
    weights = _lineage_weights(membership, np.arange(6))
    assert _onset_segments(weights, pt).tolist() == [0, 0, 1, 2, 2, 2]


def test_an_unknown_split_is_refused() -> None:
    membership, c, pt = step()
    with pytest.raises(ValueError, match="n2_split"):
        score_clonal(membership, c, pt, np.zeros(pt.size, dtype=np.int64), n2_split="quantile")  # type: ignore[arg-type]
