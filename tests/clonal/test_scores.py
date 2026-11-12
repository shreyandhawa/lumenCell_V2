import math
from collections.abc import Sequence
from dataclasses import fields

import numpy as np
import pytest
from scipy.stats import false_discovery_control

import lumencell.clonal as clonal
from lumencell.clonal.labels import CloneCodes, encode_clones
from lumencell.clonal.scores import (
    ClonalResult,
    _pseudotime_bins,
    baseline_space,
    kmeans_partition,
    score_clonal,
)

NAN = float("nan")
PT0 = np.zeros(4)
BASE0 = np.zeros(4, dtype=np.int64)


def own(rows: list[list[int]]) -> np.ndarray:
    return np.array(rows, dtype=bool)


def codes(
    labels: Sequence[object],
    donors: Sequence[object] | None = None,
    strata: Sequence[Sequence[object]] = (),
) -> CloneCodes:
    return encode_clones(
        np.array(labels, dtype=object),
        None if donors is None else np.array(donors, dtype=object),
        [np.array(s, dtype=object) for s in strata],
    )


def biased(
    n_clones: int = 60,
    size: int = 3,
    n_lineages: int = 2,
    bias: float = 0.8,
    seed: int = 0,
    by_depth: bool = False,
) -> tuple[np.ndarray, CloneCodes, np.ndarray, np.ndarray]:
    """Clones that favour one lineage, with pseudotime grouped by clone and two strata.

    ``by_depth``: each cell's lineage is its pseudotime quantile instead, so clones sit at one
    depth across lineages. Returns membership, codes (two donors, one strata column), pseudotime
    and each cell's lineage.
    """
    rng = np.random.default_rng(seed)
    clone = np.repeat(np.arange(n_clones), size)
    n = clone.size
    lineage = np.where(
        rng.random(n) < bias, rng.integers(n_lineages, size=n_clones)[clone], rng.integers(n_lineages, size=n)
    )
    pseudotime = rng.random(n_clones)[clone] + 0.2 * rng.random(n)
    if by_depth:
        lineage = np.digitize(pseudotime, np.quantile(pseudotime, np.linspace(0, 1, n_lineages + 1)[1:-1]))
    membership = np.eye(n_lineages, dtype=bool)[lineage].T
    c = codes(
        [f"c{k}" for k in clone],
        [f"d{k % 2}" for k in clone],
        [rng.choice(["t1", "t2"], size=n).tolist()],
    )
    return membership, c, pseudotime, lineage


def floats(r: ClonalResult) -> dict[str, float]:
    return {f.name: getattr(r, f.name) for f in fields(r) if isinstance(getattr(r, f.name), float)}


def test_entropy_purity_and_diagnostics_by_hand() -> None:
    # L=2, one donor. Clone A: two cells in lineage 0. Clone B: one cell in each lineage.
    r = score_clonal(
        own([[1, 1, 1, 0], [0, 0, 0, 1]]),
        codes(["A", "A", "B", "B"]),
        PT0,
        BASE0,
        min_clones=1,
        n_permutations=99,
    )
    assert r.mean_entropy == pytest.approx(0.5) and r.purity == pytest.approx(0.5)
    assert r.clone_weighted_purity == pytest.approx(0.5)
    assert r.purity_without_largest == pytest.approx(0.0)  # tie on size: the purer clone (A) is left out
    assert r.effective_clones == pytest.approx(2.0) and (r.n_cells, r.n_clones) == (4, 2)
    assert (r.status, r.mode, r.n_lineages, r.n_donors, r.n_strata) == ("ok", "own_cells", 2, 1, None)
    assert (r.min_clone_size, r.n_permutations, r.seed) == (2, 99, 0)


def test_purity_diagnostics_weigh_clones_equally_and_leave_out_the_largest() -> None:
    # Clone A: three cells in lineage 0 (H=0). Clone B: one cell in each lineage (H=1).
    r = score_clonal(
        own([[1, 1, 1, 1, 0], [0, 0, 0, 0, 1]]),
        codes(["A", "A", "A", "B", "B"]),
        np.zeros(5),
        np.zeros(5, dtype=np.int64),
        min_clones=1,
        n_permutations=9,
    )
    assert r.purity == pytest.approx(1 - 2 / 5)
    assert r.clone_weighted_purity == pytest.approx(1 - (0 + 1) / 2)  # each clone weighs the same
    assert r.purity_without_largest == pytest.approx(0.0)
    assert r.effective_clones == pytest.approx(1 / ((3 / 5) ** 2 + (2 / 5) ** 2))


def test_oracle_cells_use_fractional_counts() -> None:
    # cells: c0 in both method lineages, c1 in lineage 0, c2 in none, c3 in lineage 1;
    # clones X={c0,c1}, Y={c2,c3}
    r = score_clonal(
        np.array([[1, 1, 0, 0], [1, 0, 0, 1]], bool),
        codes(["X", "X", "Y", "Y"]),
        PT0,
        BASE0,
        oracle_cells=np.ones(4, bool),
        min_clones=1,
        n_permutations=99,
    )
    assert r.mode == "oracle_cells" and r.mean_entropy == pytest.approx(0.811278, abs=1e-6)
    assert r.n_cells == 4 and r.purity == pytest.approx(1 - 0.811278, abs=1e-6)


def test_oracle_cells_outside_the_oracle_set_are_not_scored() -> None:
    r = score_clonal(
        own([[1, 1, 1, 0], [0, 0, 0, 1]]),
        codes(["A", "A", "B", "B"]),
        PT0,
        BASE0,
        oracle_cells=np.array([True, True, True, False]),
        min_clones=1,
        n_permutations=9,
    )
    # B keeps one oracle cell, so only A (two cells, pure) is scored.
    assert (r.n_cells, r.n_clones) == (2, 1) and r.purity == pytest.approx(1.0)


def test_own_cells_mode_scores_only_single_lineage_cells() -> None:
    # Clone A: c0, c1 in lineage 0 and c2 in both. Clone B: c3, c4 in lineage 1 and c5 in none.
    membership = own([[1, 1, 1, 0, 0, 0], [0, 0, 1, 1, 1, 0]])
    c = codes(["A", "A", "A", "B", "B", "B"])
    r = score_clonal(membership, c, np.zeros(6), np.zeros(6, np.int64), min_clones=1, n_permutations=9)
    assert (r.n_cells, r.n_clones) == (4, 2) and r.purity == pytest.approx(1.0)
    # m counts the candidate cells of a clone: with m = 3 neither clone has enough.
    r3 = score_clonal(
        membership, c, np.zeros(6), np.zeros(6, np.int64), min_clone_size=3, min_clones=1, n_permutations=9
    )
    assert (r3.n_cells, r3.n_clones) == (0, 0)


def test_null_mean_and_p_match_enumeration() -> None:
    # 4 cells (L0,L0,L1,L1), clones {0,1},{2,3}: 2 of the 6 equally likely pairings give H=0, 4 give H=1.
    r = score_clonal(
        own([[1, 1, 0, 0], [0, 0, 1, 1]]),
        codes(["A", "A", "B", "B"]),
        PT0,
        BASE0,
        min_clones=1,
        n_permutations=9999,
    )
    n1 = r.nulls[0]
    assert n1.null == "N1" and n1.mean_entropy_mean == pytest.approx(2 / 3, abs=0.01)
    assert r.eta == 1.0 and n1.p == pytest.approx(1 / 3, abs=0.02)
    # H_b is 0 or 1 with P(1) = 2/3, so its sd is sqrt(2/9).
    assert n1.mean_entropy_sd == pytest.approx(math.sqrt(2 / 9), abs=0.01)


def test_p_is_never_below_one_over_b_plus_one() -> None:
    # Every clone confined to one lineage: almost no shuffle reaches H = 0, so p sits at its floor.
    membership, c, pseudotime, lineage = biased(n_clones=20, size=4, bias=1.0)
    b = 99
    r = score_clonal(membership, c, pseudotime, lineage[::-1].copy(), n_permutations=b)
    assert r.status == "ok" and r.mean_entropy == 0.0
    assert [n.null for n in r.nulls] == ["N1", "N2", "N3"]
    assert all(n.p >= 1 / (b + 1) for n in r.nulls)
    assert r.nulls[0].p == pytest.approx(1 / (b + 1))
    assert r.p >= 1 / (b + 1) and r.p_delta >= 1 / (b + 1)
    assert all(cp.p >= 1 / (b + 1) for cp in r.coupling)


def test_headline_p_is_the_largest_null_p_and_eta_is_against_n1() -> None:
    # Lineage follows depth, and clones sit at one depth: N1 and N3 call it significant, N2 does not.
    membership, c, pseudotime, lineage = biased(by_depth=True)
    r = score_clonal(membership, c, pseudotime, lineage, n_permutations=199)
    n1, n2, n3 = r.nulls
    assert n1.p == n3.p == pytest.approx(1 / 200) and n2.p == 1.0
    assert r.p == max(n.p for n in r.nulls) == n2.p
    assert r.eta == n1.eta and r.eta > 0.5
    assert n2.eta == pytest.approx(0.0, abs=1e-12) and n3.eta != n1.eta


def test_n3_exists_only_with_strata() -> None:
    membership, c, pseudotime, lineage = biased()
    plain = CloneCodes(c.clone, c.donor, None, c.n_missing_donor)
    r = score_clonal(membership, plain, pseudotime, lineage, n_permutations=19)
    assert [n.null for n in r.nulls] == ["N1", "N2"] and r.n_strata is None
    r = score_clonal(membership, c, pseudotime, lineage, n_permutations=19)
    assert [n.null for n in r.nulls] == ["N1", "N2", "N3"]
    assert r.n_strata == 4  # two donors x two strata values


def test_baseline_equal_to_the_lineages_gives_zero_delta() -> None:
    membership, c, pseudotime, lineage = biased()
    r = score_clonal(membership, c, pseudotime, lineage, n_permutations=199)
    assert math.isfinite(r.eta)
    assert r.delta_eta == pytest.approx(0.0, abs=1e-12) and r.eta_base == pytest.approx(r.eta)
    # Every paired shuffle gives the same difference, 0, so none is below the observed one.
    assert r.p_delta == 1.0


def test_a_constant_baseline_leaves_the_baseline_comparison_undefined() -> None:
    membership, c, pseudotime, _ = biased()
    r = score_clonal(membership, c, pseudotime, np.zeros(membership.shape[1], np.int64), n_permutations=99)
    assert math.isfinite(r.eta) and math.isfinite(r.p)
    assert math.isnan(r.eta_base) and math.isnan(r.delta_eta) and math.isnan(r.p_delta)


def test_a_partition_beyond_the_lineages_is_reported_against_the_same_shuffles() -> None:
    membership, c, pseudotime, lineage = biased()
    mixed = np.random.default_rng(1).integers(2, size=lineage.size)
    r = score_clonal(membership, c, pseudotime, mixed, n_permutations=499)
    assert r.eta > 0.1 and abs(r.eta_base) < 0.1
    assert r.delta_eta == pytest.approx(r.eta - r.eta_base)
    assert r.p_delta == pytest.approx(1 / 500)


def test_one_lineage_is_single_lineage_with_nan_values() -> None:
    r = score_clonal(own([[1, 1, 1, 1]]), codes(["A", "A", "B", "B"]), PT0, BASE0, n_permutations=9)
    assert r.status == "single_lineage" and math.isnan(r.purity) and r.nulls == ()
    assert r.coupling == ()
    assert (r.n_lineages, r.n_cells, r.n_clones, r.n_donors) == (1, 4, 2, 1)
    assert all(math.isnan(v) for v in floats(r).values())


def test_too_few_clones_still_reports_values() -> None:
    membership, c, pseudotime, lineage = biased(n_clones=5, size=4)
    r = score_clonal(membership, c, pseudotime, lineage, min_clones=20, n_permutations=19)
    assert r.status == "too_few_clones" and r.n_clones == 5
    assert math.isfinite(r.purity) and math.isfinite(r.p) and len(r.nulls) == 3
    assert len(r.coupling) == 1


def test_no_scored_clone_leaves_every_value_undefined() -> None:
    r = score_clonal(own([[1, 1, 0, 0], [0, 0, 1, 1]]), codes(["A", "B", "C", None]), PT0, BASE0)
    assert r.status == "too_few_clones" and (r.n_cells, r.n_clones, r.n_donors) == (0, 0, 0)
    assert r.nulls == () and r.coupling == ()
    assert all(math.isnan(v) for v in floats(r).values())


def test_one_clone_has_no_purity_without_the_largest() -> None:
    r = score_clonal(own([[1, 1, 0], [0, 0, 1]]), codes(["A", "A", "A"]), np.zeros(3), np.zeros(3, np.int64))
    assert r.n_clones == 1 and math.isfinite(r.purity)
    assert math.isnan(r.purity_without_largest)
    # One clone: every shuffle is the observation, so the null has no spread.
    n1 = r.nulls[0]
    assert n1.mean_entropy_mean == pytest.approx(r.mean_entropy) and n1.eta == pytest.approx(0.0)
    assert n1.p == 1.0
    (cp,) = r.coupling
    assert cp.observed == pytest.approx(cp.null_mean) and math.isnan(cp.z) and cp.p == 1.0


def test_coupling_by_hand() -> None:
    # same data as the enumeration test: X = [[2,0],[0,2]]; null cosine is 1 in 4/6 pairings, 0 in 2/6
    r = score_clonal(
        own([[1, 1, 0, 0], [0, 0, 1, 1]]),
        codes(["A", "A", "B", "B"]),
        PT0,
        BASE0,
        min_clones=1,
        n_permutations=9999,
    )
    (c,) = r.coupling
    assert (
        (c.a, c.b, c.observed) == (0, 1, 0.0) and c.null_mean == pytest.approx(2 / 3, abs=0.01) and c.p == 1.0
    )
    # null sd is sqrt(2/9), so z = (0 - 2/3) / sqrt(2/9) = -sqrt(2)
    assert c.z == pytest.approx(-math.sqrt(2), abs=0.05) and c.p_adj == c.p


def test_a_lineage_without_scored_cells_has_undefined_coupling() -> None:
    r = score_clonal(
        own([[1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 0, 0]]),
        codes(["A", "A", "B", "B"]),
        PT0,
        BASE0,
        min_clones=1,
        n_permutations=99,
    )
    pairs = {(cp.a, cp.b): cp for cp in r.coupling}
    assert list(pairs) == [(0, 1), (0, 2), (1, 2)]
    for key in [(0, 2), (1, 2)]:
        cp = pairs[key]
        assert all(math.isnan(v) for v in (cp.observed, cp.null_mean, cp.z, cp.p, cp.p_adj))
    # The adjustment runs over the one defined pair, so it leaves that p unchanged.
    assert pairs[(0, 1)].p_adj == pairs[(0, 1)].p == 1.0


def test_coupling_p_adj_is_benjamini_hochberg_over_the_pairs() -> None:
    membership, c, pseudotime, lineage = biased(n_lineages=3, bias=0.5)
    r = score_clonal(membership, c, pseudotime, lineage, n_permutations=199)
    assert [(cp.a, cp.b) for cp in r.coupling] == [(0, 1), (0, 2), (1, 2)]
    p = np.array([cp.p for cp in r.coupling])
    assert [cp.p_adj for cp in r.coupling] == pytest.approx(false_discovery_control(p, method="bh").tolist())


def test_cells_without_pseudotime_form_their_own_n2_stratum() -> None:
    # c0 (A, L0) and c1 (B, L1) have no pseudotime; c2 (A, L0) and c3 (B, L1) sit in deciles 0 and 5.
    # N2 can only swap c0 and c1, so H is 1 or 0 with even odds. N1 shuffles all four (mean 2/3).
    r = score_clonal(
        own([[1, 0, 1, 0], [0, 1, 0, 1]]),
        codes(["A", "B", "A", "B"]),
        np.array([NAN, NAN, 0.0, 1.0]),
        BASE0,
        min_clones=1,
        n_permutations=9999,
    )
    n1, n2 = r.nulls
    assert n1.mean_entropy_mean == pytest.approx(2 / 3, abs=0.02)
    assert n2.null == "N2" and n2.mean_entropy_mean == pytest.approx(0.5, abs=0.02)
    assert n2.p == pytest.approx(0.5, abs=0.02)


def test_pseudotime_deciles_use_average_ranks_over_finite_values() -> None:
    pt = np.array([0.3, 0.1, NAN, 0.3, 0.9])
    # average ranks of the four finite values: 2.5, 1, 2.5, 4
    assert _pseudotime_bins(pt, 10).tolist() == [3, 0, 10, 3, 7]
    assert _pseudotime_bins(pt, 2).tolist() == [0, 0, 2, 0, 1]


def test_n2_shuffles_within_donor_and_decile() -> None:
    # Each clone sits in its own decile (pseudotime is constant within a clone), so N2 cannot move
    # any clone id: every N2 shuffle equals the observation.
    membership, c, _, lineage = biased(n_clones=10, size=4, bias=0.5)
    pseudotime = np.repeat(np.arange(10.0), 4)
    r = score_clonal(membership, c, pseudotime, lineage, min_clones=1, n_permutations=49)
    n2 = r.nulls[1]
    assert n2.mean_entropy_mean == pytest.approx(r.mean_entropy) and n2.mean_entropy_sd == pytest.approx(0)
    assert n2.p == 1.0 and r.nulls[0].mean_entropy_sd > 0


def test_n1_shuffles_only_within_donor() -> None:
    # Donor d1's clones are all in lineage 0 and d2's in lineage 1. Shuffling within donor keeps
    # every clone pure; shuffling across donors would not.
    membership = own([[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]])
    c = codes(list("AABBCCDD"), ["d1"] * 4 + ["d2"] * 4)
    r = score_clonal(membership, c, np.zeros(8), np.zeros(8, np.int64), min_clones=1, n_permutations=99)
    n1 = r.nulls[0]
    assert n1.mean_entropy_mean == 0.0 and math.isnan(n1.eta) and n1.p == 1.0
    assert math.isnan(r.eta) and r.n_donors == 2


def test_scores_depend_on_the_seed_only_through_the_nulls() -> None:
    membership, c, pseudotime, lineage = biased(bias=0.5)
    a = score_clonal(membership, c, pseudotime, lineage, n_permutations=99, seed=0)
    b = score_clonal(membership, c, pseudotime, lineage, n_permutations=99, seed=0)
    other = score_clonal(membership, c, pseudotime, lineage, n_permutations=99, seed=1)
    assert a.nulls == b.nulls and a.coupling == b.coupling and a.p_delta == b.p_delta
    assert other.seed == 1 and other.purity == a.purity
    assert other.nulls[0].mean_entropy_mean != a.nulls[0].mean_entropy_mean


def test_invalid_inputs_are_refused() -> None:
    m, c = own([[1, 1, 0, 0], [0, 0, 1, 1]]), codes(["A", "A", "B", "B"])
    with pytest.raises(ValueError, match="min_clone_size"):
        score_clonal(m, c, PT0, BASE0, min_clone_size=1)
    with pytest.raises(ValueError, match="n_permutations"):
        score_clonal(m, c, PT0, BASE0, n_permutations=0)
    with pytest.raises(ValueError, match="pseudotime"):
        score_clonal(m, c, np.zeros(3), BASE0)
    with pytest.raises(ValueError, match="membership"):
        score_clonal(m[:, :3], c, PT0, BASE0)


def test_kmeans_partition_is_seeded() -> None:
    space = np.random.default_rng(0).normal(size=(200, 5))
    a = kmeans_partition(space, 3, seed=0)
    b = kmeans_partition(space, 3, seed=0)
    assert a.dtype == np.int64 and a.shape == (200,)
    np.testing.assert_array_equal(a, b)
    assert set(a.tolist()) <= {0, 1, 2}


def test_baseline_space_takes_the_leading_pcs() -> None:
    pca = np.arange(12, dtype=np.float32).reshape(3, 4)
    space = baseline_space(pca, 2)
    assert space.dtype == np.float64 and space.shape == (3, 2)
    np.testing.assert_array_equal(space, pca[:, :2])
    assert baseline_space(pca, 50).shape == (3, 4)


def test_sources_are_cited() -> None:
    doc = clonal.__doc__ or ""
    for doi in (
        "10.1038/s41586-018-0694-x",  # STARTRAC
        "10.1371/journal.pcbi.1011300",  # LRT
        "10.1038/s41587-022-01209-1",  # CoSpar
        "10.2202/1544-6115.1585",  # Phipson & Smyth
        "10.1093/bioinformatics/btaa611",  # scirpy
        "10.1038/s41587-021-00989-2",  # CoNGA
    ):
        assert doi in doc, doi
    assert "Vinh" in doc and "11" in doc and "2837" in doc
