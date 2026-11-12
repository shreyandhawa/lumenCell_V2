"""Properties of the clonal scores (spec §4.5), over random small datasets.

Relabelling clones or lineages and adding single-cell clones leave every field unchanged, exactly.
Reordering cells changes which cell each random key lands on, so it leaves the deterministic
fields unchanged and the null-based fields within Monte Carlo error (see HANDOFF, open questions).
"""

import math
from dataclasses import dataclass, fields, is_dataclass, replace

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from lumencell.clonal.labels import encode_clones
from lumencell.clonal.scores import ClonalResult, score_clonal

B = 499
ALPHABET = [f"clone{i}" for i in range(12)]
STRATA = ["t1", "t2", None]
DETERMINISTIC = (
    "n_cells",
    "n_clones",
    "n_donors",
    "n_strata",
    "effective_clones",
    "purity",
    "mean_entropy",
    "clone_weighted_purity",
    "purity_without_largest",
)


@dataclass(frozen=True)
class Case:
    membership: np.ndarray  # (L, n) bool
    labels: np.ndarray  # (n,) object, None where a cell has no receptor
    donors: np.ndarray  # (n,) object, None where the donor is missing
    pseudotime: np.ndarray  # (n,) float, some NaN
    baseline: np.ndarray  # (n,) int in range(L)
    strata: np.ndarray | None  # (n,) object, or no strata
    oracle_cells: np.ndarray | None  # (n,) bool, or own-cells mode

    @property
    def n_lineages(self) -> int:
        return int(self.membership.shape[0])

    def take(self, order: np.ndarray) -> "Case":
        return Case(
            self.membership[:, order],
            self.labels[order],
            self.donors[order],
            self.pseudotime[order],
            self.baseline[order],
            None if self.strata is None else self.strata[order],
            None if self.oracle_cells is None else self.oracle_cells[order],
        )


def placement(place: np.ndarray, n_lineages: int) -> np.ndarray:
    """0..L-1: that lineage alone; L: every lineage; L+1: none; L+2: lineages 0 and 1."""
    membership = np.zeros((n_lineages, place.size), dtype=bool)
    single = place < n_lineages
    membership[place[single], np.flatnonzero(single)] = True
    membership[:, place == n_lineages] = True
    membership[:2, place == n_lineages + 2] = True
    return membership


def objects(values: list[object]) -> np.ndarray:
    return np.array(values, dtype=object)


@st.composite
def cases(draw: st.DrawFn) -> Case:
    n_lineages = draw(st.integers(2, 4))
    n = draw(st.integers(20, 120))
    donor_names = [f"d{i}" for i in range(draw(st.integers(1, 3)))]
    labels = draw(st.lists(st.sampled_from([*ALPHABET, None]), min_size=n, max_size=n))
    donors = draw(st.lists(st.sampled_from([*donor_names * 4, None]), min_size=n, max_size=n))
    place = draw(arrays(np.int64, n, elements=st.integers(0, n_lineages + 2)))
    pseudotime = draw(arrays(np.float64, n, elements=st.one_of(st.floats(0, 1), st.just(math.nan))))
    baseline = draw(arrays(np.int64, n, elements=st.integers(0, n_lineages - 1)))
    strata = draw(st.one_of(st.none(), st.lists(st.sampled_from(STRATA), min_size=n, max_size=n)))
    oracle = draw(st.one_of(st.none(), arrays(np.bool_, n)))
    return Case(
        placement(place, n_lineages),
        objects(labels),
        objects(donors),
        pseudotime,
        baseline,
        None if strata is None else objects(strata),
        oracle,
    )


def score(case: Case) -> ClonalResult:
    codes = encode_clones(case.labels, case.donors, () if case.strata is None else [case.strata])
    return score_clonal(
        case.membership,
        codes,
        case.pseudotime,
        case.baseline,
        oracle_cells=case.oracle_cells,
        min_clones=1,
        n_permutations=B,
    )


def differences(a: object, b: object, path: str = "result") -> list[str]:
    """Where two results differ: ints and strings exactly, floats within 1e-12, NaN equal to NaN."""
    if is_dataclass(a) and not isinstance(a, type):
        if type(a) is not type(b):
            return [path]
        return [
            d
            for f in fields(a)
            for d in differences(getattr(a, f.name), getattr(b, f.name), f"{path}.{f.name}")
        ]
    if isinstance(a, tuple):
        if not isinstance(b, tuple) or len(a) != len(b):
            return [path]
        return [
            d for i, (x, y) in enumerate(zip(a, b, strict=True)) for d in differences(x, y, f"{path}[{i}]")
        ]
    if isinstance(a, float) or isinstance(b, float):
        if not (isinstance(a, float) and isinstance(b, float)):
            return [path]
        if math.isnan(a) or math.isnan(b):
            return [] if math.isnan(a) and math.isnan(b) else [f"{path}: {a} != {b}"]
        return [] if math.isclose(a, b, rel_tol=0.0, abs_tol=1e-12) else [f"{path}: {a} != {b}"]
    return [] if type(a) is type(b) and a == b else [f"{path}: {a!r} != {b!r}"]


def same(a: ClonalResult, b: ClonalResult) -> bool:
    diff = differences(a, b)
    assert not diff, diff
    return True


@settings(max_examples=30, deadline=None)
@given(cases(), st.data())
def test_relabelling_clones_or_lineages_changes_nothing(case: Case, data: st.DataObject) -> None:
    rename = dict(zip(ALPHABET, data.draw(st.permutations(ALPHABET)), strict=True))
    order = np.array(data.draw(st.permutations(range(case.n_lineages))))
    relabelled = replace(
        case,
        labels=objects([None if v is None else rename[v] for v in case.labels]),
        membership=case.membership[order],
    )
    before, after = score(case), score(relabelled)
    # Old lineage order[j] is new lineage j, so the old pair (a, b) is the new pair (new[a], new[b]).
    new = np.argsort(order)
    mapped = sorted(
        (replace(c, a=int(min(new[c.a], new[c.b])), b=int(max(new[c.a], new[c.b]))) for c in before.coupling),
        key=lambda c: (c.a, c.b),
    )
    assert same(replace(before, coupling=tuple(mapped)), after)


@settings(max_examples=30, deadline=None)
@given(cases(), st.data())
def test_reordering_cells_keeps_deterministic_fields_and_null_fields_within_monte_carlo_error(
    case: Case, data: st.DataObject
) -> None:
    order = np.array(data.draw(st.permutations(range(case.labels.size))))
    a, b = score(case), score(case.take(order))
    for name in DETERMINISTIC:
        assert not differences(getattr(a, name), getattr(b, name), name)
    assert [n.null for n in a.nulls] == [n.null for n in b.nulls]
    for x, y in zip(a.nulls, b.nulls, strict=True):
        assert abs(x.p - y.p) <= 5 * math.sqrt(0.5 / B) + 2 / (B + 1)
        if math.isnan(x.eta) or math.isnan(y.eta):
            assert math.isnan(x.eta) and math.isnan(y.eta)
            continue
        bound = max(
            5 * math.sqrt(2 / B) * a.mean_entropy * n.mean_entropy_sd / n.mean_entropy_mean**2 for n in (x, y)
        )
        assert abs(x.eta - y.eta) <= bound + 1e-9


@settings(max_examples=30, deadline=None)
@given(cases(), st.data())
def test_singleton_clones_never_change_a_score(case: Case, data: st.DataObject) -> None:
    k = data.draw(st.integers(1, 10))
    n_lineages = case.n_lineages
    extra = Case(
        data.draw(arrays(np.bool_, (n_lineages, k))),
        objects([f"single{i}" for i in range(k)]),
        objects(data.draw(st.lists(st.sampled_from(["d0", "d1", "dX", None]), min_size=k, max_size=k))),
        data.draw(arrays(np.float64, k, elements=st.one_of(st.floats(-5, 5), st.just(math.nan)))),
        data.draw(arrays(np.int64, k, elements=st.integers(0, 9))),
        None
        if case.strata is None
        else objects(data.draw(st.lists(st.sampled_from(["t1", "t3", None]), min_size=k, max_size=k))),
        None if case.oracle_cells is None else data.draw(arrays(np.bool_, k)),
    )

    def join(x: np.ndarray | None, y: np.ndarray | None, axis: int = 0) -> np.ndarray | None:
        return None if x is None or y is None else np.concatenate([x, y], axis=axis)

    grown = Case(
        np.concatenate([case.membership, extra.membership], axis=1),
        np.concatenate([case.labels, extra.labels]),
        np.concatenate([case.donors, extra.donors]),
        np.concatenate([case.pseudotime, extra.pseudotime]),
        np.concatenate([case.baseline, extra.baseline]),
        join(case.strata, extra.strata),
        join(case.oracle_cells, extra.oracle_cells),
    )
    assert same(score(case), score(grown))


@settings(max_examples=30, deadline=None)
@given(cases())
def test_purities_lie_in_the_unit_interval(case: Case) -> None:
    r = score(case)
    for value in (r.purity, r.clone_weighted_purity, r.purity_without_largest, r.mean_entropy):
        assert math.isnan(value) or 0.0 <= value <= 1.0
    for n in r.nulls:
        assert 1 / (B + 1) <= n.p <= 1.0


@settings(max_examples=30, deadline=None)
@given(cases())
def test_splitting_clones_across_donors_never_merges_them(case: Case) -> None:
    codes = encode_clones(case.labels, case.donors)
    pairs = [
        (d, lab) for d, lab in zip(case.donors, case.labels, strict=True) if d is not None and lab is not None
    ]
    ids = codes.clone[codes.clone >= 0]
    assert len(set(ids.tolist())) == len(set(pairs))
    # One id per (donor, label) pair, and one pair per id.
    assert len(set(zip(pairs, ids.tolist(), strict=True))) == len(set(pairs))


@settings(max_examples=30, deadline=None)
@given(st.lists(st.tuples(st.integers(2, 5), st.integers(2, 5)), min_size=1, max_size=5), cases())
def test_ignoring_donors_changes_the_score_when_labels_collide(
    sizes: list[tuple[int, int]], other: Case
) -> None:
    # Each shared label: donor A's copies all in lineage 0, donor B's all in lineage 1. The cells of
    # another case join in donor A, under labels no other donor uses; own-cells mode.
    labels: list[object] = []
    donors: list[object] = []
    place: list[int] = []
    for i, (in_a, in_b) in enumerate(sizes):
        labels += [f"shared{i}"] * (in_a + in_b)
        donors += ["A"] * in_a + ["B"] * in_b
        place += [0] * in_a + [1] * in_b
    membership = np.concatenate([placement(np.array(place), other.n_lineages), other.membership], axis=1)
    all_labels = objects([*labels, *(None if v is None else f"A-{v}" for v in other.labels)])
    all_donors = objects([*donors, *["A"] * other.labels.size])
    n = all_labels.size

    def purity_and_cells(donor_column: np.ndarray | None) -> tuple[float, int]:
        r = score_clonal(
            membership,
            encode_clones(all_labels, donor_column),
            np.zeros(n),
            np.zeros(n, np.int64),
            min_clones=1,
            n_permutations=B,
        )
        return r.purity, r.n_cells

    (with_donors, cells), (without, cells_without) = purity_and_cells(all_donors), purity_and_cells(None)
    assert cells == cells_without
    assert with_donors > without
