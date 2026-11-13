import dataclasses
import hashlib
import itertools
import json
from collections.abc import Callable

import anndata as ad
import numpy as np
import pytest
from scipy.stats import zipf

from lumencell.bench.gates import load_thresholds
from lumencell.bench.synth import generate, truth_from_adata
from lumencell.bench.synth_clones import (
    DEFAULT_LAYER,
    DONOR_CONFOUND_SHARE,
    FATE_BIAS,
    LABEL_RATE,
    N_DONORS,
    SINGLETON_SHARE,
    ZIPF_A,
    ClonalLayer,
    add_clones,
    clone_truth,
    draw_clone_size,
)
from lumencell.bench.truth import Truth
from lumencell.clonal.labels import encode_clones
from lumencell.clonal.scores import score_clonal


def _digest(adata: ad.AnnData) -> str:
    h = hashlib.sha256()
    for a in (adata.X.data, adata.X.indices, adata.X.indptr):
        h.update(np.ascontiguousarray(a).tobytes())
    h.update(adata.obs.to_csv().encode())
    h.update(json.dumps(adata.uns["lumencell_truth"], sort_keys=True, default=str).encode())
    return h.hexdigest()


def _data(topology: str = "bifurcation", n_cells: int = 3000, seed: int = 0) -> ad.AnnData:
    return generate(topology, "realistic", n_cells=n_cells, n_genes=50, seed=seed)


def _leaves(adata: ad.AnnData) -> list[str]:
    segments = adata.uns["lumencell_truth"]["segments"]
    parents = {v["parent"] for v in segments.values()}
    return sorted(s for s in segments if s not in parents)


def _clone_keys(adata: ad.AnnData) -> list[tuple[str, str]]:
    labels = adata.obs["clone_id"].to_numpy(dtype=object)
    donors = adata.obs["donor"].astype(str).to_numpy()
    return [(d, c) for d, c in zip(donors, labels, strict=True) if isinstance(c, str)]


def _sizes(adata: ad.AnnData) -> np.ndarray:
    _, counts = np.unique([f"{d}|{c}" for d, c in _clone_keys(adata)], return_counts=True)
    return counts


def test_generator_output_is_unchanged_byte_for_byte() -> None:
    # Computed at 142a788. If only Linux CI disagrees, stop and report; never edit this value.
    a = generate("binary2", "realistic", n_cells=400, n_genes=120, seed=7)
    assert _digest(a) == "87d55724becfbc9790e0b8697e42f39077d4c14feed972da02a1f41d45e78e3e"
    before = _digest(a)
    add_clones(a, DEFAULT_LAYER, seed=7)
    assert _digest(a) == before


def test_the_layer_adds_two_columns_and_its_spec_only() -> None:
    a = _data(n_cells=500)
    b = add_clones(a, DEFAULT_LAYER, seed=1)
    assert (b.X != a.X).nnz == 0
    assert b.var.equals(a.var)
    assert list(b.obs.columns) == [*a.obs.columns, "clone_id", "donor"]
    assert b.obs[list(a.obs.columns)].equals(a.obs)
    spec = b.uns["lumencell_truth"]["clones"]
    assert spec["fate_bias"] == 0.7
    assert spec["seed"] == 1
    assert spec["stream"] == 7100501
    assert "clones" not in a.uns["lumencell_truth"]


def test_same_seed_same_clones() -> None:
    a = _data(n_cells=500)
    x, y, z = (add_clones(a, DEFAULT_LAYER, seed=s) for s in (3, 3, 4))
    assert x.obs["clone_id"].equals(y.obs["clone_id"]) and x.obs["donor"].equals(y.obs["donor"])
    assert not x.obs["clone_id"].equals(z.obs["clone_id"])


def test_labelled_share_and_sizes_do_not_depend_on_fate_bias() -> None:
    a = _data()
    shares, singles, p75, p90 = [], [], [], []
    for bias in (0.0, 0.35, 0.7, 1.0):
        b = add_clones(a, ClonalLayer(fate_bias=bias), seed=0)
        shares.append(b.obs["clone_id"].notna().mean())
        sizes = _sizes(b)
        singles.append((sizes == 1).mean())
        p75.append(np.percentile(sizes, 75))
        p90.append(np.percentile(sizes, 90))
    assert all(0.83 <= s <= 0.87 for s in shares)
    assert max(singles) - min(singles) <= 0.06
    assert max(p75) - min(p75) <= 1
    assert max(p90) - min(p90) <= 1


@pytest.mark.parametrize("topology", ["bifurcation", "trifurcation", "binary2"])
def test_full_fate_bias_confines_every_leaf_cells_clone_to_its_leaf(topology: str) -> None:
    b = add_clones(_data(topology), ClonalLayer(fate_bias=1.0), seed=0)
    leaves = set(_leaves(b))
    seg = b.obs["true_segment"].astype(str).to_numpy()
    keys = np.array([f"{d}|{c}" for d, c in zip(b.obs["donor"].astype(str), b.obs["clone_id"], strict=True)])
    labelled = b.obs["clone_id"].notna().to_numpy()
    checked = 0
    for key in np.unique(keys[labelled]):
        segs = set(seg[labelled & (keys == key)])
        if segs & leaves:
            assert len(segs) == 1, (key, segs)
            checked += 1
    assert checked > 50


def test_label_strings_collide_across_donors() -> None:
    b = add_clones(_data(n_cells=1000), DEFAULT_LAYER, seed=0)
    for d in b.obs["donor"].cat.categories:
        assert "c0" in set(b.obs.loc[b.obs["donor"] == d, "clone_id"].dropna())
    assert list(b.obs["donor"].cat.categories) == ["d0", "d1", "d2"]


def _clone_keys_array(b: ad.AnnData) -> np.ndarray:
    donors = b.obs["donor"].astype(str)
    labels = b.obs["clone_id"].astype(str)
    return np.array([f"{d}|{c}" for d, c in zip(donors, labels, strict=True)])


def _leaf_pure_share(b: ad.AnnData) -> float:
    """Share of clones of size >= 2 holding a leaf-segment cell whose cells all sit on one leaf."""
    leaves = set(_leaves(b))
    seg = b.obs["true_segment"].astype(str).to_numpy()
    keys = _clone_keys_array(b)
    labelled = b.obs["clone_id"].notna().to_numpy()
    total = pure = 0
    for key in np.unique(keys[labelled]):
        cells = labelled & (keys == key)
        if cells.sum() < 2 or not set(seg[cells]) & leaves:
            continue
        total += 1
        pure += len(set(seg[cells])) == 1
    return pure / total


def _clone_windows(b: ad.AnnData, min_size: int) -> list[tuple[str, float, float]]:
    """(donor, min, max) of true pseudotime per clone of at least ``min_size``, by donor and clone number."""
    frame = b.obs.loc[b.obs["clone_id"].notna()]
    rows = []
    for (donor, clone), group in frame.groupby(["donor", "clone_id"], observed=True):
        if len(group) >= min_size:
            t = group["true_pseudotime"]
            rows.append((str(donor), int(str(clone)[1:]), float(t.min()), float(t.max())))
    rows.sort()
    return [(donor, lo, hi) for donor, _, lo, hi in rows]


def test_depth_only_clones_sit_in_a_pseudotime_window_and_run_in_order() -> None:
    a = _data()
    depth = _clone_windows(add_clones(a, ClonalLayer(control="depth_only"), seed=0), 3)
    mixed = _clone_windows(add_clones(a, DEFAULT_LAYER, seed=0), 3)
    assert len(depth) > 100 and len(mixed) > 100  # clones of size >= 3: singletons have span 0
    assert np.median([hi - lo for _, lo, hi in depth]) < 0.05
    assert np.median([hi - lo for _, lo, hi in mixed]) > 0.2
    for donor in {d for d, _, _ in depth}:  # clone number rises with pseudotime within a donor
        mine = [(lo, hi) for d, lo, hi in depth if d == donor]
        assert all(prev_hi <= lo for (_, prev_hi), (lo, _) in itertools.pairwise(mine))


@pytest.mark.parametrize("topology", ["bifurcation", "trifurcation", "binary2"])
def test_fate_bias_sets_the_leaf_pure_share_of_clones(topology: str) -> None:
    a = _data(topology)
    share = {
        bias: _leaf_pure_share(add_clones(a, ClonalLayer(fate_bias=bias), seed=0))
        for bias in (0, 0.35, 0.7, 1)
    }
    assert share[0] < 0.3
    assert 0.3 < share[0.35] < 0.5
    assert 0.5 < share[0.7] < 0.72
    assert share[1] == 1.0


@pytest.mark.parametrize("topology", ["bifurcation", "trifurcation", "binary2"])
def test_donor_shares_follow_the_control(topology: str) -> None:
    a = _data(topology)
    for control in ("none", "depth_only", "donor_confound"):
        b = add_clones(a, ClonalLayer(control=control), seed=0)
        truth = truth_from_adata(b)
        names = sorted(truth.names)  # lineage j: sorted leaf order
        exclusive = truth.exclusive()
        donor = b.obs["donor"].cat.codes.to_numpy()
        for j, name in enumerate(names):
            cells = exclusive == truth.names.index(name)
            shares = np.bincount(donor[cells], minlength=3) / cells.sum()
            if control == "donor_confound":
                assert abs(shares[j % 3] - 0.8) < 0.05, (name, shares)
                assert np.all(np.delete(shares, j % 3) > 0.05), (name, shares)
            else:
                assert np.all(np.abs(shares - 1 / 3) < 0.05), (control, name, shares)
        if control == "donor_confound":
            assert _leaf_pure_share(b) < 0.5  # lineage-mixed, like "none" at fate bias 0


def test_the_stream_is_seeded_by_seed_and_clone_stream_and_drawn_in_order() -> None:
    a = _data(n_cells=600)
    b = add_clones(a, DEFAULT_LAYER, seed=5)
    r = np.random.default_rng([5, 7100501])
    labelled = r.random(600) < LABEL_RATE
    donors = r.integers(N_DONORS, size=600)
    assert (labelled == b.obs["clone_id"].notna().to_numpy()).all()
    assert (donors == b.obs["donor"].cat.codes.to_numpy()).all()


def test_clone_columns_are_pinned_byte_for_byte() -> None:
    a = generate("binary2", "realistic", n_cells=400, n_genes=120, seed=7)
    b = add_clones(a, DEFAULT_LAYER, seed=7)
    text = b.obs[["clone_id", "donor"]].astype(str).to_csv().encode()
    assert (
        hashlib.sha256(text).hexdigest() == "521259033b7e26701edf5aca09b256db602c656639a51b13248b1c49c6854af8"
    )


def test_clone_sizes_follow_the_declared_model() -> None:
    rng = np.random.default_rng(123)
    n = 20000
    sizes = np.array([draw_clone_size(rng) for _ in range(n)])

    def check(p: float, observed: float) -> None:  # within 5 binomial standard deviations
        assert abs(observed - p) < 5 * np.sqrt(p * (1 - p) / n), (p, observed)

    check(0.5, (sizes == 1).mean())
    check(0.5 * zipf.pmf(1, 2.2), (sizes == 2).mean())  # 0.5 / zeta(2.2), about 0.335
    check(0.5 * zipf.sf(4, 2.2), (sizes >= 6).mean())  # upper tail: size = 1 + Zipf >= 6
    check(0.5 * zipf.sf(19, 2.2), (sizes >= 21).mean())


def _without(a: ad.AnnData, column: str) -> ad.AnnData:
    b = a.copy()
    del b.obs[column]
    return b


def _edited(a: ad.AnnData, edit: Callable[[ad.AnnData], None]) -> ad.AnnData:
    b = a.copy()
    edit(b)
    return b


def _unknown_segment(b: ad.AnnData) -> None:
    b.obs["true_segment"] = b.obs["true_segment"].astype(str)
    b.obs.iloc[0, b.obs.columns.get_loc("true_segment")] = "zzz"


def _no_segments(b: ad.AnnData) -> None:
    del b.uns["lumencell_truth"]["segments"]


def _cycle(b: ad.AnnData) -> None:  # every segment is someone's parent: no leaf
    segments = b.uns["lumencell_truth"]["segments"]
    names = list(segments)
    for name, parent in zip(names, [*names[1:], names[0]], strict=True):
        segments[name]["parent"] = parent


def _empty_leaf(b: ad.AnnData) -> None:
    b.uns["lumencell_truth"]["segments"]["extra"] = {"parent": "trunk", "length": 1.0}


def _nan_time(b: ad.AnnData) -> None:
    b.obs["true_pseudotime"] = b.obs["true_pseudotime"].astype(float)
    b.obs.iloc[0, b.obs.columns.get_loc("true_pseudotime")] = np.nan


def _add(column: str) -> Callable[[ad.AnnData], None]:
    def edit(b: ad.AnnData) -> None:
        b.obs[column] = "x"

    return edit


@pytest.mark.parametrize(
    ("layer", "edit", "message"),
    [
        (ClonalLayer(), lambda a: _without(a, "true_segment"), "true_segment"),
        (ClonalLayer(), lambda a: _without(a, "true_pseudotime"), "true_pseudotime"),
        (ClonalLayer(), lambda a: _edited(a, _no_segments), "segments"),
        (ClonalLayer(), lambda a: _edited(a, _unknown_segment), "not in segments"),
        (ClonalLayer(), lambda a: _edited(a, _cycle), "no leaf"),
        (ClonalLayer(), lambda a: _edited(a, _empty_leaf), "no cells"),
        (ClonalLayer(), lambda a: _edited(a, _nan_time), "NaN"),
        (ClonalLayer(), lambda a: _edited(a, _add("clone_id")), "already has"),
        (ClonalLayer(), lambda a: _edited(a, _add("donor")), "already has"),
    ],
)
def test_add_clones_refuses_invalid_inputs(
    layer: ClonalLayer, edit: Callable[[ad.AnnData], ad.AnnData], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        add_clones(edit(_data(n_cells=300)), layer, seed=0)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"control": "bogus"},
        {"fate_bias": -0.1},
        {"fate_bias": 1.1},
        {"label_rate": -0.1},
        {"label_rate": 1.5},
        {"n_donors": 0},
        {"n_donors": 1, "control": "donor_confound"},
    ],
)
def test_clonal_layer_refuses_invalid_values(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        ClonalLayer(**kwargs)  # type: ignore[arg-type]


def test_clone_truth_names_a_missing_column() -> None:
    a = _data(n_cells=300)
    with pytest.raises(ValueError, match="clone_id"):
        clone_truth(truth_from_adata(a), a)
    b = add_clones(a, DEFAULT_LAYER, seed=0)
    with pytest.raises(ValueError, match="donor"):
        clone_truth(truth_from_adata(b), _without(b, "donor"))


def test_clone_truth_fills_and_checks_the_truth() -> None:
    b = add_clones(_data(n_cells=300), DEFAULT_LAYER, seed=0)
    truth = clone_truth(truth_from_adata(b), b)
    assert truth.clone_ids is not None and truth.donors is not None
    assert truth.clone_ids.shape == truth.donors.shape == (300,)
    assert list(truth.donors) == list(b.obs["donor"].astype(str))
    labelled = b.obs["clone_id"].notna().to_numpy()
    assert list(truth.clone_ids[labelled]) == list(b.obs["clone_id"].astype(str)[labelled])
    base = truth_from_adata(b)
    with pytest.raises(ValueError):
        dataclasses.replace(base, clone_ids=truth.clone_ids[:10], donors=truth.donors[:10])
    with pytest.raises(ValueError):
        Truth(base.order, base.lineages, base.root, clone_ids=truth.clone_ids)
    with pytest.raises(ValueError):
        Truth(base.order, base.lineages, base.root, donors=truth.donors)


def test_draw_clone_size_is_one_or_more_and_half_singletons() -> None:
    rng = np.random.default_rng(0)
    sizes = np.array([draw_clone_size(rng) for _ in range(4000)])
    assert sizes.min() == 1
    assert abs((sizes == 1).mean() - 0.5) < 0.04


def test_layer_defaults_equal_the_thresholds_file() -> None:
    cfg = load_thresholds()["clonal"]["synthetic"]
    assert cfg["fate_bias"] == FATE_BIAS
    assert cfg["n_donors"] == N_DONORS
    assert cfg["label_rate"] == LABEL_RATE
    assert cfg["singleton_share"] == SINGLETON_SHARE
    assert cfg["zipf_a"] == ZIPF_A
    assert cfg["donor_confound_share"] == DONOR_CONFOUND_SHARE
    assert ClonalLayer(FATE_BIAS, N_DONORS, LABEL_RATE, "none") == DEFAULT_LAYER


def test_oracle_eta_rises_from_no_bias_to_full_bias() -> None:
    a = _data(n_cells=1500)

    def eta(bias: float) -> float:
        b = add_clones(a, ClonalLayer(fate_bias=bias), seed=0)
        truth = clone_truth(truth_from_adata(b), b)
        assert truth.clone_ids is not None
        codes = encode_clones(truth.clone_ids, truth.donors)
        result = score_clonal(
            truth.membership(),
            codes,
            truth.order.astype(np.float64),
            np.zeros(truth.n_obs, dtype=np.int64),
            oracle_cells=truth.exclusive() >= 0,
            n_permutations=199,
        )
        assert result.status == "ok"
        return result.eta

    assert eta(1.0) > eta(0.0)
