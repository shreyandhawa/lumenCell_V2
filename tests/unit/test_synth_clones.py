import dataclasses
import hashlib
import json

import anndata as ad
import numpy as np
import pytest

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


def test_depth_only_clones_sit_in_a_pseudotime_window() -> None:
    b = add_clones(_data(), ClonalLayer(control="depth_only"), seed=0)
    frame = b.obs.assign(key=[f"{d}|{c}" for d, c in zip(b.obs["donor"], b.obs["clone_id"], strict=True)])
    frame = frame[b.obs["clone_id"].notna().to_numpy()]
    spans = frame.groupby("key")["true_pseudotime"].agg(lambda t: float(t.max() - t.min()))
    assert spans.median() < 0.05


@pytest.mark.parametrize("topology", ["bifurcation", "binary2"])
def test_donor_confound_skews_donor_composition(topology: str) -> None:
    a = _data(topology)
    b = add_clones(a, ClonalLayer(control="donor_confound"), seed=0)
    truth = truth_from_adata(b)
    sorted_names = sorted(truth.names)
    exclusive = truth.exclusive()
    donor = b.obs["donor"].cat.codes.to_numpy()
    labelled = b.obs["clone_id"].notna().to_numpy()
    for name in sorted_names:
        j = sorted_names.index(name)  # sorted leaf order, not Truth.names order
        cells = labelled & (exclusive == truth.names.index(name))
        assert cells.sum() > 50
        assert (donor[cells] == j % 3).mean() >= 0.6
    if topology == "binary2":
        assert truth.names != sorted_names  # a topology-order mix-up would fail above


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
