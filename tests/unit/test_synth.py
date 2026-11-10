import time
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import scanpy as sc

from lumencell.bench.synth import TOPOLOGIES, generate, lineage_paths, truth_from_adata

FORK = {"bifurcation": 0.4, "trifurcation": 0.4, "binary2": 0.3}


@pytest.fixture(scope="module")
def bif() -> ad.AnnData:
    return generate("bifurcation", "realistic", n_cells=3000, n_genes=2000, seed=0)


# ---- structure ------------------------------------------------------------------


def test_counts_only_no_embedding(bif: ad.AnnData) -> None:
    assert bif.shape == (3000, 2000)
    assert len(bif.obsm) == 0, "the embedding must come from the real pipeline"
    assert bif.X.dtype == np.float32
    assert np.all(bif.X.data == np.round(bif.X.data)), "X must hold integer counts"
    assert "iroot" not in bif.uns


def test_library_sizes(bif: ad.AnnData) -> None:
    libs = np.asarray(bif.X.sum(axis=1)).ravel()
    assert np.median(libs) == pytest.approx(5000, rel=0.1)


def test_deterministic_for_a_seed() -> None:
    a = generate("bifurcation", "realistic", n_cells=300, n_genes=200, seed=3)
    b = generate("bifurcation", "realistic", n_cells=300, n_genes=200, seed=3)
    c = generate("bifurcation", "realistic", n_cells=300, n_genes=200, seed=4)
    assert (a.X != b.X).nnz == 0
    np.testing.assert_array_equal(a.obs["true_pseudotime"], b.obs["true_pseudotime"])
    assert (a.X != c.X).nnz > 0


@pytest.mark.parametrize("topology", sorted(TOPOLOGIES))
def test_truth_is_consistent(topology: str) -> None:
    adata = generate(topology, "realistic", n_cells=2000, n_genes=100, seed=1)
    t = adata.obs["true_pseudotime"].to_numpy()
    assert t.min() >= 0.0 and t.max() <= 1.0

    truth = truth_from_adata(adata)
    assert set(truth.names) == set(lineage_paths(TOPOLOGIES[topology]))
    for name, mask in truth.lineages.items():
        assert t[mask].max() > 0.95, f"lineage {name} should reach its tip"
    assert t[truth.root].max() < 0.25, "root cells are the first half of the trunk"

    shared = adata.obs["true_lineage"].astype(str).to_numpy() == "shared"
    on_many = truth.membership().sum(axis=0) > 1
    np.testing.assert_array_equal(shared, on_many)


def test_lineage_paths() -> None:
    assert lineage_paths(TOPOLOGIES["bifurcation"]) == {"A": ["trunk", "A"], "B": ["trunk", "B"]}
    assert lineage_paths(TOPOLOGIES["binary2"]) == {
        "B": ["trunk", "B"],
        "A1": ["trunk", "A", "A1"],
        "A2": ["trunk", "A", "A2"],
    }


def test_store_root_is_opt_in() -> None:
    adata = generate("bifurcation", "easy", n_cells=200, n_genes=50, seed=0, store_root=True)
    assert adata.uns["iroot"] == int(np.argmin(adata.obs["true_pseudotime"]))


def test_h5ad_round_trip_keeps_the_truth(tmp_path: Path, bif: ad.AnnData) -> None:
    path = tmp_path / "synth.h5ad"
    bif.write_h5ad(path)
    back = ad.read_h5ad(path)
    t1, t2 = truth_from_adata(bif), truth_from_adata(back)
    assert t1.names == t2.names
    np.testing.assert_array_equal(t1.membership(), t2.membership())
    np.testing.assert_array_equal(t1.root, t2.root)


def test_unknown_topology() -> None:
    with pytest.raises(ValueError, match="unknown topology"):
        generate("spiral")


def test_generation_is_fast() -> None:
    start = time.perf_counter()
    generate("binary2", "realistic", n_cells=3000, n_genes=2000, seed=0)
    assert time.perf_counter() - start < 5.0


# ---- calibration: "realistic" must be learnable but not trivial ------------------


def _calibration(adata: ad.AnnData, fork: float) -> tuple[float, float]:
    b = adata.copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    sc.pp.highly_variable_genes(b, n_top_genes=1000, flavor="seurat")
    sc.pp.pca(b, n_comps=30, mask_var="highly_variable")
    sc.pp.neighbors(b, n_neighbors=15, n_pcs=30)
    knn = b.obsp["distances"].tocsr()
    t = adata.obs["true_pseudotime"].to_numpy()
    lineage = adata.obs["true_lineage"].astype(str).to_numpy()
    deep = (lineage != "shared") & (t > fork + 0.15)
    purity = np.mean([np.mean(lineage[knn[i].indices] == lineage[i]) for i in np.flatnonzero(deep)])
    return float(b.uns["pca"]["variance_ratio"][0]), float(purity)


def _check_calibration(topology: str) -> None:
    pc1, purity = _calibration(generate(topology, "realistic", seed=0), FORK[topology])
    assert 0.03 <= pc1 <= 0.15, f"PC1 share {pc1:.3f} outside the real-data range"
    assert 0.85 <= purity < 0.99, f"branch purity {purity:.3f}: too hard or trivially separable"


def test_realistic_bifurcation_is_calibrated() -> None:
    _check_calibration("bifurcation")


@pytest.mark.slow
@pytest.mark.parametrize("topology", ["trifurcation", "binary2"])
def test_realistic_preset_is_calibrated(topology: str) -> None:
    _check_calibration(topology)
