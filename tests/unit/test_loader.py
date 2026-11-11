import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import scanpy as sc

from lumencell.bench.synth import generate
from lumencell.data.loader import DatasetError, classify_matrix, load

# The scaled-data fixtures densify on purpose.
pytestmark = pytest.mark.filterwarnings("ignore:zero-centering a sparse:UserWarning")


@pytest.fixture(scope="module")
def counts() -> ad.AnnData:
    return generate("bifurcation", "realistic", n_cells=400, n_genes=300, seed=0)


def _write(adata: ad.AnnData, tmp_path: Path, name: str = "d.h5ad") -> Path:
    path = tmp_path / name
    adata.write_h5ad(path)
    return path


def _lognorm(adata: ad.AnnData) -> ad.AnnData:
    out = adata.copy()
    sc.pp.normalize_total(out, target_sum=1e4)
    sc.pp.log1p(out)
    return out


def test_classify_matrix(counts: ad.AnnData) -> None:
    assert classify_matrix(counts.X) == "counts"
    assert classify_matrix(_lognorm(counts).X) == "lognorm"
    scaled = _lognorm(counts)
    sc.pp.scale(scaled)
    assert classify_matrix(scaled.X) == "scaled"
    cpm = counts.copy()
    sc.pp.normalize_total(cpm, target_sum=1e6)
    assert classify_matrix(cpm.X) is None


def test_counts_in_x(counts: ad.AnnData, tmp_path: Path) -> None:
    ds = load(_write(counts, tmp_path))
    assert (ds.x_kind, ds.source) == ("counts", "X")
    assert ds.matrix.shape == counts.shape
    assert ds.gene_symbol_column is None


def test_lognorm_x_without_counts(counts: ad.AnnData, tmp_path: Path) -> None:
    ds = load(_write(_lognorm(counts), tmp_path))
    assert (ds.x_kind, ds.source) == ("lognorm", "X")


def test_scaled_x_uses_counts_layer(counts: ad.AnnData, tmp_path: Path) -> None:
    scaled = _lognorm(counts)
    scaled.layers["counts"] = counts.X.copy()
    sc.pp.scale(scaled)
    ds = load(_write(scaled, tmp_path))
    assert (ds.x_kind, ds.source) == ("counts", "layers/counts")
    assert ds.original_x_kind == "scaled"


def test_hvg_subset_x_uses_full_raw_counts(counts: ad.AnnData, tmp_path: Path) -> None:
    full = counts.copy()
    full.raw = full  # raw keeps all 300 genes as counts
    subset = _lognorm(full)[:, :100].copy()
    sc.pp.scale(subset)
    ds = load(_write(subset, tmp_path))
    assert ds.source == "raw"
    assert ds.n_vars == 300


def test_scaled_without_counts_is_rejected(counts: ad.AnnData, tmp_path: Path) -> None:
    scaled = _lognorm(counts)
    sc.pp.scale(scaled)
    with pytest.raises(DatasetError, match="scaled"):
        load(_write(scaled, tmp_path))


def test_normalised_but_not_logged_is_rejected_with_a_hint(counts: ad.AnnData, tmp_path: Path) -> None:
    cpm = counts.copy()
    sc.pp.normalize_total(cpm, target_sum=1e6)
    with pytest.raises(DatasetError, match="--x-is"):
        load(_write(cpm, tmp_path))


def test_ensembl_ids_use_the_symbol_column(counts: ad.AnnData, tmp_path: Path) -> None:
    adata = counts.copy()
    adata.var["feature_name"] = [f"GENE{i}" for i in range(adata.n_vars)]
    adata.var_names = [f"ENSG{i:011d}" for i in range(adata.n_vars)]
    ds = load(_write(adata, tmp_path))
    assert ds.gene_symbol_column == "feature_name"
    assert ds.gene_symbols[0] == "GENE0"
    assert ds.var_names[0] == "ENSG00000000000"


def test_display_embeddings_are_found(counts: ad.AnnData, tmp_path: Path) -> None:
    adata = counts.copy()
    adata.obsm["X_umap"] = np.random.default_rng(0).random((adata.n_obs, 2))
    adata.obsm["X_pca"] = np.random.default_rng(0).random((adata.n_obs, 50))
    adata.obsm["one_col"] = np.zeros((adata.n_obs, 1))
    ds = load(_write(adata, tmp_path))
    assert set(ds.embeddings) == {"umap", "pca"}
    assert ds.embeddings["pca"].shape == (adata.n_obs, 2)
    assert ds.embeddings["umap"].dtype == np.float32


def test_uns_root_is_read_but_validated(counts: ad.AnnData, tmp_path: Path) -> None:
    adata = counts.copy()
    adata.uns["iroot"] = 7
    assert load(_write(adata, tmp_path, "a.h5ad")).uns_root == 7
    adata.uns["iroot"] = 10_000
    assert load(_write(adata, tmp_path, "b.h5ad")).uns_root is None


def test_fingerprint_is_stable_and_content_sensitive(counts: ad.AnnData, tmp_path: Path) -> None:
    path = _write(counts, tmp_path)
    assert load(path).fingerprint == load(path).fingerprint
    other = counts.copy()
    other.obs_names = [f"x{i}" for i in range(other.n_obs)]
    assert load(_write(other, tmp_path, "o.h5ad")).fingerprint != load(path).fingerprint


def test_fingerprint_depends_on_how_the_matrix_is_read(counts: ad.AnnData, tmp_path: Path) -> None:
    # The same file read as counts and as log-normalised values is preprocessed
    # differently, so the two readings must never share a cache.
    cpm = counts.copy()
    sc.pp.normalize_total(cpm, target_sum=1e6)
    path = _write(cpm, tmp_path)
    as_counts, as_lognorm = load(path, "counts"), load(path, "lognorm")
    assert (as_counts.x_kind, as_lognorm.x_kind) == ("counts", "lognorm")
    assert as_counts.fingerprint != as_lognorm.fingerprint


def test_unknown_x_is_value_is_rejected(counts: ad.AnnData, tmp_path: Path) -> None:
    with pytest.raises(DatasetError, match="--x-is"):
        load(_write(counts, tmp_path), "logged")  # type: ignore[arg-type]


def test_loading_never_modifies_the_file(counts: ad.AnnData, tmp_path: Path) -> None:
    path = _write(counts, tmp_path)
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    load(path)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before


def test_missing_file(tmp_path: Path) -> None:
    with pytest.raises(DatasetError, match="does not exist"):
        load(tmp_path / "nope.h5ad")


def test_the_server_can_load_a_dataset_without_its_matrix(counts: ad.AnnData, tmp_path: Path) -> None:
    from lumencell.prep.pipeline import run_pipeline

    path = tmp_path / "d.h5ad"
    counts.write_h5ad(path)
    full, light = load(path), load(path, keep_matrix=False)
    assert light.matrix is None
    assert (light.n_obs, light.n_vars) == (full.n_obs, full.n_vars) == counts.shape
    np.testing.assert_array_equal(light.genes_detected, (np.asarray(counts.X.todense()) != 0).sum(axis=1))
    assert light.fingerprint == full.fingerprint, "the same cache serves both"
    with pytest.raises(DatasetError, match="without its expression matrix"):
        run_pipeline(light)
