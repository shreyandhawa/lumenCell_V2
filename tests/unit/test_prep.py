import hashlib
from pathlib import Path

import numpy as np
import pytest
import scanpy as sc
from scipy import sparse

from lumencell.bench.synth import generate
from lumencell.data.loader import Dataset, load
from lumencell.prep.pipeline import PrepResult, run_pipeline
from lumencell.prep.stages import STAGES, PrepParams, multiscale_space


def _dataset(tmp_path: Path, n_cells: int = 1500, umap: bool = False, lognorm: bool = False) -> Dataset:
    adata = generate("bifurcation", "realistic", n_cells=n_cells, n_genes=1000, seed=0)
    if lognorm:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    if umap:
        adata.obsm["X_umap"] = np.random.default_rng(0).random((adata.n_obs, 2)).astype(np.float32)
    path = tmp_path / "d.h5ad"
    adata.write_h5ad(path)
    return load(path)


@pytest.fixture(scope="module")
def result(tmp_path_factory: pytest.TempPathFactory) -> tuple[Dataset, PrepResult]:
    ds = _dataset(tmp_path_factory.mktemp("prep"))
    return ds, run_pipeline(ds)


def test_stage_graph_is_consistent() -> None:
    produced: set[str] = set()
    names = [s.name for s in STAGES]
    assert len(set(names)) == len(names)
    for stage in STAGES:
        assert set(stage.inputs) <= produced, f"{stage.name} reads something produced later"
        assert not (set(stage.outputs) & produced), f"{stage.name} overwrites an artifact"
        assert set(stage.params) <= set(PrepParams().to_dict())
        produced |= set(stage.outputs)


def test_outputs_and_shapes(result: tuple[Dataset, PrepResult]) -> None:
    ds, res = result
    art = res.artifacts
    n = ds.n_obs
    assert sparse.issparse(art["lognorm"]) and art["lognorm"].shape == ds.matrix.shape
    assert art["hvg"].dtype == bool and art["hvg"].sum() == 1000  # n_top_genes clamps to n_vars
    assert art["pca"].shape == (n, 50) and art["pca"].dtype == np.float32
    assert art["connectivities"].shape == (n, n)
    assert art["diffmap"].shape == (n, 15)
    assert art["diffmap_evals"][0] == pytest.approx(1.0)
    assert np.all(np.diff(art["diffmap_evals"]) <= 1e-9), "eigenvalues are sorted"
    assert art["diffmap_ms"].shape[1] >= 2
    k = int(art["leiden"].max()) + 1
    assert k >= 3, "a bifurcation should split into several clusters"
    assert art["paga_connectivities"].shape == (k, k)
    assert art["paga_pos"].shape == (k, 2)
    assert art["umap"].shape == (n, 2)


def test_leiden_labels_are_ordered_by_size(result: tuple[Dataset, PrepResult]) -> None:
    sizes = np.bincount(result[1].artifacts["leiden"])
    assert np.all(np.diff(sizes) <= 0)


def test_every_stage_is_recorded(result: tuple[Dataset, PrepResult]) -> None:
    records = result[1].stages
    assert [r.name for r in records] == [s.name for s in STAGES]
    assert all(r.state == "done" for r in records)


def test_dataset_is_not_modified(tmp_path: Path) -> None:
    ds = _dataset(tmp_path, n_cells=300)
    before_matrix = ds.matrix.copy()
    before_file = hashlib.sha256(ds.path.read_bytes()).hexdigest()
    run_pipeline(ds)
    assert (ds.matrix != before_matrix).nnz == 0
    assert hashlib.sha256(ds.path.read_bytes()).hexdigest() == before_file


def test_deterministic(tmp_path: Path) -> None:
    ds = _dataset(tmp_path, n_cells=600)
    a, b = run_pipeline(ds).artifacts, run_pipeline(ds).artifacts
    np.testing.assert_array_equal(a["leiden"], b["leiden"])
    np.testing.assert_allclose(a["pca"], b["pca"], rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(a["diffmap"], b["diffmap"], rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(a["umap"], b["umap"], rtol=1e-4, atol=1e-4)


def test_existing_embedding_skips_umap(tmp_path: Path) -> None:
    ds = _dataset(tmp_path, n_cells=300, umap=True)
    res = run_pipeline(ds)
    assert "umap" not in res.artifacts
    assert res.stages[-1].state == "skipped"
    assert "umap" in run_pipeline(ds, PrepParams(recompute_umap=True)).artifacts


@pytest.mark.filterwarnings("ignore:`n_top_genes`:UserWarning")
def test_lognorm_input(tmp_path: Path) -> None:
    ds = _dataset(tmp_path, n_cells=300, lognorm=True)
    assert ds.x_kind == "lognorm"
    res = run_pipeline(ds)
    assert (res.artifacts["lognorm"] != ds.matrix).nnz == 0, "lognorm input is used as is"
    assert res.artifacts["hvg"].any()


def test_progress_reports_each_stage_in_order(tmp_path: Path) -> None:
    ds = _dataset(tmp_path, n_cells=300)
    seen: list[str] = []
    run_pipeline(ds, progress=lambda stage, _i, _n: seen.append(stage))
    assert seen == [s.name for s in STAGES] + ["done"]


def test_multiscale_space_matches_palantirs_rule() -> None:
    evals = np.array([1.0, 0.99, 0.98, 0.5, 0.4])
    evecs = np.ones((4, 5))
    ms = multiscale_space(evecs, evals)
    # Largest gap is after index 2 (0.98 -> 0.5), so n_eigs = 3 and components 1..2
    # survive, each scaled by lambda / (1 - lambda).
    np.testing.assert_allclose(ms[0], [0.99 / 0.01, 0.98 / 0.02], rtol=1e-4)


def test_multiscale_space_keeps_at_least_two_components() -> None:
    # The largest gap is right after the trivial eigenvalue: fall back to the second largest.
    evals = np.array([1.0, 0.5, 0.45, 0.2, 0.19])
    assert multiscale_space(np.ones((3, 5)), evals).shape[1] >= 2
