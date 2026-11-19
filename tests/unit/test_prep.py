import hashlib
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pytest
import scanpy as sc
from scipy import sparse

from lumencell.bench.synth import generate
from lumencell.clonal.isolation import receptor_gene_mask
from lumencell.data.loader import Dataset, load
from lumencell.prep.cell_cycle_genes import load_gene_set
from lumencell.prep.pipeline import PrepResult, run_pipeline
from lumencell.prep.stages import STAGES, PrepParams, multiscale_space
from lumencell.trust.constants import CELL_CYCLE_SEED


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
    assert {r.name: r.state for r in res.stages}["umap"] == "skipped"
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


RECEPTOR_SYMBOLS = ["TRAV1-2", "TRAV12-1", "TRAJ33", "TRAC", "TRBV5-1", "TRBV20-1", "TRBJ2-7", "TRBC2"]


def _with_receptor_genes(tmp_path: Path) -> tuple[Dataset, np.ndarray]:
    """1,000 synthetic genes plus eight receptor genes that mark the true lineages."""
    adata = generate("bifurcation", "realistic", n_cells=1500, n_genes=1000, seed=0)
    rng = np.random.default_rng(0)
    lineage = adata.obs["true_lineage"].cat.codes.to_numpy()
    cols = []
    for i in range(len(RECEPTOR_SYMBOLS)):
        lam = np.where(lineage == i % len(adata.obs["true_lineage"].cat.categories), 30.0, 0.5)
        cols.append(rng.poisson(lam).astype(np.float32))
    x = sparse.hstack([adata.X, sparse.csr_matrix(np.column_stack(cols))], format="csr")
    adata = ad.AnnData(X=x, obs=adata.obs.copy())
    adata.var_names = [*(f"G{i:05d}" for i in range(1000)), *RECEPTOR_SYMBOLS]
    path = tmp_path / "tr.h5ad"
    adata.write_h5ad(path)
    ds = load(path)
    return ds, receptor_gene_mask(ds.gene_symbols)


@pytest.fixture(scope="module")
def receptor_ds(tmp_path_factory: pytest.TempPathFactory) -> tuple[Dataset, np.ndarray]:
    return _with_receptor_genes(tmp_path_factory.mktemp("tr"))


def test_receptor_genes_are_never_hvgs(receptor_ds: tuple[Dataset, np.ndarray]) -> None:
    ds, tr = receptor_ds
    assert tr.sum() == 8
    art = run_pipeline(ds, PrepParams(n_top_genes=200)).artifacts
    assert art["hvg"].sum() == 200
    assert art["hvg"][tr].sum() == 0
    info = art["hvg_info"]
    assert info["receptor_genes"] == 8 and info["excluded"] is True
    assert info["receptor_hvg_without_exclusion"] >= 1


def test_receptor_genes_can_be_kept_on_request(receptor_ds: tuple[Dataset, np.ndarray]) -> None:
    ds, tr = receptor_ds
    art = run_pipeline(ds, PrepParams(n_top_genes=200, exclude_receptor_genes=False)).artifacts
    info = art["hvg_info"]
    assert art["hvg"].sum() == 200
    assert art["hvg"][tr].sum() == info["receptor_hvg_without_exclusion"] >= 1
    assert info["excluded"] is False


def test_synthetic_hvgs_are_unchanged_by_the_exclusion(tmp_path: Path) -> None:
    ds = _dataset(tmp_path)
    on = run_pipeline(ds, PrepParams(n_top_genes=200)).artifacts
    off = run_pipeline(ds, PrepParams(n_top_genes=200, exclude_receptor_genes=False)).artifacts
    np.testing.assert_array_equal(on["hvg"], off["hvg"])
    assert on["hvg"].sum() == 200
    assert on["hvg_info"]["receptor_genes"] == 0


def test_the_pattern_covers_human_and_mouse_segments() -> None:
    hits = ["TRAV1-2", "TRAV14DV4", "Trbc1", "trgv9", "TRDC", "TRAC", "TRBD1", "TRBJ1-1"]
    misses = ["TRADD", "TRAT1", "TRAPPC1", "TRDN", "TRDMT1", "TRBC1-AS1"]
    assert receptor_gene_mask(np.array(hits)).all()
    assert not receptor_gene_mask(np.array(misses)).any()


def test_n_top_genes_above_the_gene_count_clamps_to_the_non_receptor_genes(
    receptor_ds: tuple[Dataset, np.ndarray], monkeypatch: pytest.MonkeyPatch
) -> None:
    ds, tr = receptor_ds
    asked: list[tuple[int, int]] = []
    real = sc.pp.highly_variable_genes

    def spy(adata: ad.AnnData, **kw: Any) -> Any:
        asked.append((adata.n_vars, kw["n_top_genes"]))
        return real(adata, **kw)

    monkeypatch.setattr(sc.pp, "highly_variable_genes", spy)
    art = run_pipeline(ds, PrepParams(n_top_genes=5000)).artifacts
    # scanpy clamps by itself, so check what the stage asked for, not only the result
    assert all(n_top <= n_genes for n_genes, n_top in asked), asked
    assert art["hvg"].sum() == ds.n_vars - tr.sum()
    assert art["hvg"][tr].sum() == 0


def test_receptor_genes_are_excluded_from_lognorm_data_too(tmp_path: Path) -> None:
    ds, tr = _with_receptor_genes(tmp_path)
    adata = ad.read_h5ad(ds.path)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    path = tmp_path / "ln.h5ad"
    adata.write_h5ad(path)
    ln = load(path)
    assert ln.x_kind == "lognorm"
    art = run_pipeline(ln, PrepParams(n_top_genes=200)).artifacts
    assert art["hvg"].sum() == 200
    assert art["hvg"][tr].sum() == 0
    assert art["hvg_info"]["receptor_genes"] == 8


def _ensembl_dataset(tmp_path: Path, symbol_column: str | None) -> tuple[Dataset, np.ndarray]:
    ds, tr = _with_receptor_genes(tmp_path)
    adata = ad.read_h5ad(ds.path)
    symbols = np.array(adata.var_names)
    adata.var_names = [f"ENSG{i:011d}" for i in range(adata.n_vars)]
    if symbol_column:
        adata.var[symbol_column] = symbols
    path = tmp_path / "ens.h5ad"
    adata.write_h5ad(path)
    return load(path), tr


def test_receptor_genes_are_matched_on_the_symbol_column(tmp_path: Path) -> None:
    ds, tr = _ensembl_dataset(tmp_path, "feature_name")
    assert receptor_gene_mask(ds.var_names.to_numpy()).sum() == 0
    art = run_pipeline(ds, PrepParams(n_top_genes=200)).artifacts
    assert art["hvg"][tr].sum() == 0
    assert art["hvg_info"]["receptor_genes"] == 8
    assert art["hvg_info"]["symbol_source"] == "feature_name"
    assert art["hvg_info"]["symbols_look_like_ensembl"] is False


def test_ensembl_ids_without_symbols_warn_that_exclusion_cannot_match(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    ds, _tr = _ensembl_dataset(tmp_path, None)
    with caplog.at_level("WARNING", logger="lumencell.prep.stages"):
        art = run_pipeline(ds, PrepParams(n_top_genes=200)).artifacts
    info = art["hvg_info"]
    assert info["receptor_genes"] == 0 and info["symbol_source"] == "var_names"
    assert info["symbols_look_like_ensembl"] is True
    assert "clone identity may reach the PCA" in caplog.text


def test_plain_symbols_do_not_warn(
    receptor_ds: tuple[Dataset, np.ndarray], caplog: pytest.LogCaptureFixture
) -> None:
    ds, _ = receptor_ds
    with caplog.at_level("WARNING", logger="lumencell.prep.stages"):
        art = run_pipeline(ds, PrepParams(n_top_genes=200)).artifacts
    assert art["hvg_info"]["symbol_source"] == "var_names"
    assert not caplog.records


def test_a_dataset_made_only_of_receptor_genes_is_refused(tmp_path: Path) -> None:
    from lumencell.data.loader import DatasetError

    rng = np.random.default_rng(0)
    x = sparse.csr_matrix(rng.poisson(3.0, (200, 8)).astype(np.float32))
    adata = ad.AnnData(X=x)
    adata.var_names = RECEPTOR_SYMBOLS
    path = tmp_path / "all.h5ad"
    adata.write_h5ad(path)
    with pytest.raises(DatasetError, match="all 8 genes"):
        run_pipeline(load(path), PrepParams(n_top_genes=200))


def test_paga_with_no_edge_between_clusters_is_an_all_zero_graph(tmp_path: Path) -> None:
    from lumencell.bench.synth_controls import generate_control

    path = tmp_path / "two.h5ad"
    generate_control("two_groups", n_cells=1000, n_genes=500, seed=0).write_h5ad(path)
    art = run_pipeline(load(path)).artifacts
    k = int(np.asarray(art["leiden"]).max()) + 1
    assert k >= 2
    assert art["paga_connectivities"].shape == (k, k)
    assert art["paga_connectivities"].nnz == 0
    assert art["paga_tree"].shape == (k, k)
    assert art["paga_tree"].nnz == 0
    assert art["paga_pos"].shape == (k, 2)
    assert np.isfinite(art["paga_pos"]).all()


# ---- the report-only cell_cycle stage -------------------------------------------------


def test_no_inference_stage_reads_a_report_only_output() -> None:
    report = {o for s in STAGES if s.report_only for o in s.outputs}
    assert report
    for stage in STAGES:
        if not stage.report_only:
            assert not (set(stage.inputs) & report), stage.name


def test_cell_cycle_is_not_scored_on_synthetic_genes(result: tuple[Dataset, PrepResult]) -> None:
    info = result[1].artifacts["cc_info"]
    assert info["status"] == "not_scored"
    assert info["species"] is None
    assert "0 of 97" in info["reason"]
    assert "species" not in info["reason"]


def _cc_dataset(
    tmp_path: Path, symbols: list[str], boost: bool = False, n_cells: int = 600, name: str = "cc.h5ad"
) -> Dataset:
    adata = generate("bifurcation", "realistic", n_cells=n_cells, n_genes=400, seed=0)
    var = list(adata.var_names)
    var[: len(symbols)] = symbols
    adata.var_names = var
    if boost:
        x = adata.X.toarray() if sparse.issparse(adata.X) else np.array(adata.X)
        x[0:150, 0:43] *= 8
        x[150:300, 43:97] *= 8
        adata.X = sparse.csr_matrix(x)
    path = tmp_path / name
    adata.write_h5ad(path)
    return load(path)


def _human() -> list[str]:
    genes = load_gene_set("tirosh2016", "human")
    return [g[0].upper() for g in genes.s + genes.g2m]


def _mouse() -> list[str]:
    genes = load_gene_set("tirosh2016", "mouse")
    return [g[0] for g in genes.s + genes.g2m]


def test_an_unresolved_species_is_reported(tmp_path: Path) -> None:
    human = _human()
    symbols = [h if i % 2 == 0 else h.capitalize() for i, h in enumerate(human[:30])]
    ds = _cc_dataset(tmp_path, symbols + [f"X{i}" for i in range(67)])
    art = run_pipeline(ds).artifacts
    assert art["cc_info"]["status"] == "not_scored"
    assert "species could not be told" in art["cc_info"]["reason"]
    assert art["cc_s"].dtype == np.float32 and art["cc_g2m"].dtype == np.float32
    assert art["cc_s"].shape == (ds.n_obs,) and np.isnan(art["cc_s"]).all() and np.isnan(art["cc_g2m"]).all()
    assert (art["cc_phase"] == "unscored").all()


def test_cell_cycle_scores_human_symbols(tmp_path: Path) -> None:
    art = run_pipeline(_cc_dataset(tmp_path, _published_human(), boost=True)).artifacts
    info = art["cc_info"]
    assert info["status"] == "scored", info
    assert info["species"] == "human"
    assert info["genes_found"] == {"s": 43, "g2m": 54}
    assert art["cc_s"].dtype == np.float32 and art["cc_g2m"].dtype == np.float32
    phase = art["cc_phase"]
    assert (phase[:150] == "S").mean() >= 0.8
    assert (phase[150:300] == "G2M").mean() >= 0.8
    assert (phase[300:] == "G1").mean() > (phase[:300] == "G1").mean()


def _published_human() -> list[str]:
    # The symbols as the paper printed them: for a renamed gene, the 2016 one.
    genes = load_gene_set("tirosh2016", "human")
    return [g[0] for g in genes.s + genes.g2m]


def test_cell_cycle_resolves_mouse_from_symbol_case(tmp_path: Path) -> None:
    art = run_pipeline(_cc_dataset(tmp_path, _mouse(), boost=True)).artifacts
    info = art["cc_info"]
    assert info["status"] == "scored", info
    assert info["species"] == "mouse"
    assert info["genes_found"] == {"s": 43, "g2m": 54}


def test_duplicate_symbols_are_scored_once(tmp_path: Path) -> None:
    symbols = _published_human()
    symbols.append(symbols[0])  # the first S gene appears twice
    art = run_pipeline(_cc_dataset(tmp_path, symbols, boost=True)).artifacts
    info = art["cc_info"]
    assert info["status"] == "scored", info
    assert info["genes_found"]["s"] == 43


@pytest.mark.parametrize("error", [RuntimeError, ValueError, KeyError, ZeroDivisionError])
def test_a_scoring_error_is_reported_and_prep_still_completes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    error: type[Exception],
) -> None:
    def boom(*_a: Any, **_k: Any) -> None:
        raise error("boom")

    monkeypatch.setattr(sc.tl, "score_genes_cell_cycle", boom)
    res = run_pipeline(_cc_dataset(tmp_path, _published_human()))
    info = res.artifacts["cc_info"]
    assert info["status"] == "not_scored"
    assert f"{error.__name__}:" in info["reason"] and "boom" in info["reason"]
    assert all(r.state in ("done", "skipped") for r in res.stages)
    assert "boom" in capsys.readouterr().err


def test_a_broken_gene_list_does_not_fail_prep(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def bad(*_a: Any, **_k: Any) -> None:
        raise ValueError("bad tsv")

    monkeypatch.setattr("lumencell.prep.stages.load_gene_set", bad)
    ds = _cc_dataset(tmp_path, _published_human())
    art = run_pipeline(ds).artifacts
    info = art["cc_info"]
    assert info["status"] == "not_scored"
    assert "ValueError: bad tsv" in info["reason"]
    assert art["cc_s"].dtype == np.float32 and art["cc_g2m"].dtype == np.float32
    assert art["cc_s"].shape == (ds.n_obs,) and np.isnan(art["cc_s"]).all() and np.isnan(art["cc_g2m"]).all()
    assert (art["cc_phase"] == "unscored").all()
    assert len(capsys.readouterr().err.strip().splitlines()) == 1


def _first_genes(n_s: int, n_g2m: int) -> list[str]:
    genes = load_gene_set("tirosh2016", "human")
    return [g[0] for g in genes.s[:n_s] + genes.g2m[:n_g2m]]


def test_ten_s_and_ten_g2m_genes_are_enough_and_nineteen_are_not(tmp_path: Path) -> None:
    on = run_pipeline(_cc_dataset(tmp_path, _first_genes(10, 10), name="a.h5ad")).artifacts["cc_info"]
    assert on["status"] == "scored", on
    assert on["genes_found"] == {"s": 10, "g2m": 10}
    off = run_pipeline(_cc_dataset(tmp_path, _first_genes(9, 10), name="b.h5ad")).artifacts["cc_info"]
    assert off["status"] == "not_scored"
    assert "only 19 of 97" in off["reason"]


def test_an_empty_phase_list_is_a_gene_count_reason_not_a_failure(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    info = run_pipeline(_cc_dataset(tmp_path, _first_genes(25, 0))).artifacts["cc_info"]
    assert info["status"] == "not_scored"
    assert "only 25 of 97" in info["reason"]
    assert "scoring failed" not in info["reason"]
    assert capsys.readouterr().err == ""


def test_scoring_is_seeded_with_the_pinned_seed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    seen: dict[str, Any] = {}
    real = sc.tl.score_genes_cell_cycle

    def spy(*a: Any, **k: Any) -> Any:
        seen.update(k)
        return real(*a, **k)

    monkeypatch.setattr(sc.tl, "score_genes_cell_cycle", spy)
    run_pipeline(_cc_dataset(tmp_path, _published_human()))
    assert seen["random_state"] == CELL_CYCLE_SEED


def test_unknown_cell_cycle_params_are_refused() -> None:
    with pytest.raises(ValueError):
        PrepParams(cc_species="dog")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        PrepParams(cc_gene_set="other")  # type: ignore[arg-type]
