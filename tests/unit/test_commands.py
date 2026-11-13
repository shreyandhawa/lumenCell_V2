"""The terminal commands: prepare, roots, run, export."""

from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy import sparse
from typer.testing import CliRunner

from lumencell.bench.synth import generate
from lumencell.cli import app
from lumencell.data.loader import load
from lumencell.engine.store import RunStore
from lumencell.prep.cache import dataset_dir, prep_key
from lumencell.prep.stages import PrepParams

runner = CliRunner()


@pytest.fixture
def data(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.setenv("LUMENCELL_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("LUMENCELL_DATA_DIR", str(tmp_path / "state"))
    path = tmp_path / "d.h5ad"
    adata = generate("bifurcation", "realistic", n_cells=600, n_genes=300, seed=0, store_root=True)
    adata.obs["clone_id"] = adata.obs["true_lineage"].copy()  # scoring-only: no run may use it
    adata.write_h5ad(path)
    return path


def test_prepare_then_cached(data: Path) -> None:
    first = runner.invoke(app, ["prepare", str(data)])
    assert first.exit_code == 0, first.output
    assert "pca" in first.output and "done" in first.output
    second = runner.invoke(app, ["prepare", str(data)])
    assert "cached" in second.output


def test_roots_prints_suggestions(data: Path) -> None:
    result = runner.invoke(app, ["roots", str(data)])
    assert result.exit_code == 0, result.output
    assert "file" in result.output and "[gene_count]" in result.output


def test_run_and_export(data: Path, tmp_path: Path) -> None:
    result = runner.invoke(app, ["run", str(data), "--root", "cluster:true_label=trunk_early"])
    assert result.exit_code == 0, result.output
    run_id = result.output.split("run ")[1].split(":")[0]
    ds = load(data)
    graph = RunStore(dataset_dir(ds.fingerprint)).graph(run_id)
    assert graph["prep_key"] == prep_key(ds, PrepParams())
    assert graph["obs_keys"] == ["lumencell_leiden", "true_label"]
    out = tmp_path / "exported.h5ad"
    exported = runner.invoke(app, ["export", str(data), run_id, str(out)])
    assert exported.exit_code == 0, exported.output
    adata = ad.read_h5ad(out)
    pt = adata.obs["lumencell_pseudotime"].to_numpy()
    assert np.isfinite(pt).mean() > 0.9
    assert "shared" in adata.obs["lumencell_lineage"].cat.categories
    assert adata.uns["lumencell_run"]["root"] == "cluster:true_label=trunk_early"
    assert "X_umap" in adata.obsm
    assert "lumencell_pseudotime" not in ad.read_h5ad(data).obs, "the input is untouched"


def test_export_refuses_to_overwrite_the_input(data: Path) -> None:
    result = runner.invoke(app, ["run", str(data), "--root", "file"])
    run_id = result.output.split("run ")[1].split(":")[0]
    refused = runner.invoke(app, ["export", str(data), run_id, str(data)])
    assert refused.exit_code != 0
    assert "refusing to overwrite" in refused.output


def test_run_refuses_a_clone_column_as_the_root(data: Path) -> None:
    result = runner.invoke(app, ["run", str(data), "--root", "cluster:clone_id=A"])
    assert result.exit_code != 0
    message = " ".join(result.output.replace("│", " ").split())
    assert "'clone_id' holds clone or receptor data, which lumenCell uses only to score runs" in message


def test_run_rejects_a_bad_root(data: Path) -> None:
    result = runner.invoke(app, ["run", str(data), "--root", "auto"])
    assert result.exit_code != 0
    assert "cannot parse root" in result.output


def test_unusable_file_is_reported_clearly(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_CACHE_DIR", str(tmp_path / "cache"))
    import scanpy as sc

    adata = generate("linear", "easy", n_cells=100, n_genes=50, seed=0)
    sc.pp.normalize_total(adata, target_sum=1e6)
    path = tmp_path / "cpm.h5ad"
    adata.write_h5ad(path)
    result = runner.invoke(app, ["prepare", str(path)])
    assert result.exit_code != 0
    assert "--x-is" in result.output


def _receptor_file(tmp_path: Path, with_receptors: bool = True) -> Path:
    adata = generate("bifurcation", "realistic", n_cells=600, n_genes=300, seed=0)
    symbols = ["TRAV1-2", "TRAV12-1", "TRAJ33", "TRAC", "TRBV5-1", "TRBV20-1", "TRBJ2-7", "TRBC2"]
    if with_receptors:
        extra = np.random.default_rng(0).poisson(2.0, (adata.n_obs, len(symbols))).astype(np.float32)
        x = sparse.hstack([adata.X, sparse.csr_matrix(extra)], format="csr")
        adata = ad.AnnData(X=x, obs=adata.obs.copy())
        adata.var_names = [*(f"G{i:05d}" for i in range(300)), *symbols]
    path = tmp_path / "tr.h5ad"
    adata.write_h5ad(path)
    return path


@pytest.fixture
def cache_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_CACHE_DIR", str(tmp_path / "cache"))


def test_prepare_reports_receptor_genes(tmp_path: Path, cache_env: None) -> None:
    result = runner.invoke(app, ["prepare", str(_receptor_file(tmp_path))])
    assert result.exit_code == 0, result.output
    assert "receptor genes: 8 matched" in result.output
    assert "(symbols from var_names); excluded from HVG selection (" in result.output
    assert "would be HVGs without the exclusion)" in result.output


def test_prepare_reports_when_receptor_genes_are_kept(
    tmp_path: Path, cache_env: None, monkeypatch: pytest.MonkeyPatch
) -> None:
    from lumencell.prep.stages import PrepParams

    monkeypatch.setattr("lumencell.cli._prep_params", lambda r: PrepParams(exclude_receptor_genes=False))
    result = runner.invoke(app, ["prepare", str(_receptor_file(tmp_path))])
    assert result.exit_code == 0, result.output
    assert "kept in HVG selection (" in result.output
    assert "excluded from" not in result.output


def test_prepare_says_so_when_no_receptor_genes_matched(tmp_path: Path, cache_env: None) -> None:
    result = runner.invoke(app, ["prepare", str(_receptor_file(tmp_path, with_receptors=False))])
    assert result.exit_code == 0, result.output
    assert "receptor genes: 0 matched (symbols from var_names)" in result.output
    assert "warning" not in result.output


def test_prepare_warns_when_ids_cannot_be_matched(tmp_path: Path, cache_env: None) -> None:
    adata = ad.read_h5ad(_receptor_file(tmp_path))
    adata.var_names = [f"ENSG{i:011d}" for i in range(adata.n_vars)]
    path = tmp_path / "ens.h5ad"
    adata.write_h5ad(path)
    result = runner.invoke(app, ["prepare", str(path)])
    assert result.exit_code == 0, result.output
    assert "receptor genes: 0 matched (symbols from var_names)" in result.output
    assert "receptor genes could not be matched" in result.output
    assert "clone identity may reach the PCA" in result.output


# ---- clonal agreement ----------------------------------------------------------------------


@pytest.fixture
def clonal_data(data: Path) -> Path:
    from lumencell.bench.synth_clones import ClonalLayer, add_clones

    base = ad.read_h5ad(data)
    del base.obs["clone_id"]
    adata = add_clones(base, ClonalLayer(fate_bias=0.8), seed=0)
    adata.write_h5ad(data)
    return data


CLONAL = ["--clone-key", "clone_id", "--donor-key", "donor"]
ROOT = ["--root", "cluster:true_label=trunk_early"]


def test_run_prints_clonal_agreement(clonal_data: Path) -> None:
    result = runner.invoke(app, ["run", str(clonal_data), *ROOT, *CLONAL])
    assert result.exit_code == 0, result.output
    assert "clonal agreement: ok (own cells, m = 2)" in result.output
    assert "m = 3:" in result.output and "delta" in result.output
    assert "headline p" in result.output and "the largest over N1" in result.output
    assert "nan" not in result.output and "None" not in result.output
    rid = result.output.split("run ")[1].split(":")[0]
    from lumencell.api.models import ClonalSettings
    from lumencell.trust.clonal import trust_name

    store = RunStore(dataset_dir(load(clonal_data).fingerprint))
    name = trust_name(ClonalSettings(clone_key="clone_id", donor_key="donor"))
    assert store.trust_path(rid, name).exists()


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        (["--clone-key", "clone_id", "--donor-key", "donor", "--single-donor"], "not both"),
        (["--clone-key", "nope", "--donor-key", "donor"], "nope"),
        (["--donor-key", "donor"], "--clone-key"),
        (["--single-donor"], "--clone-key"),
        (["--clone-key", "clone_id"], "single_donor"),
    ],
)
def test_run_refuses_invalid_clonal_options(clonal_data: Path, extra: list[str], message: str) -> None:
    result = runner.invoke(app, ["run", str(clonal_data), *ROOT, *extra])
    assert result.exit_code == 2, result.output
    assert message in result.output


def test_run_refuses_a_clone_column_as_root(clonal_data: Path) -> None:
    result = runner.invoke(app, ["run", str(clonal_data), "--root", "cluster:clone_id=c0", *CLONAL])
    assert result.exit_code == 2, result.output
    assert "only to score" in result.output


def test_run_id_with_settings_matches_the_api_and_differs_without(clonal_data: Path) -> None:
    from lumencell.api import models
    from lumencell.api.convert import to_api_spec
    from lumencell.api.session import Session
    from lumencell.engine.roots import parse_root
    from lumencell.jobs.manager import JobManager

    plain = runner.invoke(app, ["run", str(clonal_data), *ROOT])
    scored = runner.invoke(app, ["run", str(clonal_data), *ROOT, *CLONAL])
    rid_plain = plain.output.split("run ")[1].split(":")[0]
    rid = scored.output.split("run ")[1].split(":")[0]
    assert rid != rid_plain

    jobs = JobManager()  # not started: the job is only queued
    session = Session(clonal_data, jobs)
    session.set_clonal(models.ClonalSettings(clone_key="clone_id", donor_key="donor"))
    created, _ = session.runs.create_run(
        models.TrajectoryRequest(method="paga_dpt", root=to_api_spec(parse_root(ROOT[1])))
    )
    assert created.run_id == rid
    jobs.shutdown()


def test_clonal_lines_print_missing_statistics_as_a_dash_and_a_refusal_with_its_note() -> None:
    from lumencell.api.models import ClonalAgreement, ClonalSettings
    from lumencell.cli import _clonal_lines

    refused = ClonalAgreement(
        status="refused",
        settings=ClonalSettings(clone_key="c", single_donor=True),
        primary=None,
        sensitivity_m3=None,
        coupling=[],
        notes=["The run used 'x', which lumenCell keeps only for scoring; run it again without it."],
    )
    text = "\n".join(_clonal_lines(refused))
    assert text.startswith("clonal agreement: refused") and "only for scoring" in text
