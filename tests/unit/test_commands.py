"""The terminal commands: prepare, roots, run, export."""

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
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
        (["--donor-key", "donor"], "--donor-key needs --clone-key"),
        (["--single-donor"], "--single-donor needs --clone-key"),
        (["--strata", "donor"], "--strata needs --clone-key"),
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


@pytest.fixture
def stratified(clonal_data: Path) -> Path:
    adata = ad.read_h5ad(clonal_data)
    n = adata.n_obs
    adata.obs["s1"] = pd.Categorical(["a", "b"] * (n // 2))
    adata.obs["s2"] = pd.Categorical(["x", "y", "z"] * (n // 3))
    adata.obs["lineage_tag"] = adata.obs["clone_id"].copy()
    adata.write_h5ad(clonal_data)
    return clonal_data


def test_strata_are_split_stripped_stored_and_stay_visible(stratified: Path) -> None:
    from lumencell.api import models
    from lumencell.api.convert import to_api_spec
    from lumencell.api.session import Session
    from lumencell.engine.roots import parse_root
    from lumencell.jobs.manager import JobManager
    from lumencell.trust.clonal import trust_name

    result = runner.invoke(
        app, ["run", str(stratified), "--root", "cluster:s1=a", *CLONAL, "--strata", "s1, s2,"]
    )
    assert result.exit_code == 0, result.output
    rid = result.output.split("run ")[1].split(":")[0]
    settings = models.ClonalSettings(clone_key="clone_id", donor_key="donor", strata_keys=["s1", "s2"])
    store = RunStore(dataset_dir(load(stratified).fingerprint))
    stored = json.loads(store.trust_path(rid, trust_name(settings)).read_text())
    assert stored["settings"]["strata_keys"] == ["s1", "s2"]
    jobs = JobManager()
    session = Session(stratified, jobs)
    session.set_clonal(settings)
    created, _ = session.runs.create_run(
        models.TrajectoryRequest(method="paga_dpt", root=to_api_spec(parse_root("cluster:s1=a")))
    )
    jobs.shutdown()
    assert created.run_id == rid


def test_strata_naming_a_missing_column_is_refused(stratified: Path) -> None:
    result = runner.invoke(app, ["run", str(stratified), *ROOT, *CLONAL, "--strata", "s1,nope"])
    assert result.exit_code == 2 and "nope" in result.output


def test_a_neutrally_named_clone_or_donor_column_is_scoring_only(stratified: Path) -> None:
    tag = ["--clone-key", "lineage_tag", "--donor-key", "donor"]
    for root in ("cluster:lineage_tag=c0", "cluster:donor=d0"):
        result = runner.invoke(app, ["run", str(stratified), "--root", root, *tag])
        assert result.exit_code == 2, result.output
        assert "only to score" in result.output
    seen: list[set[str]] = []
    import lumencell.engine.registry as registry

    real = registry.run_method

    def spy(ctx, *a, **k):  # type: ignore[no-untyped-def]
        seen.append(set(ctx.obs.columns))
        return real(ctx, *a, **k)

    import pytest as _p

    mp = _p.MonkeyPatch()
    mp.setattr(registry, "run_method", spy)
    try:
        ok = runner.invoke(app, ["run", str(stratified), *ROOT, *tag, "--strata", "s1"])
    finally:
        mp.undo()
    assert ok.exit_code == 0, ok.output
    assert "lineage_tag" not in seen[0] and "donor" not in seen[0] and "s1" in seen[0]


def test_agreement_gets_the_prep_pcs_and_the_run_obs_keys(
    clonal_data: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import lumencell.trust.clonal as clonal

    calls: list[tuple[object, object]] = []
    real = clonal.agreement_for_run

    def spy(result, obs, settings, pca, neighbors_pcs, obs_keys):  # type: ignore[no-untyped-def]
        calls.append((neighbors_pcs, list(obs_keys)))
        return real(result, obs, settings, pca, neighbors_pcs, obs_keys)

    monkeypatch.setattr(clonal, "agreement_for_run", spy)
    result = runner.invoke(app, ["run", str(clonal_data), *ROOT, *CLONAL])
    assert result.exit_code == 0, result.output
    assert calls == [(PrepParams().neighbors_pcs, ["lumencell_leiden", "true_label"])]


def test_a_trust_write_failure_keeps_the_run_and_exits_one(
    clonal_data: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def boom(self, rid, name, payload):  # type: ignore[no-untyped-def]
        raise OSError("disk full")

    monkeypatch.setattr(RunStore, "write_trust", boom)
    result = runner.invoke(app, ["run", str(clonal_data), *ROOT, *CLONAL])
    assert result.exit_code == 1
    rid = result.output.split("run ")[1].split(":")[0]
    assert RunStore(dataset_dir(load(clonal_data).fingerprint)).exists(rid)
    assert f"run {rid} was saved" in result.output and "disk full" in result.output
    assert "clonal agreement: ok" in result.output


def test_a_scoring_failure_keeps_the_run_and_exits_one(
    clonal_data: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import lumencell.trust.clonal as clonal

    def boom(*a, **k):  # type: ignore[no-untyped-def]
        raise RuntimeError("bad scores")

    monkeypatch.setattr(clonal, "agreement_for_run", boom)
    result = runner.invoke(app, ["run", str(clonal_data), *ROOT, *CLONAL])
    assert result.exit_code == 1
    assert "was saved" in result.output and "bad scores" in result.output


def _agreement(**over):  # type: ignore[no-untyped-def]
    from lumencell.api.models import ClonalAgreement, ClonalScore, ClonalSettings, NullComparison

    def score(m: int, **kw):  # type: ignore[no-untyped-def]
        base = dict(
            status="ok", mode="own_cells", n_lineages=3, min_clone_size=m, n_cells=1234, n_clones=567,
            effective_clones=45.67, n_donors=3, n_strata=None, n_permutations=100, seed=1, purity=0.12345,
            mean_entropy=0.5, nulls=[NullComparison(null=n, mean_entropy_mean=0.1, mean_entropy_sd=0.1,
            eta=0.2, p=0.01) for n in ("N1", "N2")], eta=0.2, p=0.01234, eta_base=0.3, delta_eta=-0.1,
            p_delta=0.5, clone_weighted_purity=0.1, purity_without_largest=0.1,
        )  # fmt: skip
        return ClonalScore(**{**base, **kw})

    primary = score(2, **over.get("primary", {}))
    m3 = score(3, **{"eta": 0.4, "p": 0.2, **over.get("m3", {})})
    return ClonalAgreement(
        status=over.get("status", "ok"),
        settings=ClonalSettings(clone_key="c", single_donor=True),
        primary=primary,
        sensitivity_m3=m3,
        coupling=[],
        notes=["a note"],
    )


def test_clonal_lines_are_exact() -> None:
    from lumencell.cli import _clonal_lines

    assert _clonal_lines(_agreement(status="too_few_clones")) == [
        "clonal agreement: too_few_clones (own cells, m = 2)",
        "  cells 1,234, clones 567 (effective 45.7), donors 3",
        "  purity 0.123",
        "  eta against N1 0.200",
        "  headline p 0.0123 (the largest over N1, N2)",
        "  eta_base 0.300, delta -0.100, p for delta 0.5000",
        "  m = 3: eta 0.400, headline p 0.2000",
        "  note: a note",
    ]


def test_clonal_lines_print_every_missing_value_as_a_dash() -> None:
    from lumencell.cli import _clonal_lines

    none = dict(eta=None, p=None, eta_base=None, delta_eta=None, p_delta=None, effective_clones=None,
                purity=None, nulls=[])  # fmt: skip
    text = _clonal_lines(_agreement(primary=none, m3={"eta": None, "p": None}))
    assert text[1:7] == [
        "  cells 1,234, clones 567 (effective –), donors 3",
        "  purity –",
        "  eta against N1 –",
        "  headline p –",
        "  eta_base –, delta –, p for delta –",
        "  m = 3: eta –, headline p –",
    ]
    assert not any("nan" in t or "None" in t for t in text)
