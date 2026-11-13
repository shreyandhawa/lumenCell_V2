"""Clone labels score runs and never steer them: through obs or through expression."""

import dataclasses
import math
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from lumencell.bench.synth import generate
from lumencell.clonal.isolation import (
    ISOLATION_VERSION,
    RECEPTOR_COLUMN,
    scoring_only_columns,
    scoring_only_message,
)
from lumencell.data.loader import Dataset, load
from lumencell.engine import registry
from lumencell.engine.context import (
    InferenceContext,
    ScoringOnlyColumnError,
    check_inference_column,
    inference_obs,
)
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from lumencell.prep.pipeline import PrepResult, run_pipeline
from lumencell.prep.stages import PrepParams

# The eight receptor genes of Task 4's tests (tests/unit/test_prep.py).
RECEPTOR_SYMBOLS = ["TRAV1-2", "TRAV12-1", "TRAJ33", "TRAC", "TRBV5-1", "TRBV20-1", "TRBJ2-7", "TRBC2"]
ROOT = ClusterRoot("true_label", "trunk_early")


@pytest.fixture(scope="module")
def plain(tmp_path_factory: pytest.TempPathFactory) -> Dataset:
    """A synthetic bifurcation plus eight receptor genes whose counts know nothing of lineage."""
    adata = generate("bifurcation", "realistic", n_cells=1500, n_genes=1000, seed=0)
    rng = np.random.default_rng(0)
    receptor = rng.poisson(4, size=(adata.n_obs, len(RECEPTOR_SYMBOLS))).astype(np.float32)
    x = sparse.hstack([adata.X, sparse.csr_matrix(receptor)], format="csr")
    out = ad.AnnData(X=x, obs=adata.obs.copy())
    out.var_names = [*adata.var_names, *RECEPTOR_SYMBOLS]
    path = tmp_path_factory.mktemp("isolation") / "plain.h5ad"
    out.write_h5ad(path)
    return load(path)


@pytest.fixture(scope="module")
def clonal(plain: Dataset) -> Dataset:
    """The same cells, now with a clone column and receptor genes that encode the clone.

    Each cell's eight receptor counts are the same multiset reordered so that its largest
    count sits at index (lineage code * 3) % 8, so library sizes do not change.
    """
    matrix = plain.require_matrix()
    n_receptor = len(RECEPTOR_SYMBOLS)
    n_other = plain.n_vars - n_receptor
    assert [str(s) for s in plain.gene_symbols[n_other:]] == RECEPTOR_SYMBOLS
    codes = plain.obs["true_lineage"].cat.codes.to_numpy()
    block = matrix[:, n_other:].toarray()
    for i, row in enumerate(block):
        j, k = int(np.argmax(row)), int(codes[i] * 3) % n_receptor
        row[j], row[k] = row[k], row[j]
    reordered = sparse.hstack([matrix[:, :n_other], sparse.csr_matrix(block)], format="csr", dtype=np.float32)
    obs = plain.obs.copy()
    obs["clone_id"] = obs["true_lineage"].copy()
    return dataclasses.replace(plain, obs=obs, matrix=reordered)


@pytest.fixture(scope="module")
def preps(plain: Dataset, clonal: Dataset) -> dict[str, PrepResult]:
    return {"plain": run_pipeline(plain, PrepParams()), "clonal": run_pipeline(clonal, PrepParams())}


def _ctx(ds: Dataset, prep: PrepResult, clonal_columns: tuple[str, ...] = ()) -> InferenceContext:
    return InferenceContext.from_prep(ds, prep, PrepParams(), clonal_columns=clonal_columns)


def _run(ds: Dataset, prep: PrepResult, method: str, clonal_columns: tuple[str, ...]) -> TrajectoryResult:
    return run_method(_ctx(ds, prep, clonal_columns), method, ROOT)


def _same(a: Any, b: Any, where: str) -> None:
    """Equal value for value; NaN equals NaN; arrays keep their dtype."""
    assert type(a) is type(b), f"{where}: {type(a).__name__} vs {type(b).__name__}"
    if isinstance(a, np.ndarray):
        assert a.dtype == b.dtype and a.shape == b.shape, where
        np.testing.assert_array_equal(a, b, err_msg=where)
    elif sparse.issparse(a):
        assert a.shape == b.shape, where
        np.testing.assert_array_equal(a.toarray(), b.toarray(), err_msg=where)
    elif dataclasses.is_dataclass(a) and not isinstance(a, type):
        for f in dataclasses.fields(a):
            _same(getattr(a, f.name), getattr(b, f.name), f"{where}.{f.name}")
    elif isinstance(a, dict):
        assert a.keys() == b.keys(), where
        for k in a:
            _same(a[k], b[k], f"{where}[{k!r}]")
    elif isinstance(a, list | tuple):
        assert len(a) == len(b), where
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _same(x, y, f"{where}[{i}]")
    elif isinstance(a, float) and math.isnan(a):
        assert math.isnan(b), where
    else:
        assert a == b, f"{where}: {a!r} != {b!r}"


def _assert_identical(a: TrajectoryResult, b: TrajectoryResult) -> None:
    for f in dataclasses.fields(TrajectoryResult):
        if f.name != "runtime_seconds":
            _same(getattr(a, f.name), getattr(b, f.name), f.name)


# ---- the result does not move -------------------------------------------------------------


@pytest.mark.parametrize(
    "method",
    [
        "paga_dpt",
        pytest.param("palantir", marks=pytest.mark.slow),
        pytest.param("scfates_ppt", marks=pytest.mark.slow),
    ],
)
def test_clone_columns_and_receptor_genes_change_no_result(
    method: str, plain: Dataset, clonal: Dataset, preps: dict[str, PrepResult]
) -> None:
    # The plain file has no clone_id, and naming a column the file lacks is refused.
    a = _run(plain, preps["plain"], method, ())
    b = _run(clonal, preps["clonal"], method, ("clone_id",))
    _assert_identical(a, b)


def test_the_test_has_teeth(plain: Dataset, clonal: Dataset) -> None:
    off = PrepParams(exclude_receptor_genes=False)
    assert not np.array_equal(
        run_pipeline(plain, off).artifacts["pca"], run_pipeline(clonal, off).artifacts["pca"]
    )


def test_the_clone_column_really_is_there(clonal: Dataset, preps: dict[str, PrepResult]) -> None:
    # Otherwise the identical results above would prove nothing about obs.
    assert "clone_id" in clonal.obs.columns
    assert "clone_id" not in _ctx(clonal, preps["clonal"]).obs.columns


# ---- what methods see ---------------------------------------------------------------------


@pytest.mark.parametrize(
    "column",
    [
        "clone_id",
        "clone_id_size",
        "clonotype",
        "raw_clonotype_id",
        "IR_VJ_1_junction_aa",
        "airr:clone_id",
        "receptor_type",
        "chain_pairing",
        "cdr3_nt",
    ],
)
def test_receptor_columns_never_reach_methods(
    column: str, plain: Dataset, preps: dict[str, PrepResult]
) -> None:
    obs = plain.obs.copy()
    labels = obs["true_label"].astype(str).to_numpy()
    obs[column] = labels
    obs["cell_type"] = labels
    obs["donor_age"] = 40
    obs["tissue"] = "blood"
    ctx = _ctx(dataclasses.replace(plain, obs=obs), preps["plain"])
    assert column not in ctx.obs.columns
    assert column in ctx.scoring_only
    assert {"cell_type", "donor_age", "tissue"} <= set(ctx.obs.columns)
    assert not ({"cell_type", "donor_age", "tissue"} & ctx.scoring_only)


def test_named_columns_are_scoring_only(plain: Dataset, preps: dict[str, PrepResult]) -> None:
    obs = plain.obs.copy()
    obs["donor"] = np.where(np.arange(plain.n_obs) % 2 == 0, "d1", "d2")
    ctx = _ctx(dataclasses.replace(plain, obs=obs), preps["plain"], ("donor",))
    assert "donor" not in ctx.obs.columns
    assert "donor" in ctx.scoring_only
    # Without the name, a donor column is an ordinary annotation.
    assert "donor" in _ctx(dataclasses.replace(plain, obs=obs), preps["plain"]).obs.columns


def test_the_label_key_is_never_detected_on_a_scoring_only_column(
    plain: Dataset, preps: dict[str, PrepResult]
) -> None:
    obs = plain.obs.copy()
    obs["clusters"] = obs["true_label"].copy()  # a label candidate, named as a clone column
    assert _ctx(dataclasses.replace(plain, obs=obs), preps["plain"]).label_key == "clusters"
    ctx = _ctx(dataclasses.replace(plain, obs=obs), preps["plain"], ("clusters",))
    assert ctx.label_key == "lumencell_leiden"


def test_inference_obs_copies_and_reports_the_set() -> None:
    obs = pd.DataFrame({"clone_id": ["a", "b"], "donor": ["x", "y"], "cell_type": ["t", "u"], 7: [1, 2]})
    kept, dropped = inference_obs(obs, ("donor",))
    assert list(kept.columns) == ["cell_type", 7]
    assert dropped == frozenset({"clone_id", "donor"})
    assert list(obs.columns) == ["clone_id", "donor", "cell_type", 7], "the input is never changed"


def test_scoring_only_columns_match_by_name() -> None:
    columns = ["clone_id", "Clonotype", "CDR3", "IR_VDJ_1_c_call", "chain_pairing_x", "cell_type", 3]
    assert scoring_only_columns(columns) == frozenset({"clone_id", "Clonotype", "CDR3", "IR_VDJ_1_c_call"})
    assert scoring_only_columns(columns, ("donor",)) == frozenset(
        {"clone_id", "Clonotype", "CDR3", "IR_VDJ_1_c_call", "donor"}
    )
    for ordinary in ("cell_type", "tissue", "iroot", "chain_pairing_x", "receptor"):
        assert RECEPTOR_COLUMN.search(ordinary) is None, ordinary
    assert ISOLATION_VERSION == 1


# ---- refusals ----------------------------------------------------------------------------


def test_both_wordings_of_the_refusal() -> None:
    assert scoring_only_message("clone_id", "root column") == (
        "'clone_id' holds clone or receptor data, which lumenCell uses only to score runs; "
        "it cannot be the root column"
    )
    assert scoring_only_message("donor", "cluster column") == (
        "'donor' is set as a clone or donor column, which lumenCell uses only to score runs; "
        "it cannot be the cluster column"
    )


def test_scoring_only_cluster_key_and_label_key_are_refused(
    clonal: Dataset, preps: dict[str, PrepResult], monkeypatch: pytest.MonkeyPatch
) -> None:
    ctx = _ctx(clonal, preps["clonal"])

    def no_work(*_a: object, **_k: object) -> None:
        raise AssertionError("a refused column must fail before the root is resolved")

    monkeypatch.setattr(registry, "resolve_root", no_work)
    with pytest.raises(ScoringOnlyColumnError, match=r"holds clone or receptor data.*the cluster column"):
        run_method(ctx, "paga_dpt", ROOT, params={"cluster_key": "clone_id"})
    with pytest.raises(ScoringOnlyColumnError, match=r"holds clone or receptor data.*the label column"):
        run_method(ctx, "paga_dpt", ROOT, label_key="clone_id")
    # With terminals named, lineages are named in the cluster column: that is refused once.
    with pytest.raises(ScoringOnlyColumnError, match="the cluster column"):
        run_method(ctx, "palantir", ROOT, terminals=["A"], params={"cluster_key": "clone_id"})


def test_a_named_column_is_refused_in_its_own_words(plain: Dataset, preps: dict[str, PrepResult]) -> None:
    obs = plain.obs.copy()
    obs["donor"] = obs["true_label"].copy()
    ctx = _ctx(dataclasses.replace(plain, obs=obs), preps["plain"], ("donor",))
    with pytest.raises(ScoringOnlyColumnError, match="'donor' is set as a clone or donor column"):
        run_method(ctx, "paga_dpt", ROOT, params={"cluster_key": "donor"})
    with pytest.raises(ScoringOnlyColumnError, match="'donor' is set as a clone or donor column"):
        check_inference_column(ctx, "donor", "label column")
    check_inference_column(ctx, "true_label", "label column")  # an ordinary column passes


# ---- provenance --------------------------------------------------------------------------


def test_runs_record_what_steered_them(clonal: Dataset, preps: dict[str, PrepResult]) -> None:
    ctx = _ctx(clonal, preps["clonal"])
    assert ctx.label_key == "lumencell_leiden"  # the detected label key
    result = run_method(ctx, "paga_dpt", ROOT)
    assert result.obs_keys == sorted({"true_label", "lumencell_leiden"})
    named = run_method(ctx, "paga_dpt", ROOT, params={"cluster_key": "true_segment"}, label_key="true_label")
    assert named.obs_keys == ["true_label", "true_segment"]


# ---- the context guards itself -------------------------------------------------------------


def test_a_context_built_any_other_way_still_refuses_clone_columns(
    clonal: Dataset, preps: dict[str, PrepResult]
) -> None:
    ctx = _ctx(clonal, preps["clonal"])
    with pytest.raises(ScoringOnlyColumnError, match="'clone_id'"):
        dataclasses.replace(ctx, obs=clonal.obs.copy())  # a receptor-pattern column
    obs = ctx.obs.copy()
    obs["donor"] = "d1"
    with pytest.raises(ScoringOnlyColumnError, match="'donor'"):
        dataclasses.replace(ctx, obs=obs, scoring_only=frozenset({"donor"}))  # a named column
    dataclasses.replace(ctx, obs=obs)  # an ordinary column is fine


def test_a_bare_string_is_not_a_list_of_columns(plain: Dataset, preps: dict[str, PrepResult]) -> None:
    obs = plain.obs.copy()
    obs["donor"] = "d1"
    with pytest.raises(TypeError, match="not the string 'donor'"):
        scoring_only_columns(obs.columns, "donor")
    with pytest.raises(TypeError, match="not the string 'donor'"):
        inference_obs(obs, "donor")
    with pytest.raises(TypeError, match="not the string 'donor'"):
        InferenceContext.from_prep(dataclasses.replace(plain, obs=obs), preps["plain"], PrepParams(), "donor")


def test_a_named_column_the_file_lacks_is_refused(plain: Dataset) -> None:
    obs = plain.obs.copy()
    obs["donor"] = "d1"
    with pytest.raises(ValueError, match="'Donor' are not in the file's obs columns"):
        inference_obs(obs, ("Donor",))
