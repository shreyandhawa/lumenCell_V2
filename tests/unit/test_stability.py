import math
import statistics
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

import lumencell.engine.stability as stability
from lumencell.bench.synth import generate
from lumencell.data.loader import Dataset, load
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import run_method
from lumencell.engine.roots import CellRoot, ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from lumencell.engine.stability import (
    ReplicateRow,
    ReplicateScores,
    RunSpec,
    StabilityConfig,
    StabilityError,
    map_leiden_terminals,
    pad_result,
    replicate_plan,
    reseed_prep,
    restrict,
    run_stability,
    score_replicate,
    subsample_cells,
    subsample_prep,
    summarise,
)
from lumencell.prep.pipeline import PrepResult, run_pipeline
from lumencell.prep.stages import STAGE_BY_NAME, PrepParams, run_stage

ROOT = ClusterRoot("true_label", "trunk_early")


class Bif:
    def __init__(self, ds: Dataset, prep: PrepResult, params: PrepParams, reference: TrajectoryResult):
        self.ds, self.prep, self.params, self.reference = ds, prep, params, reference


@pytest.fixture(scope="module")
def bif(tmp_path_factory: pytest.TempPathFactory) -> Bif:
    path = tmp_path_factory.mktemp("stab") / "bif.h5ad"
    generate("bifurcation", "realistic", n_cells=600, n_genes=400, seed=0).write_h5ad(path)
    ds = load(path)
    params = PrepParams()
    prep = run_pipeline(ds, params)
    reference = run_method(InferenceContext.from_prep(ds, prep, params), "paga_dpt", ROOT)
    return Bif(ds, prep, params, reference)


class Spy:
    """Records every run_method call and forwards it to the real one (or raises on chosen calls)."""

    def __init__(self, fail_on: set[int] = frozenset()):  # type: ignore[assignment]
        self.calls: list[dict[str, Any]] = []
        self.fail_on = fail_on

    def __call__(
        self,
        ctx: InferenceContext,
        method: str,
        root: Any,
        terminals: Any = None,
        params: Any = None,
        **kw: Any,
    ) -> TrajectoryResult:
        self.calls.append({"root": root, "terminals": terminals, "params": params, "ctx": ctx, **kw})
        if len(self.calls) in self.fail_on:
            raise RuntimeError("boom")
        return run_method(ctx, method, root, terminals, params, **kw)


def test_replicate_plan_follows_the_declared_seeds() -> None:
    plan = replicate_plan(StabilityConfig(), 0, {"seed": 20})
    assert [(p.kind, p.index) for p in plan] == [
        ("reseed", 0),
        ("reseed", 1),
        ("subsample", 0),
        ("subsample", 1),
        ("subsample", 2),
    ]
    assert [(p.seed, p.method_seed) for p in plan] == [(1, 21), (2, 22), (1, 20), (2, 20), (3, 20)]
    assert all(p.method_seed is None for p in replicate_plan(StabilityConfig(), 0, {}))
    assert [p.seed for p in replicate_plan(StabilityConfig(), 5, {})][:2] == [6, 7]


def test_subsample_keeps_the_root_and_draws_the_fraction() -> None:
    a = subsample_cells(1000, 17, 0.8, 1)
    assert a.size == 800 and 17 in a
    assert np.array_equal(a, np.unique(a))
    assert np.array_equal(a, subsample_cells(1000, 17, 0.8, 1))
    assert not np.array_equal(a, subsample_cells(1000, 17, 0.8, 2))


def test_restrict_slices_the_cells_and_maps_the_stored_root(bif: Bif) -> None:
    ds = replace(bif.ds, uns_root=5, embeddings={"e": np.arange(1200, dtype=np.float32).reshape(600, 2)})
    cells = np.array([1, 5, 9], dtype=np.int64)
    sub = restrict(ds, cells)
    assert sub.n_obs == 3 and sub.n_vars == ds.n_vars and sub.matrix is None
    assert sub.obs.index.equals(ds.obs.index[cells])
    np.testing.assert_array_equal(sub.genes_detected, ds.genes_detected[cells])
    np.testing.assert_array_equal(sub.embeddings["e"], ds.embeddings["e"][cells])
    assert sub.uns_root == 1
    assert restrict(ds, np.array([0, 2], dtype=np.int64)).uns_root is None
    assert sub.fingerprint != ds.fingerprint
    assert restrict(ds, np.array([1, 6, 9], dtype=np.int64)).fingerprint != sub.fingerprint


def test_reseed_reruns_only_leiden_and_paga(bif: Bif) -> None:
    out, params = reseed_prep(bif.ds, bif.prep, bif.params, 1)
    assert params == replace(bif.params, seed=1)
    for key in ("pca", "connectivities", "diffmap", "hvg"):
        assert out.artifacts[key] is bif.prep.artifacts[key]
    expected, _ = run_stage(STAGE_BY_NAME["leiden"], bif.ds, bif.prep.artifacts, params)
    np.testing.assert_array_equal(out.artifacts["leiden"], expected["leiden"])
    assert out.artifacts["paga_tree"] is not bif.prep.artifacts["paga_tree"]


def test_subsample_prep_reuses_pca_and_recomputes_the_graph(bif: Bif) -> None:
    cells = subsample_cells(bif.ds.n_obs, 0, 0.8, 1)
    k = cells.size
    sub_ds, out = subsample_prep(bif.ds, bif.prep, bif.params, cells)
    np.testing.assert_array_equal(out.artifacts["pca"], bif.prep.artifacts["pca"][cells])
    np.testing.assert_array_equal(out.artifacts["hvg"], bif.prep.artifacts["hvg"])
    assert out.artifacts["connectivities"].shape == (k, k)
    assert out.artifacts["leiden"].shape == (k,)
    assert out.artifacts["diffmap"].shape[0] == k
    assert "lognorm" not in out.artifacts
    assert sub_ds.obs.index.equals(bif.ds.obs.index[cells])


def test_leiden_terminals_map_to_the_largest_overlap() -> None:
    ref = np.array([0, 0, 0, 1, 1, 2])
    rep = np.array([4, 4, 3, 3, 3])
    cells = np.arange(5)
    assert map_leiden_terminals(["0", "1"], ref, rep, cells) == ["4", "3"]
    assert map_leiden_terminals(["2"], ref, rep, cells) == []
    assert map_leiden_terminals(["0", "0"], ref, rep, cells) == ["4"]
    # both reference clusters land in replicate cluster 7: they collapse
    assert map_leiden_terminals(["0", "1"], ref, np.array([7, 7, 7, 7, 7]), cells) == ["7"]
    # a tie goes to the lowest label
    assert map_leiden_terminals(["1"], ref, np.array([0, 0, 0, 9, 2]), cells) == ["2"]


def test_pad_result_validates_and_matches_on_the_subset(bif: Bif) -> None:
    cells = subsample_cells(bif.ds.n_obs, bif.reference.root.cell_index, 0.8, 1)
    sub_ds, sub_prep = subsample_prep(bif.ds, bif.prep, bif.params, cells)
    ctx = InferenceContext.from_prep(sub_ds, sub_prep, bif.params)
    rep = run_method(ctx, "paga_dpt", CellRoot(bif.reference.root.obs_name))
    padded = pad_result(rep, cells, bif.ds.n_obs)
    padded.validate()
    assert padded.n_obs == bif.ds.n_obs
    assert padded.root.cell_index == bif.reference.root.cell_index
    outside = np.setdiff1d(np.arange(bif.ds.n_obs), cells)
    assert np.all(np.isnan(padded.pseudotime[outside]))
    assert not padded.membership[:, outside].any()
    assert np.all(padded.segment[outside] == -1)
    np.testing.assert_array_equal(padded.pseudotime[cells], rep.pseudotime)
    np.testing.assert_array_equal(padded.membership[:, cells], rep.membership)
    assert padded.milestone_weights.shape == (bif.ds.n_obs, len(rep.milestones))


def test_a_replicate_identical_to_the_reference_scores_one(bif: Bif) -> None:
    r = bif.reference
    n_finite = int(np.isfinite(r.pseudotime).sum())
    got = score_replicate(r, r, np.ones(r.n_obs, dtype=bool), 0.6)
    # spearmanr of a vector with itself is 1 to rounding (0.9999999999999999 here)
    assert replace(got, spearman=1.0) == ReplicateScores(n_finite, len(r.lineages), 1.0, 1.0, 1.0)
    assert got.spearman == pytest.approx(1.0, abs=1e-12)


def test_undefined_scores_are_nan(bif: Bif) -> None:
    r = bif.reference
    shared = np.zeros(r.n_obs, dtype=bool)
    shared[np.flatnonzero(np.isfinite(r.pseudotime))[:2]] = True
    got = score_replicate(r, r, shared, 0.6)
    assert got.n_compared == 2 and math.isnan(got.spearman)


def _spec(**kw: Any) -> RunSpec:
    return RunSpec(**{"method": "paga_dpt", "params": {}, "terminals": None, "label_key": None, **kw})


def test_every_replicate_is_rooted_at_the_reference_root_cell(
    bif: Bif, monkeypatch: pytest.MonkeyPatch
) -> None:
    spy = Spy()
    monkeypatch.setattr(stability, "run_method", spy)
    seen: list[tuple[str, int, int]] = []
    run_stability(bif.ds, bif.prep, bif.params, bif.reference, _spec(), progress=lambda *a: seen.append(a))
    assert [c["root"] for c in spy.calls] == [CellRoot(bif.reference.root.obs_name)] * 5
    assert seen == [
        ("reseed 1", 0, 5),
        ("reseed 2", 1, 5),
        ("subsample 1", 2, 5),
        ("subsample 2", 3, 5),
        ("subsample 3", 4, 5),
    ]


def test_user_column_terminals_pass_through(bif: Bif, monkeypatch: pytest.MonkeyPatch) -> None:
    spy = Spy()
    monkeypatch.setattr(stability, "run_method", spy)
    spec = _spec(params={"cluster_key": "true_label"}, terminals=["A_late", "B_late"])
    res = run_stability(bif.ds, bif.prep, bif.params, bif.reference, spec)
    assert [c["terminals"] for c in spy.calls] == [["A_late", "B_late"]] * 5
    assert all(row.terminals == ("A_late", "B_late") for row in (*res.reseeds.rows, *res.subsamples.rows))


def test_leiden_terminals_are_mapped_per_replicate(bif: Bif, monkeypatch: pytest.MonkeyPatch) -> None:
    spy = Spy()
    monkeypatch.setattr(stability, "run_method", spy)
    ref_leiden = np.asarray(bif.prep.artifacts["leiden"])
    terminal = str(ref_leiden[bif.ds.obs["true_label"].astype(str).to_numpy() == "A_late"][0])
    run_stability(bif.ds, bif.prep, bif.params, bif.reference, _spec(terminals=[terminal]))
    for call in spy.calls:
        ctx = call["ctx"]
        names = ctx.obs_names
        cells = np.flatnonzero(bif.ds.obs.index.astype(str).isin(names))
        expected = map_leiden_terminals([terminal], ref_leiden, ctx.leiden, cells)
        assert call["terminals"] == expected


def test_clonal_columns_reach_every_replicate(bif: Bif, monkeypatch: pytest.MonkeyPatch) -> None:
    obs = bif.ds.obs.copy()
    obs["donor"] = "d1"
    ds = replace(bif.ds, obs=obs)
    spy = Spy()
    monkeypatch.setattr(stability, "run_method", spy)
    run_stability(ds, bif.prep, bif.params, bif.reference, _spec(clonal_columns=("donor",)))
    assert len(spy.calls) == 5
    assert all("donor" in c["ctx"].scoring_only and "donor" not in c["ctx"].obs.columns for c in spy.calls)


def test_a_failing_replicate_counts_as_no_agreement(bif: Bif, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(stability, "run_method", Spy(fail_on={2}))
    res = run_stability(bif.ds, bif.prep, bif.params, bif.reference, _spec())
    bad = res.reseeds.rows[1]
    assert bad.error == "RuntimeError: boom"
    assert res.n_failed == 1
    assert math.isnan(bad.share_matched) and math.isnan(bad.spearman) and bad.n_lineages == 0
    assert res.reseeds.share_matched_min == 0.0
    assert res.reseeds.n_lineages_min == 0
    shares = [
        0.0 if math.isnan(r.share_matched) else r.share_matched
        for r in (*res.reseeds.rows, *res.subsamples.rows)
    ]
    assert len(shares) == 5
    assert res.lineage_stability == statistics.median(shares)
    assert all(r.error is None for r in (res.reseeds.rows[0], *res.subsamples.rows))


def test_summaries_count_undefined_values_as_zero() -> None:
    def row(share: float, n_lin: int = 2, error: str | None = None) -> ReplicateRow:
        return ReplicateRow("subsample", 0, 1, None, 10, 10, n_lin, share, share, share, None, error)

    s = summarise("subsample", 0.8, [row(1.0), row(math.nan, 0, "RuntimeError: x"), row(1.0)])
    assert s.share_matched_median == 1.0 and s.share_matched_min == 0.0
    assert s.spearman_min == 0.0 and s.mean_jaccard_min == 0.0
    assert s.n_lineages_median == 2.0 and s.n_lineages_min == 0
    assert s.r == 3 and s.fraction == 0.8 and s.seeds == (1, 1, 1)


def test_a_repeated_root_name_is_refused(bif: Bif) -> None:
    names = bif.ds.obs.index.astype(str).to_numpy().copy()
    root = bif.reference.root.cell_index
    other = 0 if root != 0 else 1
    names[other] = names[root]
    obs = bif.ds.obs.copy()
    obs.index = pd.Index(names)
    with pytest.raises(StabilityError, match="more than once"):
        run_stability(replace(bif.ds, obs=obs), bif.prep, bif.params, bif.reference, _spec())


def test_stability_runs_two_reseeds_and_three_subsamples(bif: Bif) -> None:
    res = run_stability(bif.ds, bif.prep, bif.params, bif.reference, _spec())
    assert res.reseeds.r == 2 and res.subsamples.r == 3
    assert res.subsamples.fraction == 0.8 and res.reseeds.fraction == 1.0
    assert res.reseeds.seeds == (1, 2) and res.subsamples.seeds == (1, 2, 3)
    assert all(row.n_cells == 480 for row in res.subsamples.rows)
    assert all(row.n_cells == 600 for row in res.reseeds.rows)
    assert all(row.error is None for row in (*res.reseeds.rows, *res.subsamples.rows))
    assert res.n_obs == 600 and res.reference_lineages == len(bif.reference.lineages)
    assert 0 <= res.lineage_stability <= 1
    assert -1 <= res.pseudotime_stability <= 1


def test_replicates_never_touch_the_disk(bif: Bif, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    run_stability(
        bif.ds, bif.prep, bif.params, bif.reference, _spec(), StabilityConfig(n_reseeds=1, n_subsamples=1)
    )
    assert list(tmp_path.iterdir()) == []


def test_a_method_default_seed_is_reseeded(bif: Bif, monkeypatch: pytest.MonkeyPatch) -> None:
    # palantir's seed defaults to 20; a run that left it unset still has its method reseeded.
    calls: list[dict[str, Any]] = []

    def fake(
        ctx: InferenceContext, method: str, root: Any, terminals: Any = None, params: Any = None, **kw: Any
    ) -> TrajectoryResult:
        calls.append(dict(params or {}))
        return run_method(ctx, "paga_dpt", root, terminals, {}, **kw)

    monkeypatch.setattr(stability, "run_method", fake)
    res = run_stability(bif.ds, bif.prep, bif.params, bif.reference, _spec(method="palantir"))
    assert [c.get("seed") for c in calls] == [21, 22, 20, 20, 20]
    assert [r.method_seed for r in res.reseeds.rows] == [21, 22]


def test_paga_dpt_reseeds_vary_only_prep(bif: Bif, monkeypatch: pytest.MonkeyPatch) -> None:
    spy = Spy()
    monkeypatch.setattr(stability, "run_method", spy)
    res = run_stability(bif.ds, bif.prep, bif.params, bif.reference, _spec(), StabilityConfig(n_subsamples=0))
    assert all("seed" not in (c["params"] or {}) for c in spy.calls)
    assert [r.method_seed for r in res.reseeds.rows] == [None, None]
    assert [c["ctx"].prep_params.seed for c in spy.calls] == [1, 2]


def test_a_run_where_every_replicate_fails_counts_them_all(bif: Bif, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(stability, "run_method", Spy(fail_on={1, 2, 3, 4, 5}))
    res = run_stability(bif.ds, bif.prep, bif.params, bif.reference, _spec())
    assert res.n_failed == 5
    assert res.lineage_stability == 0.0 and res.pseudotime_stability == 0.0
    assert all(r.error is not None for r in (*res.reseeds.rows, *res.subsamples.rows))


def test_the_leiden_terminal_test_reads_the_methods_default_cluster_key(
    bif: Bif, monkeypatch: pytest.MonkeyPatch
) -> None:
    from lumencell.engine.context import LEIDEN_KEY

    for m in ("paga_dpt", "palantir", "scfates_ppt"):
        assert stability._terminal_key_default(m) == LEIDEN_KEY
    monkeypatch.setattr(stability, "_terminal_key_default", lambda _m: "true_label")
    spy = Spy()
    monkeypatch.setattr(stability, "run_method", spy)
    run_stability(bif.ds, bif.prep, bif.params, bif.reference, _spec(terminals=["A_late"]))
    assert [c["terminals"] for c in spy.calls] == [["A_late"]] * 5  # not mapped as Leiden labels


def test_subsample_cells_use_seed_plus_one_plus_j() -> None:
    # root 2 of 10, half kept; default_rng(4) and default_rng(5) give these draws
    assert subsample_cells(10, 2, 0.5, 4).tolist() == [2, 5, 7, 8, 9]
    assert subsample_cells(10, 2, 0.5, 5).tolist() == [0, 2, 5, 6, 8]


def test_each_subsample_replicate_draws_with_seed_plus_one_plus_j(
    bif: Bif, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[np.ndarray] = []
    real = stability.subsample_prep

    def spy(ds: Dataset, prep: PrepResult, params: PrepParams, cells: np.ndarray) -> Any:
        seen.append(cells)
        return real(ds, prep, params, cells)

    monkeypatch.setattr(stability, "subsample_prep", spy)
    config = StabilityConfig(n_reseeds=0, n_subsamples=2, seed=3)
    run_stability(bif.ds, bif.prep, bif.params, bif.reference, _spec(), config)
    root = bif.reference.root.cell_index
    others = np.delete(np.arange(bif.ds.n_obs), root)
    for j, cells in enumerate(seen):
        drawn = np.random.default_rng(3 + 1 + j).choice(
            others, size=round(0.8 * bif.ds.n_obs) - 1, replace=False
        )
        assert cells.tolist() == sorted([*drawn.tolist(), root])
    assert len(seen) == 2


class _R:
    def __init__(self, n_lineages: int) -> None:
        self.lineages = list(range(n_lineages))


def _fake_compare(jaccards: list[float]) -> Any:
    from lumencell.engine.compare import LineagePair, RunComparison

    pairs = [LineagePair(i, i, j) for i, j in enumerate(jaccards)]
    return lambda a, b, shared: RunComparison(10, 0.5, pairs, np.zeros(1, dtype=np.float32))


def test_share_matched_divides_by_the_reference_lineage_count(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(stability, "compare", _fake_compare([0.9, 0.1]))
    got = score_replicate(_R(2), _R(3), np.ones(10, dtype=bool), 0.6)  # type: ignore[arg-type]
    assert got.share_matched == 0.5 and got.n_lineages == 3


def test_a_jaccard_equal_to_the_threshold_is_matched(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(stability, "compare", _fake_compare([0.6, 0.5]))
    got = score_replicate(_R(2), _R(2), np.ones(10, dtype=bool), 0.6)  # type: ignore[arg-type]
    assert got.share_matched == 0.5


def test_pad_result_leaves_every_outside_array_unassigned(bif: Bif) -> None:
    cells = subsample_cells(bif.ds.n_obs, bif.reference.root.cell_index, 0.8, 1)
    sub_ds, sub_prep = subsample_prep(bif.ds, bif.prep, bif.params, cells)
    rep = run_method(
        InferenceContext.from_prep(sub_ds, sub_prep, bif.params),
        "paga_dpt",
        CellRoot(bif.reference.root.obs_name),
    )
    padded = pad_result(rep, cells, bif.ds.n_obs)
    outside = np.setdiff1d(np.arange(bif.ds.n_obs), cells)
    assert outside.size
    assert np.all(np.isnan(padded.progress[..., outside]))
    assert not padded.membership[:, outside].any()
    assert np.all(padded.segment[outside] == -1)
    if padded.fate is not None:
        assert np.all(np.isnan(padded.fate[..., outside]))
    for v in padded.extras.values():
        assert np.all(np.isnan(v[..., outside]))
