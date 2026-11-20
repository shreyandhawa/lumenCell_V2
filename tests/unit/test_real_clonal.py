import json
import math
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any, get_args

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

import lumencell.bench.real_clonal as real_clonal_mod
from lumencell.bench import report, runner
from lumencell.bench.datasets import DATASETS, RealClonal, RealDataset, truth_for
from lumencell.bench.download import Source, sha256_of
from lumencell.bench.gates import REPO, check_gates, load_thresholds
from lumencell.bench.metrics import Scores
from lumencell.bench.real_clonal import adjusted_rand_index, clone_recurrence, reported_metrics
from lumencell.bench.runner import CaseResult, ClonalRow, ClonalRowStatus, run_real_case
from lumencell.bench.synth import generate
from lumencell.bench.synth_clones import DEFAULT_LAYER, add_clones
from lumencell.bench.truth import Truth
from lumencell.cli import app
from lumencell.engine import stability as stability_mod
from lumencell.engine.context import InferenceContext

BASE_GATES = load_thresholds()["real"]["daniel2022"]["gates"]
RESULTS_PATH = REPO / "benchmarks/results/real-test.json"
cli = CliRunner()


def _raise(*a: Any, **k: Any) -> None:
    raise AssertionError("nothing may be scored")


@pytest.fixture
def fake(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> str:
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path / "datasets"))
    monkeypatch.setenv("LUMENCELL_CACHE_DIR", str(tmp_path / "cache"))
    (tmp_path / "datasets").mkdir()
    base = generate("bifurcation", "realistic", n_cells=1500, n_genes=1000, seed=0)
    adata = add_clones(base, DEFAULT_LAYER, seed=0)
    obs = adata.obs
    obs["state"] = obs["true_label"].astype(str)
    obs["capture"] = obs["donor"].astype(str)
    clone = obs["clone_id"].to_numpy(dtype=object)
    trb = clone.copy()
    trb[::5] = np.nan
    obs["clone_id_trb"] = trb
    has = ~pd.isna(clone)
    raw = np.full(len(obs), np.nan, dtype=object)
    raw[has] = obs["capture"].to_numpy(dtype=object)[has] + ":" + clone[has].astype(str)
    obs["raw_clonotype_id"] = raw
    obs["receptor_n_tra"] = 1
    obs["receptor_n_trb"] = 1
    t = obs["true_pseudotime"].to_numpy()
    obs["day"] = np.digitize(t, np.quantile(t, [1 / 3, 2 / 3])).choose([5, 8, 21])
    path = tmp_path / "datasets" / "fake_tcell.h5ad"
    adata.write_h5ad(path)
    paths = base.uns["lumencell_truth"]["lineages"]
    lineages = {k: [f"{s}_{p}" for s in v for p in ("early", "late")] for k, v in paths.items()}
    ds = RealDataset(
        id="fake_tcell",
        title="fake",
        citation="none",
        sources=(
            Source("https://example.org/f.h5ad", "fake_tcell.h5ad", path.stat().st_size, sha256_of(path)),
        ),
        licence="test",
        label_key="state",
        root_label="trunk_early",
        lineages=lineages,
        gated=list(lineages),
        clonal=RealClonal(
            "clone_id", "donor", "clone_id_trb", "raw_clonotype_id", "capture", True, "test", "2026-01-02"
        ),
        report_time_key="day",
    )
    monkeypatch.setitem(DATASETS, "fake_tcell", ds)
    return "fake_tcell"


def test_a_real_dataset_with_clones_is_scored_in_oracle_cells_mode(fake: str) -> None:
    row = run_real_case(fake, "paga_dpt").row()
    assert row["clonal_status"] in get_args(ClonalRowStatus) and row["clonal_status"] != "not_applicable"
    assert math.isfinite(row["clonal_eta_oracle"])
    assert {
        "reported_day_spearman",
        "reported_clone_recurrence",
        "reported_raw_clonotype_ari",
        "reported_labelled_share_of_raw",
        "reported_trb_status",
    } <= row.keys()
    assert {
        "clonal_n_cells",
        "clonal_effective_clones",
        "clonal_eta_m3",
        "clonal_p_m3",
    } <= row.keys()  # phase 1, reported only


def test_every_clone_column_of_a_built_dataset_is_scoring_only(
    fake: str, monkeypatch: pytest.MonkeyPatch
) -> None:  # Review Focus 4
    seen: dict[str, InferenceContext] = {}
    real = runner.run_method
    monkeypatch.setattr(
        runner, "run_method", lambda ctx, *a, **k: (seen.setdefault("ctx", ctx), real(ctx, *a, **k))[1]
    )
    run_real_case(fake, "paga_dpt")
    cols = {"clone_id", "clone_id_trb", "raw_clonotype_id", "receptor_n_tra", "receptor_n_trb", "donor"}
    assert cols <= seen["ctx"].scoring_only and not cols & set(seen["ctx"].obs.columns)
    assert {"state", "day", "capture"} <= set(seen["ctx"].obs.columns)


@pytest.mark.slow
def test_stability_replicates_keep_clone_columns_out_too(fake: str, monkeypatch: pytest.MonkeyPatch) -> None:
    ctxs: list[InferenceContext] = []
    real = stability_mod.run_method
    monkeypatch.setattr(
        stability_mod, "run_method", lambda ctx, *a, **k: (ctxs.append(ctx), real(ctx, *a, **k))[1]
    )
    run_real_case(fake, "paga_dpt", stability=True)
    assert ctxs and all(
        {"clone_id_trb", "donor"} <= c.scoring_only and "donor" not in c.obs.columns for c in ctxs
    )


def test_adjusted_rand_index_by_hand() -> None:
    assert adjusted_rand_index(np.array([0, 0, 1, 1]), np.array([0, 0, 1, 2])) == pytest.approx(4 / 7)
    assert adjusted_rand_index(np.array([0, 0, 1, 1]), np.array([5, 5, 7, 7])) == 1.0
    assert math.isnan(adjusted_rand_index(np.array([0]), np.array([0])))


def test_clone_recurrence_by_hand() -> None:
    share, n = clone_recurrence(
        np.array(["X", "X", "Y", "Z", "Z", None], object), np.array(["A", "B", "A", "B", "B", "A"], object)
    )
    assert share == pytest.approx(1 / 3) and n == 3


def test_a_declared_not_assessable_result_is_shown_not_omitted(
    fake: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    ds = DATASETS[fake]
    monkeypatch.setitem(
        DATASETS,
        fake,
        replace(ds, clonal=replace(ds.clonal, assessable=False, reason="r", declared="2026-01-02")),
    )
    monkeypatch.setattr(real_clonal_mod, "clonal_row", _raise)  # nothing may be scored
    row = run_real_case(fake, "paga_dpt").row()
    assert (row["clonal_status"], row["clonal_reason"]) == ("not_assessable", "declared 2026-01-02, r")
    assert row["reported_trb_status"] == "not_assessable"
    assert '"clonal_reason": "declared 2026-01-02, r"' in json.dumps(row)  # the results JSON
    assert not [f for f in check_gates(row, BASE_GATES) if f.metric.startswith("clonal")]  # not gated
    monkeypatch.setattr(
        report,
        "latest_real",
        lambda: (RESULTS_PATH, {"commit": "x", "methods": {"paga_dpt": {"daniel2022": row}}}),
    )
    assert "clonal: not assessable (declared 2026-01-02, r)" in report.render()  # the README table


def test_day_spearman_is_nan_with_one_day(fake: str, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    ds = DATASETS[fake]
    obs = pd.DataFrame({"day": [8] * 10, "clone_id": [None] * 10})
    result: Any = type("R", (), {"pseudotime": np.linspace(0, 1, 10)})
    truth: Any = None
    out = real_clonal_mod.reported_metrics(replace(ds, clonal=None), obs, result, truth, np.zeros((10, 2)), 2)
    assert math.isnan(out["reported_day_spearman"])


def _daniel_row(**kw: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "case": "daniel2022",
        "method": "paga_dpt",
        "pt_spearman_lineage": math.nan,
        "branch_accuracy": math.nan,
        "lineage_recall": math.nan,
        "terminal_share": math.nan,
        "pt_spearman_global": math.nan,
        "clonal_status": "not_assessable",
        "clonal_reason": "declared 2026-11-20, because",
        "reported_day_spearman": 0.31,
        "reported_clone_recurrence": 0.1,
        "reported_n_clone_strings": 40,
        "reported_trb_status": "not_assessable",
    }
    row.update(kw)
    return row


def _patch_results(monkeypatch: pytest.MonkeyPatch, row: dict[str, Any]) -> None:
    monkeypatch.setattr(
        report,
        "latest_real",
        lambda: (RESULTS_PATH, {"commit": "x", "methods": {"paga_dpt": {"daniel2022": row}}}),
    )


def test_the_readme_shows_reported_fields(monkeypatch: pytest.MonkeyPatch) -> None:
    row = _daniel_row(clonal_status="ok", clonal_excess_ratio=0.8, clonal_p=0.001, reported_trb_status="ok")
    _patch_results(monkeypatch, row)
    text = report.render()
    assert "Reported, not gated" in text and "clonal: 0.80" in text and "η at m = 3" in text


NOTE = "lineage gates: not measured (no gated lineage, declared 2026-11-20)"
_FAIL_ROW = {"pt_spearman_lineage": 0.1}  # a real gate failure (needs >= 0.70)


def _daniel_cell(monkeypatch: pytest.MonkeyPatch, **kw: Any) -> str:
    _patch_results(monkeypatch, _daniel_row(**kw))
    return next(x for x in report.render().splitlines() if x.startswith("| PAGA + DPT |") and "Spearman" in x)


def _run_cli(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, row: dict[str, Any]) -> str:
    monkeypatch.setattr("lumencell.bench.datasets.is_cached", lambda d: True)
    monkeypatch.setattr(
        "lumencell.bench.runner.run_real_case",
        lambda *a, **k: type(
            "C", (), {"row": lambda self: dict(row, runtime_seconds=1.0, tree_support="not_measured")}
        )(),
    )
    res = cli.invoke(
        app,
        [
            "bench",
            "real",
            "--datasets",
            "daniel2022",
            "--methods",
            "paga_dpt",
            "--out",
            str(tmp_path / "r.json"),
        ],
    )
    return res.output


def test_p14_readme_cell_without_gated_lineage_says_not_measured_and_never_pass(
    monkeypatch: pytest.MonkeyPatch,
) -> None:  # controller ruling P14
    ds = DATASETS["daniel2022"]
    assert ds.gated == [] and ds.clonal is not None and ds.clonal.declared == "2026-11-20"
    cell = _daniel_cell(monkeypatch)
    assert NOTE in cell and "clonal: not assessable (declared 2026-11-20, because)" in cell
    assert "PASS" not in cell and "pass" not in cell


def test_p14_cli_line_without_gated_lineage_says_not_measured_and_never_pass(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    output = _run_cli(monkeypatch, tmp_path, _daniel_row())
    line = next(x for x in output.splitlines() if x.startswith("daniel2022"))
    assert NOTE in line and "clonal=not_assessable" in line
    assert "PASS" not in line and "pass" not in line


def test_p14_never_hides_a_real_failure_in_the_readme_cell(monkeypatch: pytest.MonkeyPatch) -> None:
    cell = _daniel_cell(monkeypatch, **_FAIL_ROW)
    assert "fail: pt_spearman_lineage" in cell
    assert NOTE not in cell


def test_p14_never_hides_a_real_failure_in_the_cli_line(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from lumencell.bench.gates import GateFailure

    monkeypatch.setattr(
        "lumencell.bench.gates.check_gates",
        lambda row, gates: [GateFailure("daniel2022", "pt_spearman_lineage", 0.1, ">= 0.7")],
    )
    output = _run_cli(monkeypatch, tmp_path, _daniel_row())
    line = next(x for x in output.splitlines() if x.startswith("daniel2022"))
    assert "FAIL" in line and NOTE not in line
    assert "gate     paga_dpt: daniel2022: pt_spearman_lineage" in output


def test_p14_without_a_clonal_declaration_there_is_no_date_to_cite(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(DATASETS, "daniel2022", replace(DATASETS["daniel2022"], clonal=None))
    cell = _daniel_cell(monkeypatch)
    assert "lineage gates: not measured (no gated lineage)" in cell and "pass" not in cell


def test_p14_a_dataset_with_gated_lineages_has_no_note() -> None:
    pancreas = DATASETS["pancreas"]
    assert pancreas.gated and report.lineage_gates_note(pancreas) is None


# ---- Fix round 1 -----------------------------------------------------------------------


def _mini_ds(**kw: Any) -> RealDataset:
    base = RealDataset(
        id="mini",
        title="mini",
        citation="none",
        sources=(Source("https://example.org/m.h5ad", "m.h5ad", 1, "0" * 64),),
        licence="test",
        label_key="state",
        root_label="a",
        lineages={"a": ["a"], "b": ["b"]},
        gated=["a", "b"],
        clonal=RealClonal("clone", "donor", "trb", "raw", "capture", True, "t", "2026-01-02"),
    )
    return replace(base, **kw)


def _mini_truth(n: int) -> Truth:
    mask = np.arange(n) % 2 == 0
    return Truth(order=np.arange(n, dtype=float), lineages={"a": mask, "b": ~mask}, root=mask)


def _metrics(ds: RealDataset, obs: pd.DataFrame, pt: Any = None) -> dict[str, Any]:
    n = len(obs)
    pt = np.linspace(0, 1, n) if pt is None else np.asarray(pt, dtype=float)
    return reported_metrics(ds, obs, SimpleNamespace(pseudotime=pt), _mini_truth(n), None, 2)  # type: ignore[arg-type]


def _obs(**cols: Any) -> pd.DataFrame:
    n = len(next(iter(cols.values())))
    base: dict[str, Any] = {k: [None] * n for k in ("clone", "donor", "trb", "raw", "capture")}
    base.update(cols)
    return pd.DataFrame({k: pd.Series(v, dtype=object) for k, v in base.items()})


@pytest.fixture
def seen_trb(monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    seen: dict[str, Any] = {}

    def fake_row(ds: Any, truth: Truth, *a: Any) -> ClonalRow:
        seen["truth"] = truth
        return ClonalRow("ok", excess_ratio=0.5, p=0.02)

    monkeypatch.setattr(real_clonal_mod, "real_clonal_row", fake_row)
    return seen


# day parsing (ruling: numeric columns are used directly)


def test_a_numeric_day_column_is_used_as_numbers_not_parsed_from_text() -> None:
    from scipy.stats import spearmanr

    day = [1e-05, 3e-05, 0.5, np.nan, 0.7]
    pt = [0.1, 0.2, 0.3, 0.4, 0.5]
    ds = _mini_ds(clonal=None, report_time_key="day")
    out = _metrics(ds, pd.DataFrame({"day": day}), pt)
    assert out["reported_day_spearman"] == pytest.approx(
        spearmanr([0.1, 0.2, 0.3, 0.5], [1e-05, 3e-05, 0.5, 0.7]).statistic
    )
    assert out["reported_day_spearman"] == pytest.approx(1.0)


def test_a_string_day_column_reads_the_first_number_and_skips_unparsable_cells() -> None:
    from scipy.stats import spearmanr

    ds = _mini_ds(clonal=None, report_time_key="day")
    obs = pd.DataFrame({"day": pd.Series(["D5", "D8", "D21", None], dtype=object)})
    out = _metrics(ds, obs, [0.1, 0.5, 0.9, 0.4])
    assert out["reported_day_spearman"] == pytest.approx(spearmanr([0.1, 0.5, 0.9], [5, 8, 21]).statistic)
    assert out["reported_day_spearman"] > 0


def test_an_all_unparsable_day_column_is_nan() -> None:
    ds = _mini_ds(clonal=None, report_time_key="day")
    out = _metrics(ds, pd.DataFrame({"day": pd.Series(["a", "b", "c"], dtype=object)}))
    assert math.isnan(out["reported_day_spearman"])


# P16: the declared labelled share


def test_the_labelled_share_is_clone_cells_over_raw_cells_as_declared_and_can_exceed_one(
    seen_trb: dict[str, Any],
) -> None:
    obs = _obs(
        clone=["x", "x", None, "y", "w", "z"],
        raw=["r1", "r1", "r2", None, None, "r3"],
        capture=["A"] * 6,
        donor=["d"] * 6,
        trb=[None] * 6,
    )
    out = _metrics(_mini_ds(), obs)
    assert out["reported_labelled_share_of_raw"] == pytest.approx(5 / 4)  # 5 clone cells, 4 raw cells


def test_the_labelled_share_is_nan_without_raw_cells(seen_trb: dict[str, Any]) -> None:
    obs = _obs(clone=["x", "y"], raw=[None, None], capture=["A", "A"], donor=["d", "d"])
    assert math.isnan(_metrics(_mini_ds(), obs)["reported_labelled_share_of_raw"])


# raw clonotype ARI


def test_the_raw_ari_is_keyed_by_capture_and_restricted_to_cells_with_both_ids(
    seen_trb: dict[str, Any],
) -> None:
    # "c" is the clone string in both captures; the raw ids differ by capture. Keyed by capture the
    # partitions agree; un-keyed, clone_id is one cluster against two raw ids. The last two cells
    # have only one of the two ids and would break the agreement if they were counted.
    obs = _obs(
        clone=["c", "c", "c", "c", "c", None],
        raw=["A1", "A1", "B1", "B1", None, "A1"],
        capture=["A", "A", "B", "B", "B", "A"],
        donor=["d"] * 6,
    )
    assert _metrics(_mini_ds(), obs)["reported_raw_clonotype_ari"] == 1.0


# TRB-only


def test_trb_only_scores_the_trb_column_with_the_donor_column(seen_trb: dict[str, Any]) -> None:
    obs = _obs(
        clone=["x", "y", "z", "w"],
        trb=["t1", None, "t3", "t4"],
        donor=["m1", "m1", "m2", "m2"],
        capture=["c1", "c2", "c3", "c4"],
        raw=["r"] * 4,
    )
    out = _metrics(_mini_ds(), obs)
    truth = seen_trb["truth"]
    assert list(truth.clone_ids) == ["t1", None, "t3", "t4"]
    assert list(truth.donors) == ["m1", "m1", "m2", "m2"]
    assert (out["reported_trb_status"], out["reported_trb_excess_ratio"], out["reported_trb_p"]) == (
        "ok",
        0.5,
        0.02,
    )


def test_recurrence_without_any_clone_is_nan_and_zero(seen_trb: dict[str, Any]) -> None:
    obs = _obs(clone=[None, None], capture=["A", "B"], raw=[None, None], donor=["d", "d"])
    out = _metrics(_mini_ds(), obs)
    assert math.isnan(out["reported_clone_recurrence"]) and out["reported_n_clone_strings"] == 0


def test_clone_recurrence_drops_cells_whose_capture_is_missing() -> None:
    share, n = clone_recurrence(np.array(["X", "X", "Y"], object), np.array(["A", np.nan, "A"], object))
    assert share == 0.0 and n == 2  # X is in one capture, not in "A" and a capture called "nan"


# ARI edge cases


def test_ari_of_two_single_clusters_and_of_two_all_singletons_is_one() -> None:
    assert adjusted_rand_index(np.zeros(4, np.int64), np.full(4, 3, np.int64)) == 1.0
    assert adjusted_rand_index(np.arange(4), np.arange(4)[::-1].copy()) == 1.0


def test_ari_of_one_cluster_against_all_singletons_is_zero() -> None:
    assert adjusted_rand_index(np.zeros(4, np.int64), np.arange(4)) == 0.0


def test_ari_of_a_partial_overlap_is_by_hand() -> None:
    # a = {0,1},{2}; b = {0},{1,2}: index 0, rows 1, cols 1, total 3, expected 1/3, max 1.
    assert adjusted_rand_index(np.array([0, 0, 1]), np.array([0, 1, 1])) == pytest.approx(-0.5)


def test_ari_of_two_cells_identical_up_to_relabelling_is_one() -> None:
    assert adjusted_rand_index(np.array([0, 1]), np.array([5, 7])) == 1.0
    assert adjusted_rand_index(np.array([0, 0]), np.array([4, 4])) == 1.0


def test_ari_with_a_zero_denominator_always_means_identical_partitions() -> None:
    # expected == maximum only when both are one cluster or both are all singletons, so the
    # "not identical" outcome of the zero-denominator branch cannot occur.
    for n in range(2, 7):
        for a_kind in ("one", "singletons"):
            a = np.zeros(n, np.int64) if a_kind == "one" else np.arange(n)
            for b_kind in ("one", "singletons"):
                b = np.zeros(n, np.int64) if b_kind == "one" else np.arange(n)
                value = adjusted_rand_index(a, b)
                assert value == (1.0 if a_kind == b_kind else 0.0)


# CaseResult.row


def test_row_rounds_reported_floats_to_four_places_and_keeps_nan_strings_and_ints() -> None:
    fields = {f: 0.0 for f in Scores.__dataclass_fields__}
    fields.update(n_spurious=0)
    reported = {
        "reported_a": 0.123456,
        "reported_nan": math.nan,
        "reported_s": "ok",
        "reported_n": 40,
    }
    row = CaseResult("c", "m", Scores(**fields), 3, 1.0, reported=reported).row()
    assert row["reported_a"] == 0.1235
    assert math.isnan(row["reported_nan"])
    assert row["reported_s"] == "ok" and row["reported_n"] == 40 and isinstance(row["reported_n"], int)


# truth_for


def test_truth_for_takes_clone_ids_and_donors_from_their_own_columns() -> None:
    ds = _mini_ds()
    obs = pd.DataFrame(
        {
            "state": ["a", "b", "a", "b"],
            "clone": ["x", "y", None, "z"],
            "donor": ["m1", "m1", "m2", "m2"],
            "capture": ["c1", "c2", "c3", "c4"],
        }
    )
    t = truth_for(ds, obs)
    assert t.clone_ids is not None and t.donors is not None
    assert list(t.clone_ids[:2]) == ["x", "y"] and pd.isna(t.clone_ids[2]) and t.clone_ids[3] == "z"
    assert list(t.donors) == ["m1", "m1", "m2", "m2"]


def test_truth_for_has_no_clones_when_not_assessable_or_without_clonal() -> None:
    obs = pd.DataFrame({"state": ["a", "b"], "clone": ["x", "y"], "donor": ["m", "m"]})
    declared = _mini_ds()
    na = replace(declared, clonal=replace(declared.clonal, assessable=False))
    for ds in (na, _mini_ds(clonal=None)):
        t = truth_for(ds, obs)
        assert t.clone_ids is None and t.donors is None


# report variants


def test_a_clonal_status_other_than_ok_is_shown_as_the_status(monkeypatch: pytest.MonkeyPatch) -> None:
    assert "clonal: too_few_clones" in _daniel_cell(monkeypatch, clonal_status="too_few_clones")


def test_a_not_applicable_clonal_status_adds_nothing(monkeypatch: pytest.MonkeyPatch) -> None:
    assert "clonal" not in _daniel_cell(monkeypatch, clonal_status="not_applicable")


def test_a_not_assessable_result_without_a_reason_says_so(monkeypatch: pytest.MonkeyPatch) -> None:
    cell = _daniel_cell(monkeypatch, clonal_reason="")
    assert "clonal: not assessable (no reason recorded)" in cell


def test_the_reported_table_has_its_exact_header_and_shows_counts_as_integers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _patch_results(monkeypatch, _daniel_row())
    lines = report.render().splitlines()
    assert "Reported, not gated (daniel2022):" in lines
    assert (
        "| Method | clone recurrence | day spearman | n clone strings | trb status | cells | "
        "effective clones | η at m = 3 | p at m = 3 |"
    ) in lines
    row = next(x for x in lines if x.startswith("| PAGA + DPT |") and "0.10" in x)
    assert " 40 |" in row and "40.00" not in row


def test_a_plain_dataset_gets_no_clonal_columns_and_no_table_without_reported_keys(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    plain = {k: v for k, v in _daniel_row(case="pancreas").items() if not k.startswith("reported_")}
    plain["clonal_status"] = "not_applicable"
    monkeypatch.setattr(
        report,
        "latest_real",
        lambda: (RESULTS_PATH, {"commit": "x", "methods": {"paga_dpt": {"pancreas": plain}}}),
    )
    real = report.render().split("**Real datasets**")[1]
    assert "Reported, not gated" not in real and "clonal" not in real


def test_a_plain_dataset_with_reported_keys_gets_no_clonal_columns(monkeypatch: pytest.MonkeyPatch) -> None:
    row = _daniel_row(case="pancreas", clonal_status="not_applicable")
    monkeypatch.setattr(
        report,
        "latest_real",
        lambda: (RESULTS_PATH, {"commit": "x", "methods": {"paga_dpt": {"pancreas": row}}}),
    )
    text = report.render()
    assert "Reported, not gated (pancreas):" in text
    assert "η at m = 3" not in text and "effective clones" not in text
