import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any, get_args

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

import lumencell.bench.real_clonal as real_clonal_mod
from lumencell.bench import report, runner
from lumencell.bench.datasets import DATASETS, RealClonal, RealDataset
from lumencell.bench.download import Source, sha256_of
from lumencell.bench.gates import REPO, check_gates, load_thresholds
from lumencell.bench.real_clonal import adjusted_rand_index, clone_recurrence
from lumencell.bench.runner import ClonalRowStatus, run_real_case
from lumencell.bench.synth import generate
from lumencell.bench.synth_clones import DEFAULT_LAYER, add_clones
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


def test_a_dataset_with_no_gated_lineage_is_shown_as_not_measured(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:  # controller ruling P14
    note = "lineage gates: not measured (no gated lineage, declared 2026-11-20)"
    ds = DATASETS["daniel2022"]
    assert ds.gated == [] and ds.clonal is not None and ds.clonal.declared == "2026-11-20"
    row = _daniel_row()
    _patch_results(monkeypatch, row)
    cell = next(x for x in report.render().splitlines() if x.startswith("| PAGA + DPT |") and "Spearman" in x)
    assert note in cell and "clonal: not assessable (declared 2026-11-20, because)" in cell
    assert "PASS" not in cell and "pass" not in cell

    # The CLI line: the dataset runs, its gates score NaN and are skipped, and nothing says PASS.
    monkeypatch.setattr("lumencell.bench.datasets.is_cached", lambda d: True)
    monkeypatch.setattr(
        "lumencell.bench.runner.run_real_case",
        lambda *a, **k: type(
            "C", (), {"row": lambda self: dict(row, runtime_seconds=1.0, tree_support="not_measured")}
        )(),
    )
    out = tmp_path / "r.json"
    res = cli.invoke(
        app, ["bench", "real", "--datasets", "daniel2022", "--methods", "paga_dpt", "--out", str(out)]
    )
    line = next(x for x in res.output.splitlines() if x.startswith("daniel2022"))
    assert res.exit_code == 0 and note in line and "clonal=not_assessable" in line
    assert "PASS" not in line and "pass" not in line

    # No clonal declaration: no date to cite, said plainly. Gated datasets render as before.
    monkeypatch.setitem(DATASETS, "daniel2022", replace(ds, clonal=None))
    cell = next(x for x in report.render().splitlines() if x.startswith("| PAGA + DPT |") and "Spearman" in x)
    assert "lineage gates: not measured (no gated lineage)" in cell and "pass" not in cell
    pancreas = DATASETS["pancreas"]
    assert pancreas.gated and report.lineage_gates_note(pancreas) is None
