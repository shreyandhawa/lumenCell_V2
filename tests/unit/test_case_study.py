"""The case-study numbers (phase 3 Task 12): each patient alone, the pooled run as a diagnostic.

Synthetic data only. A ``timepoint`` column by true-pseudotime tertile stands in for Sheih's time
point: it is a stratum, visible to methods, and it defines the root.
"""

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from lumencell.api.models import ClonalSettings
from lumencell.bench import case_study
from lumencell.bench.case_study import (
    CaseStudyDeclaration,
    load_declaration,
    render_markdown,
    run_case_study,
    run_units,
    split_units,
    write_report,
)
from lumencell.bench.report import replace_section
from lumencell.bench.synth import generate
from lumencell.bench.synth_clones import DEFAULT_LAYER, add_clones
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.stability import StabilityError
from lumencell.engine.tree_support import (
    ClonalEvidence,
    StabilityEvidence,
    calibrated_thresholds,
    verdict,
)
from lumencell.trust.clonal import settings_columns
from lumencell.trust.constants import VERDICT_VERIFIED_ON

SETTINGS = ClonalSettings(clone_key="clone_id", donor_key="donor", strata_keys=["timepoint"])
ROOT = ClusterRoot("timepoint", "t0")


def _write(path: Path, n_cells: int) -> Path:
    adata = add_clones(
        generate("bifurcation", "realistic", n_cells=n_cells, n_genes=1000, seed=0), DEFAULT_LAYER, 0
    )
    t = adata.obs["true_pseudotime"].to_numpy()
    cuts = np.quantile(t, [1 / 3, 2 / 3])
    adata.obs["timepoint"] = pd.Categorical(
        np.array(["t0", "t1", "t2"])[np.searchsorted(cuts, t, side="right")], categories=["t0", "t1", "t2"]
    )
    adata.write_h5ad(path)
    return path


def _run(path: Path, cache: Path, methods: list[str], include_stability: bool = False) -> dict[str, Any]:
    return run_case_study(
        path, SETTINGS, ROOT, methods, include_stability=include_stability, cache_root=cache
    )


@pytest.fixture(scope="module")
def data600(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return _write(tmp_path_factory.mktemp("cs600") / "d.h5ad", 600)


@pytest.fixture(scope="module")
def report(data600: Path, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    return _run(data600, tmp_path_factory.mktemp("cache600"), ["paga_dpt"])


@pytest.fixture(scope="module")
def units(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    base = tmp_path_factory.mktemp("cs900")
    path = _write(base / "d.h5ad", 900)
    decl = CaseStudyDeclaration(
        settings=SETTINGS, root=ROOT, methods=["paga_dpt"], unit_key="donor", pooled_diagnostic=True
    )
    return run_units(path, decl, base / "work", include_stability=False)


def _assert_capped_verdict(m: dict[str, Any]) -> None:  # (amended 2026-11-19, phase 2 X-1)
    ts = m["tree_support"]
    state = {c["id"]: c["state"] for c in ts["checks"]}
    # supported is capped while VERDICT_VERIFIED is False
    assert ts["level"] in {"provisional", "weak", "not_supported"}
    assert state["C2"] == state["C3"] == state["C4"] == "excluded"  # every calibrated tau is None
    if ts["level"] == "provisional":  # C1 passed and C5 did not fail: the cap is the reason
        assert state["C1"] == "pass" and state["C5"] != "fail"
        assert any("failed its verification" in r for r in ts["reasons"])
    if ts["level"] == "weak":
        assert state["C1"] == "pass" and state["C5"] == "fail"
    if ts["level"] == "not_supported":
        assert state["C1"] == "fail"


def _raise(name: str) -> None:
    raise ValueError(f"not strict JSON: {name}")


def _drop_volatile(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _drop_volatile(v) for k, v in obj.items() if k not in {"generated", "runtime_seconds"}}
    if isinstance(obj, list):
        return [_drop_volatile(v) for v in obj]
    return obj


def _table_after(text: str, heading: str) -> list[list[str]]:
    """The body rows of the first table after ``heading``, as cells."""
    lines = text[text.index(heading) :].splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("|"))
    rows = []
    for line in lines[start + 2 :]:
        if not line.startswith("|"):
            break
        rows.append([c.strip() for c in line.strip("|").split("|")])
    return rows


# ---- the plan's tests -----------------------------------------------------------------


def test_every_block_is_reported(report: dict[str, Any]) -> None:
    m = report["methods"]["paga_dpt"]
    assert [n["null"] for n in m["clonal"]["primary"]["nulls"]] == ["N1", "N2", "N3"]
    assert m["clonal"]["sensitivity_m3"]["min_clone_size"] == 3 and "eta_base" in m["clonal"]["primary"]
    assert {"coupling", "notes"} <= m["clonal"].keys() and {"tree_support", "cell_cycle"} <= m.keys()


def test_the_time_point_can_be_the_root(report: dict[str, Any]) -> None:
    assert report["methods"]["paga_dpt"]["clonal"]["status"] != "refused"


def test_each_patient_runs_alone_and_the_pooled_run_is_a_labelled_diagnostic(units: dict[str, Any]) -> None:
    assert set(units["units"]) == {"d0", "d1", "d2"}
    assert all(
        u["methods"]["paga_dpt"]["clonal"]["primary"]["n_donors"] == 1 for u in units["units"].values()
    )
    assert units["pooled_diagnostic"]["methods"]["paga_dpt"]["clonal"]["primary"]["n_donors"] == 3
    text = render_markdown(units)
    assert text.index("Per patient (primary)") < text.index("Pooled diagnostic (not the headline")


def test_unit_files_are_named_by_hash(tmp_path: Path) -> None:
    src = tmp_path / "in"
    src.mkdir()
    path = src / "d.h5ad"
    obs = pd.DataFrame(
        {"u": pd.Categorical(["../x", "a", "../x", None, "a", "../x"]), "k": list("pqrstu")},
        index=[f"c{i}" for i in range(6)],
    )
    var = pd.DataFrame(index=["g0", "g1", "g2"])
    x = sparse.csr_matrix(np.arange(18, dtype=np.float32).reshape(6, 3))
    ad.AnnData(X=x, obs=obs, var=var).write_h5ad(path)
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    work = tmp_path / "work"

    files = split_units(path, "u", work)

    assert set(files) == {"../x", "a"}
    for value, f in files.items():
        assert f == work / f"unit-{hashlib.sha256(value.encode()).hexdigest()[:12]}.h5ad"
        assert f.exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["in", "work"]
    assert sorted(p.name for p in src.iterdir()) == ["d.h5ad"]
    assert sorted(p.name for p in work.iterdir()) == sorted(f.name for f in files.values())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
    x_unit = ad.read_h5ad(files["../x"])
    assert list(x_unit.obs_names) == ["c0", "c2", "c5"]
    assert list(x_unit.var_names) == ["g0", "g1", "g2"]
    assert list(x_unit.obs["k"]) == ["p", "r", "u"]
    np.testing.assert_array_equal(x_unit.X.toarray(), x.toarray()[[0, 2, 5]])


# ---- the listed extra tests -----------------------------------------------------------


def test_the_report_is_strict_json_and_reproducible(
    report: dict[str, Any], data600: Path, tmp_path: Path
) -> None:
    again = _run(data600, tmp_path / "cache", ["paga_dpt"])
    texts = []
    for name, r in (("a.json", report), ("b.json", again)):
        write_report(r, tmp_path / name, "abc1234")
        texts.append((tmp_path / name).read_text(encoding="utf-8"))
    a, b = (json.loads(t, parse_constant=_raise) for t in texts)
    assert a["commit"] == "abc1234" and "generated" in a
    assert _drop_volatile(a) == _drop_volatile(b)
    assert not list(tmp_path.glob("*.tmp*"))


def test_write_report_replaces_the_file_atomically(report: dict[str, Any], tmp_path: Path) -> None:
    target = tmp_path / "r.json"
    target.write_text("old", encoding="utf-8")
    write_report(report, target, "abc1234")
    assert json.loads(target.read_text(encoding="utf-8"))["commit"] == "abc1234"
    bad = {**report, "x": math.nan}
    with pytest.raises(ValueError):
        write_report(bad, target, "abc1234")
    assert json.loads(target.read_text(encoding="utf-8"))["commit"] == "abc1234"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["r.json"]


def test_render_has_one_row_per_patient_and_method_and_no_nan(units: dict[str, Any]) -> None:
    text = render_markdown(units)
    rows = _table_after(text, "Per patient (primary)")
    assert sorted((r[0], r[1]) for r in rows) == [("d0", "paga_dpt"), ("d1", "paga_dpt"), ("d2", "paga_dpt")]
    pooled = _table_after(text, "Pooled diagnostic (not the headline")
    assert [(r[0], r[1]) for r in pooled] == [("pooled", "paga_dpt")]
    assert not re.search(r"\bnan\b", text, re.IGNORECASE)
    assert "None" not in text


def test_render_shows_a_missing_value_as_a_dash(units: dict[str, Any]) -> None:
    broken = json.loads(json.dumps(units))
    broken["units"]["d0"]["methods"]["paga_dpt"]["clonal"]["primary"]["eta_base"] = None
    rows = _table_after(render_markdown(broken), "Per patient (primary)")
    assert "–" in next(r for r in rows if r[0] == "d0")


def test_render_holds_only_tables_and_the_fixed_footnote(units: dict[str, Any]) -> None:
    text = render_markdown(units)
    for line in text.splitlines():
        if not line.strip() or line.startswith("|") or line.startswith("#"):
            continue
        assert line.startswith("Verdict: the rule failed its verification on"), line


def test_render_escapes_pipes_in_names(units: dict[str, Any]) -> None:
    renamed = json.loads(json.dumps(units))
    renamed["units"] = {"a|b": renamed["units"]["d0"]}
    rows = _table_after(render_markdown(renamed), "Per patient (primary)")
    assert rows[0][0] == "a\\" and rows[0][1].startswith("b")  # split on the bare pipe only
    assert "a\\|b" in render_markdown(renamed)


def test_without_a_pooled_run_there_is_no_pooled_heading(units: dict[str, Any]) -> None:
    text = render_markdown({**units, "pooled_diagnostic": None})
    assert "Pooled diagnostic" not in text and "Per patient (primary)" in text


def test_stability_gets_the_reference_params(
    monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path
) -> None:
    specs: list[Any] = []
    results: list[Any] = []
    real_run_method = case_study.run_method

    def spy_method(*args: Any, **kwargs: Any) -> Any:
        results.append(real_run_method(*args, **kwargs))
        return results[-1]

    def spy_stability(ds: Any, prep: Any, params: Any, reference: Any, spec: Any, *a: Any, **k: Any) -> Any:
        specs.append(spec)
        raise StabilityError("recorded")

    monkeypatch.setattr(case_study, "run_method", spy_method)
    monkeypatch.setattr(case_study, "run_stability", spy_stability)
    out = _run(data600, tmp_path, ["paga_dpt"], include_stability=True)
    (spec,), (result,) = specs, results
    assert spec.params == result.params
    assert spec.clonal_columns == settings_columns(SETTINGS)
    assert spec.method == "paga_dpt" and spec.terminals is None
    stab = out["methods"]["paga_dpt"]["stability"]
    assert stab["status"] == "failed" and stab["notes"] == ["recorded"]


def test_without_stability_the_block_is_none_and_the_evidence_pending(report: dict[str, Any]) -> None:
    m = report["methods"]["paga_dpt"]
    assert m["stability"] is None
    reasons = " ".join(m["tree_support"]["reasons"])
    assert "waiting for the stability run" in reasons


def test_a_nan_tree_likeness_is_written_as_null_with_a_reason(
    monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path
) -> None:
    monkeypatch.setattr(case_study, "tree_likeness", lambda *a: math.nan)
    out = _run(data600, tmp_path / "cache", ["paga_dpt"])
    write_report(out, tmp_path / "r.json", "abc1234")
    text = (tmp_path / "r.json").read_text(encoding="utf-8")
    assert '"tree_likeness": null' in text
    loaded = json.loads(text)
    assert loaded["dataset"]["tree_likeness"] is None and loaded["dataset"]["tree_likeness_note"]


def test_a_finite_tree_likeness_has_no_note(report: dict[str, Any]) -> None:
    assert 0 <= report["dataset"]["tree_likeness"] <= 1
    assert report["dataset"]["tree_likeness_note"] is None


def test_the_report_names_its_dataset_settings_root_and_versions(report: dict[str, Any]) -> None:
    d = report["dataset"]
    assert d["n_cells"] == 600 and d["n_genes"] == 1000 and d["fingerprint"] and d["name"] == "d.h5ad"
    assert report["settings"] == SETTINGS.model_dump(mode="json")
    assert report["root"] == "cluster:timepoint=t0"
    assert set(report["versions"]) == {"lumencell", "scanpy", "numpy"}
    m = report["methods"]["paga_dpt"]
    assert m["n_lineages"] >= 2 and 0 <= m["coverage"] <= 1 and m["runtime_seconds"] > 0
    assert m["cell_cycle"]["status"] in {"scored", "not_scored"}


def test_settings_that_do_not_fit_raise(data600: Path, tmp_path: Path) -> None:
    bad = ClonalSettings(clone_key="no_such_column", donor_key="donor")
    with pytest.raises(ValueError, match="no_such_column"):
        run_case_study(data600, bad, ROOT, ["paga_dpt"], include_stability=False, cache_root=tmp_path)


def test_the_input_file_is_not_modified(data600: Path, report: dict[str, Any]) -> None:
    # the fixture's run (and every other run in this module) read data600 without writing it
    before = hashlib.sha256(data600.read_bytes()).hexdigest()
    _run(data600, data600.parent / "cache-again", ["paga_dpt"])
    assert hashlib.sha256(data600.read_bytes()).hexdigest() == before


@pytest.mark.slow
def test_reseeds_use_new_method_seeds(monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path) -> None:
    results: list[Any] = []
    real_run_method = case_study.run_method

    def spy_method(*args: Any, **kwargs: Any) -> Any:
        results.append(real_run_method(*args, **kwargs))
        return results[-1]

    monkeypatch.setattr(case_study, "run_method", spy_method)
    out = _run(data600, tmp_path, ["palantir"], include_stability=True)
    (result,) = results
    seeds = [r["method_seed"] for r in out["methods"]["palantir"]["stability"]["reseeds"]["replicates"]]
    assert seeds and all(s is not None and s != result.params["seed"] for s in seeds)


def test_the_declaration_is_read_from_the_thresholds_file() -> None:
    decl = load_declaration("sheih2020")
    assert decl.unit_key == "patient" and decl.pooled_diagnostic is True
    assert decl.settings == ClonalSettings(
        clone_key="clone_id", donor_key="patient", strata_keys=["timepoint"]
    )
    assert decl.root == ClusterRoot("timepoint", "IP")
    assert decl.methods == ["paga_dpt", "palantir", "scfates_ppt"]


def test_a_declaration_without_a_field_is_refused() -> None:
    with pytest.raises(KeyError):
        load_declaration("sheih2020", {"case_study": {"sheih2020": {"clone_key": "clone_id"}}})
    with pytest.raises(KeyError):
        load_declaration("nope", {"case_study": {}})


@pytest.mark.slow
def test_stability_and_the_verdict_are_reported(data600: Path, tmp_path: Path) -> None:
    out = _run(data600, tmp_path, ["paga_dpt"], include_stability=True)
    m = out["methods"]["paga_dpt"]
    assert m["stability"] is not None and m["stability"]["status"] == "ok"
    assert m["stability"]["reseeds"]["r"] == 2 and m["stability"]["subsamples"]["r"] == 3
    _assert_capped_verdict(m)


def test_the_verdict_is_capped_and_the_report_says_so(
    report: dict[str, Any],
    units: dict[str, Any],
    data600: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _assert_capped_verdict(report["methods"]["paga_dpt"])
    text = render_markdown(units)
    assert "failed its verification" in text and VERDICT_VERIFIED_ON in text

    from lumencell.engine import tree_support
    from lumencell.trust import constants

    monkeypatch.setattr(tree_support, "VERDICT_VERIFIED", True)
    monkeypatch.setattr(constants, "VERDICT_VERIFIED", True)
    assert "failed its verification" not in render_markdown(units)
    # supported is possible once verified: the cap, not the data or the test, excludes it
    evidence = (
        1.0,
        0.9,
        StabilityEvidence("pending"),
        ClonalEvidence("not_configured", 2),
        calibrated_thresholds(),
    )
    assert verdict(*evidence).level == "supported"
    verified = _run(data600, tmp_path, ["paga_dpt"])["methods"]["paga_dpt"]["tree_support"]
    before = report["methods"]["paga_dpt"]["tree_support"]["level"]
    assert verified["level"] == ("supported" if before == "provisional" else before)


def test_replace_section_rewrites_only_between_the_markers(tmp_path: Path) -> None:
    doc = tmp_path / "doc.md"
    doc.write_text("head\n<!-- s -->\nold\n<!-- e -->\ntail\n", encoding="utf-8")
    replace_section(doc, "new\n", "<!-- s -->", "<!-- e -->")
    assert doc.read_text(encoding="utf-8") == "head\n<!-- s -->\nnew\n<!-- e -->\ntail\n"
    with pytest.raises(ValueError):
        replace_section(doc, "x\n", "<!-- case-study:start -->", "<!-- case-study:end -->")


# ---- fix round 1: the wiring, the tables and the failure path --------------------------


def _load_script() -> Any:
    import importlib.util

    path = Path(case_study.__file__).parents[3] / "benchmarks" / "case_studies" / "sheih2020.py"
    spec = importlib.util.spec_from_file_location("sheih2020_script", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _header_after(text: str, heading: str) -> list[str]:
    lines = text[text.index(heading) :].splitlines()
    line = next(line for line in lines if line.startswith("|"))
    return [c.strip() for c in line.strip("|").split("|")]


def _copy(obj: Any) -> Any:
    return json.loads(json.dumps(obj))


def test_the_context_gets_the_clonal_columns_and_hides_them_from_methods(
    monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path
) -> None:
    seen: dict[str, Any] = {"columns": [], "ctx": []}
    real_from_prep = case_study.InferenceContext.from_prep
    real_run_method = case_study.run_method

    def spy_from_prep(ds: Any, prep: Any, params: Any, clonal_columns: Any = ()) -> Any:
        seen["columns"].append(list(clonal_columns))
        return real_from_prep(ds, prep, params, clonal_columns=clonal_columns)

    def spy_method(ctx: Any, *a: Any, **k: Any) -> Any:
        seen["ctx"].append(ctx)
        return real_run_method(ctx, *a, **k)

    monkeypatch.setattr(case_study.InferenceContext, "from_prep", staticmethod(spy_from_prep))
    monkeypatch.setattr(case_study, "run_method", spy_method)
    _run(data600, tmp_path, ["paga_dpt"])
    assert seen["columns"] == [list(settings_columns(SETTINGS))] and seen["columns"][0]
    (ctx,) = seen["ctx"]
    assert "clone_id" not in ctx.obs.columns and "donor" not in ctx.obs.columns


def test_settings_only_check_settings_rejects_raise_before_any_method_runs(
    monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path
) -> None:
    def no_method(*a: Any, **k: Any) -> Any:
        raise AssertionError("a method ran")

    monkeypatch.setattr(case_study, "run_method", no_method)
    same = ClonalSettings(clone_key="clone_id", donor_key="clone_id")
    with pytest.raises(ValueError, match="must be different"):
        run_case_study(data600, same, ROOT, ["paga_dpt"], include_stability=False, cache_root=tmp_path)
    adata = ad.read_h5ad(data600)
    adata.obs["clone_id"] = pd.Categorical([None] * adata.n_obs, categories=["c"])
    empty = tmp_path / "empty.h5ad"
    adata.write_h5ad(empty)
    with pytest.raises(ValueError, match="no clone label"):
        run_case_study(empty, SETTINGS, ROOT, ["paga_dpt"], include_stability=False, cache_root=tmp_path)


def test_the_verdict_gets_the_lineage_count_and_tree_likeness_the_paga_artifacts_in_order(
    monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path
) -> None:
    calls: dict[str, Any] = {"evidence": [], "tree": [], "results": [], "preps": []}
    real_run_method = case_study.run_method
    real_evidence = case_study.clonal_evidence
    real_tree = case_study.tree_likeness
    real_cache_run = case_study.PrepCache.run

    def spy_method(*a: Any, **k: Any) -> Any:
        calls["results"].append(real_run_method(*a, **k))
        return calls["results"][-1]

    def spy_evidence(agreement: Any, n: Any) -> Any:
        calls["evidence"].append(n)
        return real_evidence(agreement, n)

    def spy_tree(*args: Any) -> Any:
        calls["tree"].append(args)
        return real_tree(*args)

    def spy_cache_run(self: Any, params: Any) -> Any:
        calls["preps"].append(real_cache_run(self, params))
        return calls["preps"][-1]

    monkeypatch.setattr(case_study, "run_method", spy_method)
    monkeypatch.setattr(case_study, "clonal_evidence", spy_evidence)
    monkeypatch.setattr(case_study, "tree_likeness", spy_tree)
    monkeypatch.setattr(case_study.PrepCache, "run", spy_cache_run)
    _run(data600, tmp_path, ["paga_dpt"])
    (result,), (prep,) = calls["results"], calls["preps"]
    assert calls["evidence"] == [len(result.lineages)] and calls["evidence"][0] >= 2
    (args,) = calls["tree"]
    assert args[0] is prep.artifacts["paga_connectivities"]
    assert args[1] is prep.artifacts["paga_tree"]


def test_a_stability_error_that_is_not_a_stability_failure_propagates(
    monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path
) -> None:
    def boom(*a: Any, **k: Any) -> Any:
        raise RuntimeError("not a StabilityError")

    monkeypatch.setattr(case_study, "run_stability", boom)
    with pytest.raises(RuntimeError, match="not a StabilityError"):
        _run(data600, tmp_path, ["paga_dpt"], include_stability=True)


@pytest.fixture(scope="module")
def failed_report(data600: Path, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """paga_dpt runs; palantir raises, and the others still report."""
    real_run_method = case_study.run_method

    def run_method_failing(ctx: Any, method: str, root: Any) -> Any:
        if method == "palantir":
            raise RuntimeError("boom /private/secret/path")
        return real_run_method(ctx, method, root)

    mp = pytest.MonkeyPatch()
    mp.setattr(case_study, "run_method", run_method_failing)
    try:
        out = _run(data600, tmp_path_factory.mktemp("cachefail"), ["paga_dpt", "palantir"])
    finally:
        mp.undo()
    return {
        "declaration": {},
        "units": {"d0": _copy(out), "d1": _copy(out)},
        "pooled_diagnostic": _copy(out),
    }


def test_a_failed_method_is_recorded_and_the_others_still_report(failed_report: dict[str, Any]) -> None:
    methods = failed_report["units"]["d0"]["methods"]
    assert methods["palantir"]["error"] == "RuntimeError: boom /private/secret/path"
    assert methods["palantir"]["clonal"] is None and methods["palantir"]["stability"] is None
    assert methods["paga_dpt"]["clonal"] is not None and "error" not in methods["paga_dpt"]
    assert case_study.failed_methods(failed_report) == [
        ("d0", "palantir", "RuntimeError"),
        ("d1", "palantir", "RuntimeError"),
        ("pooled", "palantir", "RuntimeError"),
    ]


def test_a_failed_method_is_shown_in_the_tables_without_its_message(failed_report: dict[str, Any]) -> None:
    text = render_markdown(failed_report)
    per = _table_after(text, "Per patient (primary)")
    pooled = _table_after(text, "Pooled diagnostic (not the headline")
    assert {len(r) for r in per} == {len(case_study.ROW_HEADER)}
    assert {len(r) for r in pooled} == {len(case_study.ROW_HEADER) + 1}
    ok, failed = (next(r for r in per if r[0] == "d0" and r[1] == m) for m in ("paga_dpt", "palantir"))
    assert failed[3] == "failed" and ok[3] != "failed"
    assert [c for c in failed[5:] if c != "–"] == []
    failed_rows = _table_after(text, "#### Failed methods")
    assert failed_rows == [["d0", "palantir", "RuntimeError"], ["d1", "palantir", "RuntimeError"]]
    assert text.count("#### Failed methods") == 2  # per patient and pooled, each its own
    assert "boom" not in text and "secret" not in text
    coupling = _table_after(text, "#### Lineage-pair coupling")
    assert coupling and all(r[1] == "paga_dpt" for r in coupling)
    stability = _table_after(text, "#### Stability")
    assert next(r for r in stability if r[:2] == ["d0", "palantir"])[:3] == ["d0", "palantir", "not measured"]


def test_the_script_writes_the_report_then_exits_non_zero_when_a_method_failed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, failed_report: dict[str, Any]
) -> None:
    script = _load_script()
    written: list[Path] = []
    clean = {"declaration": {}, "units": {"d0": {"methods": {"paga_dpt": {}}}}, "pooled_diagnostic": None}
    monkeypatch.setattr(script, "is_cached", lambda item: True)
    monkeypatch.setattr(script, "local_path", lambda item: tmp_path / "x.h5ad")
    monkeypatch.setattr(script, "write_report", lambda report, out, commit: written.append(out))
    monkeypatch.setattr(
        "sys.argv", ["sheih2020.py", "--out", str(tmp_path / "r.json"), "--workdir", str(tmp_path)]
    )

    monkeypatch.setattr(script, "run_units", lambda *a, **k: failed_report)
    assert script.main() == 1
    assert written == [tmp_path / "r.json"]  # written first, then the failure is reported

    monkeypatch.setattr(script, "run_units", lambda *a, **k: clean)
    assert script.main() == 0


def _stub_run_case_study(monkeypatch: pytest.MonkeyPatch) -> list[tuple[Path, Any]]:
    calls: list[tuple[Path, Any]] = []

    def stub(
        path: Path, settings: Any, root: Any, methods: Any, include_stability: bool, *, cache_root: Any = None
    ) -> Any:
        calls.append((path, cache_root))
        return {"methods": {}}

    monkeypatch.setattr(case_study, "run_case_study", stub)
    return calls


def test_the_pooled_run_happens_only_when_declared_and_caches_live_under_the_workdir(
    monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path
) -> None:
    calls = _stub_run_case_study(monkeypatch)
    decl = CaseStudyDeclaration(SETTINGS, ROOT, ["paga_dpt"], "donor", pooled_diagnostic=False)
    out = run_units(data600, decl, tmp_path / "w1")
    assert out["pooled_diagnostic"] is None and len(calls) == 3
    assert all(c != data600 for c, _ in calls)
    assert {cache for _, cache in calls} == {tmp_path / "w1" / "cache"}

    calls.clear()
    pooled = CaseStudyDeclaration(SETTINGS, ROOT, ["paga_dpt"], "donor", pooled_diagnostic=True)
    out = run_units(data600, pooled, tmp_path / "w2")
    assert out["pooled_diagnostic"] == {"methods": {}} and len(calls) == 4
    assert calls[-1] == (data600, tmp_path / "w2" / "cache")


def test_the_report_carries_the_declaration_as_written(
    monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path
) -> None:
    _stub_run_case_study(monkeypatch)
    section = {
        "declared": "2026-11-19",
        "methods": ["paga_dpt"],
        "root": "cluster:timepoint=t0",
        "clone_key": "clone_id",
        "donor_key": "donor",
        "strata_keys": ["timepoint"],
        "unit_key": "donor",
        "pooled_diagnostic": False,
        "chains": "paired",
    }
    decl = load_declaration("t", {"case_study": {"t": section}})
    assert run_units(data600, decl, tmp_path / "w")["declaration"] == section
    built = CaseStudyDeclaration(SETTINGS, ROOT, ["paga_dpt"], "donor", pooled_diagnostic=False)
    record = run_units(data600, built, tmp_path / "w3")["declaration"]
    assert record == {
        **SETTINGS.model_dump(mode="json"),
        "root": "cluster:timepoint=t0",
        "methods": ["paga_dpt"],
        "unit_key": "donor",
        "pooled_diagnostic": False,
    }


def test_the_sheih_declaration_is_a_single_donor_false_namespace() -> None:
    assert load_declaration("sheih2020").settings.single_donor is False


def test_a_failed_replace_leaves_the_old_report_and_no_temp_file(
    monkeypatch: pytest.MonkeyPatch, report: dict[str, Any], tmp_path: Path
) -> None:
    target = tmp_path / "r.json"
    target.write_text("old", encoding="utf-8")

    def broken(*a: Any, **k: Any) -> None:
        raise OSError("disk")

    monkeypatch.setattr(case_study.os, "replace", broken)
    with pytest.raises(OSError, match="disk"):
        write_report(report, target, "abc1234")
    assert target.read_text(encoding="utf-8") == "old"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["r.json"]


def test_a_failed_replace_leaves_no_temp_unit_file(
    monkeypatch: pytest.MonkeyPatch, data600: Path, tmp_path: Path
) -> None:
    def broken(*a: Any, **k: Any) -> None:
        raise OSError("disk")

    monkeypatch.setattr(case_study.os, "replace", broken)
    with pytest.raises(OSError, match="disk"):
        split_units(data600, "donor", tmp_path / "w")
    assert list((tmp_path / "w").iterdir()) == []


@pytest.mark.parametrize("only", ["<!-- s -->", "<!-- e -->"])
def test_replace_section_needs_both_markers_and_leaves_the_file_alone(tmp_path: Path, only: str) -> None:
    doc = tmp_path / "doc.md"
    doc.write_text(f"head\n{only}\nold\ntail\n", encoding="utf-8")
    with pytest.raises(ValueError):
        replace_section(doc, "new\n", "<!-- s -->", "<!-- e -->")
    assert doc.read_text(encoding="utf-8") == f"head\n{only}\nold\ntail\n"


@pytest.mark.parametrize("v", [None, math.nan, math.inf, -math.inf, np.float64("nan")])
def test_a_missing_or_non_finite_value_is_a_dash(v: Any) -> None:
    assert case_study._cell(v) == "–" and case_study._p(v) == "–"
    assert case_study._cell(0.5) == "0.500" and case_study._p(0.0123456) == "0.0123"
    assert case_study._cell(True) == "yes" and case_study._cell(7) == "7"


def test_the_pooled_table_has_a_donors_column_and_every_row_matches_its_header(units: dict[str, Any]) -> None:
    text = render_markdown(units)
    per_head = _header_after(text, "Per patient (primary)")
    pooled_head = _header_after(text, "Pooled diagnostic (not the headline")
    assert "Donors scored" not in per_head and pooled_head[8] == "Donors scored"
    assert len(pooled_head) == len(per_head) + 1
    (pooled_row,) = _table_after(text, "Pooled diagnostic (not the headline")
    assert len(pooled_row) == len(pooled_head) and pooled_row[8] == "3"
    assert {len(r) for r in _table_after(text, "Per patient (primary)")} == {len(per_head)}


def test_two_methods_render_the_whole_patient_by_method_product(report: dict[str, Any]) -> None:
    m = report["methods"]["paga_dpt"]
    both = {"methods": {"a": _copy(m), "b": _copy(m)}}
    text = render_markdown({"units": {"p1": both, "p2": _copy(both)}, "pooled_diagnostic": None})
    rows = _table_after(text, "Per patient (primary)")
    assert [(r[0], r[1]) for r in rows] == [("p1", "a"), ("p1", "b"), ("p2", "a"), ("p2", "b")]


def test_coupling_stability_and_cell_cycle_tables_show_the_given_numbers(report: dict[str, Any]) -> None:
    base = _copy(report["methods"]["paga_dpt"])
    base["clonal"]["coupling"] = [{"a": 0, "b": 2, "observed": 0.25, "z": -1.5, "p_adj": 0.5}]

    def rep(**k: Any) -> dict[str, Any]:
        return {
            "replicates": [],
            "r": 3,
            "share_matched_median": 0.5,
            "share_matched_min": 0.25,
            "spearman_median": 0.75,
            "spearman_min": 0.125,
            "mean_jaccard_median": 0.625,
            "mean_jaccard_min": 0.375,
            "n_lineages_min": 4,
            **k,
        }

    runs = {"a": base, "b": _copy(base), "c": _copy(base)}
    runs["a"]["stability"] = {
        "status": "ok",
        "reseeds": rep(replicates=[{"error": None}, {"error": "x"}, {"error": None}]),
        "subsamples": None,
    }
    runs["b"]["stability"] = None
    runs["c"]["cell_cycle"] = {
        "status": "scored",
        "dataset_share": 0.125,
        "warnings": [{"message": "one"}, {"message": "two"}],
    }
    runs["c"]["stability"] = None
    text = render_markdown({"units": {"p": {"methods": runs}}, "pooled_diagnostic": None})

    coupling = _table_after(text, "#### Lineage-pair coupling")
    assert coupling[0] == ["p", "a", "0", "2", "0.250", "-1.500", "0.5"] and len(coupling) == 3
    stability = _table_after(text, "#### Stability")
    assert stability[0] == [
        "p",
        "a",
        "reseeds",
        "3",
        "1",
        "0.500",
        "0.250",
        "0.750",
        "0.125",
        "0.625",
        "0.375",
        "4",
    ]
    assert stability[1] == ["p", "a", "subsamples: ok"] + ["–"] * 9
    assert stability[2] == ["p", "b", "not measured"] + ["–"] * 9
    cycle = _table_after(text, "#### Cell cycle")
    assert cycle[2] == ["p", "c", "scored", "0.125", "one; two"]


def test_the_verdict_footnote_comes_once_after_the_patients_and_before_the_pooled_run(
    units: dict[str, Any],
) -> None:
    text = render_markdown(units)
    note = "Verdict: the rule failed its verification on"
    assert text.count(note) == 1
    assert text.index("Per patient (primary)") < text.index(note) < text.index("Pooled diagnostic")
    assert text.index(note) < text.index("#### Lineage-pair coupling")
