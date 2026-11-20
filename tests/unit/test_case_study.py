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
