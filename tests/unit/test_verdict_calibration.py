"""The tau calibration rule ([tree_support.calibration]) and its recorded run.

The rule is pure and tested on hand-made records. ``test_calibration_on_the_declared_seeds`` is the
slow, run-once measurement: ``LUMENCELL_VERDICT_CALIBRATION_OUT=<path>`` writes its payload before
any assertion, so an honest result is always recorded.
"""

import json
import math
import os
from collections.abc import Sequence
from pathlib import Path

import pytest

from lumencell.bench.gates import load_thresholds
from lumencell.bench.verdict_calibration import (
    Calibration,
    CheckCalibration,
    Group,
    RunRecord,
    Verification,
    calibration_payload,
    collect_records,
    derive,
    expected_runs,
    expected_verification_runs,
    records_from_payload,
    verification_payload,
    verify,
    write_json,
)
from lumencell.engine.tree_support import Thresholds, calibrated_thresholds

CFG = load_thresholds()
CAL = CFG["tree_support"]["calibration"]
FACTORS = tuple(CAL["c2_resolution_factors"])


def _rec(
    group: Group,
    dataset: str,
    method: str = "paga_dpt",
    *,
    lin: float = 1.0,
    pt: float = 1.0,
    tree: float | Sequence[float] = 1.0,
    error: str | None = None,
) -> RunRecord:
    trees = [tree] * len(FACTORS) if isinstance(tree, float) else list(tree)
    return RunRecord(
        group=group,
        dataset=dataset,
        seed=0,
        method=method,
        coverage=1.0,
        tree_likeness=tuple(zip(FACTORS, trees, strict=True)),
        lineage_stability=lin,
        pseudotime_stability=pt,
        n_lineages=2,
        error=error,
    )


def _check(cal: Calibration, check: str, factor: float | None = None) -> CheckCalibration:
    (c,) = [c for c in cal.checks if c.check == check and c.factor == factor]
    return c


def _spread(matrix: Sequence[float], control: Sequence[float], field: str = "lin") -> list[RunRecord]:
    out = [_rec("matrix", f"m/seed{i}", **{field: v}) for i, v in enumerate(matrix)]  # type: ignore[arg-type]
    out += [_rec("control", f"c/seed{i}", **{field: v}) for i, v in enumerate(control)]  # type: ignore[arg-type]
    return out


def test_tau_is_the_midpoint_of_the_declared_quantiles() -> None:
    cal = derive(_spread([0.6, 0.7, 0.8, 0.9, 1.0], [0.1, 0.2, 0.3, 0.4, 0.5]), CAL)
    c3 = _check(cal, "C3")
    assert c3.matrix_quantile == pytest.approx(0.64)
    assert c3.control_quantile == pytest.approx(0.46)
    assert c3.tau == pytest.approx(0.55)
    assert c3.discriminating
    assert (c3.n_matrix, c3.n_control) == (5, 5)
    assert cal.tau_lin == c3.tau
    assert "C3" not in cal.non_discriminating


@pytest.mark.parametrize(
    ("matrix", "control"),
    [
        ([0.1, 0.2, 0.3, 0.4, 0.5], [0.6, 0.7, 0.8, 0.9, 1.0]),  # overlap the wrong way
        ([0.5, 0.5, 0.5], [0.5, 0.5, 0.5]),  # touch: equal quantiles
    ],
)
def test_ranges_that_overlap_or_touch_are_non_discriminating(
    matrix: list[float], control: list[float]
) -> None:
    cal = derive(_spread(matrix, control, "pt"), CAL)
    c4 = _check(cal, "C4")
    assert c4.discriminating is False
    assert c4.tau is None
    assert cal.tau_pt is None
    assert "C4" in cal.non_discriminating
    assert "C3" in cal.non_discriminating  # lin is 1.0 everywhere, so it touches too


def _c2_records(control_at: dict[float, float]) -> list[RunRecord]:
    """Matrix tree-likeness 0.9 at every factor; controls at ``control_at[factor]``."""
    recs = [_rec("matrix", f"m/seed{i}", tree=0.9, lin=0.9) for i in range(3)]
    recs += [_rec("control", f"c/seed{i}", tree=[control_at[f] for f in FACTORS], lin=0.1) for i in range(3)]
    return recs


def test_c2_is_kept_only_if_it_discriminates_at_all_three_resolutions() -> None:
    overlapping = derive(_c2_records({0.5: 0.5, 1.0: 0.5, 2.0: 0.95}), CAL)
    assert _check(overlapping, "C2", 0.5).discriminating
    assert _check(overlapping, "C2", 1.0).discriminating
    assert not _check(overlapping, "C2", 2.0).discriminating
    assert overlapping.tau_tree is None
    assert "C2" in overlapping.non_discriminating

    kept = derive(_c2_records({0.5: 0.3, 1.0: 0.5, 2.0: 0.7}), CAL)
    assert kept.tau_tree == _check(kept, "C2", 1.0).tau == pytest.approx(0.7)
    assert "C2" not in kept.non_discriminating
    assert kept.non_discriminating == ("C4",)  # pt is 1.0 everywhere


def test_c2_counts_each_dataset_once() -> None:
    recs = [
        _rec("matrix", f"m/seed{i}", method, tree=0.1 * (i + 1))
        for i in range(4)
        for method in ("paga_dpt", "palantir", "scfates_ppt")
    ]
    recs += [_rec("control", "c/seed0", method) for method in ("paga_dpt", "palantir")]
    cal = derive(recs, CAL)
    c2 = _check(cal, "C2", 1.0)
    assert (c2.n_matrix, c2.n_control) == (4, 1)
    assert c2.matrix_quantile == pytest.approx(0.13)  # quantile of 0.1, 0.2, 0.3, 0.4
    assert _check(cal, "C3").n_matrix == 12


def test_failed_runs_are_left_out_and_counted() -> None:
    recs = _spread([0.6, 0.7, 0.8, 0.9, 1.0], [0.1, 0.2, 0.3, 0.4, 0.5])
    nan = math.nan
    recs += [
        _rec("matrix", "m/seed9", lin=nan, pt=nan, tree=nan, error="RuntimeError: boom"),
        _rec("control", "c/seed9", lin=0.99, pt=0.99, tree=0.99, error="ValueError: no root"),
    ]
    cal = derive(recs, CAL)
    c3 = _check(cal, "C3")
    assert (c3.n_matrix, c3.n_control) == (5, 5)
    assert c3.tau == pytest.approx(0.55)
    assert _check(cal, "C2", 1.0).n_control == 5
    assert cal.n_failed == 2


def test_an_undefined_tree_likeness_counts_as_zero() -> None:
    recs = [_rec("matrix", f"m/seed{i}", tree=0.8) for i in range(2)]
    recs += [_rec("control", "c/seed0", tree=math.nan), _rec("control", "c/seed1", tree=0.4)]
    c2 = _check(derive(recs, CAL), "C2", 1.0)
    assert c2.n_control == 2
    assert c2.control_quantile == pytest.approx(0.36)  # 0.9 quantile of 0.0 and 0.4


def test_an_empty_group_is_non_discriminating() -> None:
    cal = derive([_rec("matrix", "m/seed0")], CAL)
    c3 = _check(cal, "C3")
    assert c3.n_control == 0
    assert math.isnan(c3.control_quantile)
    assert not c3.discriminating
    assert cal.non_discriminating == ("C2", "C3", "C4")


def test_the_payload_is_strict_json(tmp_path: Path) -> None:
    recs = _spread([0.6, 0.7], [0.1, 0.2])
    recs.append(_rec("matrix", "m/seed9", lin=math.nan, pt=math.nan, tree=math.nan, error="E: x"))
    cal = derive(recs, CAL)
    payload = calibration_payload(recs, cal, CFG, "abc123")
    path = tmp_path / "out.json"
    write_json(path, payload)
    doc = json.loads(path.read_text())  # strict: NaN would not be valid here
    assert "NaN" not in path.read_text()
    assert set(doc) >= {
        "generated",
        "commit",
        "config",
        "records",
        "checks",
        "taus",
        "non_discriminating",
        "n_failed",
        "complete",
    }
    assert doc["commit"] == "abc123"
    assert doc["n_failed"] == 1
    assert doc["complete"] is False
    assert doc["records"][-1]["lineage_stability"] is None
    assert doc["records"][-1]["tree_likeness"] == [[f, None] for f in FACTORS]
    assert doc["config"]["calibration"] == CAL
    assert set(doc["taus"]) == {"tau_tree", "tau_lin", "tau_pt"}


def test_records_read_back_from_a_written_payload_derive_the_same(tmp_path: Path) -> None:
    recs = _spread([0.6, 0.7, 0.9], [0.1, 0.2, 0.3])
    recs.append(_rec("control", "c/seed8", tree=math.nan))
    recs.append(_rec("matrix", "m/seed9", lin=math.nan, pt=math.nan, tree=math.nan, error="E: x"))
    path = tmp_path / "out.json"
    write_json(path, calibration_payload(recs, derive(recs, CAL), CFG, "abc"))
    back = records_from_payload(json.loads(path.read_text()))
    assert derive(back, CAL) == derive(recs, CAL)
    assert back[:-2] == recs[:-2]


def test_write_json_refuses_infinity(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        write_json(tmp_path / "x.json", {"v": math.inf})


def test_expected_runs_are_every_dataset_and_method() -> None:
    runs = expected_runs(CFG)
    topologies, kinds = CFG["synthetic"]["topologies"], CFG["synthetic"]["controls"]["kinds"]
    n_datasets = len(topologies) * len(CAL["matrix_seeds"]) + len(kinds) * len(CAL["control_seeds"])
    assert len(runs) == n_datasets * len(CAL["methods"]) == 72
    assert ("matrix", "bifurcation/seed0", "palantir") in runs
    assert ("control", "loop/seed2", "scfates_ppt") in runs


def test_a_payload_with_every_run_is_complete() -> None:
    recs = [_rec(g, d, m) for g, d, m in sorted(expected_runs(CFG))]
    payload = calibration_payload(recs, derive(recs, CAL), CFG, "abc")
    assert payload["complete"] is True


@pytest.mark.slow
def test_calibration_on_the_declared_seeds(tmp_path_factory: pytest.TempPathFactory) -> None:
    from lumencell.bench.calibration import _commit

    seeds: dict[Group, list[int]] = {"matrix": CAL["matrix_seeds"], "control": CAL["control_seeds"]}
    records = collect_records(CFG, seeds, CAL["methods"], tmp_path_factory.mktemp("verdict-cal"))
    cal = derive(records, CAL)
    payload = calibration_payload(records, cal, CFG, _commit())
    out = os.environ.get("LUMENCELL_VERDICT_CALIBRATION_OUT")
    if out:
        write_json(Path(out), payload)
    for c in cal.checks:
        print(c)
    print(f"taus: {payload['taus']}; non_discriminating: {list(cal.non_discriminating)}")
    print(f"n_failed: {cal.n_failed}")
    assert payload["complete"]
    assert len(records) == 72


REPO = Path(__file__).resolve().parents[2]


def test_the_committed_taus_follow_from_the_recorded_calibration() -> None:
    newest = sorted((REPO / "benchmarks" / "results").glob("verdict-calibration-*.json"))[-1]
    doc = json.loads(newest.read_text())
    assert doc["complete"] is True
    assert not doc["commit"].endswith("-dirty")
    cal = derive(records_from_payload(doc), CAL)
    expected: dict[str, object] = {"non_discriminating": list(cal.non_discriminating)}
    for key, tau in (("tau_tree", cal.tau_tree), ("tau_lin", cal.tau_lin), ("tau_pt", cal.tau_pt)):
        if tau is not None:
            expected[key] = tau
    assert CFG["tree_support"]["calibrated"] == expected
    assert doc["taus"] == {"tau_tree": cal.tau_tree, "tau_lin": cal.tau_lin, "tau_pt": cal.tau_pt}


# ---- Task 9: verification -------------------------------------------------------------------

VER = CFG["tree_support"]["verification"]
LOOSE = Thresholds(tau_tree=None, tau_lin=None, tau_pt=None)  # only C1 decides


def _cov(rec: RunRecord, coverage: float) -> RunRecord:
    return RunRecord(**{**rec.__dict__, "coverage": coverage})


def test_a_supported_control_fails_verification() -> None:
    recs = [_rec("matrix", "m/seed3"), _rec("control", "loop/seed3")]
    ver = verify(recs, LOOSE, 0)
    assert ver.controls_supported == ("loop/seed3 paga_dpt",)
    assert ver.passed is False
    assert verify(recs, LOOSE, 1).passed is True


def test_a_control_that_is_not_supported_passes() -> None:
    recs = [_rec("matrix", "m/seed3"), _cov(_rec("control", "two_groups/seed3"), 0.5)]
    ver = verify(recs, LOOSE, 0)
    assert ver.controls_supported == ()
    assert ver.passed is True
    assert ("two_groups/seed3", "paga_dpt", "not_supported") in ver.verdicts


def test_the_matrix_share_is_reported_not_gated() -> None:
    # Every matrix record is weak (C3 fails against tau 0.9), no control is supported.
    tau = Thresholds(tau_tree=None, tau_lin=0.9, tau_pt=None)
    recs = [_rec("matrix", f"m/seed{i}", lin=0.1) for i in range(3)] + [
        _cov(_rec("control", "loop/seed3"), 0.5)
    ]
    ver = verify(recs, tau, 0)
    assert {lvl for d, _, lvl in ver.verdicts if d.startswith("m/")} == {"weak"}
    assert ver.matrix_supported_share == 0.0
    assert ver.passed is True


def test_the_matrix_share_counts_supported_matrix_runs() -> None:
    recs = [_rec("matrix", "m/seed3"), _cov(_rec("matrix", "m/seed4"), 0.5)]
    assert verify(recs, LOOSE, 0).matrix_supported_share == 0.5


def test_failed_runs_get_no_verdict_and_are_counted() -> None:
    recs = [_rec("matrix", "m/seed3"), _rec("control", "loop/seed3", error="boom")]
    ver = verify(recs, LOOSE, 0)
    assert ver.n_failed == 1
    assert [d for d, _, _ in ver.verdicts] == ["m/seed3"]
    assert ver.passed is True


def test_the_verification_payload_is_strict_json_and_has_passed(tmp_path: Path) -> None:
    recs = [_rec("matrix", "m/seed3"), _rec("control", "loop/seed3", error="x", tree=math.nan)]
    ver = verify(recs, LOOSE, 0)
    payload = verification_payload(recs, ver, LOOSE, CFG, "abc")
    assert payload["passed"] is True
    assert payload["complete"] is False
    assert payload["level_counts"]["matrix"] == {"supported": 1}
    out = tmp_path / "v.json"
    write_json(out, payload)
    doc = json.loads(out.read_text(), parse_constant=lambda c: pytest.fail(c))
    assert doc["controls_supported"] == []
    again = verify(records_from_payload(doc), LOOSE, doc["max_controls_supported"])
    assert again == ver


def test_expected_verification_runs_use_the_fresh_seeds() -> None:
    runs = expected_verification_runs(CFG)
    assert {int(d.rsplit("seed", 1)[1]) for _, d, _ in runs} == set(VER["matrix_seeds"]) | set(
        VER["control_seeds"]
    )
    assert not {0, 1, 2} & {int(d.rsplit("seed", 1)[1]) for _, d, _ in runs}


@pytest.mark.slow
def test_verification_on_fresh_seeds(tmp_path_factory: pytest.TempPathFactory) -> None:
    from lumencell.bench.calibration import _commit

    seeds: dict[Group, list[int]] = {"matrix": VER["matrix_seeds"], "control": VER["control_seeds"]}
    records = collect_records(CFG, seeds, CAL["methods"], tmp_path_factory.mktemp("verdict-ver"))
    thresholds = calibrated_thresholds()
    ver: Verification = verify(records, thresholds, VER["max_controls_supported"])
    payload = verification_payload(records, ver, thresholds, CFG, _commit())
    out = os.environ.get("LUMENCELL_VERDICT_VERIFICATION_OUT")
    if out:
        write_json(Path(out), payload)
    print(f"level counts: {payload['level_counts']}")
    print(f"controls supported: {payload['controls_supported_by_kind']} of {len(ver.controls_supported)}")
    print(f"matrix supported share: {ver.matrix_supported_share:.3f}; n_failed: {ver.n_failed}")
    assert payload["complete"]
    assert ver.passed, f"controls verdicted supported: {ver.controls_supported}"


def test_the_recorded_verification_follows_from_its_records() -> None:
    found = sorted((REPO / "benchmarks" / "results").glob("verdict-verification-*.json"))
    if not found:
        pytest.skip("verification not recorded yet")
    doc = json.loads(found[-1].read_text())
    assert doc["complete"] is True
    assert not doc["commit"].endswith("-dirty")
    ver = verify(records_from_payload(doc), calibrated_thresholds(), VER["max_controls_supported"])
    assert doc["passed"] == ver.passed
    assert doc["controls_supported"] == list(ver.controls_supported)
