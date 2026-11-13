"""The counting rules of the calibration report, on made-up p-values (fast; the run itself is slow)."""

import json
import math
from pathlib import Path

import pytest
from scipy.stats import binom

from lumencell.bench.calibration import (
    CalibrationReport,
    MonotonicCase,
    RejectionCount,
    count_rejections,
    rejection_bound,
    write_report,
)
from lumencell.bench.gates import REPO, load_thresholds

NAN = math.nan


def test_a_nan_p_is_counted_apart_and_the_bound_uses_the_scored_seeds() -> None:
    count = count_rejections("null/x", "headline", True, [0.01, 0.2, NAN, 0.05, 0.9, NAN], 0.05, 0.99)
    assert (count.n, count.n_nan, count.n_scored, count.rejections) == (6, 2, 4, 2)  # 0.05 rejects
    assert count.bound == int(binom.ppf(0.99, 4, 0.05))
    assert count.p_values == (0.01, 0.2, None, 0.05, 0.9, None)


def test_a_nan_p_is_neither_a_rejection_nor_a_pass() -> None:
    all_nan = count_rejections("null/x", "headline", True, [NAN] * 5, 0.05, 0.99)
    assert (all_nan.n_scored, all_nan.rejections, all_nan.bound) == (0, 0, 0)
    # Five NaN seeds out of 200 lower the bound only through n_scored: 195 scored seeds, not 200.
    ps = [0.5] * 195 + [NAN] * 5
    assert count_rejections("c", "N1", True, ps, 0.05, 0.99).bound == int(binom.ppf(0.99, 195, 0.05))


@pytest.mark.parametrize(
    ("expect_calibrated", "rejections", "passed"),
    [(True, 18, True), (True, 19, False), (False, 18, False), (False, 19, True)],
)
def test_passed_compares_rejections_with_the_bound(
    expect_calibrated: bool, rejections: int, passed: bool
) -> None:
    count = RejectionCount("c", "N1", 200, 0, rejections, 18, expect_calibrated, ())
    assert count.passed is passed


def test_the_report_is_strict_json_with_null_for_nan(tmp_path: Path) -> None:
    check = count_rejections("null/x", "headline", True, [0.5, NAN], 0.05, 0.99)
    mono = MonotonicCase("bifurcation", 0, (0.0, 0.7), (NAN, 0.5), NAN)
    assert not mono.passed
    report = CalibrationReport((check,), (mono,), {"methods": ["x"]}, "abc")
    path = tmp_path / "out" / "r.json"
    write_report(report, path)

    def refuse(token: str) -> None:
        raise AssertionError(f"non-strict JSON constant {token}")

    doc = json.loads(path.read_text(), parse_constant=refuse)
    assert doc["passed"] is False
    assert doc["checks"][0]["p_values"] == [0.5, None]
    assert doc["checks"][0]["n_nan"] == 1
    assert doc["monotonic"][0]["eta"] == [None, 0.5]
    assert doc["monotonic"][0]["p_at_zero"] is None


def test_monotonic_needs_strictly_rising_eta_and_p_above_one_percent() -> None:
    assert MonotonicCase("t", 0, (0.0, 0.5), (0.0, 0.3), 0.2).passed
    assert not MonotonicCase("t", 0, (0.0, 0.5), (0.3, 0.3), 0.2).passed
    assert not MonotonicCase("t", 0, (0.0, 0.5), (0.0, 0.3), 0.01).passed


def _mono(eta: tuple[float, ...], p0: float, bias: tuple[float, ...] | None = None) -> MonotonicCase:
    return MonotonicCase("t", 0, bias or tuple(range(len(eta))), eta, p0)  # type: ignore[arg-type]


def test_monotonic_rule_in_each_part() -> None:
    assert _mono((0.0, 0.3, 0.6, 1.0), 0.2).passed
    assert not _mono((0.0, NAN, 0.5), 0.2).passed  # NaN eta alone, p fine
    assert not _mono((0.0, 0.3, NAN), 0.2).passed  # NaN as the last eta
    assert not _mono((0.0, 0.3, 0.2), 0.2).passed  # a fall after the first pair
    assert not _mono((0.0, 0.3), NAN).passed  # NaN p at zero
    assert _mono((0.0, 0.3), 0.0100001).passed
    assert not _mono((0.0, 0.3), 0.01).passed


def _report(
    *counts: RejectionCount, mono_ok: bool = True, errors: list[str] | None = None
) -> CalibrationReport:
    mono = (_mono((0.0, 0.5), 0.2 if mono_ok else 0.001),)
    return CalibrationReport(counts, mono, {"errors": errors or []}, "abc")


GOOD = RejectionCount("a", "N1", 200, 0, 18, 18, True, ())
BAD = RejectionCount("b", "N1", 200, 0, 19, 18, True, ())
TEETH = RejectionCount("c", "N1", 200, 0, 19, 18, False, ())


@pytest.mark.parametrize(
    ("report", "passed"),
    [
        (_report(GOOD, TEETH), True),
        (_report(BAD), False),
        (_report(GOOD, BAD), False),
        (_report(BAD, GOOD), False),
        (_report(GOOD, mono_ok=False), False),
        (_report(GOOD, errors=["seed 100 paga_dpt: boom"]), False),
    ],
)
def test_report_passes_only_if_every_part_passes(
    report: CalibrationReport, passed: bool, tmp_path: Path
) -> None:
    assert report.passed is passed
    write_report(report, tmp_path / "r.json")
    assert json.loads((tmp_path / "r.json").read_text())["passed"] is passed


@pytest.mark.parametrize("expect_calibrated", [True, False])
def test_a_check_with_no_scored_seed_never_passes(expect_calibrated: bool) -> None:
    count = count_rejections("c", "N1", expect_calibrated, [NAN] * 200, 0.05, 0.99)
    assert (count.n_scored, count.rejections, count.bound) == (0, 0, 0)
    assert not count.passed


def test_the_json_carries_every_field(tmp_path: Path) -> None:
    check = count_rejections("null/x", "headline", True, [0.5, NAN], 0.05, 0.99)
    mono = MonotonicCase("bifurcation", 0, (0.0, 0.7), (0.1, 0.5), 0.2)
    write_report(CalibrationReport((check,), (mono,), {}, "abc"), tmp_path / "r.json")
    doc = json.loads((tmp_path / "r.json").read_text())
    assert doc["commit"] == "abc"
    row = doc["checks"][0]
    assert (row["name"], row["null"], row["n"], row["n_scored"], row["bound"]) == (
        "null/x",
        "headline",
        2,
        1,
        1,
    )
    assert row["passed"] is True and row["expect_calibrated"] is True
    assert doc["monotonic"][0]["passed"] is True
    assert doc["monotonic"][0]["topology"] == "bifurcation"


def test_count_rejections_keeps_labels_and_expectation() -> None:
    c = count_rejections("depth_only/N1", "N1", False, [0.01] * 30 + [0.9] * 170, 0.05, 0.99)
    assert (c.name, c.null, c.expect_calibrated) == ("depth_only/N1", "N1", False)
    assert c.rejections == 30 and c.passed


def test_alpha_and_quantile_are_used_as_given() -> None:
    c = count_rejections("c", "N1", True, [0.07] * 10 + [0.5] * 90, 0.10, 0.95)
    assert (c.rejections, c.bound) == (10, int(binom.ppf(0.95, 100, 0.10)))
    assert c.bound == 15
    edge = count_rejections("c", "N1", True, [0.10, 0.1000001], 0.10, 0.95)
    assert edge.rejections == 1


def test_the_recorded_run_matches_the_declared_configuration() -> None:
    cfg = load_thresholds()["clonal"]["calibration"]
    doc = json.loads((REPO / "benchmarks" / "results" / "clonal-calibration-2026-11-13.json").read_text())
    assert doc["config"]["calibration"] == cfg
    assert doc["config"]["errors"] == []
    teeth = {"donor_confound/N1_donors_ignored", "depth_only/N1"}
    names = {c["name"] for c in doc["checks"]}
    assert names == teeth | {"donor_confound/N1", "depth_only/N2"} | {
        f"null/{m}" for m in doc["config"]["methods"]
    }
    for c in doc["checks"]:
        scored = [p for p in c["p_values"] if p is not None]
        assert c["n"] == cfg["n_seeds"] == len(c["p_values"]) == 200
        assert c["n_nan"] == c["n"] - len(scored)
        assert c["bound"] == int(binom.ppf(cfg["binomial_quantile"], len(scored), cfg["alpha"])) == 18
        assert c["rejections"] == sum(p <= cfg["alpha"] for p in scored)
        assert c["expect_calibrated"] is (c["name"] not in teeth)
    # Recorded outcome: depth_only/N2 is the one failed check.
    assert {c["name"] for c in doc["checks"] if not c["passed"]} == {"depth_only/N2"}
    assert doc["passed"] is False


def test_the_declared_bound_is_18() -> None:
    cfg = load_thresholds()["clonal"]["calibration"]
    assert rejection_bound(cfg["n_seeds"], cfg["alpha"], cfg["binomial_quantile"]) == 18
