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
    write_report,
)

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
