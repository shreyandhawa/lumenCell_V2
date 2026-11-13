"""Calibration of the clonal nulls: 200 pre-declared seeds, run once and recorded (slow).

``LUMENCELL_CALIBRATION_OUT=<path>`` writes the report before any assertion, so an honest failure is
still recorded. Definitions: ``[clonal.calibration]`` in benchmarks/thresholds.toml.
"""

import os
from pathlib import Path

import pytest
from scipy.stats import binom

from lumencell.bench.calibration import (
    CalibrationReport,
    RejectionCount,
    rejection_bound,
    run_calibration,
    write_report,
)
from lumencell.bench.gates import load_thresholds
from lumencell.engine.registry import METHODS

pytestmark = pytest.mark.slow

CFG = load_thresholds()
CAL = CFG["clonal"]["calibration"]


@pytest.fixture(scope="module")
def report(tmp_path_factory: pytest.TempPathFactory) -> CalibrationReport:
    result = run_calibration(CFG, sorted(METHODS), tmp_path_factory.mktemp("calibration"))
    out = os.environ.get("LUMENCELL_CALIBRATION_OUT")
    if out:
        write_report(result, Path(out))
    return result


def _check(report: CalibrationReport, name: str) -> RejectionCount:
    (count,) = [c for c in report.checks if c.name == name]
    print(
        f"{name}: n={count.n} n_nan={count.n_nan} rejections={count.rejections} "
        f"bound={count.bound} expect_calibrated={count.expect_calibrated}"
    )
    return count


def test_binomial_bound_is_18() -> None:
    # Reads the declared numbers, so this is not a check of SciPy alone.
    assert CAL["n_seeds"] == 200
    assert rejection_bound(CAL["n_seeds"], CAL["alpha"], CAL["binomial_quantile"]) == 18
    assert int(binom.ppf(CAL["binomial_quantile"], CAL["n_seeds"], CAL["alpha"])) == 18


@pytest.mark.parametrize("method", sorted(METHODS))
def test_the_null_is_calibrated_on_each_methods_own_masks(report: CalibrationReport, method: str) -> None:
    assert _check(report, f"null/{method}").passed


def test_donor_confound_is_calibrated_under_n1(report: CalibrationReport) -> None:
    assert _check(report, "donor_confound/N1").passed


def test_donor_confound_with_donors_ignored_is_not(report: CalibrationReport) -> None:
    assert _check(report, "donor_confound/N1_donors_ignored").passed


def test_depth_only_is_calibrated_under_n2(report: CalibrationReport) -> None:
    assert _check(report, "depth_only/N2").passed


def test_depth_only_is_not_calibrated_under_n1(report: CalibrationReport) -> None:
    assert _check(report, "depth_only/N1").passed


def test_oracle_eta_rises_with_fate_bias(report: CalibrationReport) -> None:
    for case in report.monotonic:
        print(f"{case.topology} seed {case.seed}: eta {case.eta} at fate_bias {case.fate_bias}")
    bad = [
        (c.topology, c.seed, c.eta)
        for c in report.monotonic
        if not all(a < b for a, b in zip(c.eta, c.eta[1:], strict=False))
    ]
    assert not bad, bad


def test_oracle_p_at_no_bias_is_above_0_01(report: CalibrationReport) -> None:
    for case in report.monotonic:
        print(f"{case.topology} seed {case.seed}: p at fate_bias 0 = {case.p_at_zero}")
    bad = [(c.topology, c.seed, c.p_at_zero) for c in report.monotonic if not c.p_at_zero > 0.01]
    assert not bad, bad
