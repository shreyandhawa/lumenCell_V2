"""Calibration of the N2a amendment: 200 fresh pre-declared seeds, run once and recorded (slow).

``LUMENCELL_CALIBRATION_N2_OUT=<path>`` writes the report before any assertion, so an honest failure
is still recorded. Definitions: ``[clonal]`` (the amendment and its adoption rule) and
``[clonal.calibration_n2_amended]`` in benchmarks/thresholds.toml. Every check uses the amended N2.
"""

import os
from pathlib import Path

import pytest

from lumencell.bench.calibration import CalibrationReport, RejectionCount, run_calibration, write_report
from lumencell.bench.gates import load_thresholds
from lumencell.engine.registry import METHODS

pytestmark = pytest.mark.slow

CFG = load_thresholds()
TABLE = "calibration_n2_amended"


@pytest.fixture(scope="module")
def report(tmp_path_factory: pytest.TempPathFactory) -> CalibrationReport:
    result = run_calibration(
        CFG, sorted(METHODS), tmp_path_factory.mktemp("calibration-n2"), table=TABLE, n2_split="lineage_onset"
    )
    out = os.environ.get("LUMENCELL_CALIBRATION_N2_OUT")
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


def test_the_run_is_the_declared_one(report: CalibrationReport) -> None:
    assert report.config["calibration"] == CFG["clonal"][TABLE]
    assert (report.config["table"], report.config["n2_split"]) == (TABLE, "lineage_onset")
    assert report.config["errors"] == []


@pytest.mark.parametrize("method", sorted(METHODS))
def test_the_null_is_calibrated_on_each_methods_own_masks(report: CalibrationReport, method: str) -> None:
    assert _check(report, f"null/{method}").passed


def test_donor_confound_is_calibrated_under_n1(report: CalibrationReport) -> None:
    assert _check(report, "donor_confound/N1").passed


def test_donor_confound_with_donors_ignored_is_not(report: CalibrationReport) -> None:
    assert _check(report, "donor_confound/N1_donors_ignored").passed


def test_depth_only_is_calibrated_under_the_amended_n2(report: CalibrationReport) -> None:
    assert _check(report, "depth_only/N2").passed


def test_depth_only_is_not_calibrated_under_n1(report: CalibrationReport) -> None:
    assert _check(report, "depth_only/N1").passed


def test_oracle_eta_rises_and_p_at_no_bias_is_above_0_01(report: CalibrationReport) -> None:
    for case in report.monotonic:
        print(f"{case.topology} seed {case.seed}: eta {case.eta}, p at fate_bias 0 = {case.p_at_zero}")
    assert all(case.passed for case in report.monotonic)


def test_the_whole_report_passes(report: CalibrationReport) -> None:
    assert report.passed
