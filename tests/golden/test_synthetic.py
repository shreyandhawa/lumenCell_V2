"""Golden accuracy benchmark on synthetic data with known trajectories.

Fast (runs in make ci): the bifurcation case, which has the short-trunk,
long-branch shape that fooled V1's automatic root, must stay within
``regression.max_drop`` of benchmarks/baseline.json.

Slow (make py-test-all / make bench-synth): the full matrix against the
baseline, and a check that each method's ``experimental`` flag matches whether
it misses any gate.
"""

from pathlib import Path

import pytest

from lumencell.bench.gates import check_gates, check_regression, load_baseline, load_thresholds
from lumencell.bench.runner import run_synthetic_case
from lumencell.engine.registry import METHODS
from lumencell.prep.stages import PrepParams

CFG = load_thresholds()
SYN = CFG["synthetic"]
BASELINE = load_baseline()

pytestmark = pytest.mark.golden


def _run(tmp_path: Path, method: str, topology: str, seed: int, **kw: object) -> dict[str, object]:
    res = run_synthetic_case(
        tmp_path, topology, method=method, preset=SYN["preset"], seed=seed,
        n_cells=SYN["n_cells"], n_genes=SYN["n_genes"], **kw,  # type: ignore[arg-type]
    )  # fmt: skip
    return res.row()


def test_v1_trap_bifurcation_holds_its_baseline(tmp_path: Path) -> None:
    row = _run(tmp_path, "paga_dpt", "bifurcation", 0)
    base = BASELINE["paga_dpt"][str(row["case"])]
    assert not check_regression(row, base, CFG["regression"]["max_drop"])
    # The V1 lesson, stated directly: with an explicit root, pseudotime on the
    # bifurcation follows the truth (V1's default view scored -0.21 to -0.33).
    assert float(row["pt_spearman_global"]) >= 0.90


@pytest.mark.slow
@pytest.mark.parametrize("method", sorted(METHODS))
def test_full_matrix_holds_baseline_and_experimental_flag_is_honest(tmp_path: Path, method: str) -> None:
    failures, regressions = [], []
    for topology in SYN["topologies"]:
        for seed in SYN["seeds"]:
            row = _run(tmp_path, method, topology, seed)
            failures += check_gates(row, SYN["gates"])
            base = BASELINE.get(method, {}).get(str(row["case"]))
            assert base is not None, f"no baseline for {method} {row['case']}: run make bench-baseline"
            regressions += check_regression(row, base, CFG["regression"]["max_drop"])
    assert not regressions, "\n".join(map(str, regressions))
    assert METHODS[method].experimental == bool(failures), (
        f"{method}.experimental must be {bool(failures)}; gate failures:\n" + "\n".join(map(str, failures))
    )


@pytest.mark.slow
@pytest.mark.parametrize("resolution", [0.5, 2.0])
def test_bifurcation_is_found_across_leiden_resolutions(tmp_path: Path, resolution: float) -> None:
    row = _run(tmp_path, "paga_dpt", "bifurcation", 0, prep_params=PrepParams(leiden_resolution=resolution))
    assert row["lineage_recall"] == 1.0
    assert float(row["pt_spearman_lineage"]) >= 0.85
