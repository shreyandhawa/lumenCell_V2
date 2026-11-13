"""The bench targets must cover every method, or a re-measure silently drops rows."""

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("target", ["bench-synth", "bench-baseline"])
def test_bench_targets_cover_every_method(target: str) -> None:
    out = subprocess.run(
        ["make", "-n", target], capture_output=True, text=True, cwd=REPO_ROOT, check=True
    ).stdout
    assert "--methods paga_dpt,palantir,scfates_ppt" in out
