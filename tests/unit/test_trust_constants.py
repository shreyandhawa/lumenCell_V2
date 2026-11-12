from functools import reduce
from typing import Any

from lumencell.bench.gates import load_thresholds
from lumencell.trust import constants

# constant name -> key path in benchmarks/thresholds.toml. Phase 2 appends entries.
MIRRORED: dict[str, tuple[str, ...]] = {
    "CLONAL_METRIC_VERSION": ("clonal", "metric_version"),
    "MIN_CLONE_SIZE": ("clonal", "min_clone_size"),
    "SENSITIVITY_CLONE_SIZE": ("clonal", "sensitivity_clone_size"),
    "MIN_CLONES": ("clonal", "min_clones"),
    "N_PERMUTATIONS": ("clonal", "n_permutations"),
    "CLONAL_SEED": ("clonal", "seed"),
    "PSEUDOTIME_BINS": ("clonal", "pseudotime_bins"),
    "CLONAL_P_MAX": ("clonal", "p_max"),
    "ORACLE_ETA_MIN": ("clonal", "synthetic", "oracle_eta_min"),
    "ORACLE_P_MAX": ("clonal", "synthetic", "oracle_p_max"),
}


def test_runtime_constants_equal_the_thresholds_file() -> None:
    thresholds = load_thresholds()
    for name, path in MIRRORED.items():
        value: Any = reduce(lambda node, key: node[key], path, thresholds)
        assert getattr(constants, name) == value, name
    declared = {n for n in vars(constants) if n.isupper()} - {"MAX_DONORS"}
    assert set(MIRRORED) == declared
