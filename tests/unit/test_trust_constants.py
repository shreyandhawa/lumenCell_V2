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
    "CELL_CYCLE_MIN_GENES": ("cell_cycle", "min_genes"),
    "CELL_CYCLE_SEED": ("cell_cycle", "seed"),
    "CELL_CYCLE_RATIO_MIN": ("cell_cycle", "lineage_ratio_min"),
    "CELL_CYCLE_SHARE_MIN": ("cell_cycle", "lineage_share_min"),
    "STABILITY_METRIC_VERSION": ("stability", "metric_version"),
    "N_RESEEDS": ("stability", "n_reseeds"),
    "N_SUBSAMPLES": ("stability", "n_subsamples"),
    "SUBSAMPLE_FRACTION": ("stability", "subsample_fraction"),
    "STABILITY_SEED": ("stability", "seed"),
    "MATCH_JACCARD": ("stability", "match_jaccard"),
    "COVERAGE_MIN": ("tree_support", "coverage_min"),
}


def test_runtime_constants_equal_the_thresholds_file() -> None:
    thresholds = load_thresholds()
    for name, path in MIRRORED.items():
        value: Any = reduce(lambda node, key: node[key], path, thresholds)
        assert getattr(constants, name) == value, name
    declared = {n for n in vars(constants) if n.isupper()} - {"MAX_DONORS"}
    assert set(MIRRORED) | set(CALIBRATED) == declared


# tau constant -> (its check, its key in [tree_support.calibrated]). None when non-discriminating.
CALIBRATED: dict[str, tuple[str, str]] = {
    "TAU_TREE": ("C2", "tau_tree"),
    "TAU_LIN": ("C3", "tau_lin"),
    "TAU_PT": ("C4", "tau_pt"),
}


def test_calibrated_taus_equal_the_thresholds_file() -> None:
    calibrated = load_thresholds()["tree_support"]["calibrated"]
    for name, (check, key) in CALIBRATED.items():
        value = getattr(constants, name)
        if check in calibrated["non_discriminating"]:
            assert value is None, name
            assert key not in calibrated, name
        else:
            assert value == calibrated[key], name


def test_every_scalar_clonal_key_is_mirrored() -> None:
    clonal = load_thresholds()["clonal"]
    scalars = {k for k, v in clonal.items() if not isinstance(v, dict)}
    assert scalars == {path[1] for path in MIRRORED.values() if len(path) == 2 and path[0] == "clonal"}


# Robustness scalars the worker or server never reads at run time, so they have no constant.
NOT_MIRRORED = {("cell_cycle", "gene_set")}


def test_every_scalar_robustness_key_is_mirrored() -> None:
    cfg = load_thresholds()
    scalars = {
        (table, k)
        for table in ("cell_cycle", "stability", "tree_support")
        for k, v in cfg[table].items()
        if not isinstance(v, dict)
    }
    mirrored = {path for path in MIRRORED.values() if len(path) == 2 and path[0] != "clonal"}
    assert scalars - NOT_MIRRORED == mirrored


def test_the_clonal_p_gate_equals_the_runtime_p_max() -> None:
    gates = load_thresholds()["synthetic"]["gates"]
    assert gates["clonal_p"]["max"] == constants.CLONAL_P_MAX
