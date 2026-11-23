"""Runtime values of the clonal and robustness definitions.

bench/policy/thresholds.toml is not in the wheel, so the worker reads its values from here.
tests/unit/test_trust_constants.py keeps the two equal: change the TOML first, with a dated
changelog entry, then this file.
"""

CLONAL_METRIC_VERSION = 2  # 2: N2 with the lineage-onset split (N2a, 2026-11-17)
MIN_CLONE_SIZE = 2
SENSITIVITY_CLONE_SIZE = 3
MIN_CLONES = 20
N_PERMUTATIONS = 1000
CLONAL_SEED = 0
PSEUDOTIME_BINS = 10
CLONAL_P_MAX = 0.01
ORACLE_ETA_MIN = 0.05
ORACLE_P_MAX = 0.01

CELL_CYCLE_MIN_GENES = 20
CELL_CYCLE_SEED = 0
CELL_CYCLE_RATIO_MIN = 2.0
CELL_CYCLE_SHARE_MIN = 0.5
STABILITY_METRIC_VERSION = 1
N_RESEEDS = 2
N_SUBSAMPLES = 3
SUBSAMPLE_FRACTION = 0.8
STABILITY_SEED = 0
MATCH_JACCARD = 0.6
COVERAGE_MIN = 0.95
# [tree_support.calibrated] (2026-11-17): None means non-discriminating, so the check is left out.
TAU_TREE: float | None = None
TAU_LIN: float | None = None
TAU_PT: float | None = None
# [tree_support.verified] (2026-11-17): the recorded verification on fresh seeds. While it is False
# the served supported verdict is capped at provisional (ruling F14, amended by H-1).
VERDICT_VERIFIED = False
VERDICT_VERIFIED_ON = "2026-11-17"

# A validation limit (spec 4.1), not a threshold: not mirrored in the TOML.
MAX_DONORS = 1000
