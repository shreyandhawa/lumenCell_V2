"""Runtime values of the clonal definitions.

benchmarks/thresholds.toml is not in the wheel, so the worker reads its values from here.
tests/unit/test_trust_constants.py keeps the two equal: change the TOML first, with a dated
changelog entry, then this file.
"""

CLONAL_METRIC_VERSION = 1
MIN_CLONE_SIZE = 2
SENSITIVITY_CLONE_SIZE = 3
MIN_CLONES = 20
N_PERMUTATIONS = 1000
CLONAL_SEED = 0
PSEUDOTIME_BINS = 10
CLONAL_P_MAX = 0.01
ORACLE_ETA_MIN = 0.05
ORACLE_P_MAX = 0.01

# A validation limit (spec 4.1), not a threshold: not mirrored in the TOML.
MAX_DONORS = 1000
