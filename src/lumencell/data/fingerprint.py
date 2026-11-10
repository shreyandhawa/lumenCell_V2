"""A cheap identity for a dataset file, used to key the on-disk cache.

Hashing a multi-gigabyte .h5ad on every launch would be slow, so the fingerprint
combines the file's size and modification time with a hash of the cell and
gene names. Copying or touching the file gives a new fingerprint, which only
costs a cache miss; two different files never share one in practice.
"""

import hashlib
from collections.abc import Iterable
from pathlib import Path

FINGERPRINT_VERSION = "1"


def _names_digest(names: Iterable[str]) -> str:
    h = hashlib.sha256()
    for name in names:
        h.update(name.encode("utf-8"))
        h.update(b"\x00")
    return h.hexdigest()


def fingerprint(
    path: Path, obs_names: Iterable[str], var_names: Iterable[str], n_obs: int, n_vars: int
) -> str:
    stat = path.stat()
    h = hashlib.sha256()
    for part in (
        FINGERPRINT_VERSION,
        str(stat.st_size),
        str(stat.st_mtime_ns),
        str(n_obs),
        str(n_vars),
        _names_digest(obs_names),
        _names_digest(var_names),
    ):
        h.update(part.encode("utf-8"))
        h.update(b"|")
    return h.hexdigest()[:24]
