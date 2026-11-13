"""A synthetic clone layer on top of a synthetic dataset, and two controls.

Clone labels and donors are drawn from the true lineage structure with a fate bias: with
probability ``fate_bias`` a labelled cell on a leaf segment draws its clone from that leaf's pool,
otherwise from the shared trunk pool. Everything is drawn from its own stream,
``np.random.default_rng([seed, CLONE_STREAM])``, so the layer never touches the generator's draws.

"Lineage j" means the index of the true lineage in SORTED leaf-name order (a leaf is a segment
that is no one's parent). ``Truth.names`` follows topology order, which differs for binary2; the
two are never mixed here.

Draw order (fixed): (1) ``labelled = rng.random(n) < label_rate``; (2) donors: uniform
``rng.integers(n_donors, size=n)``; for ``donor_confound`` that draw is followed by
``rng.random(n)`` (keep the home donor?) and ``rng.integers(n_donors - 1, size=n)`` (which other
donor), used only for cells exclusive to one lineage; (3) ``join_leaf = rng.random(n) < fate_bias``
(control ``none`` only); (4) per donor in ascending order, pools in the order trunk then leaves by
sorted name: each pool is permuted with ``rng.permutation`` and cut into consecutive clones of
``draw_clone_size`` cells, only the last truncated. Labels are ``c0, c1, ...``, counted across a
donor's pools and restarting in each donor.
"""

import dataclasses
from dataclasses import dataclass
from typing import Literal

import anndata as ad
import numpy as np
import numpy.typing as npt
import pandas as pd

from lumencell.bench.truth import Truth

FATE_BIAS = 0.7
N_DONORS = 3
LABEL_RATE = 0.85
SINGLETON_SHARE = 0.5
ZIPF_A = 2.2
DONOR_CONFOUND_SHARE = 0.8
CLONE_STREAM = 7100501

Control = Literal["none", "depth_only", "donor_confound"]
IntArray = npt.NDArray[np.int64]


@dataclass(frozen=True)
class ClonalLayer:
    fate_bias: float = FATE_BIAS
    n_donors: int = N_DONORS
    label_rate: float = LABEL_RATE
    control: Control = "none"


DEFAULT_LAYER = ClonalLayer()


def draw_clone_size(rng: np.random.Generator) -> int:
    """1 with probability ``SINGLETON_SHARE``, else ``1 + Zipf(ZIPF_A)``."""
    if rng.random() < SINGLETON_SHARE:
        return 1
    return 1 + int(rng.zipf(ZIPF_A))


def add_clones(adata: ad.AnnData, layer: ClonalLayer, seed: int) -> ad.AnnData:
    """A copy of ``adata`` with obs ``clone_id`` (NaN where unlabelled) and ``donor`` (every cell).

    Needs obs ``true_segment`` and ``true_pseudotime`` and ``uns["lumencell_truth"]["segments"]``.
    The layer, seed and stream are recorded in ``uns["lumencell_truth"]["clones"]``.
    """
    n = adata.n_obs
    segments: dict[str, dict[str, str]] = dict(adata.uns["lumencell_truth"]["segments"])
    seg = adata.obs["true_segment"].astype(str).to_numpy()
    t = adata.obs["true_pseudotime"].to_numpy(dtype=np.float64)
    leaves = sorted(set(segments) - {str(v["parent"]) for v in segments.values()})
    rng = np.random.default_rng([seed, CLONE_STREAM])

    labelled = rng.random(n) < layer.label_rate
    donor = _draw_donors(rng, layer, seg, segments, leaves)
    donor_ids = np.full(n, -1, dtype=np.int64)
    clone_ids = np.full(n, None, dtype=object)

    if layer.control == "none":
        join_leaf = rng.random(n) < layer.fate_bias
    for d in range(layer.n_donors):
        mine = labelled & (donor == d)
        if layer.control == "depth_only":
            pools = [np.flatnonzero(mine)[np.argsort(t[mine], kind="stable")]]
            permute = False
        elif layer.control == "donor_confound":
            pools, permute = [np.flatnonzero(mine)], True
        else:
            on_leaf = mine & join_leaf
            pools = [np.flatnonzero(mine & ~(on_leaf & np.isin(seg, leaves)))]
            pools += [np.flatnonzero(on_leaf & (seg == leaf)) for leaf in leaves]
            permute = True
        count = 0
        for pool in pools:
            cells = rng.permutation(pool) if permute else pool
            pos = 0
            while pos < len(cells):
                size = draw_clone_size(rng)
                clone_ids[cells[pos : pos + size]] = f"c{count}"
                pos += size
                count += 1
        donor_ids[donor == d] = d

    out = adata.copy()
    out.obs["clone_id"] = pd.Categorical(clone_ids)
    out.obs["donor"] = pd.Categorical(
        [f"d{d}" for d in donor_ids], categories=[f"d{d}" for d in range(layer.n_donors)]
    )
    out.uns["lumencell_truth"]["clones"] = {**dataclasses.asdict(layer), "seed": seed, "stream": CLONE_STREAM}
    return out


def clone_truth(truth: Truth, adata: ad.AnnData) -> Truth:
    """``truth`` with the clone labels and donors of ``adata`` (made by :func:`add_clones`)."""
    return dataclasses.replace(
        truth,
        clone_ids=adata.obs["clone_id"].to_numpy(dtype=object),
        donors=adata.obs["donor"].to_numpy(dtype=object),
    )


def _draw_donors(
    rng: np.random.Generator,
    layer: ClonalLayer,
    seg: npt.NDArray[np.str_],
    segments: dict[str, dict[str, str]],
    leaves: list[str],
) -> IntArray:
    n = len(seg)
    donor = rng.integers(layer.n_donors, size=n)
    if layer.control != "donor_confound":
        return donor
    keep_home = rng.random(n) < DONOR_CONFOUND_SHARE
    other = rng.integers(max(layer.n_donors - 1, 1), size=n)
    home = np.full(n, -1, dtype=np.int64)  # j mod D for cells on exactly one lineage
    paths = [_path(segments, leaf) for leaf in leaves]
    for name in segments:
        owners = [j for j, path in enumerate(paths) if name in path]
        if len(owners) == 1:
            home[seg == name] = owners[0] % layer.n_donors
    away = (home + 1 + other) % layer.n_donors
    return np.where(home < 0, donor, np.where(keep_home, home, away)).astype(np.int64)


def _path(segments: dict[str, dict[str, str]], leaf: str) -> set[str]:
    path, cur = set(), leaf
    while cur:
        path.add(cur)
        cur = str(segments[cur]["parent"])
    return path
