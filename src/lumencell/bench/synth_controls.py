"""Negative controls for the tree verdict (spec §5.3).

Four datasets with no branching tree in them: a blob (no structure), two unconnected groups, a
closed loop and a two-dimensional sheet. They are built here, never by the tree generator, and sit
outside the gated matrix (``[synthetic.controls]`` in ``benchmarks/thresholds.toml``): no gate,
baseline or experimental flag reads them. Counts use the same noise model as the tree generator.

Like the tree generator they write raw counts only, with no embedding and no ground-truth tree.
"""

from typing import Literal, get_args

import anndata as ad
import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import sparse

from lumencell.bench.synth import PRESETS, sample_counts

ControlKind = Literal["blob", "two_groups", "loop", "sheet"]
KINDS: tuple[ControlKind, ...] = get_args(ControlKind)

# Mirrored from [synthetic.controls]; tests/unit/test_synth_controls.py asserts they are equal.
PROGRAMME_SHARE = 0.10
LOOP_PROGRAMME_SHARE = 0.20
LFC_SCALE = 1.0
GROUP_LFC_SCALE = 2.0

CONTROL_STREAM = 7100502

_TWO_PI = 2.0 * np.pi


def _lfc(rng: np.random.Generator, scale: float, size: int) -> npt.NDArray[np.float64]:
    """Signed log fold changes with median magnitude ``scale`` (the tree generator's convention)."""
    magnitude = rng.lognormal(mean=np.log(scale), sigma=0.5, size=size)
    sign = np.where(rng.random(size) < 0.5, -1.0, 1.0)
    out: npt.NDArray[np.float64] = sign * magnitude
    return out


def generate_control(
    kind: ControlKind,
    n_cells: int = 3000,
    n_genes: int = 2000,
    seed: int = 0,
    preset: str = "realistic",
) -> ad.AnnData:
    """Raw counts for one negative control. Same arguments give identical counts."""
    if kind not in KINDS:
        raise ValueError(f"unknown control {kind!r}; choose from {list(KINDS)}")
    rng = np.random.default_rng([CONTROL_STREAM, KINDS.index(kind), seed])
    base = rng.gamma(shape=0.6, scale=1 / 0.3, size=n_genes)
    log_base = np.log(np.maximum(base, 1e-4))
    log_mu = np.broadcast_to(log_base, (n_cells, n_genes)).copy()
    extra: dict[str, npt.NDArray[np.float64] | npt.NDArray[np.int64]] = {}
    programme = np.zeros(n_genes, dtype=bool)
    root = 0

    if kind == "two_groups":
        group = (rng.random(n_cells) < 0.5).astype(np.int64)
        programme = rng.random(n_genes) < PROGRAMME_SHARE
        lfc = np.where(programme, _lfc(rng, GROUP_LFC_SCALE, n_genes), 0.0)
        log_mu += group[:, None] * lfc[None, :]
        extra["control_group"] = group
        root = int(np.flatnonzero(group == 0)[0])
    elif kind == "loop":
        theta = rng.uniform(0.0, _TWO_PI, n_cells)
        programme = rng.random(n_genes) < LOOP_PROGRAMME_SHARE
        magnitude = np.where(programme, np.abs(_lfc(rng, LFC_SCALE, n_genes)), 0.0)
        phase = rng.uniform(0.0, _TWO_PI, n_genes)
        log_mu += magnitude[None, :] * np.cos(theta[:, None] - phase[None, :])
        extra["control_angle"] = theta
        root = int(np.argmin(np.minimum(theta, _TWO_PI - theta)))
    elif kind == "sheet":
        x = rng.random(n_cells)
        y = rng.random(n_cells)
        draw = rng.random(n_genes)  # one draw per gene keeps the two sets disjoint
        x_genes = draw < PROGRAMME_SHARE
        y_genes = (draw >= PROGRAMME_SHARE) & (draw < 2 * PROGRAMME_SHARE)
        programme = x_genes | y_genes
        lfc_x = np.where(x_genes, _lfc(rng, LFC_SCALE, n_genes), 0.0)
        lfc_y = np.where(y_genes, _lfc(rng, LFC_SCALE, n_genes), 0.0)
        log_mu += (2 * x - 1)[:, None] * lfc_x[None, :] + (2 * y - 1)[:, None] * lfc_y[None, :]
        extra["control_x"] = x
        extra["control_y"] = y
        root = int(np.argmin(x * x + y * y))

    counts = sample_counts(log_mu, PRESETS[preset], rng)
    obs = pd.DataFrame(extra, index=pd.Index([f"cell{i:06d}" for i in range(n_cells)]))
    var = pd.DataFrame(index=pd.Index([f"G{i:05d}" for i in range(n_genes)]))
    adata = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
    adata.uns["lumencell_control"] = {
        "kind": kind,
        "seed": int(seed),
        "preset": preset,
        "root_cell": str(obs.index[root]),
        "programme_genes": int(programme.sum()),
    }
    return adata


def control_root_cell(adata: ad.AnnData) -> str:
    """The obs name a run on this control is rooted at (as ``cell:<name>``)."""
    return str(adata.uns["lumencell_control"]["root_cell"])
