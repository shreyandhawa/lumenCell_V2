"""Synthetic single-cell data with a known trajectory.

Modelled on Splatter's ``splatPaths`` (Zappia et al. 2017): cells sit on the
segments of a known tree, each gene's mean expression changes along the path
from the root, and counts are drawn from an overdispersed gamma-Poisson model.

The generator writes **raw counts only, with no embedding**. V1's demo data
placed cells directly on a drawn tree, which made trajectory inference trivially
easy; here the embedding must always come from lumenCell's real preprocessing.

Expression model, per gene ``g`` and cell ``c`` on segment ``s`` at fraction ``u``:

* base mean ``mu_g ~ Gamma(0.6, rate 0.3)`` (Splatter's defaults);
* each segment changes a random subset of genes (probability ``de_prob``) by a
  log fold change of median magnitude ``de_scale``, ramped in smoothly along the
  segment; a child segment inherits its ancestors' full changes, so expression
  is continuous across forks;
* a fraction of genes also get a transient bump in pseudotime;
* library size ``L_c ~ LogNormal(log lib_size, lib_sd)``;
* counts ``~ Poisson(Gamma(1/bcv^2, m * bcv^2))`` with
  ``bcv = bcv_common + 1/sqrt(m)`` for cell mean ``m`` (Splatter's BCV trend).
"""

from collections.abc import Iterator
from dataclasses import asdict, dataclass, replace
from typing import Any

import anndata as ad
import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import sparse

from lumencell.bench.truth import Truth

ROOT_SEGMENT = "trunk"


@dataclass(frozen=True)
class Segment:
    name: str
    parent: str | None
    length: float
    weight: float = 1.0


# Every topology's root-to-tip paths have length 1, so true pseudotime is in [0, 1].
# The bifurcation has a short trunk and long branches: its tree diameter runs from
# tip to tip, which is exactly the shape that fooled V1's automatic root.
TOPOLOGIES: dict[str, tuple[Segment, ...]] = {
    "linear": (Segment("trunk", None, 0.5), Segment("tail", "trunk", 0.5)),
    "bifurcation": (Segment("trunk", None, 0.4), Segment("A", "trunk", 0.6), Segment("B", "trunk", 0.6)),
    "trifurcation": (
        Segment("trunk", None, 0.4),
        Segment("A", "trunk", 0.6),
        Segment("B", "trunk", 0.6),
        Segment("C", "trunk", 0.6),
    ),
    "binary2": (
        Segment("trunk", None, 0.3),
        Segment("A", "trunk", 0.35),
        Segment("B", "trunk", 0.7),
        Segment("A1", "A", 0.35),
        Segment("A2", "A", 0.35),
    ),
}


@dataclass(frozen=True)
class SynthParams:
    de_prob: float
    de_scale: float
    bcv_common: float = 0.2
    transient_frac: float = 0.05
    lib_size: float = 5000.0
    lib_sd: float = 0.25

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


# Calibrated (before any inference method existed) against two properties after a
# standard log-normalise + HVG + PCA + kNN pipeline, 3k cells x 2k genes:
#   PC1 variance share: real scRNA-seq is typically 3-15%.
#   kNN branch purity well past a fork: how separable the branches are.
#                 PC1 share      branch purity
#   easy          5-12%          0.91-0.99
#   realistic     3.1-6.5%       0.89-0.97   (gated benchmark)
#   hard          1.1-2.4%       0.82-0.90   (reported only)
PRESETS: dict[str, SynthParams] = {
    "easy": SynthParams(de_prob=0.2, de_scale=1.0),
    "realistic": SynthParams(de_prob=0.10, de_scale=1.0, bcv_common=0.3),
    "hard": SynthParams(de_prob=0.08, de_scale=0.6, bcv_common=0.3),
}

_CHUNK = 5000
_BCV_MAX = 10.0  # caps 1/sqrt(m) for near-zero means so the gamma shape stays positive


def _children(topology: tuple[Segment, ...]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {s.name: [] for s in topology}
    for s in topology:
        if s.parent is not None:
            out[s.parent].append(s.name)
    return out


def _ancestors(topology: tuple[Segment, ...], name: str) -> list[str]:
    """Segments from the root down to, but excluding, ``name``."""
    by_name = {s.name: s for s in topology}
    chain: list[str] = []
    cur = by_name[name].parent
    while cur is not None:
        chain.append(cur)
        cur = by_name[cur].parent
    return chain[::-1]


def lineage_paths(topology: tuple[Segment, ...]) -> dict[str, list[str]]:
    """Root-to-leaf segment paths, keyed by the leaf segment's name."""
    kids = _children(topology)
    leaves = [s.name for s in topology if not kids[s.name]]
    return {leaf: [*_ancestors(topology, leaf), leaf] for leaf in leaves}


def _segment_starts(topology: tuple[Segment, ...]) -> dict[str, float]:
    by_name = {s.name: s for s in topology}
    return {s.name: sum(by_name[a].length for a in _ancestors(topology, s.name)) for s in topology}


def _smoothstep(u: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    out: npt.NDArray[np.float64] = u * u * (3.0 - 2.0 * u)
    return out


@dataclass(frozen=True)
class _GeneProgram:
    log_base: npt.NDArray[np.float64]  # (G,)
    seg_lfc: dict[str, npt.NDArray[np.float64]]  # segment -> (G,) log fold change
    bump_amp: npt.NDArray[np.float64]  # (G,), 0 for non-transient genes
    bump_centre: npt.NDArray[np.float64]
    bump_width: npt.NDArray[np.float64]


def _gene_program(
    topology: tuple[Segment, ...], n_genes: int, params: SynthParams, rng: np.random.Generator
) -> _GeneProgram:
    base = rng.gamma(shape=0.6, scale=1 / 0.3, size=n_genes)
    log_base = np.log(np.maximum(base, 1e-4))
    seg_lfc: dict[str, npt.NDArray[np.float64]] = {}
    for seg in topology:
        de = rng.random(n_genes) < params.de_prob
        magnitude = rng.lognormal(mean=np.log(params.de_scale), sigma=0.5, size=n_genes)
        sign = np.where(rng.random(n_genes) < 0.5, -1.0, 1.0)
        seg_lfc[seg.name] = np.where(de, sign * magnitude, 0.0)
    transient = rng.random(n_genes) < params.transient_frac
    bump_amp = np.where(
        transient, rng.lognormal(mean=np.log(params.de_scale * 1.5), sigma=0.3, size=n_genes), 0.0
    )
    return _GeneProgram(
        log_base=log_base,
        seg_lfc=seg_lfc,
        bump_amp=bump_amp,
        bump_centre=rng.random(n_genes),
        bump_width=rng.uniform(0.05, 0.15, n_genes),
    )


def _cells_to_segments(
    topology: tuple[Segment, ...], n_cells: int, rng: np.random.Generator
) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.float64]]:
    weights = np.array([s.length * s.weight for s in topology])
    seg_idx = rng.choice(len(topology), size=n_cells, p=weights / weights.sum())
    return seg_idx, rng.random(n_cells)


def _log_means(
    topology: tuple[Segment, ...],
    program: _GeneProgram,
    seg_idx: npt.NDArray[np.int64],
    u: npt.NDArray[np.float64],
    t: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """(cells, genes) log of relative expression, before library-size scaling."""
    out = np.broadcast_to(program.log_base, (seg_idx.size, program.log_base.size)).copy()
    for i, seg in enumerate(topology):
        sel = seg_idx == i
        if not sel.any():
            continue
        inherited = sum((program.seg_lfc[a] for a in _ancestors(topology, seg.name)), np.zeros_like(out[0]))
        ramp = _smoothstep(u[sel])[:, None]
        out[sel] += inherited[None, :] + ramp * program.seg_lfc[seg.name][None, :]
    bumps = program.bump_amp > 0
    if bumps.any():
        z = (t[:, None] - program.bump_centre[None, bumps]) / program.bump_width[None, bumps]
        out[:, bumps] += program.bump_amp[None, bumps] * np.exp(-0.5 * z * z)
    return out


def _sample_counts(
    log_mu: npt.NDArray[np.float64], params: SynthParams, rng: np.random.Generator
) -> npt.NDArray[np.float32]:
    rel = np.exp(log_mu - log_mu.max(axis=1, keepdims=True))
    rel /= rel.sum(axis=1, keepdims=True)
    lib = rng.lognormal(mean=np.log(params.lib_size), sigma=params.lib_sd, size=(log_mu.shape[0], 1))
    m = rel * lib
    bcv = params.bcv_common + np.minimum(1.0 / np.sqrt(np.maximum(m, 1e-12)), _BCV_MAX)
    shape = 1.0 / (bcv * bcv)
    lam = rng.gamma(shape=shape, scale=m / shape)
    return rng.poisson(lam).astype(np.float32)


def _chunks(n: int) -> Iterator[slice]:
    for start in range(0, n, _CHUNK):
        yield slice(start, min(start + _CHUNK, n))


def generate(
    topology: str = "bifurcation",
    preset: str | SynthParams = "realistic",
    n_cells: int = 3000,
    n_genes: int = 2000,
    seed: int = 0,
    store_root: bool = False,
) -> ad.AnnData:
    """Generate a dataset. Same arguments and seed give identical counts.

    ``store_root=True`` also writes ``uns["iroot"]`` (the earliest cell), so the
    "root stored in the file" path can be exercised; it is off by default because
    lumenCell must never rely on a root it was not given.
    """
    if topology not in TOPOLOGIES:
        raise ValueError(f"unknown topology {topology!r}; choose from {sorted(TOPOLOGIES)}")
    params = PRESETS[preset] if isinstance(preset, str) else preset
    preset_name = preset if isinstance(preset, str) else "custom"
    segs = TOPOLOGIES[topology]
    rng = np.random.default_rng(seed)

    program = _gene_program(segs, n_genes, params, rng)
    seg_idx, u = _cells_to_segments(segs, n_cells, rng)
    starts = _segment_starts(segs)
    seg_start = np.array([starts[s.name] for s in segs])[seg_idx]
    seg_len = np.array([s.length for s in segs])[seg_idx]
    t = seg_start + u * seg_len

    blocks = []
    for sl in _chunks(n_cells):
        log_mu = _log_means(segs, program, seg_idx[sl], u[sl], t[sl])
        blocks.append(sparse.csr_matrix(_sample_counts(log_mu, params, rng)))
    X = sparse.vstack(blocks, format="csr", dtype=np.float32)

    names = [s.name for s in segs]
    paths = lineage_paths(segs)
    seg_names = np.array(names)[seg_idx]
    on_lineages = {leaf: np.isin(seg_names, path) for leaf, path in paths.items()}
    n_on = np.sum(np.stack(list(on_lineages.values())), axis=0)
    true_lineage = np.where(
        n_on > 1, "shared", np.array(list(paths))[np.argmax(np.stack(list(on_lineages.values())), axis=0)]
    )
    stage = np.where(u < 0.5, "early", "late")
    labels = np.char.add(np.char.add(seg_names, "_"), stage)
    label_order = [f"{s}_{p}" for s in names for p in ("early", "late")]

    obs = pd.DataFrame(
        {
            "true_pseudotime": t.astype(np.float32),
            "true_segment": pd.Categorical(seg_names, categories=names),
            "true_lineage": pd.Categorical(true_lineage, categories=["shared", *paths]),
            "true_label": pd.Categorical(labels, categories=label_order),
        },
        index=pd.Index([f"cell{i:06d}" for i in range(n_cells)]),
    )
    var = pd.DataFrame(index=pd.Index([f"G{i:05d}" for i in range(n_genes)]))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["lumencell_truth"] = {
        "topology": topology,
        "preset": preset_name,
        "params": params.to_dict(),
        "seed": seed,
        "root_label": f"{ROOT_SEGMENT}_early",
        "segments": {s.name: {"parent": s.parent or "", "length": s.length} for s in segs},
        "lineages": {leaf: list(path) for leaf, path in paths.items()},
        "terminal_labels": {leaf: f"{leaf}_late" for leaf in paths},
    }
    if store_root:
        adata.uns["iroot"] = int(np.argmin(t))
    return adata


def truth_from_adata(adata: ad.AnnData) -> Truth:
    """Build the scoring :class:`Truth` from a dataset made by :func:`generate`."""
    meta: dict[str, Any] = dict(adata.uns["lumencell_truth"])
    seg = adata.obs["true_segment"].astype(str).to_numpy()
    lineages = {leaf: np.isin(seg, [str(s) for s in path]) for leaf, path in dict(meta["lineages"]).items()}
    root = (adata.obs["true_label"].astype(str) == meta["root_label"]).to_numpy()
    return Truth(
        order=adata.obs["true_pseudotime"].to_numpy(dtype=np.float64),
        lineages=lineages,
        root=root,
        meta=meta,
    )


def with_params(preset: str, **overrides: float) -> SynthParams:
    """A preset with some parameters changed (used by calibration sweeps)."""
    return replace(PRESETS[preset], **overrides)
