"""Draw an inferred trajectory on a display embedding (2-D, or the 3-D diffusion view).

Inference happens in PCA / diffusion space; UMAP and friends are only for
looking at it. Everything drawn here is therefore a *projection* of the
inferred trajectory, and says so:

* milestone positions are the weighted mean display position of the cells
  associated with each milestone (the method's ``milestone_weights``);
* each lineage's curve is either the method's own fitted curve (a principal
  tree: each node placed at the weighted mean of its cells), or, for methods
  without one, the mean display position of its cells in 40 pseudotime quantile
  bins, lightly smoothed, starting at the root milestone;
* ``fidelity`` is the Spearman correlation between a cell's pseudotime and the
  arc-length position of its nearest point on the drawn curve. A low value means
  the embedding distorts that lineage, and the UI should say so.
"""

import re
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
from scipy import sparse
from scipy.stats import spearmanr

from lumencell.engine.schema import Notice, TrajectoryResult

N_BINS = 40
FIDELITY_WARN = 0.8
# The diffusion views: DCi is column i of the diffusion map (column 0 is the constant one).
_DC_ALIASES = {"dc12": (1, 2), "dc123": (1, 2, 3)}
_DC_NAME = re.compile(r"dc:([1-9]\d*),([1-9]\d*)(?:,([1-9]\d*))?")


def diffusion_axes(name: str) -> tuple[int, ...] | None:
    """The diffusion components a view name asks for, in order: ``dc12``, ``dc123`` or ``dc:i,j[,k]``.

    None for any other name. Whether the dataset has those components is the caller's check.
    """
    if name in _DC_ALIASES:
        return _DC_ALIASES[name]
    m = _DC_NAME.fullmatch(name)
    return tuple(int(g) for g in m.groups() if g is not None) if m else None


@dataclass(frozen=True)
class LineageCurve:
    lineage: int
    points: list[tuple[float, ...]]
    fidelity: float | None


@dataclass(frozen=True)
class Projection:
    milestones: dict[int, tuple[float, ...]]
    curves: list[LineageCurve]
    curve_source: str = "binned"
    warnings: list[Notice] = field(default_factory=list)


def _smooth(points: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    if len(points) < 3:
        return points
    padded = np.vstack([points[:1], points, points[-1:]])
    out: npt.NDArray[np.float64] = (padded[:-2] + padded[1:-1] + padded[2:]) / 3.0
    out[0], out[-1] = points[0], points[-1]
    return out


def _arc_position(cells: npt.NDArray[np.float64], curve: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Arc-length position (0..1) of each cell's nearest point on a polyline."""
    seg_start, seg_end = curve[:-1], curve[1:]
    seg_vec = seg_end - seg_start
    seg_len = np.linalg.norm(seg_vec, axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg_len)])
    total = cum[-1] if cum[-1] > 0 else 1.0
    best_d = np.full(len(cells), np.inf)
    best_pos = np.zeros(len(cells))
    for k in range(len(seg_vec)):
        denom = float(seg_vec[k] @ seg_vec[k])
        t = (
            np.clip(((cells - seg_start[k]) @ seg_vec[k]) / denom, 0, 1)
            if denom > 0
            else np.zeros(len(cells))
        )
        proj = seg_start[k] + t[:, None] * seg_vec[k]
        d = np.linalg.norm(cells - proj, axis=1)
        better = d < best_d
        best_d = np.where(better, d, best_d)
        best_pos = np.where(better, cum[k] + t * seg_len[k], best_pos)
    out: npt.NDArray[np.float64] = best_pos / total
    return out


def project(result: TrajectoryResult, coords: npt.NDArray[np.floating], embedding: str) -> Projection:
    xy = np.asarray(coords, dtype=np.float64)[:, :3]
    w = result.milestone_weights.tocsc()
    col_sum = np.asarray(w.sum(axis=0)).ravel()
    milestone_xy = (w.T @ xy) / np.maximum(col_sum, 1e-12)[:, None]
    milestones = {
        m.id: tuple(float(v) for v in milestone_xy[i])
        for i, m in enumerate(result.milestones)
        if col_sum[i] > 0
    }
    root_id = next(m.id for m in result.milestones if m.kind == "root")

    axes = diffusion_axes(embedding)
    if axes:
        *rest, last = [f"DC{i}" for i in axes]
        where, hint = f"{', '.join(rest)} and {last}", "Try other diffusion components."
    else:
        where, hint = f"the {embedding} embedding", "Try the diffusion (DC) view."

    pt = result.pseudotime.astype(np.float64)
    curves: list[LineageCurve] = []
    warnings: list[Notice] = []
    fitted = result.curve_nodes is not None and result.curve_paths is not None
    node_xy = None
    if fitted:
        cn = sparse.csc_matrix(result.curve_nodes)
        node_sum = np.asarray(cn.sum(axis=0)).ravel()
        node_xy = (cn.T @ xy) / np.maximum(node_sum, 1e-12)[:, None]
    for lin in result.lineages:
        members = result.membership[lin.id] & np.isfinite(pt)
        idx = np.flatnonzero(members)
        if node_xy is not None and result.curve_paths is not None:
            curve = np.asarray(node_xy[result.curve_paths[lin.id]], dtype=np.float64)
            if len(curve) < 2 or idx.size < 3:
                curves.append(LineageCurve(lin.id, [], None))
                continue
        else:
            if members.sum() < N_BINS:
                curves.append(LineageCurve(lin.id, [], None))
                continue
            edges = np.quantile(pt[idx], np.linspace(0, 1, N_BINS + 1))
            bins = np.clip(np.searchsorted(edges, pt[idx], side="right") - 1, 0, N_BINS - 1)
            means = np.array([xy[idx[bins == b]].mean(axis=0) for b in range(N_BINS) if np.any(bins == b)])
            if root_id in milestones:
                means = np.vstack([np.array(milestones[root_id])[None, :], means])
            curve = _smooth(means)
        pos = _arc_position(xy[idx], curve)
        fid = float(spearmanr(pt[idx], pos).statistic) if np.ptp(pos) > 0 else None
        if fid is not None and fid < FIDELITY_WARN:
            warnings.append(
                Notice(
                    "embedding_distorts_lineage",
                    "warning",
                    f"On {where}, lineage {lin.name!r} is drawn with fidelity {fid:.2f}: "
                    f"positions along the curve do not follow pseudotime well. {hint}",
                )
            )
        curves.append(LineageCurve(lin.id, [tuple(float(v) for v in p) for p in curve], fid))
    return Projection(
        milestones=milestones, curves=curves, curve_source="fitted" if fitted else "binned", warnings=warnings
    )
