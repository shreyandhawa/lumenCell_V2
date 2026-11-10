"""Draw an inferred trajectory on a 2-D (or 3-D diffusion) display embedding.

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

from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt
from scipy import sparse
from scipy.stats import spearmanr

from lumencell.engine.schema import Notice, TrajectoryResult

N_BINS = 40
FIDELITY_WARN = 0.8


Coords = tuple[float, float] | tuple[float, float, float]


def _coords(row: npt.NDArray[np.float64]) -> Coords:
    if len(row) == 3:
        return (float(row[0]), float(row[1]), float(row[2]))
    return (float(row[0]), float(row[1]))


@dataclass(frozen=True)
class LineageCurve:
    lineage: int
    points: list[Coords]
    fidelity: float | None


@dataclass(frozen=True)
class Projection:
    milestones: dict[int, Coords]
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
    """Project onto a 2-D or 3-D embedding; extra columns beyond the third are ignored."""
    xy = np.asarray(coords, dtype=np.float64)[:, :3]
    w = result.milestone_weights.tocsc()
    col_sum = np.asarray(w.sum(axis=0)).ravel()
    milestone_xy = (w.T @ xy) / np.maximum(col_sum, 1e-12)[:, None]
    milestones = {m.id: _coords(milestone_xy[i]) for i, m in enumerate(result.milestones) if col_sum[i] > 0}
    root_id = next(m.id for m in result.milestones if m.kind == "root")

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
                    f"On the {embedding} embedding, lineage {lin.name!r} is drawn with fidelity {fid:.2f}: "
                    "positions along the curve do not follow pseudotime well. Try the diffusion (DC) view.",
                )
            )
        curves.append(LineageCurve(lin.id, [_coords(p) for p in curve], fid))
    return Projection(
        milestones=milestones, curves=curves, curve_source="fitted" if fitted else "binned", warnings=warnings
    )
