"""Where a trajectory starts.

lumenCell never guesses the root. V1 picked one end of the tree's diameter,
which on a Y-shaped tree is a branch tip, and silently reversed pseudotime for
the default view. Here the root is always an explicit choice:

* ``cell:<obs_name>``: one cell;
* ``cluster:<obs key>=<value>``: a cluster, resolved deterministically to the
  cell in it lying farthest, in diffusion space, from the cells outside it
  (the cluster's outer end);
* ``file``: the root stored in the .h5ad (``uns['iroot']`` or a boolean obs
  ``root`` / ``is_root`` column), used only when asked for.

:func:`suggest_roots` offers candidates with their evidence. Nothing here applies
a suggestion; the caller (the user) decides.
"""

import math
import re
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt
from scipy import sparse, stats
from scipy.sparse.csgraph import connected_components

from lumencell.engine.context import LEIDEN_KEY, InferenceContext
from lumencell.engine.schema import Notice, ResolvedRoot

ROOT_COLUMNS = ("root", "is_root")
TIME_COLUMN = re.compile(r"(^|_)(day|days|time|timepoint|time_point|stage|age|hpf|dpf|hours?)($|_)", re.I)
_NUMBER = re.compile(r"[-+]?\d*\.?\d+")
SMALL_COMPONENT = 0.5  # warn when the root's connected group holds less than this share of cells


class RootError(ValueError):
    """The root spec cannot be resolved; the message says why and what is available."""


@dataclass(frozen=True)
class CellRoot:
    obs_name: str


@dataclass(frozen=True)
class ClusterRoot:
    key: str
    value: str


@dataclass(frozen=True)
class FileRoot:
    pass


RootSpec = CellRoot | ClusterRoot | FileRoot


def parse_root(text: str) -> RootSpec:
    """Parse the CLI form: ``cell:NAME``, ``cluster:KEY=VALUE`` or ``file``."""
    text = text.strip()
    if text == "file":
        return FileRoot()
    kind, sep, rest = text.partition(":")
    if sep and kind == "cell" and rest:
        return CellRoot(rest)
    if sep and kind == "cluster":
        key, eq, value = rest.partition("=")
        if eq and key and value:
            return ClusterRoot(key, value)
    raise RootError(f"cannot parse root {text!r}; use cell:<name>, cluster:<column>=<value> or file")


def format_root(spec: RootSpec) -> str:
    if isinstance(spec, CellRoot):
        return f"cell:{spec.obs_name}"
    if isinstance(spec, ClusterRoot):
        return f"cluster:{spec.key}={spec.value}"
    return "file"


# ---- resolution ---------------------------------------------------------------------


def _extreme_cell(ctx: InferenceContext, members: npt.NDArray[np.bool_]) -> int:
    """The member farthest in diffusion space from the centroid of the non-members."""
    space = ctx.diffusion_space()
    idx = np.flatnonzero(members)
    centre = space[~members].mean(axis=0) if (~members).any() else space.mean(axis=0)
    dist = np.linalg.norm(space[idx] - centre[None, :], axis=1)
    return int(idx[int(np.argmax(dist))])  # argmax returns the first maximum: ties -> lowest index


def _file_root_cells(ctx: InferenceContext) -> tuple[npt.NDArray[np.bool_] | None, str]:
    if ctx.uns_root is not None:
        mask = np.zeros(ctx.n_obs, dtype=bool)
        mask[ctx.uns_root] = True
        return mask, f"uns['iroot'] = {ctx.uns_root}"
    for col in ROOT_COLUMNS:
        if col in ctx.obs.columns:
            values = ctx.obs[col]
            # {0, 1} also matches True/False, since True == 1 and False == 0.
            if values.dtype == bool or set(values.dropna().unique()) <= {0, 1}:
                mask = values.fillna(False).astype(bool).to_numpy()
                if mask.any():
                    return mask, f"obs['{col}'] ({int(mask.sum())} cells)"
    return None, ""


def _component_warnings(ctx: InferenceContext, root: int) -> list[Notice]:
    graph = sparse.csr_matrix(ctx.artifacts["connectivities"])
    _n, labels = connected_components(graph, directed=False)
    in_root = labels == labels[root]
    share = float(in_root.mean())
    out: list[Notice] = []
    n_out = int((~in_root).sum())
    if share < SMALL_COMPONENT:
        out.append(
            Notice(
                "root_in_small_group",
                "warning",
                f"The root is in a group holding only {share:.0%} of cells; the rest are not "
                "connected to it in the neighbour graph and will get no pseudotime.",
            )
        )
    elif n_out:
        out.append(
            Notice(
                "disconnected_cells",
                "info",
                f"{n_out} cells ({1 - share:.1%}) are not connected to the root and get no pseudotime.",
            )
        )
    return out


def resolve_root(ctx: InferenceContext, spec: RootSpec) -> ResolvedRoot:
    if isinstance(spec, CellRoot):
        matches = np.flatnonzero(ctx.obs_names == spec.obs_name)
        if matches.size == 0:
            raise RootError(f"no cell is named {spec.obs_name!r}")
        cell = int(matches[0])
        how = "the cell you chose"
        key: str | None = None
        label: str | None = None
    elif isinstance(spec, ClusterRoot):
        if spec.key not in ctx.obs.columns:
            raise RootError(f"there is no obs column {spec.key!r}")
        labels = ctx.labels(spec.key)
        members = labels == spec.value
        if not members.any():
            available = sorted(set(labels.tolist()))
            shown = ", ".join(available[:20]) + (" …" if len(available) > 20 else "")
            raise RootError(f"no cell has {spec.key} = {spec.value!r}; values are: {shown}")
        if members.all():
            raise RootError(f"every cell has {spec.key} = {spec.value!r}; a root cluster must be a subset")
        cell = _extreme_cell(ctx, members)
        how = (
            f"the cell of {spec.key} = {spec.value} farthest, in diffusion space, "
            "from the cells outside that cluster"
        )
        key, label = spec.key, spec.value
    else:
        mask, source = _file_root_cells(ctx)
        if mask is None:
            raise RootError("the file stores no root (no uns['iroot'] and no boolean obs 'root' column)")
        cell = int(np.flatnonzero(mask)[0]) if mask.sum() == 1 else _extreme_cell(ctx, mask)
        how = f"the root stored in the file: {source}"
        key, label = None, None

    if key is None:
        key, label = LEIDEN_KEY, str(ctx.labels(LEIDEN_KEY)[cell])
    return ResolvedRoot(
        spec=format_root(spec),
        cell_index=cell,
        obs_name=str(ctx.obs_names[cell]),
        cluster_key=key,
        cluster_label=label,
        resolved_by=how,
        warnings=tuple(_component_warnings(ctx, cell)),
    )


def interior_root_warning(
    ctx: InferenceContext, root: ResolvedRoot, pseudotime: npt.NDArray[np.floating]
) -> Notice | None:
    """Warn when much of the root's own cluster lies far along: the root may be mid-trajectory."""
    leiden = ctx.leiden
    same = leiden == leiden[root.cell_index]
    pt = pseudotime[same]
    pt = pt[np.isfinite(pt)]
    if pt.size == 0:
        return None
    share_far = float(np.mean(pt > 0.3))
    if share_far > 0.3:
        return Notice(
            "root_interior",
            "warning",
            f"{share_far:.0%} of the root's cluster lies past pseudotime 0.3. The root may be in "
            "the middle of the trajectory rather than at its start.",
        )
    return None


# ---- suggestions ----------------------------------------------------------------------


@dataclass(frozen=True)
class RootSuggestion:
    spec: ClusterRoot | FileRoot
    source: Literal["file", "time_column", "gene_count"]
    score: float
    evidence: str
    is_extremity: bool


def _paga_degree(ctx: InferenceContext) -> npt.NDArray[np.int64]:
    tree = sparse.csr_matrix(ctx.artifacts["paga_tree"])
    sym = ((tree + tree.T) > 0).astype(np.int64)
    return np.asarray(sym.sum(axis=1)).ravel()


def _is_extremity(ctx: InferenceContext, members: npt.NDArray[np.bool_]) -> bool:
    leiden = ctx.leiden[members]
    if leiden.size == 0:
        return False
    main = int(np.bincount(leiden).argmax())
    return bool(_paga_degree(ctx)[main] <= 1)


def _parse_times(values: npt.NDArray[np.str_]) -> npt.NDArray[np.float64] | None:
    out = np.full(values.size, np.nan)
    for i, v in enumerate(values):
        m = _NUMBER.search(v)
        if m:
            out[i] = float(m.group())
    if np.mean(np.isfinite(out)) < 0.8 or np.unique(out[np.isfinite(out)]).size < 2:
        return None
    return out


def _time_suggestion(ctx: InferenceContext, key: str) -> RootSuggestion | None:
    clusters = ctx.labels(key)
    for col in ctx.obs.columns:
        if col == key or not TIME_COLUMN.search(str(col)):
            continue
        times = _parse_times(ctx.obs[col].astype(str).to_numpy())
        if times is None:
            continue
        earliest = times == np.nanmin(times)
        n_total, n_early = ctx.n_obs, int(earliest.sum())
        best: tuple[float, str, float] | None = None
        for c in np.unique(clusters):
            members = clusters == c
            k = int((members & earliest).sum())
            n = int(members.sum())
            if k == 0 or k / n <= n_early / n_total:
                continue
            p = float(stats.hypergeom.sf(k - 1, n_total, n_early, n))
            if best is None or p < best[0]:
                best = (p, str(c), k / n)
        if best is None:
            continue
        p, c, share = best
        earliest_value = ctx.obs[col].astype(str).to_numpy()[earliest][0]
        return RootSuggestion(
            spec=ClusterRoot(key, c),
            source="time_column",
            score=min(1.0, -math.log10(max(p, 1e-300)) / 50.0),
            evidence=(
                f"{share:.0%} of its cells come from the earliest {col} ({earliest_value}); "
                f"enrichment p = {p:.1e}"
            ),
            is_extremity=_is_extremity(ctx, clusters == c),
        )
    return None


def _gene_count_suggestion(ctx: InferenceContext, key: str) -> RootSuggestion | None:
    if ctx.counts is None:
        return None
    detected = np.diff(ctx.counts.indptr).astype(np.float64)
    clusters = ctx.labels(key)
    medians = {str(c): float(np.median(detected[clusters == c])) for c in np.unique(clusters)}
    best = max(medians, key=lambda c: (medians[c], c))
    overall = float(np.median(detected))
    return RootSuggestion(
        spec=ClusterRoot(key, best),
        source="gene_count",
        score=0.3,
        evidence=(
            f"median {medians[best]:,.0f} genes detected per cell vs {overall:,.0f} overall. "
            "Heuristic (CytoTRACE, Gulati et al. 2020): less differentiated cells tend to express "
            "more genes. Unreliable for cycling or doublet-rich clusters."
        ),
        is_extremity=_is_extremity(ctx, clusters == best),
    )


def suggest_roots(ctx: InferenceContext, key: str | None = None) -> list[RootSuggestion]:
    """Candidate roots, best first. Each carries its evidence; none is ever applied."""
    key = key or ctx.label_key
    out: list[RootSuggestion] = []
    mask, source = _file_root_cells(ctx)
    if mask is not None:
        out.append(
            RootSuggestion(
                FileRoot(), "file", 1.0, f"the file stores a root: {source}", _is_extremity(ctx, mask)
            )
        )
    time_s = _time_suggestion(ctx, key)
    if time_s:
        out.append(time_s)
    gene_s = _gene_count_suggestion(ctx, key)
    if gene_s:
        out.append(gene_s)
    return out
