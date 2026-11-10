"""The PAGA cluster tree, oriented from the root. Shared by methods that need it.

PAGA (Wolf et al. 2019) summarises the kNN graph as connectivities between
clusters; its maximum-connectivity spanning tree is a coarse picture of the
trajectory. Here weak tree edges are cut and the tree is oriented by
breadth-first search from the root cell's cluster.

paga_dpt uses the whole tree for its lineages. Palantir only uses the leaves,
as suggested terminal states when the user names none.
"""

from collections import deque
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd
import scanpy as sc
from scipy import sparse

from lumencell.engine.context import LEIDEN_KEY, InferenceContext
from lumencell.prep.stages import neighbors_view


@dataclass(frozen=True)
class ClusterTree:
    cats: list[str]  # group labels, indexed by group id
    group_of: npt.NDArray[np.int64]  # group id per cell
    connectivities: npt.NDArray[np.float64]  # symmetric PAGA connectivities
    root_group: int
    parent: dict[int, int]  # -1 for the root group; unreachable groups are absent
    order: list[int]  # breadth-first order from the root group

    @property
    def children(self) -> dict[int, list[int]]:
        out: dict[int, list[int]] = {g: [] for g in self.order}
        for g in self.order[1:]:
            out[self.parent[g]].append(g)
        return out

    @property
    def unreachable(self) -> list[int]:
        return [g for g in range(len(self.cats)) if g not in self.parent]

    @property
    def leaves(self) -> list[int]:
        kids = self.children
        leaves = [g for g in self.order if not kids[g] and g != self.root_group]
        return leaves or [self.root_group]

    def path_to(self, group: int) -> list[int]:
        path = [group]
        while self.parent[path[-1]] != -1:
            path.append(self.parent[path[-1]])
        return path[::-1]


def _paga(
    ctx: InferenceContext, key: str
) -> tuple[list[str], npt.NDArray[np.int64], sparse.csr_matrix, sparse.csr_matrix]:
    labels = ctx.labels(key)
    if key == LEIDEN_KEY:
        cats = [str(i) for i in range(int(ctx.leiden.max()) + 1)]
        return (
            cats,
            ctx.leiden.astype(np.int64),
            sparse.csr_matrix(ctx.artifacts["paga_connectivities"]),
            sparse.csr_matrix(ctx.artifacts["paga_tree"]),
        )
    cats = sorted(set(labels.tolist()))
    view = neighbors_view(ctx.n_obs, ctx.artifacts, ctx.prep_params)
    view.obs["groups"] = pd.Categorical(labels, categories=cats)
    sc.tl.paga(view, groups="groups")
    lookup = {c: i for i, c in enumerate(cats)}
    return (
        cats,
        np.array([lookup[v] for v in labels], dtype=np.int64),
        sparse.csr_matrix(view.uns["paga"]["connectivities"]),
        sparse.csr_matrix(view.uns["paga"]["connectivities_tree"]),
    )


def cluster_tree(
    ctx: InferenceContext, root_cell: int, key: str = LEIDEN_KEY, min_connectivity: float = 0.05
) -> ClusterTree:
    cats, group_of, conn, tree = _paga(ctx, key)
    sym = conn.maximum(conn.T).toarray()
    adj = ((tree + tree.T) > 0).toarray() & (sym >= min_connectivity)
    root_group = int(group_of[root_cell])
    parent = {root_group: -1}
    order = [root_group]
    queue = deque([root_group])
    while queue:
        g = queue.popleft()
        for nb in (int(x) for x in np.flatnonzero(adj[g])):
            if nb not in parent:
                parent[nb] = g
                order.append(nb)
                queue.append(nb)
    return ClusterTree(cats, group_of, sym, root_group, parent, order)


def leaf_terminal_cells(ctx: InferenceContext, tree: ClusterTree, root_cell: int) -> list[int]:
    """For each leaf group, the cell lying farthest from the root in diffusion space."""
    space = ctx.diffusion_space()
    dist = np.linalg.norm(space - space[root_cell], axis=1)
    out = []
    for g in tree.leaves:
        if g == tree.root_group:
            continue
        cells = np.flatnonzero(tree.group_of == g)
        out.append(int(cells[np.argmax(dist[cells])]))
    return out
