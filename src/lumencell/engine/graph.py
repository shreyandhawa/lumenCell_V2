"""Refuse a neighbour graph that splits into parts a method cannot cross."""

import numpy as np
import numpy.typing as npt
from scipy import sparse
from scipy.sparse.csgraph import connected_components


class DisconnectedGraphError(ValueError):
    """The neighbour graph a method needs has more than one connected part."""


def require_connected(graph: sparse.spmatrix | npt.NDArray[np.float64], root_cell: int, what: str) -> None:
    """Raise ``DisconnectedGraphError`` unless ``graph`` is one connected part.

    The message names counts only (parts, cells in the root's part and in the largest other
    part), never cell names, since it reaches the API, the CLI and job logs.
    """
    # A stored explicit zero is not an edge (ruling R19): drop zeros on a copy, never the caller's matrix.
    edges = sparse.csr_matrix(graph, copy=True)
    edges.eliminate_zeros()
    # Weak connectivity (directed=False) is intended: one-way kNN edges still join cells.
    n_parts, labels = connected_components(edges, directed=False)
    if n_parts <= 1:
        return
    sizes = np.bincount(labels, minlength=n_parts)
    root_part = int(labels[root_cell])
    others = np.delete(sizes, root_part)
    raise DisconnectedGraphError(
        f"The {what} has {n_parts} disconnected parts: {int(sizes[root_part]):,} cells in the part "
        f"holding the root and {int(others.max()):,} in the largest other part. This method cannot "
        "order cells across parts. Subset the data to one part (for example the cell type you are "
        "studying), or raise the number of neighbours so the parts join."
    )
