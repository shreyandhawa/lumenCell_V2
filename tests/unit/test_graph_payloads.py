"""The cluster and neighbour graphs drawn over an embedding."""

import numpy as np
from scipy import sparse

from lumencell.api.dataset_view import knn_edges, paga_layout


def test_paga_nodes_sit_at_the_mean_of_their_cells() -> None:
    coords = np.array([[0, 0], [2, 2], [np.nan, 5], [10, 0], [np.inf, 1]], dtype=np.float32)
    labels = np.array([0, 0, 0, 1, 2])
    paga = sparse.csr_matrix(np.array([[0, 0.5, 0.2], [0.5, 0, 0.9], [0.2, 0.9, 0]], dtype=np.float32))
    nodes, edges = paga_layout(coords, labels, paga)
    # Cluster 0's NaN cell is left out of the mean but counted; cluster 2 has no finite cell, so it
    # and its edges are not drawn.
    assert [(n.id, n.x, n.y, n.z, n.n_cells) for n in nodes] == [
        (0, 1.0, 1.0, None, 3),
        (1, 10.0, 0.0, None, 1),
    ]
    assert [(e.source, e.target, e.weight) for e in edges] == [(0, 1, 0.5)]


def test_paga_on_a_3d_view_sets_z_and_a_graph_without_edges_has_none() -> None:
    coords = np.array([[0, 0, 1], [1, 1, 3]], dtype=np.float32)
    nodes, edges = paga_layout(coords, np.array([0, 1]), sparse.csr_matrix((2, 2), dtype=np.float32))
    assert [n.z for n in nodes] == [1.0, 3.0] and edges == []


def test_each_neighbour_edge_is_sent_once() -> None:
    dense = np.zeros((4, 4), dtype=np.float32)
    dense[0, 1] = dense[1, 0] = 0.5
    dense[2, 3] = 0.7  # only one direction stored: still one edge
    dense[3, 3] = 1.0  # a self-loop is not an edge
    src, dst, w, total = knn_edges(sparse.csr_matrix(dense), max_edges=10, seed=0)
    assert (src.tolist(), dst.tolist(), w.tolist(), total) == ([0, 2], [1, 3], [0.5, np.float32(0.7)], 2)
    assert [a.dtype.name for a in (src, dst, w)] == ["uint32", "uint32", "float32"]


def test_above_the_cap_a_seeded_sample_is_sent() -> None:
    conn = sparse.random(300, 300, density=0.1, random_state=1, format="csr", dtype=np.float32)
    all_src, all_dst, _, total = knn_edges(conn, max_edges=10**6, seed=0)
    src, dst, _, again = knn_edges(conn, max_edges=50, seed=0)
    assert again == total > 50 and len(src) == 50
    assert {(a, b) for a, b in zip(src, dst, strict=True)} <= set(zip(all_src, all_dst, strict=True))
    same, _, _, _ = knn_edges(conn, max_edges=50, seed=0)
    other, _, _, _ = knn_edges(conn, max_edges=50, seed=1)
    assert np.array_equal(src, same) and not np.array_equal(src, other)
