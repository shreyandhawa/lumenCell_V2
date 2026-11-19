"""A disconnected kNN graph fails early with a message a user can act on (phase 2 follow-up)."""

import numpy as np
import pytest
from scipy import sparse

from lumencell.bench.synth_controls import control_root_cell, generate_control
from lumencell.data.loader import load
from lumencell.engine.context import InferenceContext
from lumencell.engine.graph import DisconnectedGraphError, require_connected
from lumencell.engine.registry import run_method
from lumencell.engine.roots import CellRoot
from lumencell.prep.pipeline import run_pipeline
from lumencell.prep.stages import PrepParams


def test_require_connected_counts_parts_and_names_no_cells() -> None:
    blocks = sparse.block_diag([np.ones((4, 4)), np.ones((2, 2)), np.ones((3, 3))])
    with pytest.raises(DisconnectedGraphError) as exc:
        require_connected(blocks, root_cell=5, what="kNN graph")
    msg = str(exc.value)
    assert "3 disconnected parts" in msg and "2 cells in the part holding the root" in msg
    assert "4 in the largest other part" in msg and "Subset" in msg and "neighbours" in msg


def test_require_connected_accepts_one_part() -> None:
    require_connected(sparse.csr_matrix(np.ones((5, 5))), root_cell=0, what="kNN graph")


@pytest.fixture(scope="module")
def split(tmp_path_factory: pytest.TempPathFactory) -> tuple[InferenceContext, str]:
    adata = generate_control("two_groups", n_cells=1000, n_genes=500, seed=0)
    path = tmp_path_factory.mktemp("split") / "two_groups.h5ad"
    adata.write_h5ad(path)
    ds, params = load(path), PrepParams()
    return InferenceContext.from_prep(ds, run_pipeline(ds, params), params), control_root_cell(adata)


@pytest.mark.parametrize("method", ["palantir", "scfates_ppt"])
def test_methods_that_cannot_cross_parts_fail_with_a_clear_error(
    split: tuple[InferenceContext, str], method: str
) -> None:
    pytest.importorskip("palantir" if method == "palantir" else "scFates")
    ctx, root = split
    with pytest.raises(
        DisconnectedGraphError, match=r"2 disconnected parts: [\d,]+ cells in the part"
    ) as exc:
        run_method(ctx, method, CellRoot(root))
    assert root not in str(exc.value)


def test_paga_dpt_still_runs_on_the_split_graph(split: tuple[InferenceContext, str]) -> None:
    ctx, root = split
    run_method(ctx, "paga_dpt", CellRoot(root)).validate()
