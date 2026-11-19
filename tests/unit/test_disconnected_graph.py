"""A disconnected kNN graph fails early with a message a user can act on (phase 2 follow-up)."""

from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from typer.testing import CliRunner

from lumencell.bench.synth_controls import control_root_cell, generate_control
from lumencell.cli import app
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


def test_error_is_a_value_error_so_the_cli_and_api_report_it_cleanly() -> None:
    assert issubclass(DisconnectedGraphError, ValueError)


def test_stored_explicit_zero_is_not_an_edge_and_the_input_is_not_mutated() -> None:
    rows, cols = [0, 1, 2, 3, 1, 2], [1, 0, 3, 2, 2, 1]
    data = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0]  # the last two are stored zeros
    csr = sparse.csr_matrix((data, (rows, cols)), shape=(4, 4))
    assert csr.nnz == 6
    with pytest.raises(DisconnectedGraphError, match="2 disconnected parts"):
        require_connected(csr, root_cell=0, what="kNN graph")
    assert csr.nnz == 6


def test_connectivity_is_weak_so_a_one_way_chain_is_one_part() -> None:
    chain = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 2])), shape=(3, 3))
    require_connected(chain, root_cell=0, what="kNN graph")
    with_isolated = sparse.csr_matrix(([1.0, 1.0], ([0, 1], [1, 2])), shape=(4, 4))
    with pytest.raises(DisconnectedGraphError, match="2 disconnected parts"):
        require_connected(with_isolated, root_cell=0, what="kNN graph")


def test_largest_other_part_excludes_the_roots_part() -> None:
    blocks = sparse.block_diag([np.ones((4, 4)), np.ones((2, 2)), np.ones((3, 3))])
    with pytest.raises(DisconnectedGraphError) as exc:
        require_connected(blocks, root_cell=0, what="kNN graph")
    assert "4 cells in the part holding the root and 3 in the largest other part" in str(exc.value)


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


@pytest.mark.parametrize("method", ["palantir", "scfates_ppt"])
def test_the_message_reports_the_part_of_the_run_root(
    split: tuple[InferenceContext, str], method: str
) -> None:
    """Rooted in the part that does NOT hold cell 0, so a check passed index 0 reports the wrong part."""
    pytest.importorskip("palantir" if method == "palantir" else "scFates")
    ctx, _ = split
    _, labels = connected_components(sparse.csr_matrix(ctx.artifacts["connectivities"]), directed=False)
    sizes = np.bincount(labels)
    other = int(np.flatnonzero(np.arange(len(sizes)) != labels[0])[0])
    assert sizes[other] != sizes[labels[0]]  # the two parts differ in size, so the count tells them apart
    root = str(ctx.obs_names[int(np.flatnonzero(labels == other)[0])])
    with pytest.raises(DisconnectedGraphError) as exc:
        run_method(ctx, method, CellRoot(root))
    assert f"{int(sizes[other]):,} cells in the part holding the root" in str(exc.value)


def test_cli_run_on_a_disconnected_graph_is_a_clean_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("scFates")
    monkeypatch.setenv("LUMENCELL_CACHE_DIR", str(tmp_path / "cache"))
    adata = generate_control("two_groups", n_cells=1000, n_genes=500, seed=0)
    path = tmp_path / "two_groups.h5ad"
    adata.write_h5ad(path)
    result = CliRunner().invoke(
        app, ["run", str(path), "--method", "scfates_ppt", "--root", f"cell:{control_root_cell(adata)}"]
    )
    assert result.exit_code == 2
    assert "disconnected parts" in result.output
    assert isinstance(result.exception, SystemExit), "a message, not a traceback"
