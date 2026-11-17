import inspect
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from lumencell.bench.gates import load_thresholds
from lumencell.bench.synth import generate
from lumencell.bench.synth_controls import (
    CONTROL_STREAM,
    GROUP_LFC_SCALE,
    KINDS,
    LFC_SCALE,
    LOOP_PROGRAMME_SHARE,
    PROGRAMME_SHARE,
    ControlKind,
    control_root_cell,
    generate_control,
)
from lumencell.data.loader import load
from lumencell.prep.stages import STAGES, PrepParams, run_stage

N_CELLS, N_GENES = 1000, 500


def _make(kind: ControlKind, seed: int = 0) -> ad.AnnData:
    return generate_control(kind, n_cells=N_CELLS, n_genes=N_GENES, seed=seed)


def _knn_edges(adata: ad.AnnData, tmp_path: Path) -> tuple[np.ndarray, np.ndarray]:
    path = tmp_path / "c.h5ad"
    adata.write_h5ad(path)
    ds, params = load(path), PrepParams()
    artifacts: dict[str, Any] = {}
    # Only as far as the kNN graph: full prep (PAGA) fails on two clusters with no edge between them.
    for stage in STAGES[: [s.name for s in STAGES].index("leiden")]:
        artifacts.update(run_stage(stage, ds, artifacts, params)[0])
    coo = sparse.coo_matrix(artifacts["connectivities"])
    keep = coo.row < coo.col
    return coo.row[keep], coo.col[keep]


@pytest.mark.parametrize("kind", KINDS)
def test_controls_are_deterministic_and_differ_by_seed(kind: ControlKind) -> None:
    a, b, c = _make(kind, 0), _make(kind, 0), _make(kind, 1)
    assert (a.X != b.X).nnz == 0
    assert a.obs.equals(b.obs)
    assert a.uns["lumencell_control"] == b.uns["lumencell_control"]
    assert (a.X != c.X).nnz > 0


def test_controls_are_counts_with_synthetic_gene_names() -> None:
    for kind in KINDS:
        a = _make(kind)
        assert np.array_equal(a.X.data, np.round(a.X.data))
        assert a.X.data.min() >= 0
        assert a.var_names[0] == "G00000"
        assert a.obs_names[0] == "cell000000"
        assert "lumencell_truth" not in a.uns
        assert not a.obsm
        assert set(a.uns["lumencell_control"]) == {"kind", "seed", "preset", "root_cell", "programme_genes"}
    with pytest.raises(ValueError):
        generate("blob")
    with pytest.raises(ValueError):
        generate_control("tree", n_cells=10, n_genes=10)  # type: ignore[arg-type]


def test_blob_has_no_programme() -> None:
    a = _make("blob")
    assert a.uns["lumencell_control"]["programme_genes"] == 0
    assert not [c for c in a.obs.columns if str(c).startswith("control_")]


def test_two_groups_are_unconnected(tmp_path: Path) -> None:
    a = _make("two_groups")
    row, col = _knn_edges(a, tmp_path)
    g = a.obs["control_group"].to_numpy()
    assert 0.2 < g.mean() < 0.8
    assert (g[row] != g[col]).mean() < 0.001


def test_the_loop_closes(tmp_path: Path) -> None:
    a = _make("loop")
    row, col = _knn_edges(a, tmp_path)
    theta = a.obs["control_angle"].to_numpy()
    gap = np.abs(theta[row] - theta[col])
    gap = np.minimum(gap, 2 * np.pi - gap)
    assert np.median(gap) < 0.3
    low, high = theta < 0.2, theta > 2 * np.pi - 0.2
    assert (low[row] & high[col]).any() or (high[row] & low[col]).any()


def test_the_sheet_resolves_both_gradients(tmp_path: Path) -> None:
    a = _make("sheet")
    row, col = _knn_edges(a, tmp_path)
    x, y = a.obs["control_x"].to_numpy(), a.obs["control_y"].to_numpy()
    assert np.median(np.abs(x[row] - x[col])) < 0.15
    assert np.median(np.abs(y[row] - y[col])) < 0.15


@pytest.mark.parametrize("kind", KINDS)
def test_root_cells_sit_where_declared(kind: ControlKind) -> None:
    a = _make(kind)
    root = control_root_cell(a)
    assert root in a.obs_names
    i = a.obs_names.get_loc(root)
    if kind == "blob":
        assert root == "cell000000"
    elif kind == "two_groups":
        g = a.obs["control_group"].to_numpy()
        assert g[i] == 0 and i == int(np.flatnonzero(g == 0)[0])
    elif kind == "loop":
        t = a.obs["control_angle"].to_numpy()
        d = np.minimum(t, 2 * np.pi - t)
        assert d[i] == d.min()
    else:
        x, y = a.obs["control_x"].to_numpy(), a.obs["control_y"].to_numpy()
        r = x * x + y * y
        assert r[i] == r.min()


def test_constants_equal_the_thresholds_file() -> None:
    t: dict[str, Any] = load_thresholds()["synthetic"]
    c = t["controls"]
    assert tuple(c["kinds"]) == KINDS
    assert c["programme_share"] == PROGRAMME_SHARE
    assert c["loop_programme_share"] == LOOP_PROGRAMME_SHARE
    assert c["lfc_scale"] == LFC_SCALE
    assert c["group_lfc_scale"] == GROUP_LFC_SCALE
    assert (c["n_cells"], c["n_genes"]) == (3000, 2000)
    assert t["preset"] == "realistic"
    d = inspect.signature(generate_control).parameters
    assert (d["n_cells"].default, d["n_genes"].default, d["preset"].default) == (3000, 2000, t["preset"])
    assert CONTROL_STREAM == 7100502
