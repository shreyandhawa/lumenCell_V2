from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from tests.conftest import Prepared

from lumencell.engine.context import LEIDEN_KEY, InferenceContext
from lumencell.engine.roots import (
    CellRoot,
    ClusterRoot,
    FileRoot,
    RootError,
    format_root,
    interior_root_warning,
    parse_root,
    resolve_root,
    suggest_roots,
)

# ---- parsing ------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("text", "spec"),
    [
        ("cell:cell000012", CellRoot("cell000012")),
        ("cluster:clusters=Ductal", ClusterRoot("clusters", "Ductal")),
        ("cluster:stage=E7.5=early", ClusterRoot("stage", "E7.5=early")),
        ("file", FileRoot()),
    ],
)
def test_parse_and_format_round_trip(text: str, spec: object) -> None:
    assert parse_root(text) == spec
    assert format_root(parse_root(text)) == text


@pytest.mark.parametrize("text", ["", "cell:", "cluster:noequals", "cluster:=x", "tip", "auto"])
def test_unparseable_roots(text: str) -> None:
    with pytest.raises(RootError, match="cannot parse"):
        parse_root(text)


# ---- resolution ---------------------------------------------------------------------


def test_cell_root(bifurcation: Prepared) -> None:
    name = str(bifurcation.ctx.obs_names[42])
    r = resolve_root(bifurcation.ctx, CellRoot(name))
    assert (r.cell_index, r.obs_name) == (42, name)
    assert r.cluster_key == LEIDEN_KEY


def test_unknown_cell(bifurcation: Prepared) -> None:
    with pytest.raises(RootError, match="no cell is named"):
        resolve_root(bifurcation.ctx, CellRoot("nope"))


def test_cluster_root_picks_the_clusters_outer_end(bifurcation: Prepared) -> None:
    # true_label is only used here to name a known start cluster; it is the first
    # half of the trunk, true pseudotime 0..0.2.
    r = resolve_root(bifurcation.ctx, ClusterRoot("true_label", "trunk_early"))
    t = bifurcation.ds.obs["true_pseudotime"].to_numpy()
    assert t[r.cell_index] < 0.05, "the farthest cell from the rest should be at the very start"
    assert r.cluster_label == "trunk_early"
    assert "farthest" in r.resolved_by


def test_cluster_root_is_deterministic(bifurcation: Prepared) -> None:
    spec = ClusterRoot("true_label", "trunk_early")
    assert resolve_root(bifurcation.ctx, spec) == resolve_root(bifurcation.ctx, spec)


def test_unknown_cluster_value_lists_the_options(bifurcation: Prepared) -> None:
    with pytest.raises(RootError, match=r"values are: .*trunk_early"):
        resolve_root(bifurcation.ctx, ClusterRoot("true_label", "nope"))


def test_unknown_column(bifurcation: Prepared) -> None:
    with pytest.raises(RootError, match="no obs column"):
        resolve_root(bifurcation.ctx, ClusterRoot("nope", "x"))


def test_file_root(bifurcation: Prepared) -> None:
    r = resolve_root(bifurcation.ctx, FileRoot())
    assert r.cell_index == bifurcation.ds.uns_root
    assert "uns['iroot']" in r.resolved_by


def test_file_root_from_obs_column(bifurcation: Prepared) -> None:
    obs = bifurcation.ctx.obs.copy()
    obs["is_root"] = False
    obs.iloc[7, obs.columns.get_loc("is_root")] = True
    ctx = replace(bifurcation.ctx, obs=obs, uns_root=None)
    assert resolve_root(ctx, FileRoot()).cell_index == 7


def test_missing_file_root(bifurcation: Prepared) -> None:
    ctx = replace(bifurcation.ctx, uns_root=None)
    with pytest.raises(RootError, match="stores no root"):
        resolve_root(ctx, FileRoot())


def test_connected_root_has_no_warnings(bifurcation: Prepared) -> None:
    assert resolve_root(bifurcation.ctx, ClusterRoot("true_label", "trunk_early")).warnings == ()


def _disconnect(ctx: InferenceContext, cells: np.ndarray) -> InferenceContext:
    conn = sparse.lil_matrix(ctx.artifacts["connectivities"])
    conn[cells, :] = 0
    conn[:, cells] = 0
    art = dict(ctx.artifacts)
    art["connectivities"] = sparse.csr_matrix(conn)
    return replace(ctx, artifacts=art)


def test_disconnected_cells_are_reported(bifurcation: Prepared) -> None:
    ctx = _disconnect(bifurcation.ctx, np.arange(1000, 1020))
    r = resolve_root(ctx, CellRoot(str(ctx.obs_names[0])))
    assert [w.code for w in r.warnings] == ["disconnected_cells"]
    assert "20 cells" in r.warnings[0].message


def test_root_in_a_small_group_is_a_warning(bifurcation: Prepared) -> None:
    ctx = _disconnect(bifurcation.ctx, np.arange(5, 1500))
    r = resolve_root(ctx, CellRoot(str(ctx.obs_names[0])))
    assert r.warnings[0].code == "root_in_small_group"


def test_interior_root_warning(bifurcation: Prepared) -> None:
    ctx = bifurcation.ctx
    root = resolve_root(ctx, CellRoot(str(ctx.obs_names[0])))
    same = ctx.leiden == ctx.leiden[0]
    pt = np.where(same, 0.8, 0.1)
    pt[0] = 0.0
    warning = interior_root_warning(ctx, root, pt)
    assert warning is not None and warning.code == "root_interior"
    assert interior_root_warning(ctx, root, np.where(same, 0.05, 0.6)) is None


# ---- suggestions ----------------------------------------------------------------------


def test_file_root_is_suggested_first_but_not_applied(bifurcation: Prepared) -> None:
    before = bifurcation.ctx.obs.copy()
    suggestions = suggest_roots(bifurcation.ctx)
    assert suggestions[0].source == "file"
    assert isinstance(suggestions[0].spec, FileRoot)
    pd.testing.assert_frame_equal(bifurcation.ctx.obs, before)


def test_gene_count_heuristic_is_labelled_as_one(bifurcation: Prepared) -> None:
    s = [x for x in suggest_roots(bifurcation.ctx) if x.source == "gene_count"]
    assert len(s) == 1
    assert "Heuristic" in s[0].evidence and "Unreliable" in s[0].evidence


def test_time_column_suggests_the_earliest_enriched_cluster(bifurcation: Prepared) -> None:
    ctx = bifurcation.ctx
    t = bifurcation.ds.obs["true_pseudotime"].to_numpy()
    obs = ctx.obs.copy()
    obs["day"] = pd.Categorical([f"day {int(v * 4)}" for v in np.minimum(t, 0.999)])
    ctx = replace(ctx, obs=obs, uns_root=None)
    s = next(x for x in suggest_roots(ctx) if x.source == "time_column")
    assert isinstance(s.spec, ClusterRoot) and s.spec.key == LEIDEN_KEY
    assert "earliest day (day 0)" in s.evidence
    resolved = resolve_root(ctx, s.spec)
    assert t[resolved.cell_index] < 0.15
