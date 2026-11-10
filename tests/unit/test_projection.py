import numpy as np
import pytest

from lumencell.engine.projection import N_BINS, _arc_position, project
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from tests.conftest import Prepared


@pytest.fixture(scope="module")
def result(bifurcation: Prepared) -> TrajectoryResult:
    return run_method(bifurcation.ctx, "paga_dpt", ClusterRoot("true_label", "trunk_early"))


def test_arc_position_on_a_straight_line() -> None:
    curve = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    cells = np.array([[0.0, 0.5], [1.0, -0.2], [2.0, 0.1], [3.0, 0.0]])
    np.testing.assert_allclose(_arc_position(cells, curve), [0.0, 0.5, 1.0, 1.0])


def test_milestones_and_curves_on_umap(result: TrajectoryResult, bifurcation: Prepared) -> None:
    umap = bifurcation.prep.artifacts["umap"]
    proj = project(result, umap, "umap")
    assert set(proj.milestones) == {m.id for m in result.milestones}
    assert len(proj.curves) == len(result.lineages)
    root_id = next(m.id for m in result.milestones if m.kind == "root")
    for curve in proj.curves:
        assert len(curve.points) == N_BINS + 1, "root milestone + one point per bin"
        assert curve.points[0] == pytest.approx(proj.milestones[root_id])
        assert curve.fidelity is not None and curve.fidelity > 0.8
    assert proj.curve_source == "binned"


def test_curves_on_three_diffusion_components(result: TrajectoryResult, bifurcation: Prepared) -> None:
    dc123 = np.asarray(bifurcation.prep.artifacts["diffmap"])[:, 1:4]
    proj = project(result, dc123, "dc123")
    root_id = next(m.id for m in result.milestones if m.kind == "root")
    assert all(len(p) == 3 for p in proj.milestones.values())
    for curve in proj.curves:
        assert all(len(p) == 3 for p in curve.points)
        assert curve.points[0] == pytest.approx(proj.milestones[root_id])
        assert curve.fidelity is not None and curve.fidelity > 0.8


def test_scrambled_embedding_triggers_the_fidelity_warning(result: TrajectoryResult) -> None:
    noise = np.random.default_rng(0).random((result.n_obs, 2))
    proj = project(result, noise, "noise")
    assert any(w.code == "embedding_distorts_lineage" for w in proj.warnings)
    assert all(c.fidelity is None or c.fidelity < 0.8 for c in proj.curves)
