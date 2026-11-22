import numpy as np
import pytest

from lumencell.engine.projection import N_BINS, _arc_position, diffusion_axes, project
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


def test_scrambled_embedding_triggers_the_fidelity_warning(result: TrajectoryResult) -> None:
    noise = np.random.default_rng(0).random((result.n_obs, 2))
    proj = project(result, noise, "noise")
    assert any(w.code == "embedding_distorts_lineage" for w in proj.warnings)
    assert all(c.fidelity is None or c.fidelity < 0.8 for c in proj.curves)


@pytest.mark.parametrize(
    ("name", "axes"),
    [
        ("dc12", (1, 2)),
        ("dc123", (1, 2, 3)),
        ("dc:1,2,4", (1, 2, 4)),
        ("dc:3,1", (3, 1)),
        ("dc:12,13", (12, 13)),
        ("umap", None),
        ("dc:1", None),
        ("dc:1,2,3,4", None),
        ("dc:0,1", None),
        ("dc:01,2", None),
        ("dc:1, 2", None),
        ("dc:-1,2", None),
        ("dc:1,2,", None),
    ],
)
def test_diffusion_view_names(name: str, axes: tuple[int, ...] | None) -> None:
    # Repeats and components past the last are left to the caller, which knows how many there are.
    assert diffusion_axes(name) == axes


def test_the_fidelity_warning_names_the_diffusion_axes(result: TrajectoryResult) -> None:
    noise = np.random.default_rng(0).random((result.n_obs, 3))
    messages = [w.message for w in project(result, noise, "dc:1,2,4").warnings]
    assert messages and all(
        "On DC1, DC2 and DC4," in m and "Try other diffusion components." in m for m in messages
    )
    assert all("On the noise embedding," in w.message for w in project(result, noise, "noise").warnings)
