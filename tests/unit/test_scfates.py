import builtins
from pathlib import Path

import numpy as np
import pytest

from lumencell.engine.methods.scfates_adapter import INSTALL_HINT, ScFates
from lumencell.engine.projection import project
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from lumencell.engine.store import RunStore
from tests.conftest import Prepared

pytest.importorskip("scFates")

ROOT = ClusterRoot("true_label", "trunk_early")


@pytest.fixture(scope="module")
def result(bifurcation: Prepared) -> TrajectoryResult:
    return run_method(bifurcation.ctx, "scfates_ppt", ROOT)


def test_result_is_valid_with_a_fitted_tree(result: TrajectoryResult) -> None:
    result.validate()
    assert result.topology_kind == "tree"
    assert result.inference_rep == "X_diffmap_ms"
    assert result.pseudotime[result.root.cell_index] == 0.0
    assert result.curve_nodes is not None and result.curve_paths is not None
    assert len(result.curve_paths) == len(result.lineages) == 2
    kinds = sorted(m.kind for m in result.milestones)
    assert kinds.count("root") == 1 and kinds.count("tip") == 2 and kinds.count("fork") >= 1


def test_pseudotime_follows_the_truth(result: TrajectoryResult, bifurcation: Prepared) -> None:
    t = bifurcation.ds.obs["true_pseudotime"].to_numpy()
    ok = np.isfinite(result.pseudotime)
    rank = lambda v: np.argsort(np.argsort(v))  # noqa: E731
    assert np.corrcoef(rank(result.pseudotime[ok]), rank(t[ok]))[0, 1] > 0.95


def test_drawn_curve_is_the_fitted_tree(result: TrajectoryResult, bifurcation: Prepared) -> None:
    proj = project(result, bifurcation.prep.artifacts["umap"], "umap")
    assert proj.curve_source == "fitted"
    for curve, path in zip(proj.curves, result.curve_paths or [], strict=True):
        assert len(curve.points) == len(path)
        assert curve.fidelity is not None and curve.fidelity > 0.8


def test_fitted_tree_survives_the_store(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    store.save("s", result)
    back = store.load("s")
    back.validate()
    assert back.curve_paths == result.curve_paths
    assert back.curve_nodes is not None and back.curve_nodes.shape == result.curve_nodes.shape  # type: ignore[union-attr]


def test_named_terminal_keeps_one_lineage(bifurcation: Prepared) -> None:
    r = run_method(
        bifurcation.ctx, "scfates_ppt", ROOT, terminals=["A_late"], params={"cluster_key": "true_label"}
    )
    assert len(r.lineages) == 1
    assert r.lineages[0].terminal_label == "A_late"
    assert r.n_unassigned > 0


def test_missing_install_gives_a_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    real_import = builtins.__import__

    def fake(name: str, *args: object, **kwargs: object) -> object:
        if name == "scFates":
            raise ImportError("no scFates")
        return real_import(name, *args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(builtins, "__import__", fake)
    avail = ScFates().available()
    assert not avail.ok and avail.install_hint == INSTALL_HINT
