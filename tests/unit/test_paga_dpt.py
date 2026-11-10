import numpy as np
import pytest
from pydantic import ValidationError

from lumencell.engine.context import LEIDEN_KEY
from lumencell.engine.registry import METHODS, UnknownMethodError, get_method, run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from tests.conftest import Prepared

ROOT = ClusterRoot("true_label", "trunk_early")


@pytest.fixture(scope="module")
def result(bifurcation: Prepared) -> TrajectoryResult:
    return run_method(bifurcation.ctx, "paga_dpt", ROOT)


def test_result_is_valid_and_rooted(result: TrajectoryResult, bifurcation: Prepared) -> None:
    result.validate()
    assert result.pseudotime[result.root.cell_index] == 0.0
    assert result.root.spec == "cluster:true_label=trunk_early"
    assert result.inference_rep == "X_diffmap"
    t = bifurcation.ds.obs["true_pseudotime"].to_numpy()
    ok = np.isfinite(result.pseudotime)
    assert (
        np.corrcoef(np.argsort(np.argsort(result.pseudotime[ok])), np.argsort(np.argsort(t[ok])))[0, 1] > 0.9
    )


def test_two_lineages_named_after_their_tip_clusters(result: TrajectoryResult) -> None:
    assert len(result.lineages) >= 2
    for lin in result.lineages:
        assert lin.name.startswith("Cluster ") or lin.name.startswith("Lineage ")
        assert lin.terminal_label is not None
        assert lin.n_cells > 0
    assert len({lin.name for lin in result.lineages}) == len(result.lineages)


def test_divergence_is_computed_for_every_pair(result: TrajectoryResult) -> None:
    n = len(result.lineages)
    assert len(result.divergence) == n * (n - 1) // 2
    assert all(0.0 <= d.pseudotime <= 1.0 for d in result.divergence)


def test_deterministic(bifurcation: Prepared, result: TrajectoryResult) -> None:
    again = run_method(bifurcation.ctx, "paga_dpt", ROOT)
    np.testing.assert_array_equal(again.pseudotime, result.pseudotime)
    np.testing.assert_array_equal(again.membership, result.membership)


def test_named_terminals_restrict_the_lineages(bifurcation: Prepared, result: TrajectoryResult) -> None:
    first_tip = result.lineages[0].terminal_label
    assert first_tip is not None
    only = run_method(bifurcation.ctx, "paga_dpt", ROOT, terminals=[first_tip])
    assert len(only.lineages) == 1
    assert only.n_unassigned > 0, "cells on the other branch get no pseudotime"
    assert any(w.code == "cells_off_lineages" for w in only.warnings)


def test_unknown_terminal_is_an_error(bifurcation: Prepared) -> None:
    with pytest.raises(ValueError, match="terminal clusters not found"):
        run_method(bifurcation.ctx, "paga_dpt", ROOT, terminals=["no-such-cluster"])


def test_params_are_validated(bifurcation: Prepared) -> None:
    with pytest.raises(ValidationError):
        run_method(bifurcation.ctx, "paga_dpt", ROOT, params={"n_dcs": 1})
    with pytest.raises(ValidationError):
        run_method(bifurcation.ctx, "paga_dpt", ROOT, params={"unknown": 1})


def test_custom_cluster_key_and_the_sibling_chain_diagnostic(bifurcation: Prepared) -> None:
    # With whole segments as groups, PAGA's tree chains the two sister branches
    # (trunk -> A -> B): A-B connectivity beats either branch's link to the trunk.
    # The method cannot fix that, but it must say so.
    r = run_method(bifurcation.ctx, "paga_dpt", ROOT, params={"cluster_key": "true_segment"})
    assert {m.label for m in r.milestones} == {"trunk", "A", "B"}
    assert len(r.lineages) == 1
    assert any(w.code == "possible_sibling_chain" for w in r.warnings)


def test_no_sibling_warning_on_a_correct_tree(result: TrajectoryResult) -> None:
    assert not any(w.code == "possible_sibling_chain" for w in result.warnings)


def test_unknown_method() -> None:
    with pytest.raises(UnknownMethodError, match="available: paga_dpt"):
        get_method("slingshot")


def test_registry_metadata() -> None:
    m = METHODS["paga_dpt"]
    assert m.available().ok
    schema = m.Params.model_json_schema()
    assert {"n_dcs", "min_connectivity", "cluster_key"} <= set(schema["properties"])
    assert schema["properties"]["cluster_key"]["default"] == LEIDEN_KEY
