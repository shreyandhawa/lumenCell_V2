import builtins
from pathlib import Path

import numpy as np
import pytest

from lumencell.api import binary
from lumencell.engine.methods.palantir_adapter import INSTALL_HINT, Palantir
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from lumencell.engine.store import RunStore
from tests.conftest import Prepared, prepare

pytest.importorskip("palantir")

ROOT = ClusterRoot("true_label", "trunk_early")


@pytest.fixture(scope="module")
def result(bifurcation: Prepared) -> TrajectoryResult:
    return run_method(bifurcation.ctx, "palantir", ROOT)


def test_result_is_valid_with_fates(result: TrajectoryResult) -> None:
    result.validate()
    assert result.topology_kind == "fate_star"
    assert result.pseudotime[result.root.cell_index] == 0.0
    assert result.fate is not None
    np.testing.assert_allclose(result.fate.sum(axis=0), 1.0, atol=1e-3)
    assert "entropy" in result.extras
    assert len(result.lineages) >= 2


def test_root_is_used_as_given(result: TrajectoryResult, bifurcation: Prepared) -> None:
    # Palantir would otherwise swap in a diffusion-component extreme for the start.
    t = bifurcation.ds.obs["true_pseudotime"].to_numpy()
    early = t < 0.1
    assert np.median(result.pseudotime[early]) < 0.25


def test_terminal_source_is_reported(result: TrajectoryResult) -> None:
    notice = next(w for w in result.warnings if w.code == "terminals_source")
    assert "PAGA" in notice.message


def test_named_terminals(bifurcation: Prepared) -> None:
    r = run_method(
        bifurcation.ctx,
        "palantir",
        ROOT,
        terminals=["A_late", "B_late"],
        params={"cluster_key": "true_label"},
    )
    assert len(r.lineages) == 2
    assert {lin.terminal_label for lin in r.lineages} == {"A_late", "B_late"}
    assert "you named" in next(w.message for w in r.warnings if w.code == "terminals_source")


def test_unknown_terminal(bifurcation: Prepared) -> None:
    with pytest.raises(ValueError, match="terminal clusters not found"):
        run_method(
            bifurcation.ctx, "palantir", ROOT, terminals=["nope"], params={"cluster_key": "true_label"}
        )


def test_linear_trajectory_falls_back_to_one_terminal(tmp_path: Path) -> None:
    linear = prepare(tmp_path / "lin.h5ad", topology="linear", n_cells=800)
    r = run_method(linear.ctx, "palantir", ROOT, params={"terminal_detection": "palantir"})
    r.validate()
    assert len(r.lineages) >= 1


def test_fates_and_extras_are_stored(tmp_path: Path, result: TrajectoryResult) -> None:
    store = RunStore(tmp_path)
    store.save("p", result)
    arrays = binary.decode(store.cells_bytes("p")).arrays
    assert arrays["fate"].shape == (len(result.lineages), result.n_obs)
    assert "extras.entropy" in arrays
    back = store.load("p")
    assert back.fate is not None and "entropy" in back.extras


def test_missing_install_gives_a_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    real_import = builtins.__import__

    def fake(name: str, *args: object, **kwargs: object) -> object:
        if name == "palantir":
            raise ImportError("no palantir")
        return real_import(name, *args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(builtins, "__import__", fake)
    avail = Palantir().available()
    assert not avail.ok and avail.install_hint == INSTALL_HINT
