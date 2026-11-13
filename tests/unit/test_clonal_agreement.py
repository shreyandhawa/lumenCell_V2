"""Own-cells clonal agreement of a run, as the score_clonal job builds and stores it."""

import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest

from lumencell.api.models import ClonalSettings
from lumencell.bench.synth import generate
from lumencell.bench.synth_clones import DEFAULT_LAYER, add_clones
from lumencell.data.loader import load
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from lumencell.engine.store import TRUST_NAME, RunStore
from lumencell.prep.pipeline import run_pipeline
from lumencell.prep.stages import PrepParams
from lumencell.trust.clonal import (
    COUPLING_NOTE,
    agreement_for_run,
    not_scored,
    to_api_score,
    trust_hash,
    trust_name,
)
from lumencell.trust.constants import MIN_CLONES

S = ClonalSettings(clone_key="clone_id", donor_key="donor")


@dataclass(frozen=True)
class Case:
    result: TrajectoryResult
    obs: pd.DataFrame
    pca: npt.NDArray[np.floating[Any]]
    neighbors_pcs: int


@pytest.fixture(scope="module")
def case(tmp_path_factory: pytest.TempPathFactory) -> Case:
    adata = add_clones(
        generate("bifurcation", "realistic", n_cells=1500, n_genes=1000, seed=0), DEFAULT_LAYER, 0
    )
    path = tmp_path_factory.mktemp("agreement") / "d.h5ad"
    adata.write_h5ad(path)
    ds = load(path)
    pp = PrepParams()
    prep = run_pipeline(ds, pp)
    ctx = InferenceContext.from_prep(ds, prep, pp, clonal_columns=("clone_id", "donor"))
    result = run_method(ctx, "paga_dpt", ClusterRoot("true_label", "trunk_early"))
    assert result.membership.shape[0] >= 2
    return Case(result, ds.obs, np.asarray(prep.artifacts["pca"]), pp.neighbors_pcs)


def _agree(case: Case, settings: ClonalSettings = S, **kw: Any) -> Any:
    args: dict[str, Any] = {
        "result": case.result,
        "obs": case.obs,
        "settings": settings,
        "pca": case.pca,
        "neighbors_pcs": case.neighbors_pcs,
        "obs_keys": case.result.obs_keys,
    }
    args.update(kw)
    return agreement_for_run(**args)


def _exclusive(case: Case, lineage: int) -> npt.NDArray[np.intp]:
    m = case.result.membership
    return np.flatnonzero(m[lineage] & (m.sum(axis=0) == 1))


def _raise(constant: str) -> None:
    raise ValueError(f"not strict JSON: {constant}")


def test_a_run_with_clones_scores_ok(case: Case) -> None:
    got = _agree(case)
    n_lineages = case.result.membership.shape[0]
    assert got.status == "ok"
    assert got.settings == S
    assert got.primary is not None and got.sensitivity_m3 is not None
    assert got.primary.status == "ok"
    assert got.primary.mode == "own_cells"
    assert got.primary.min_clone_size == 2
    assert got.sensitivity_m3.min_clone_size == 3
    assert got.primary.n_lineages == n_lineages
    assert got.primary.n_clones >= 20
    assert [n.null for n in got.primary.nulls] == ["N1", "N2"]
    assert got.primary.eta is not None and got.primary.p is not None
    assert len(got.coupling) == n_lineages * (n_lineages - 1) // 2
    assert COUPLING_NOTE in got.notes
    assert not any("no donor" in n for n in got.notes)


def test_strata_add_the_third_null(case: Case) -> None:
    got = _agree(case, ClonalSettings(clone_key="clone_id", donor_key="donor", strata_keys=["true_segment"]))
    assert got.primary is not None
    assert [n.null for n in got.primary.nulls] == ["N1", "N2", "N3"]
    assert got.primary.n_strata is not None and got.primary.n_strata >= 3


def test_a_run_steered_by_a_clone_column_is_refused(case: Case) -> None:
    got = _agree(case, obs_keys=["true_label", "clone_id"])
    assert got.status == "refused"
    assert got.primary is None and got.sensitivity_m3 is None and got.coupling == []
    assert any("'clone_id'" in n for n in got.notes)
    assert not any("true_label" in n for n in got.notes)


def test_a_run_steered_by_the_donor_column_is_refused(case: Case) -> None:
    got = _agree(case, obs_keys=["donor"])
    assert got.status == "refused"
    assert any("'donor'" in n for n in got.notes)


def test_one_lineage_is_single_lineage(case: Case) -> None:
    r = case.result
    one = replace(r, membership=r.membership[:1], lineages=r.lineages[:1])
    got = _agree(case, result=one)
    assert got.status == "single_lineage"
    assert got.primary is None and got.coupling == []
    assert got.notes


def test_missing_donors_are_noted(case: Case) -> None:
    obs = case.obs.copy()
    labelled = np.flatnonzero(obs["clone_id"].notna().to_numpy())[:10]
    donor = obs["donor"].to_numpy(dtype=object).copy()
    donor[labelled] = None
    obs["donor"] = donor
    got = _agree(case, obs=obs)
    assert any(n.startswith("10 cells with a receptor have no donor") for n in got.notes)


def _with_clones(case: Case, labels: npt.NDArray[np.object_]) -> pd.DataFrame:
    obs = case.obs.copy()
    obs["clone_id"] = labels
    return obs


@pytest.mark.parametrize("n_clones", [0, 1])
def test_undefined_values_are_null_and_the_file_is_strict_json(
    case: Case, tmp_path: Path, n_clones: int
) -> None:
    labels = np.full(case.result.n_obs, None, dtype=object)
    if n_clones == 1:
        # One clone over two lineages: in one lineage it would be single_lineage instead.
        # One donor, as clone identity is the pair (donor, label).
        donor = case.obs["donor"].to_numpy(dtype=object)
        a, b = _exclusive(case, 0), _exclusive(case, 1)
        d = donor[a[0]]
        labels[a[donor[a] == d][:3]] = "c1"
        labels[b[donor[b] == d][:3]] = "c1"
    got = _agree(case, obs=_with_clones(case, labels))
    assert got.status == "too_few_clones"
    assert got.primary is not None and got.primary.status == "too_few_clones"
    assert got.primary.n_clones == n_clones
    assert any(
        n == f"Only {n_clones} clones were scored; at least {MIN_CLONES} are needed before the score counts."
        for n in got.notes
    )
    if n_clones == 1:
        assert got.primary.purity_without_largest is None
        assert got.primary.purity is not None
    else:
        assert got.primary.purity is None and got.primary.eta is None and got.primary.nulls == []

    store = RunStore(tmp_path)
    store.save("r", case.result)
    name = trust_name(S)
    target = store.write_trust("r", name, got.model_dump(mode="json"))
    parsed = json.loads(target.read_text(encoding="utf-8"), parse_constant=_raise)
    assert parsed["status"] == "too_few_clones"
    assert parsed["primary"]["purity_without_largest"] is None


def test_an_m3_shortfall_shows_under_an_ok_primary(case: Case) -> None:
    # Every clone has exactly two cells, all in one donor: m = 2 is ok, m = 3 scores none.
    labels = np.full(case.result.n_obs, None, dtype=object)
    donor = case.obs["donor"].to_numpy(dtype=object)
    a, b = _exclusive(case, 0), _exclusive(case, 1)
    a, b = a[donor[a] == donor[a[0]]], b[donor[b] == donor[a[0]]]
    n = min(a.size, b.size) // 2
    assert n >= MIN_CLONES
    for i in range(n):
        labels[[a[2 * i], b[2 * i + 1] if i % 2 else a[2 * i + 1]]] = f"c{i}"
    got = _agree(case, obs=_with_clones(case, labels))
    assert got.status == "ok"
    assert got.primary is not None and got.primary.n_clones == n
    assert got.sensitivity_m3 is not None
    assert got.sensitivity_m3.status == "too_few_clones"
    assert got.sensitivity_m3.n_clones == 0


def test_the_trust_name_hides_column_names(case: Case) -> None:
    sneaky = ClonalSettings(clone_key="../../etc", donor_key="donor")
    obs = case.obs.rename(columns={"clone_id": "../../etc"})
    got = _agree(case, settings=sneaky, obs=obs)
    assert got.status == "ok"
    name = trust_name(sneaky)
    assert TRUST_NAME.fullmatch(name)
    assert name == f"clonal-{trust_hash(sneaky)}"
    assert "etc" not in name and "/" not in name
    ab = ClonalSettings(clone_key="c", donor_key="d", strata_keys=["a", "b"])
    ba = ClonalSettings(clone_key="c", donor_key="d", strata_keys=["b", "a"])
    assert trust_hash(ab) == trust_hash(ba)
    assert trust_hash(ab) != trust_hash(ab.model_copy(update={"donor_key": "e"}))
    assert trust_hash(ab) != trust_hash(
        ClonalSettings(clone_key="c", single_donor=True, strata_keys=["a", "b"])
    )


def test_the_trust_hash_covers_the_scoring_parameters(monkeypatch: pytest.MonkeyPatch) -> None:
    before = trust_hash(S)
    monkeypatch.setattr("lumencell.trust.clonal.N_PERMUTATIONS", 999)
    assert trust_hash(S) != before


def test_not_scored_carries_one_note() -> None:
    got = not_scored("failed", S, "the run was prepared with other settings")
    assert got.status == "failed" and got.settings == S
    assert got.primary is None and got.sensitivity_m3 is None and got.coupling == []
    assert got.notes == ["the run was prepared with other settings"]


def test_nan_maps_to_none() -> None:
    from lumencell.clonal.scores import ClonalResult

    nan = float("nan")
    r = ClonalResult(
        status="too_few_clones", mode="own_cells", n_lineages=2, n_cells=0, n_clones=0,
        effective_clones=nan, n_donors=0, n_strata=None, min_clone_size=2, n_permutations=1000, seed=0,
        purity=nan, mean_entropy=nan, nulls=(), eta=nan, p=nan, eta_base=nan, delta_eta=nan, p_delta=nan,
        clone_weighted_purity=nan, purity_without_largest=nan, coupling=(),
    )  # fmt: skip
    got = to_api_score(r)
    assert got.effective_clones is None and got.p is None and got.delta_eta is None
    assert got.status == "too_few_clones"


def test_no_note_claims_direction() -> None:
    for note in (COUPLING_NOTE,):
        assert "it does not say which came first" in note
        for word in ("precede", "gives rise", "differentiates into", "before it", "ancestor of"):
            assert word not in note.lower()
