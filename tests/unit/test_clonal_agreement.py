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
from lumencell.clonal.labels import encode_clones
from lumencell.clonal.scores import ClonalResult, Coupling, NullComparison, score_with_baseline
from lumencell.data.loader import load
from lumencell.engine.context import InferenceContext
from lumencell.engine.registry import run_method
from lumencell.engine.roots import ClusterRoot
from lumencell.engine.schema import TrajectoryResult
from lumencell.engine.store import TRUST_NAME, RunStore
from lumencell.prep.pipeline import run_pipeline
from lumencell.prep.stages import PrepParams
from lumencell.trust import clonal
from lumencell.trust.clonal import (
    COUPLING_NOTE,
    agreement_for_run,
    not_scored,
    to_api_coupling,
    to_api_score,
    trust_hash,
    trust_name,
)
from lumencell.trust.constants import CLONAL_SEED, MIN_CLONES, N_PERMUTATIONS

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
    assert got.primary.n_donors >= 2 and got.primary.n_cells > 0
    assert got.primary.n_permutations == N_PERMUTATIONS and got.primary.seed == CLONAL_SEED
    assert got.primary.eta_base is not None and got.primary.delta_eta is not None
    assert [(c.a, c.b) for c in got.coupling] == [
        (a, b) for a in range(n_lineages) for b in range(a + 1, n_lineages)
    ]


def _expected(case: Case, obs: pd.DataFrame, min_clone_size: int, k: int | None = None) -> ClonalResult:
    """The score the agreement should hold, computed here with the brief's wiring."""
    codes = encode_clones(obs["clone_id"].to_numpy(dtype=object), obs["donor"].to_numpy(dtype=object))
    r = case.result
    return score_with_baseline(
        r.membership,
        codes,
        r.pseudotime,
        case.pca,
        case.neighbors_pcs,
        r.membership.shape[0] if k is None else k,
        min_clone_size=min_clone_size,
        seed=CLONAL_SEED,
    )


def test_the_scores_use_a_k_equals_l_baseline_and_the_clonal_seed(case: Case) -> None:
    got = _agree(case)
    primary, m3 = _expected(case, case.obs, 2), _expected(case, case.obs, 3)
    assert got.primary == to_api_score(primary)
    assert got.sensitivity_m3 == to_api_score(m3)
    assert got.coupling == [to_api_coupling(c) for c in primary.coupling]
    # k pins: another k gives another baseline.
    other = _expected(case, case.obs, 2, k=case.result.membership.shape[0] + 1)
    assert other.eta_base != primary.eta_base


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


def test_the_refusal_names_every_steered_column_once(case: Case) -> None:
    got = _agree(case, obs_keys=["clone_id", "donor", "clone_id"])
    assert got.status == "refused"
    (note,) = got.notes
    assert "'clone_id'" in note and "'donor'" in note
    assert note.count("'clone_id'") == 1


def test_settings_that_do_not_fit_the_file_fail_with_the_problems(case: Case) -> None:
    got = _agree(case, settings=ClonalSettings(clone_key="nope", donor_key="donor"))
    assert got.status == "failed"
    assert got.primary is None and got.coupling == []
    (note,) = got.notes
    assert "there is no obs column named 'nope'" in note
    obs = case.obs.copy()
    obs["clone_id"] = np.arange(case.result.n_obs, dtype=float)
    got = _agree(case, obs=obs)
    assert got.status == "failed"
    assert "must hold text labels" in got.notes[0]


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
    assert got.notes == ["The run has one lineage; clones are scored over two or more."]


def test_clones_in_one_lineage_are_single_lineage_with_a_reason(case: Case) -> None:
    labels = np.full(case.result.n_obs, None, dtype=object)
    donor = case.obs["donor"].to_numpy(dtype=object)
    a = _exclusive(case, 0)
    a = a[donor[a] == donor[a[0]]]
    for i in range(a.size // 2):
        labels[a[2 * i : 2 * i + 2]] = f"c{i}"
    got = _agree(case, obs=_with_clones(case, labels))
    assert got.status == "single_lineage"
    assert got.primary is not None and got.primary.status == "single_lineage"
    assert got.coupling == []
    assert got.notes == ["The scored clones fall in only one lineage; clones are scored over two or more."]


def test_missing_donors_are_noted(case: Case) -> None:
    obs = case.obs.copy()
    labelled = np.flatnonzero(obs["clone_id"].notna().to_numpy())[:10]
    donor = obs["donor"].to_numpy(dtype=object).copy()
    donor[labelled] = None
    obs["donor"] = donor
    got = _agree(case, obs=obs)
    assert any(n.startswith("10 cells with a receptor have no donor") for n in got.notes)


def test_one_missing_donor_is_noted(case: Case) -> None:
    obs = case.obs.copy()
    donor = obs["donor"].to_numpy(dtype=object).copy()
    donor[np.flatnonzero(obs["clone_id"].notna().to_numpy())[0]] = None
    obs["donor"] = donor
    got = _agree(case, obs=obs)
    assert "1 cells with a receptor have no donor and were not scored." in got.notes


def _with_clones(case: Case, labels: npt.NDArray[np.object_]) -> pd.DataFrame:
    obs = case.obs.copy()
    obs["clone_id"] = labels
    return obs


@pytest.mark.parametrize("n_clones", [0, 1])
def test_undefined_values_are_null_and_the_file_is_strict_json(
    case: Case, tmp_path: Path, n_clones: int
) -> None:
    labels = np.full(case.result.n_obs, None, dtype=object)
    # Singletons only: the settings are valid, but no clone reaches m = 2.
    labels[_exclusive(case, 0)[:5]] = [f"s{i}" for i in range(5)]
    if n_clones == 1:
        # One two-cell clone over two lineages: in one lineage it would be single_lineage instead.
        # One donor, as clone identity is the pair (donor, label). At m = 3 it scores no clone,
        # so the note must come from the primary.
        donor = case.obs["donor"].to_numpy(dtype=object)
        a, b = _exclusive(case, 0), _exclusive(case, 1)
        d = donor[a[0]]
        labels[a[donor[a] == d][:1]] = "c1"
        labels[b[donor[b] == d][:1]] = "c1"
    got = _agree(case, obs=_with_clones(case, labels))
    assert got.status == "too_few_clones"
    assert got.primary is not None and got.primary.status == "too_few_clones"
    assert got.primary.n_clones == n_clones
    assert any(
        n == f"Only {n_clones} clones were scored; at least {MIN_CLONES} are needed before the score counts."
        for n in got.notes
    )
    assert got.sensitivity_m3 is not None and got.sensitivity_m3.n_clones == 0
    if n_clones == 1:
        assert got.primary.purity_without_largest is None
        assert got.primary.purity is not None
        assert got.coupling and COUPLING_NOTE in got.notes
    else:
        assert got.primary.purity is None and got.primary.eta is None and got.primary.nulls == []
        assert got.coupling == [] and COUPLING_NOTE not in got.notes

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
    # Coupling comes from the primary: m = 3 has none here.
    expected = _expected(case, _with_clones(case, labels), 2)
    assert got.coupling and got.coupling == [to_api_coupling(c) for c in expected.coupling]


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


@pytest.mark.parametrize(
    "name",
    [
        "MIN_CLONE_SIZE",
        "SENSITIVITY_CLONE_SIZE",
        "MIN_CLONES",
        "N_PERMUTATIONS",
        "CLONAL_SEED",
        "PSEUDOTIME_BINS",
        "CLONAL_METRIC_VERSION",
    ],
)
def test_the_trust_hash_covers_each_scoring_parameter(monkeypatch: pytest.MonkeyPatch, name: str) -> None:
    before = trust_hash(S)
    monkeypatch.setattr(f"lumencell.trust.clonal.{name}", getattr(clonal, name) + 1)
    assert trust_hash(S) != before


def test_not_scored_carries_one_note() -> None:
    got = not_scored("failed", S, "the run was prepared with other settings")
    assert got.status == "failed" and got.settings == S
    assert got.primary is None and got.sensitivity_m3 is None and got.coupling == []
    assert got.notes == ["the run was prepared with other settings"]


NAN, INF = float("nan"), float("inf")


def _result(**kw: Any) -> ClonalResult:
    base: dict[str, Any] = {
        "status": "too_few_clones", "mode": "own_cells", "n_lineages": 2, "n_cells": 0, "n_clones": 0,
        "effective_clones": NAN, "n_donors": 0, "n_strata": None, "min_clone_size": 2,
        "n_permutations": 1000, "seed": 0, "purity": NAN, "mean_entropy": NAN, "nulls": (), "eta": NAN,
        "p": NAN, "eta_base": NAN, "delta_eta": NAN, "p_delta": NAN, "clone_weighted_purity": NAN,
        "purity_without_largest": NAN, "coupling": (),
    }  # fmt: skip
    return ClonalResult(**{**base, **kw})


@pytest.mark.parametrize("bad", [NAN, INF, -INF])
def test_non_finite_values_map_to_none(bad: float) -> None:
    null = NullComparison(null="N1", mean_entropy_mean=bad, mean_entropy_sd=bad, eta=bad, p=bad)
    fields = {f: bad for f in ("effective_clones", "purity", "mean_entropy", "eta", "p", "eta_base")}
    got = to_api_score(_result(nulls=(null,), delta_eta=bad, clone_weighted_purity=bad, **fields))
    dumped = got.model_dump()
    floats = [k for k, v in dumped.items() if k not in ("nulls", "n_strata") and not isinstance(v, int | str)]
    assert floats and all(dumped[k] is None for k in floats)
    assert got.nulls[0].model_dump() == {
        "null": "N1", "mean_entropy_mean": None, "mean_entropy_sd": None, "eta": None, "p": None
    }  # fmt: skip
    c = to_api_coupling(Coupling(a=0, b=1, observed=bad, null_mean=bad, z=bad, p=bad, p_adj=bad))
    assert c.model_dump() == {
        "a": 0,
        "b": 1,
        "observed": None,
        "null_mean": None,
        "z": None,
        "p": None,
        "p_adj": None,
    }
    json.dumps(got.model_dump(mode="json"), allow_nan=False)


def test_finite_values_map_through_unchanged() -> None:
    nulls = (
        NullComparison(null="N1", mean_entropy_mean=0.7, mean_entropy_sd=0.01, eta=0.3, p=0.002),
        NullComparison(null="N2", mean_entropy_mean=0.6, mean_entropy_sd=0.02, eta=0.2, p=0.004),
    )
    r = _result(
        status="ok", mode="oracle_cells", n_lineages=3, n_cells=91, n_clones=33, effective_clones=21.5,
        n_donors=4, n_strata=7, min_clone_size=3, n_permutations=17, seed=5, purity=0.51,
        mean_entropy=0.49, nulls=nulls, eta=0.3, p=0.004, eta_base=0.1, delta_eta=0.2, p_delta=0.03,
        clone_weighted_purity=0.55, purity_without_largest=0.47,
    )  # fmt: skip
    assert to_api_score(r).model_dump() == {
        "status": "ok", "mode": "oracle_cells", "n_lineages": 3, "min_clone_size": 3, "n_cells": 91,
        "n_clones": 33, "effective_clones": 21.5, "n_donors": 4, "n_strata": 7, "n_permutations": 17,
        "seed": 5, "purity": 0.51, "mean_entropy": 0.49,
        "nulls": [
            {"null": "N1", "mean_entropy_mean": 0.7, "mean_entropy_sd": 0.01, "eta": 0.3, "p": 0.002},
            {"null": "N2", "mean_entropy_mean": 0.6, "mean_entropy_sd": 0.02, "eta": 0.2, "p": 0.004},
        ],
        "eta": 0.3, "p": 0.004, "eta_base": 0.1, "delta_eta": 0.2, "p_delta": 0.03,
        "clone_weighted_purity": 0.55, "purity_without_largest": 0.47,
    }  # fmt: skip
    c = Coupling(a=1, b=2, observed=0.4, null_mean=0.2, z=3.1, p=0.01, p_adj=0.02)
    assert to_api_coupling(c).model_dump() == {
        "a": 1, "b": 2, "observed": 0.4, "null_mean": 0.2, "z": 3.1, "p": 0.01, "p_adj": 0.02
    }  # fmt: skip


def test_p_delta_is_none_when_delta_eta_is_undefined() -> None:
    got = to_api_score(_result(delta_eta=NAN, p_delta=0.5))
    assert got.delta_eta is None and got.p_delta is None


def test_no_note_claims_direction(case: Case) -> None:
    obs = case.obs.copy()
    donor = obs["donor"].to_numpy(dtype=object).copy()
    donor[np.flatnonzero(obs["clone_id"].notna().to_numpy())[:3]] = None
    obs["donor"] = donor
    r = case.result
    agreements = [
        _agree(case, obs=obs),
        _agree(case, obs_keys=["clone_id"]),
        _agree(case, result=replace(r, membership=r.membership[:1], lineages=r.lineages[:1])),
        _agree(case, obs=_with_clones(case, np.array([f"s{i}" for i in range(r.n_obs)], dtype=object))),
        _agree(case, settings=ClonalSettings(clone_key="nope", donor_key="donor")),
    ]
    notes = {n for a in agreements for n in a.notes}
    assert COUPLING_NOTE in notes and len(notes) >= 5
    assert "it does not say which came first" in COUPLING_NOTE
    for note in notes:
        for word in (
            "precede",
            "gives rise",
            "differentiates into",
            "before it",
            "ancestor of",
            "came first",
        ):
            assert word not in note.lower() or note == COUPLING_NOTE, note
