import math
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast, get_args

import numpy as np
import pytest

from lumencell.bench.gates import check_gates, load_thresholds
from lumencell.bench.report import render
from lumencell.bench.runner import ClonalRow, ClonalRowStatus, clonal_row, run_real_case, run_synthetic_case
from lumencell.bench.synth import generate, truth_from_adata
from lumencell.bench.synth_clones import DEFAULT_LAYER, ClonalLayer, add_clones, clone_truth
from lumencell.bench.truth import Truth
from lumencell.clonal.labels import encode_clones
from lumencell.clonal.scores import ClonalResult, kmeans_partition, score_clonal, score_with_baseline
from lumencell.trust.constants import MIN_CLONE_SIZE, ORACLE_ETA_MIN, ORACLE_P_MAX, SENSITIVITY_CLONE_SIZE

GATES = load_thresholds()["synthetic"]["gates"]


def _row(tmp_path: Path, topology: str, **kw: Any) -> dict[str, Any]:
    return run_synthetic_case(tmp_path, topology, n_cells=1500, n_genes=1000, **kw).row()


def _clone_truth(layer: ClonalLayer = DEFAULT_LAYER, topology: str = "bifurcation") -> Truth:
    adata = generate(topology, "realistic", n_cells=1500, n_genes=200, seed=0)
    return clone_truth(truth_from_adata(adata), add_clones(adata, layer, seed=0))


def _result(membership: Any, truth: Truth) -> Any:
    return SimpleNamespace(membership=membership, pseudotime=truth.order.astype(np.float32))


def _pca(n: int) -> Any:
    return np.random.default_rng(0).normal(size=(n, 10))


def _gate_failures(row: ClonalRow) -> list[Any]:
    values = {
        "case": "c",
        "clonal_status": row.status,
        "clonal_excess_ratio": row.excess_ratio,
        "clonal_p": row.p,
    }
    return check_gates(values, {k: v for k, v in GATES.items() if k.startswith("clonal")})


def test_branching_cases_carry_a_clonal_score(tmp_path: Path) -> None:
    row = _row(tmp_path, "bifurcation")
    assert row["clonal_status"] in get_args(ClonalRowStatus)
    assert row["clonal_eta_oracle"] >= 0.05
    assert row["clonal_n_clones"] >= 20
    assert row["clonal_n_cells"] >= 2 * row["clonal_n_clones"]
    assert row["clonal_effective_clones"] > 0
    assert math.isfinite(row["clonal_eta_m3"])
    assert 0 < row["clonal_p_m3"] <= 1


def test_linear_cases_are_not_applicable(tmp_path: Path) -> None:
    row = _row(tmp_path, "linear")
    assert row["clonal_status"] == "not_applicable"
    assert math.isnan(row["clonal_excess_ratio"])
    assert row["clonal_n_cells"] == 0
    assert math.isnan(row["clonal_effective_clones"])
    assert math.isnan(row["clonal_eta_m3"])
    assert not [f for f in check_gates(row, GATES) if f.metric.startswith("clonal")]


def test_existing_metrics_do_not_move_with_the_clonal_layer(tmp_path: Path) -> None:
    with_layer = _row(tmp_path, "bifurcation", clonal_layer=DEFAULT_LAYER)
    without = _row(tmp_path, "bifurcation", clonal_layer=None)
    assert without["clonal_status"] == "not_measured"
    skip = {"runtime_seconds"}
    keys = {k for k in with_layer if k not in skip and not k.startswith("clonal_")}
    assert keys
    # repr, because NaN != NaN (terminal_share is NaN on a synthetic case)
    assert repr({k: with_layer[k] for k in sorted(keys)}) == repr({k: without[k] for k in sorted(keys)})


def test_a_single_lineage_method_fails_the_gate() -> None:
    truth = _clone_truth()
    row = clonal_row(truth, _result(truth.membership()[:1], truth), _pca(truth.n_obs), 10)
    assert row.status == "single_lineage"
    assert row.eta_oracle >= 0.05 and row.n_clones >= 20
    assert math.isnan(row.excess_ratio) and math.isnan(row.eta_m3)
    assert {f.metric for f in _gate_failures(row)} == {"clonal_excess_ratio", "clonal_p"}


def test_an_oracle_below_its_floor_is_not_assessable() -> None:
    truth = _clone_truth(ClonalLayer(fate_bias=0.0))
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert row.status == "not_assessable"
    assert math.isfinite(row.eta_method) and math.isfinite(row.eta_m3)
    assert math.isnan(row.excess_ratio)
    assert _gate_failures(row)


def test_a_scored_method_gets_a_ratio() -> None:
    truth = _clone_truth()
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert row.status == "ok"
    assert row.excess_ratio == pytest.approx(row.eta_method / row.eta_oracle)
    assert row.excess_ratio == pytest.approx(1.0)  # the method is the truth
    assert row.p <= 0.01 and not _gate_failures(row)


def test_a_truth_without_clones_is_not_measured_and_fails_the_gates() -> None:
    adata = generate("bifurcation", "realistic", n_cells=600, n_genes=100, seed=0)
    truth = truth_from_adata(adata)
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert row.status == "not_measured"
    assert {f.metric for f in _gate_failures(row)} == {"clonal_excess_ratio", "clonal_p"}


def test_a_truth_with_one_lineage_is_not_applicable() -> None:
    truth = _clone_truth(topology="linear")
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert repr(row) == repr(ClonalRow("not_applicable"))


def test_an_error_is_recorded_as_failed(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def boom(*args: object, **kwargs: object) -> None:
        raise RuntimeError("no scores today")

    monkeypatch.setattr("lumencell.clonal.scores.score_clonal", boom)
    truth = _clone_truth()
    row = clonal_row(truth, _result(truth.membership(), truth), _pca(truth.n_obs), 10)
    assert row.status == "failed"
    assert row.reason == "RuntimeError: no scores today"
    assert "no scores today" in capsys.readouterr().err


def test_the_report_says_the_clonal_gate_is_not_a_branch_gate() -> None:
    assert "not a branch gate" in render()


# ---- fix round 1: the row's fields, thresholds, helper wiring and run wiring, with a scripted scorer


@pytest.fixture(scope="module")
def truth3() -> Truth:
    adata = generate("trifurcation", "realistic", n_cells=1500, n_genes=200, seed=0)
    return clone_truth(truth_from_adata(adata), add_clones(adata, DEFAULT_LAYER, seed=0))


@pytest.fixture(scope="module")
def real_result() -> ClonalResult:
    truth = _clone_truth()
    assert truth.clone_ids is not None and truth.donors is not None
    codes = encode_clones(truth.clone_ids, truth.donors)
    return score_clonal(
        truth.membership(),
        codes,
        truth.order,
        kmeans_partition(_pca(truth.n_obs), 2, 0),
        oracle_cells=truth.exclusive() >= 0,
        n_permutations=20,
    )


class _Scripted:
    """Stands in for ``score_with_baseline``: answers from a list, one per call, and records the calls."""

    def __init__(self, template: ClonalResult, *answers: dict[str, Any]) -> None:
        self.template, self.answers = template, list(answers)
        self.calls: list[tuple[tuple[Any, ...], dict[str, Any]]] = []

    def __call__(self, *args: Any, **kwargs: Any) -> ClonalResult:
        self.calls.append((args, kwargs))
        return replace(self.template, **self.answers[len(self.calls) - 1])


def _scripted(monkeypatch: pytest.MonkeyPatch, template: ClonalResult, *answers: dict[str, Any]) -> _Scripted:
    spy = _Scripted(template, *answers)
    monkeypatch.setattr("lumencell.bench.runner.score_with_baseline", spy)
    return spy


def _run_row(truth: Truth, n_lineages: int | None = None) -> ClonalRow:
    membership = truth.membership()[:n_lineages]
    return clonal_row(truth, _result(membership, truth), _pca(truth.n_obs), 10)


def test_an_ok_row_takes_its_fields_from_the_method_and_the_m3_call(
    monkeypatch: pytest.MonkeyPatch, real_result: ClonalResult
) -> None:
    spy = _scripted(
        monkeypatch,
        real_result,
        {"eta": 0.6, "p": 0.001},
        {"eta": 0.3, "p": 0.004, "delta_eta": 0.07, "status": "ok"},
        {"eta": 0.2, "p": 0.03},
    )
    row = _run_row(_clone_truth())
    assert row.status == "ok"
    assert row.excess_ratio == pytest.approx(0.5)  # method over oracle, not 2.0
    assert (row.p, row.delta_base, row.eta_oracle, row.eta_method) == (0.004, 0.07, 0.6, 0.3)
    assert (row.eta_m3, row.p_m3) == (0.2, 0.03)
    sizes = [kw.get("min_clone_size", MIN_CLONE_SIZE) for _, kw in spy.calls]
    assert sizes == [MIN_CLONE_SIZE, MIN_CLONE_SIZE, SENSITIVITY_CLONE_SIZE]


ETA_MIN, P_MAX = ORACLE_ETA_MIN, ORACLE_P_MAX


@pytest.mark.parametrize(
    ("eta", "p", "status"),
    [
        (ETA_MIN, P_MAX, "ok"),  # both boundaries are assessable
        (ETA_MIN - 1e-9, 0.001, "not_assessable"),  # the eta floor alone
        (0.5, P_MAX + 1e-9, "not_assessable"),  # the p ceiling alone
        (math.nan, 0.001, "not_assessable"),
        (0.5, math.nan, "not_assessable"),
    ],
)
def test_the_oracle_floors_are_each_a_clause_and_inclusive(
    monkeypatch: pytest.MonkeyPatch, real_result: ClonalResult, eta: float, p: float, status: str
) -> None:
    _scripted(
        monkeypatch,
        real_result,
        {"eta": eta, "p": p},
        {"eta": 0.3, "p": 0.004, "status": "ok"},
        {"eta": 0.2, "p": 0.03},
    )
    row = _run_row(_clone_truth())
    assert row.status == status
    assert math.isfinite(row.eta_method) and math.isfinite(row.eta_m3)
    assert math.isnan(row.excess_ratio) == (status != "ok")


def test_too_few_clones_is_passed_through_with_its_values(
    monkeypatch: pytest.MonkeyPatch, real_result: ClonalResult
) -> None:
    _scripted(
        monkeypatch,
        real_result,
        {"eta": 0.6, "p": 0.001},
        {"eta": 0.3, "status": "too_few_clones"},
        {"eta": 0.2, "p": 0.03},
    )
    row = _run_row(_clone_truth())
    assert row.status == "too_few_clones"
    assert (row.eta_method, row.eta_m3, row.p_m3) == (0.3, 0.2, 0.03)
    assert math.isnan(row.excess_ratio)


def test_a_method_with_mass_in_one_lineage_is_single_lineage(
    monkeypatch: pytest.MonkeyPatch, real_result: ClonalResult
) -> None:
    spy = _scripted(monkeypatch, real_result, {"eta": 0.6}, {"status": "single_lineage"})
    row = _run_row(_clone_truth())
    assert row.status == "single_lineage"
    assert row.eta_oracle == 0.6 and row.n_clones == real_result.n_clones
    assert math.isnan(row.eta_method) and math.isnan(row.eta_m3) and math.isnan(row.excess_ratio)
    assert len(spy.calls) == 2, "no m = 3 call for a method that is not scored"


@pytest.mark.parametrize(("n_lineages", "n_calls"), [(1, 1), (2, 3)])
def test_a_method_needs_two_lineages_to_be_scored(
    monkeypatch: pytest.MonkeyPatch, real_result: ClonalResult, n_lineages: int, n_calls: int
) -> None:
    spy = _scripted(monkeypatch, real_result, {"eta": 0.6}, {"status": "ok"}, {"status": "ok"})
    row = _run_row(_clone_truth(), n_lineages)
    assert len(spy.calls) == n_calls
    assert (row.status == "single_lineage") == (n_lineages == 1)


@pytest.mark.parametrize(
    ("topology", "n_method"), [("bifurcation", 2), ("trifurcation", 2), ("trifurcation", 3)]
)
def test_oracle_and_method_get_their_own_k_and_the_oracle_cells(
    monkeypatch: pytest.MonkeyPatch, real_result: ClonalResult, truth3: Truth, topology: str, n_method: int
) -> None:
    truth = truth3 if topology == "trifurcation" else _clone_truth()
    n_truth = len(truth.names)
    assert n_truth == (3 if topology == "trifurcation" else 2)
    spy = _scripted(monkeypatch, real_result, {"eta": 0.6, "p": 0.001}, {"status": "ok"}, {"status": "ok"})
    _run_row(truth, n_lineages=n_method)
    cells = truth.exclusive() >= 0
    assert len(spy.calls) == 3
    (oracle_args, oracle_kw), (method_args, method_kw), (m3_args, m3_kw) = spy.calls
    assert oracle_args[5] == n_truth and oracle_args[0].shape[0] == n_truth
    for args, kw in ((method_args, method_kw), (m3_args, m3_kw)):
        assert args[5] == n_method and args[0].shape[0] == n_method
        assert np.array_equal(kw["oracle_cells"], cells)
    assert np.array_equal(oracle_kw["oracle_cells"], cells)


def test_the_helper_passes_k_and_the_scoring_inputs_through(monkeypatch: pytest.MonkeyPatch) -> None:
    seen: dict[str, Any] = {}
    partition = np.arange(30) % 3

    def fake_kmeans(space: Any, k: int, seed: int) -> Any:
        seen.update(space=space, k=k, seed=seed)
        return partition

    def fake_score(membership: Any, codes: Any, pseudotime: Any, baseline: Any, **kw: Any) -> Any:
        seen.update(baseline=baseline, kw=kw)
        return "scored"

    monkeypatch.setattr("lumencell.clonal.scores.kmeans_partition", fake_kmeans)
    monkeypatch.setattr("lumencell.clonal.scores.score_clonal", fake_score)
    cells = np.ones(30, dtype=bool)
    out: Any = score_with_baseline(
        np.ones((3, 30), bool),
        cast(Any, None),
        np.zeros(30),
        _pca(30),
        5,
        3,
        oracle_cells=cells,
        min_clone_size=3,
        seed=4,
    )
    assert out == "scored"
    assert (seen["k"], seen["seed"], seen["space"].shape) == (3, 4, (30, 5))
    assert seen["baseline"] is partition
    assert seen["kw"]["oracle_cells"] is cells
    assert (seen["kw"]["min_clone_size"], seen["kw"]["seed"]) == (3, 4)


class _Stop(Exception):
    pass


@pytest.mark.parametrize(
    ("topology", "layer", "expected"),
    [
        ("bifurcation", DEFAULT_LAYER, ("clone_id", "donor")),
        ("bifurcation", None, ()),
        ("linear", DEFAULT_LAYER, ()),
    ],
)
def test_inference_is_told_the_clone_columns_only_when_the_layer_was_applied(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    topology: str,
    layer: ClonalLayer | None,
    expected: tuple[str, ...],
) -> None:
    seen: list[Any] = []

    def spy(ds: Any, prep: Any, params: Any, clonal_columns: Any = ()) -> None:
        seen.append(tuple(clonal_columns))
        raise _Stop

    monkeypatch.setattr("lumencell.bench.runner.InferenceContext.from_prep", spy)
    with pytest.raises(_Stop):
        run_synthetic_case(tmp_path, topology, n_cells=400, n_genes=200, clonal_layer=layer)
    assert seen == [expected]


def test_a_real_case_row_is_not_applicable(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    truth = _clone_truth()
    spec = SimpleNamespace(root_spec="r", label_key=None)
    monkeypatch.setattr("lumencell.bench.datasets.DATASETS", {"d": spec})
    monkeypatch.setattr("lumencell.bench.datasets.is_cached", lambda ds: True)
    monkeypatch.setattr("lumencell.bench.datasets.local_path", lambda ds: tmp_path)
    monkeypatch.setattr("lumencell.bench.datasets.truth_for", lambda ds, obs: truth)
    monkeypatch.setattr("lumencell.engine.roots.parse_root", lambda s: s)
    monkeypatch.setattr("lumencell.bench.runner.load", lambda p: SimpleNamespace(obs=None))
    monkeypatch.setattr(
        "lumencell.prep.cache.PrepCache",
        lambda data: SimpleNamespace(run=lambda pp: SimpleNamespace(artifacts={})),
    )
    monkeypatch.setattr(
        "lumencell.bench.runner.InferenceContext.from_prep",
        lambda *a, **k: SimpleNamespace(leiden=np.array([0, 1])),
    )
    monkeypatch.setattr(
        "lumencell.bench.runner.run_method",
        lambda *a, **k: SimpleNamespace(
            pseudotime=truth.order.astype(np.float32), membership=truth.membership(), warnings=[]
        ),
    )
    res = run_real_case("d")
    assert repr(res.clonal) == repr(ClonalRow("not_applicable"))


def test_the_report_counts_failures_per_gate() -> None:
    text = render()
    line = next(x for x in text.splitlines() if x.startswith("Cases failing each gate"))
    assert "`clonal_excess_ratio`" in line
    assert "PAGA + DPT" in line


def test_the_clonal_column_never_rounds_up_across_the_gate() -> None:
    from lumencell.bench.report import _range

    rows = [{"clonal_excess_ratio": 0.4996}, {"clonal_excess_ratio": 0.5004}]
    assert _range(rows, "clonal_excess_ratio", floor=True) == "0.499–0.500"
    assert _range([{"clonal_excess_ratio": 0.4996}], "clonal_excess_ratio", floor=True) == "0.499–0.499"
