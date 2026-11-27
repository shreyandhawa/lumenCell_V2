"""Run admission: what a request for a run becomes, or why it is refused."""

from dataclasses import replace
from typing import Any

import pytest

from lumencell.engine import registry
from lumencell.engine.admission import AdmittedRun, RunRefused, admit_run
from lumencell.engine.context import InferenceContext
from lumencell.engine.methods.base import Availability
from lumencell.engine.roots import ClusterRoot, format_root, parse_root
from lumencell.engine.store import run_id
from lumencell.prep.stages import PrepParams
from tests.conftest import Prepared

ROOT = "cluster:true_label=trunk_early"
KEY = "prep"
BAD_PARAMS = {"n_dcs": 1}  # below paga_dpt's minimum of 2


@pytest.fixture(scope="module")
def ctx(bifurcation: Prepared) -> InferenceContext:
    """The bifurcation plus a clone column, a donor column named as clonal, and a column named with '='."""
    n = bifurcation.ds.n_obs
    obs = bifurcation.ds.obs.assign(
        clone_id=[f"c{i % 7}" for i in range(n)],
        donor=["d1" if i % 2 else "d2" for i in range(n)],
        **{"batch=1": ["A " if i % 3 else "B" for i in range(n)]},
    )
    ds = replace(bifurcation.ds, obs=obs)
    return InferenceContext.from_prep(ds, bifurcation.prep, PrepParams(), clonal_columns=("donor",))


def admit(ctx: InferenceContext, **request: Any) -> AdmittedRun:
    asked = {"method": "paga_dpt", "params": None, "root": ROOT, "terminals": None, **request}
    return admit_run(asked["method"], asked["params"], asked["root"], asked["terminals"], lambda: ctx, KEY)


def refusal(ctx: InferenceContext, **request: Any) -> RunRefused:
    with pytest.raises(RunRefused) as caught:
        admit(ctx, **request)
    return caught.value


def test_an_admitted_run_is_canonical_and_carries_its_id(ctx: InferenceContext) -> None:
    escaped = 'json:{"kind": "cluster", "key": "true_label", "value": "trunk_early"}'
    admitted = admit(ctx, root=escaped, terminals=["A_late"])
    spec = admitted.spec
    assert spec.root == ROOT, "a root in any spelling is kept as its canonical text"
    assert spec.params == registry.get_method("paga_dpt").Params.model_validate({}).model_dump()
    assert spec.terminals == ["A_late"]
    assert spec.clonal_columns == ("donor",), "the clonal columns are the context's"
    assert admitted.run_id == run_id(spec, KEY)
    assert admit(ctx, terminals=["A_late"]).run_id == admitted.run_id, "one run, however it was asked for"


def test_a_root_the_plain_text_cannot_carry_comes_back_unchanged(ctx: InferenceContext) -> None:
    root = ClusterRoot("batch=1", "A ")
    assert parse_root(admit(ctx, root=format_root(root)).spec.root) == root


@pytest.mark.parametrize(
    ("request_", "code"),
    [
        # Each request also has every later problem, so only the order decides the code.
        ({"method": "nope", "root": "", "params": BAD_PARAMS}, "unknown_method"),
        ({"root": "", "params": BAD_PARAMS}, "root_required"),
        ({"root": "auto", "params": BAD_PARAMS}, "invalid_params"),
        ({"root": "auto", "params": {"cluster_key": "clone_id"}}, "invalid_root"),
        ({"root": "cluster:clone_id=nowhere"}, "scoring_only_column"),
        ({"root": "cluster:true_label=nowhere"}, "invalid_root"),
    ],
)
def test_the_first_problem_is_the_one_refused(
    ctx: InferenceContext, request_: dict[str, Any], code: str
) -> None:
    assert refusal(ctx, **request_).code == code


def test_bad_params_are_refused_with_each_error(ctx: InferenceContext) -> None:
    refused = refusal(ctx, params=BAD_PARAMS)
    assert refused.message == "The parameters for PAGA + DPT are not valid."
    assert refused.details is not None and [e["loc"] for e in refused.details["errors"]] == [["n_dcs"]]


def test_an_uninstalled_method_is_refused_with_its_install_hint(
    ctx: InferenceContext, monkeypatch: pytest.MonkeyPatch
) -> None:
    missing = Availability(ok=False, install_hint="pip install it")
    monkeypatch.setattr(registry.METHODS["paga_dpt"], "available", lambda: missing)
    refused = refusal(ctx, root="")
    assert refused.code == "method_unavailable"
    assert refused.message == "PAGA + DPT is not installed. pip install it"


def test_the_context_is_asked_for_only_once_the_request_is_sound() -> None:
    class NotReady(Exception):
        pass

    def not_ready() -> InferenceContext:
        raise NotReady  # as a dataset still being prepared does

    unsound = (("nope", None, ROOT), ("paga_dpt", None, ""), ("paga_dpt", BAD_PARAMS, ROOT))
    for method, params, root in unsound:
        with pytest.raises(RunRefused):
            admit_run(method, params, root, None, not_ready, KEY)
    with pytest.raises(NotReady):
        admit_run("paga_dpt", None, ROOT, None, not_ready, KEY)


def test_every_scoring_only_column_is_named_once(ctx: InferenceContext) -> None:
    both = refusal(ctx, root="cluster:clone_id=c1", params={"cluster_key": "donor"})
    assert both.code == "scoring_only_column"
    assert both.details == {"columns": ["clone_id", "donor"]}
    assert both.message.endswith("it cannot be the root column")

    twice = refusal(ctx, root="cluster:clone_id=c1", params={"cluster_key": "clone_id"})
    assert twice.details == {"columns": ["clone_id"]}, "each column is named once"

    donor = refusal(ctx, root="cluster:donor=d1", params={"cluster_key": "donor"})
    assert donor.details == {"columns": ["donor"]}
    assert donor.message == (
        "'donor' is set as a clone or donor column, which lumenCell uses only to score runs; "
        "it cannot be the root column"
    )
    cluster = refusal(ctx, params={"cluster_key": "donor"})
    assert cluster.message.endswith("it cannot be the cluster column")
