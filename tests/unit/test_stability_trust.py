import json
import math
from pathlib import Path

import pytest

from lumencell.api import models
from lumencell.engine.stability import KindSummary, ReplicateRow, StabilityResult, summarise
from lumencell.engine.store import TRUST_NAME, RunStore
from lumencell.trust import stability as trust_stability
from lumencell.trust.stability import STABILITY_NOTE, not_measured, stability_name, to_api_stability

nan = math.nan


def _row(kind: str, index: int, error: str | None = None, spearman: float = 0.9) -> ReplicateRow:
    failed = error is not None
    return ReplicateRow(
        kind,  # type: ignore[arg-type]
        index,
        index + 1,
        None,
        10,
        0 if failed else 10,
        0 if failed else 2,
        nan if failed else spearman,
        nan if failed else 1.0,
        nan if failed else 0.8,
        None if failed else ("a", "b"),
        error,
    )


def _result() -> StabilityResult:
    reseeds = summarise("reseed", 1.0, [_row("reseed", 0), _row("reseed", 1, "RuntimeError: boom")])
    subsamples = summarise(
        "subsample", 0.8, [_row("subsample", 0, spearman=nan), _row("subsample", 1), _row("subsample", 2)]
    )
    return StabilityResult(
        n_obs=10,
        reference_lineages=2,
        match_jaccard=0.6,
        reseeds=reseeds,
        subsamples=subsamples,
        n_failed=1,
        lineage_stability=0.5,
        pseudotime_stability=nan,
    )


def test_the_name_is_a_valid_trust_name_and_follows_the_parameters(monkeypatch: pytest.MonkeyPatch) -> None:
    name = stability_name()
    assert TRUST_NAME.fullmatch(name) and name.startswith("stability-")
    monkeypatch.setattr(trust_stability, "MATCH_JACCARD", 0.7)
    assert stability_name() != name


@pytest.mark.parametrize(
    ("constant", "value"),
    [
        ("N_RESEEDS", 7),
        ("N_SUBSAMPLES", 9),
        ("SUBSAMPLE_FRACTION", 0.5),
        ("STABILITY_SEED", 11),
        ("MATCH_JACCARD", 0.7),
    ],
)
def test_every_stability_parameter_changes_the_name(
    monkeypatch: pytest.MonkeyPatch, constant: str, value: float
) -> None:
    name = stability_name()
    monkeypatch.setattr(trust_stability, constant, value)
    assert stability_name() != name


def test_the_name_follows_the_metric_version(monkeypatch: pytest.MonkeyPatch) -> None:
    name = stability_name()
    monkeypatch.setattr(trust_stability, "STABILITY_METRIC_VERSION", 99)
    assert stability_name() != name


def _raise(constant: str) -> None:
    raise ValueError(f"not strict JSON: {constant}")


def test_undefined_values_are_null_and_the_file_is_strict_json(tmp_path: Path) -> None:
    api = to_api_stability(_result())
    assert api.status == "ok" and api.n_failed == 1
    assert api.pseudotime_stability is None and api.lineage_stability == 0.5
    assert api.subsamples is not None and api.reseeds is not None
    assert api.subsamples.replicates[0].spearman is None
    failed = api.reseeds.replicates[1]
    assert failed.error == "RuntimeError: boom"
    assert failed.spearman is None and failed.share_matched is None and failed.mean_jaccard is None
    assert failed.terminals is None

    store = _store(tmp_path)
    target = store.write_trust("r", stability_name(), api.model_dump(mode="json"))
    parsed = json.loads(target.read_text(encoding="utf-8"), parse_constant=_raise)
    assert parsed["pseudotime_stability"] is None
    assert parsed["reseeds"]["replicates"][1]["spearman"] is None
    assert models.Stability.model_validate(parsed) == api


def test_an_empty_kind_is_null_too(tmp_path: Path) -> None:
    empty = KindSummary("reseed", 0, 1.0, (), nan, nan, nan, nan, nan, nan, nan, 0, ())
    api = to_api_stability(StabilityResult(10, 2, 0.6, empty, empty, 0, nan, nan))
    assert api.reseeds is not None and api.reseeds.spearman_median is None
    _store(tmp_path).write_trust("r", stability_name(), api.model_dump(mode="json"))


def test_failed_rows_are_noted() -> None:
    notes = to_api_stability(_result()).notes
    assert notes == [STABILITY_NOTE, "reseed 2 failed: RuntimeError: boom"]


def test_not_measured_carries_its_reason(tmp_path: Path) -> None:
    got = not_measured("failed", "the run was prepared with other settings")
    assert got.status == "failed" and got.notes == ["the run was prepared with other settings"]
    assert got.reseeds is None and got.lineage_stability is None and got.n_failed == 0
    _store(tmp_path).write_trust("r", stability_name(), got.model_dump(mode="json"))


def _store(tmp_path: Path) -> RunStore:
    store = RunStore(tmp_path)
    (store.path("r")).mkdir(parents=True)
    (store.path("r") / ".complete").write_text("x")
    return store
