import hashlib
import io
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lumencell.bench import datasets as dsmod
from lumencell.bench.datasets import DATASETS, FetchError, fetch, truth_for
from lumencell.bench.metrics import score


def _pancreas_obs() -> pd.DataFrame:
    labels = (
        ["Ductal"] * 4
        + ["Ngn3 low EP"] * 2
        + ["Pre-endocrine"] * 2
        + ["Alpha"] * 3
        + ["Beta"] * 3
        + ["Delta"] * 2
    )
    return pd.DataFrame({"clusters": pd.Categorical(labels)})


def test_registry_is_complete() -> None:
    for ds in DATASETS.values():
        assert ds.sha256 and len(ds.sha256) == 64
        assert ds.size > 0 and ds.url.startswith("https://")
        assert set(ds.gated) <= set(ds.lineages)
        assert all(stages[0] == ds.root_label for stages in ds.lineages.values())


def test_pancreas_truth_uses_stage_order_and_gated_lineages() -> None:
    obs = _pancreas_obs()
    truth = truth_for(DATASETS["pancreas"], obs)
    assert truth.names == ["Alpha", "Beta"], "Delta and Epsilon are reported, not gated"
    assert list(truth.order[:4]) == [0, 0, 0, 0] and truth.order[8] == 4
    assert truth.root.sum() == 4
    assert truth.terminals is not None and truth.terminals["Alpha"].sum() == 3
    assert not truth.lineages["Alpha"][14], "Delta cells are on no gated lineage"


def test_time_order_for_gastrulation() -> None:
    obs = pd.DataFrame(
        {"celltype": ["Blood progenitors 1", "Erythroid3", "Erythroid1"], "stage": ["E7.0", "E8.5", "E7.75"]}
    )
    truth = truth_for(DATASETS["gastrulation_erythroid"], obs)
    np.testing.assert_allclose(truth.order, [7.0, 8.5, 7.75])


def test_terminal_share() -> None:
    truth = truth_for(DATASETS["pancreas"], _pancreas_obs())
    n = len(truth.order)
    pt = np.nan_to_num(truth.order, nan=0.0) / 4
    alpha_cells = truth.terminals["Alpha"]  # type: ignore[index]
    half = alpha_cells.copy()
    half[np.flatnonzero(alpha_cells)[:1]] = False  # one of three Alpha cells elsewhere
    membership = np.stack([half | (truth.order <= 3), truth.lineages["Beta"]])
    s = score(pt, membership[:, :n], truth)
    assert s.terminal_share == pytest.approx(2 / 3)


class _FakeResponse(io.BytesIO):
    def __enter__(self) -> "_FakeResponse":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()


def _fake_dataset(payload: bytes, sha: str | None) -> dsmod.RealDataset:
    return replace(DATASETS["pancreas"], filename="fake.h5ad", size=len(payload), sha256=sha)


def test_fetch_checks_the_checksum_and_cleans_up(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))
    payload = b"not really an h5ad"
    monkeypatch.setattr(dsmod.urllib.request, "urlopen", lambda *a, **k: _FakeResponse(payload))
    with pytest.raises(FetchError, match="checksum"):
        fetch(_fake_dataset(payload, "0" * 64))
    assert list(tmp_path.iterdir()) == [], "a failed download leaves nothing behind"
    good = _fake_dataset(payload, hashlib.sha256(payload).hexdigest())
    assert fetch(good).read_bytes() == payload


def test_fetch_checks_the_size(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))
    monkeypatch.setattr(dsmod.urllib.request, "urlopen", lambda *a, **k: _FakeResponse(b"short"))
    ds = replace(_fake_dataset(b"x" * 100, None), size=100)
    with pytest.raises(FetchError, match="expected 100 bytes"):
        fetch(ds)
