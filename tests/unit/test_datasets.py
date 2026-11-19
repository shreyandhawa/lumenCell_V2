import hashlib
import io
import re
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from lumencell.bench import download as dl
from lumencell.bench.datasets import (
    DATASETS,
    FetchError,
    RealDataset,
    fetch,
    is_cached,
    local_path,
    truth_for,
)
from lumencell.bench.download import Build, Source, build_dataset, fetch_source, source_path
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
        assert ds.licence.strip() and (ds.build is not None or len(ds.sources) == 1)
        for s in ds.sources:
            assert s.sha256 and len(s.sha256) == 64 and s.size > 0 and s.url.startswith("https://")
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


def test_single_file_datasets_keep_their_cache_paths(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))
    assert local_path(DATASETS["pancreas"]) == tmp_path / "pancreas_endocrinogenesis_day15.h5ad"


_GOOD: dict[str, Any] = {"url": "https://h/x", "filename": "x", "size": 1, "sha256": None}


@pytest.mark.parametrize(
    "bad",
    [{"filename": "../x"}, {"filename": "a/b"}, {"url": "http://h/x"}, {"size": 0}, {"sha256": "ABC"}],
)
def test_source_fields_are_checked(bad: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        Source(**(_GOOD | bad))


def test_source_accepts_good_fields() -> None:
    Source(**_GOOD)
    Source(**(_GOOD | {"sha256": "a" * 64}))


def _fake_single(payload: bytes, sha: str | None) -> RealDataset:
    src = Source(url="https://h/fake", filename="fake.h5ad", size=len(payload), sha256=sha)
    return RealDataset(**(_REAL_FIELDS | {"id": "fake", "sources": (src,)}))


_REAL_FIELDS: dict[str, Any] = {
    "id": "x",
    "title": "t",
    "citation": "c",
    "sources": (),
    "licence": "test licence",
    "label_key": "clusters",
    "root_label": "r",
    "lineages": {"A": ["r", "a"]},
    "gated": ["A"],
}


def test_fetch_checks_the_checksum_and_cleans_up(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))
    payload = b"not really an h5ad"
    monkeypatch.setattr(dl, "_open", lambda *a, **k: _FakeResponse(payload))
    with pytest.raises(FetchError, match="checksum"):
        fetch_source(_fake_single(payload, "0" * 64).sources[0])
    assert list(tmp_path.iterdir()) == [], "a failed download leaves nothing behind"
    good = _fake_single(payload, hashlib.sha256(payload).hexdigest())
    assert fetch(good).read_bytes() == payload


def test_fetch_checks_the_size(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))
    monkeypatch.setattr(dl, "_open", lambda *a, **k: _FakeResponse(b"short"))
    src = Source(url="https://h/fake", filename="fake.h5ad", size=100, sha256=None)
    with pytest.raises(FetchError, match="expected 100 bytes"):
        fetch_source(src)


PAYLOADS = {"alpha.txt": b"alpha", "beta.txt": b"beta" * 10}


class TwoFile:
    """A two-source dataset served by a fake urlopen, with counters."""

    def __init__(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))
        self.dir = tmp_path
        self.opens = 0
        self.builds = 0
        self.fail = False
        sources = tuple(
            Source(
                url=f"https://h/{name}",
                filename=name,
                size=len(data),
                sha256=hashlib.sha256(data).hexdigest(),
            )
            for name, data in PAYLOADS.items()
        )
        self.ds = RealDataset(
            **(
                _REAL_FIELDS
                | {
                    "id": "twofile",
                    "sources": sources,
                    "build": Build(builder=self._builder, version=1, params={"k": 1}),
                }
            )
        )

        def fake_open(request: Any) -> _FakeResponse:
            self.opens += 1
            return _FakeResponse(PAYLOADS[request.full_url.rsplit("/", 1)[1]])

        monkeypatch.setattr(dl, "_open", fake_open)

    def _builder(self, paths: Any, out: Path, params: Any) -> None:
        self.builds += 1
        a = ad.AnnData(np.ones((2, 2), np.float32))
        a.uns["seen"] = sorted(paths)
        a.write_h5ad(out)
        if self.fail:
            raise RuntimeError("builder broke")


@pytest.fixture
def two_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> TwoFile:
    return TwoFile(tmp_path, monkeypatch)


def test_fetch_downloads_every_source_then_builds(two_file: TwoFile) -> None:
    ds = two_file.ds
    assert fetch(ds) == local_path(ds)
    assert sorted(ad.read_h5ad(local_path(ds)).uns["seen"]) == sorted(PAYLOADS)
    for s in ds.sources:
        assert source_path(s).read_bytes() == PAYLOADS[s.filename]
    assert not [p for p in two_file.dir.iterdir() if ".tmp" in p.name]
    leftovers = {p.name for p in two_file.dir.iterdir()} - {local_path(ds).name, *PAYLOADS}
    assert not leftovers, "no temporary file or directory remains"
    assert is_cached(ds)


def test_downloaded_sources_are_not_fetched_again(two_file: TwoFile) -> None:
    fetch(two_file.ds)
    opens, builds = two_file.opens, two_file.builds
    fetch(two_file.ds)
    assert (two_file.opens, two_file.builds) == (opens, builds)


def test_a_failed_build_leaves_no_derived_file(two_file: TwoFile) -> None:
    two_file.fail = True
    with pytest.raises(RuntimeError, match="builder broke"):
        fetch(two_file.ds)
    assert not local_path(two_file.ds).exists()
    assert not list(two_file.dir.glob("*.tmp.h5ad"))
    assert not is_cached(two_file.ds)


def test_build_refuses_a_source_whose_checksum_changed(two_file: TwoFile) -> None:
    ds = two_file.ds
    fetch(ds)
    victim = source_path(ds.sources[0])
    victim.write_bytes(b"x" * len(PAYLOADS[victim.name]))
    local_path(ds).unlink()
    builds = two_file.builds
    with pytest.raises(
        FetchError, match=rf"checksum.*{re.escape(victim.name)}|{re.escape(victim.name)}.*checksum"
    ):
        fetch(ds)
    assert two_file.builds == builds


def test_build_needs_recorded_checksums(two_file: TwoFile) -> None:
    ds = two_file.ds
    fetch(ds)
    first = Source(ds.sources[0].url, ds.sources[0].filename, ds.sources[0].size, None)
    loose = RealDataset(
        **(_REAL_FIELDS | {"id": "twofile", "sources": (first, ds.sources[1]), "build": ds.build})
    )
    with pytest.raises(FetchError, match="no recorded checksum"):
        build_dataset(loose)


def test_the_derived_path_follows_the_build_inputs(two_file: TwoFile) -> None:
    ds = two_file.ds
    assert ds.build is not None
    base = local_path(ds)
    assert re.fullmatch(rf"{ds.id}-[0-9a-f]{{12}}\.h5ad", base.name)

    def with_build(**kw: Any) -> RealDataset:
        b = Build(**({"builder": ds.build.builder, "version": 1, "params": {"k": 1}} | kw))  # type: ignore[union-attr]
        return RealDataset(**(_REAL_FIELDS | {"id": ds.id, "sources": ds.sources, "build": b}))

    assert local_path(with_build()) == base
    assert local_path(with_build(params={"k": 2})) != base
    assert local_path(with_build(version=2)) != base
    other = Source(ds.sources[0].url, ds.sources[0].filename, ds.sources[0].size, "b" * 64)
    changed = RealDataset(
        **(_REAL_FIELDS | {"id": ds.id, "sources": (other, ds.sources[1]), "build": ds.build})
    )
    assert local_path(changed) != base


def test_dataset_rules() -> None:
    a = Source("https://h/a", "a", 1, None)
    with pytest.raises(ValueError):
        RealDataset(**(_REAL_FIELDS | {"sources": (a, Source("https://h/b", "b", 1, None))}))  # two, no build
    b = Build(builder=lambda *_: None, version=1, params={})
    with pytest.raises(ValueError):
        RealDataset(**(_REAL_FIELDS | {"sources": (a, a), "build": b}))  # duplicate filenames
    with pytest.raises(ValueError):
        RealDataset(**(_REAL_FIELDS | {"sources": (a,), "licence": " "}))


def test_each_download_gets_its_own_empty_directory(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))
    seen: list[list[str]] = []
    payload = b"abc"

    def fake_open(request: Any) -> _FakeResponse:
        (work,) = [p for p in tmp_path.iterdir() if p.is_dir()]
        seen.append([p.name for p in work.iterdir()])
        return _FakeResponse(payload)

    monkeypatch.setattr(dl, "_open", fake_open)
    for name in ("a.bin", "b.bin"):
        fetch_source(Source("https://h/" + name, name, len(payload), None))
    assert seen == [["a.bin"], ["b.bin"]], "each download's directory held only its own file"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.bin", "b.bin"]


def test_the_tls_context_verifies_certificates_with_certifi(monkeypatch: pytest.MonkeyPatch) -> None:
    import ssl

    import certifi

    calls: list[dict[str, Any]] = []
    real = ssl.create_default_context

    def spy(*a: Any, **k: Any) -> ssl.SSLContext:
        calls.append(k)
        return real(*a, **k)

    monkeypatch.setattr(dl.ssl, "create_default_context", spy)
    ctx = dl._tls_context()
    assert ctx.verify_mode == ssl.CERT_REQUIRED and ctx.check_hostname
    assert calls == [{"cafile": certifi.where()}]


def test_the_opener_has_the_redirect_guard(monkeypatch: pytest.MonkeyPatch) -> None:
    import urllib.request

    seen: list[Any] = []
    real = urllib.request.build_opener

    def spy(*handlers: Any) -> Any:
        seen.extend(handlers)
        raise OSError("stop before any network use")

    monkeypatch.setattr(dl.urllib.request, "build_opener", spy)
    with pytest.raises(OSError, match="stop before"):
        dl._open(urllib.request.Request("https://h/x"))
    assert real is not spy
    assert any(isinstance(h, dl._HttpsOnlyRedirect) for h in seen)
    assert any(isinstance(h, urllib.request.HTTPSHandler) for h in seen)


def test_a_redirect_to_a_non_https_address_is_refused() -> None:
    import urllib.request

    handler = dl._HttpsOnlyRedirect()
    req = urllib.request.Request("https://h/x")
    for target in ("http://evil/x", "ftp://evil/x", "file:///etc/passwd"):
        with pytest.raises(FetchError, match="non-https"):
            handler.redirect_request(req, None, 302, "Found", {}, target)  # type: ignore[arg-type]
    other_host = handler.redirect_request(req, None, 302, "Found", {}, "https://s3.example/x")  # type: ignore[arg-type]
    assert other_host is not None and other_host.full_url == "https://s3.example/x"


def test_fetch_refuses_a_redirect_end_to_end(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))

    def refuse(request: Any) -> Any:
        raise FetchError("refusing a redirect to a non-https address: 'http://evil/x'")

    monkeypatch.setattr(dl, "_open", refuse)
    with pytest.raises(FetchError, match="non-https"):
        fetch_source(Source("https://h/a", "a", 3, None))
    assert list(tmp_path.iterdir()) == []


def test_an_empty_build_is_refused_and_leaves_no_file(two_file: TwoFile) -> None:
    ds = two_file.ds
    assert ds.build is not None
    empty = Build(builder=lambda paths, out, params: None, version=1, params={"k": 1})
    from dataclasses import replace

    ds = replace(ds, build=empty)
    with pytest.raises(FetchError, match="empty"):
        fetch(ds)
    assert not local_path(ds).exists()
    assert not list(two_file.dir.glob("*.tmp.h5ad"))
    assert not is_cached(ds)


def test_build_dataset_names_a_source_that_is_not_downloaded(two_file: TwoFile) -> None:
    with pytest.raises(FetchError, match="not downloaded"):
        build_dataset(two_file.ds)
    assert two_file.builds == 0


def test_is_downloaded_rejects_a_wrong_size_and_fetch_downloads_again(two_file: TwoFile) -> None:
    from lumencell.bench.download import is_downloaded

    src = two_file.ds.sources[0]
    source_path(src).write_bytes(b"x" * (src.size + 3))
    assert not is_downloaded(src)
    fetch_source(src)
    assert two_file.opens == 1 and source_path(src).read_bytes() == PAYLOADS[src.filename]


def test_a_cached_check_needs_every_source_even_when_the_derived_file_exists(two_file: TwoFile) -> None:
    ds = two_file.ds
    fetch(ds)
    assert is_cached(ds)
    source_path(ds.sources[0]).unlink()
    assert local_path(ds).exists() and not is_cached(ds)
