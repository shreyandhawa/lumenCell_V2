import time
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

from lumencell.bench.synth import generate
from lumencell.data.loader import Dataset, load
from lumencell.prep import cache as cache_mod
from lumencell.prep.cache import PrepCache, stage_keys
from lumencell.prep.pipeline import run_pipeline
from lumencell.prep.stages import STAGES, PrepParams


@pytest.fixture(scope="module")
def ds(tmp_path_factory: pytest.TempPathFactory) -> Dataset:
    path = tmp_path_factory.mktemp("cache-ds") / "d.h5ad"
    generate("bifurcation", "realistic", n_cells=600, n_genes=400, seed=0).write_h5ad(path)
    return load(path)


def test_first_run_computes_and_matches_the_in_memory_pipeline(ds: Dataset, tmp_path: Path) -> None:
    res = PrepCache(ds, tmp_path).run()
    assert all(r.state == "done" for r in res.stages)
    ref = run_pipeline(ds).artifacts
    np.testing.assert_array_equal(res.artifacts["leiden"], ref["leiden"])
    np.testing.assert_allclose(res.artifacts["pca"], ref["pca"], atol=1e-5)


def test_second_run_is_a_fast_cache_hit_with_identical_artifacts(ds: Dataset, tmp_path: Path) -> None:
    first = PrepCache(ds, tmp_path).run().artifacts
    start = time.perf_counter()
    second = PrepCache(ds, tmp_path).run()
    assert time.perf_counter() - start < 2.0
    assert all(r.state == "cached" for r in second.stages)
    for name, value in first.items():
        got = second.artifacts[name]
        if hasattr(value, "toarray"):
            np.testing.assert_array_equal(got.toarray(), value.toarray())
        else:
            np.testing.assert_array_equal(np.asarray(got), np.asarray(value))


def test_changing_a_parameter_recomputes_only_downstream_stages(ds: Dataset, tmp_path: Path) -> None:
    PrepCache(ds, tmp_path).run()
    res = PrepCache(ds, tmp_path).run(PrepParams(n_top_genes=300))
    states = {r.name: r.state for r in res.stages}
    assert states["normalize"] == "cached"
    assert all(states[s] == "done" for s in ("hvg", "pca", "neighbors", "leiden", "diffmap", "paga", "umap"))


def test_keys_follow_the_dependency_graph(ds: Dataset) -> None:
    base = stage_keys(ds, PrepParams())
    leiden_changed = stage_keys(ds, PrepParams(leiden_resolution=2.0))
    assert base["neighbors"] == leiden_changed["neighbors"]
    assert base["diffmap"] == leiden_changed["diffmap"], "diffmap does not use Leiden"
    assert base["leiden"] != leiden_changed["leiden"]
    assert base["paga"] != leiden_changed["paga"]


def test_status_reports_cached_and_pending_stages(ds: Dataset, tmp_path: Path) -> None:
    cache = PrepCache(ds, tmp_path)
    assert set(cache.status(PrepParams()).values()) == {"pending"}
    cache.run()
    assert set(cache.status(PrepParams()).values()) == {"done"}
    assert cache.manifest() is not None


def test_cache_hits_keep_the_compute_timings(ds: Dataset, tmp_path: Path) -> None:
    cache = PrepCache(ds, tmp_path)
    cache.run()
    computed = cache.manifest()
    assert computed is not None
    assert all(s["state"] == "done" for s in computed["stages"].values())
    PrepCache(ds, tmp_path).run()  # every stage is now a cache hit
    assert cache.manifest() == computed, "loading from the cache must not overwrite the timings"
    assert not [p.name for p in cache.dir.iterdir() if p.suffix == ".tmp"]


def test_an_unreadable_manifest_is_ignored_and_rewritten(ds: Dataset, tmp_path: Path) -> None:
    cache = PrepCache(ds, tmp_path)
    cache.run()
    (cache.dir / "manifest.json").write_text("{broken", encoding="utf-8")
    assert cache.manifest() is None
    cache.run()
    rewritten = cache.manifest()
    assert rewritten is not None
    assert {s["seconds"] for s in rewritten["stages"].values()} == {None}, "load times are not recorded"


def test_a_failed_stage_leaves_no_partial_artifacts(ds: Dataset, tmp_path: Path) -> None:
    cache = PrepCache(ds, tmp_path)
    boom = RuntimeError("disk full")
    with mock.patch.object(cache_mod, "_save_artifact", side_effect=boom), pytest.raises(RuntimeError):
        cache.run()
    assert not any(p.name.startswith(".") for p in cache.stages_dir.iterdir()), "temp dirs are removed"
    assert set(cache.status(PrepParams()).values()) == {"pending"}


def test_incomplete_stage_directory_is_not_trusted(ds: Dataset, tmp_path: Path) -> None:
    cache = PrepCache(ds, tmp_path)
    cache.run()
    keys = stage_keys(ds, PrepParams())
    marker = cache.stages_dir / f"pca-{keys['pca']}" / ".complete"
    marker.unlink()
    assert cache.status(PrepParams())["pca"] == "pending"


def test_cache_root_can_be_overridden(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("LUMENCELL_CACHE_DIR", str(tmp_path))
    assert cache_mod.default_cache_root() == tmp_path


def test_every_stage_has_a_key(ds: Dataset) -> None:
    assert set(stage_keys(ds, PrepParams())) == {s.name for s in STAGES}
