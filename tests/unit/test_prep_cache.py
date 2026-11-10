import shutil
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


# ---- damaged stages ----------------------------------------------------------------------


def _artifact(cache: PrepCache, stage: str, filename: str) -> Path:
    key = stage_keys(cache.ds, PrepParams())[stage]
    return cache.stages_dir / f"{stage}-{key}" / filename


def _assert_same_artifacts(got: dict[str, object], want: dict[str, object]) -> None:
    assert set(got) == set(want)
    for name, value in want.items():
        g, w = (got[name].toarray(), value.toarray()) if hasattr(value, "toarray") else (got[name], value)  # type: ignore[attr-defined]
        g, w = np.asarray(g), np.asarray(w)
        if np.issubdtype(w.dtype, np.floating):
            np.testing.assert_allclose(g, w, atol=1e-5, err_msg=name)
        else:
            np.testing.assert_array_equal(g, w, err_msg=name)


def test_a_stage_with_a_missing_artifact_counts_as_not_cached(
    ds: Dataset, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    clean = PrepCache(ds, tmp_path).run().artifacts
    artifact = _artifact(PrepCache(ds, tmp_path), "pca", "pca.npy")
    artifact.unlink()
    cache = PrepCache(ds, tmp_path)
    assert cache.status(PrepParams())["pca"] == "pending"
    capsys.readouterr()
    healed = cache.run()
    err = capsys.readouterr().err
    assert "pca" in err and str(artifact.parent) in err, "replacing the incomplete copy is reported"
    assert not [p.name for p in cache.stages_dir.iterdir() if p.name.startswith(".")]
    states = {r.name: r.state for r in healed.stages}
    assert states.pop("pca") == "done"
    assert set(states.values()) == {"cached"}, "only the damaged stage is recomputed"
    _assert_same_artifacts(healed.artifacts, clean)
    again = PrepCache(ds, tmp_path).run()
    assert {r.state for r in again.stages} == {"cached"}, "the healed stage is on disk"


@pytest.mark.parametrize(
    ("stage", "filename", "keep"),
    [
        ("pca", "pca.npy", 0.5),
        ("pca", "pca.npy", 0.0),
        ("normalize", "lognorm.npz", 0.5),
        ("neighbors", "connectivities.npz", 0.5),
        ("neighbors", "connectivities.npz", 0.0),
    ],
)
def test_a_truncated_artifact_is_recomputed(
    ds: Dataset, tmp_path: Path, capsys: pytest.CaptureFixture[str], stage: str, filename: str, keep: float
) -> None:
    clean = PrepCache(ds, tmp_path).run().artifacts
    artifact = _artifact(PrepCache(ds, tmp_path), stage, filename)
    artifact.write_bytes(artifact.read_bytes()[: int(artifact.stat().st_size * keep)])
    capsys.readouterr()
    healed = PrepCache(ds, tmp_path).run()
    states = {r.name: r.state for r in healed.stages}
    assert states.pop(stage) == "done"
    assert set(states.values()) == {"cached"}, "only the damaged stage is recomputed"
    _assert_same_artifacts(healed.artifacts, clean)
    err = capsys.readouterr().err
    assert stage in err and str(artifact.parent) in err, "the damage is reported, not hidden"
    assert not [p.name for p in PrepCache(ds, tmp_path).stages_dir.iterdir() if p.name.startswith(".")]
    again = PrepCache(ds, tmp_path).run()
    assert {r.state for r in again.stages} == {"cached"}, "the healed stage is on disk"


def test_load_all_never_computes(ds: Dataset, tmp_path: Path) -> None:
    cache = PrepCache(ds, tmp_path)
    with pytest.raises(cache_mod.CacheMiss):
        cache.load_all(PrepParams())
    assert not cache.stages_dir.exists() or not any(cache.stages_dir.iterdir())


def test_load_all_on_a_damaged_stage_takes_it_out_and_writes_nothing(ds: Dataset, tmp_path: Path) -> None:
    cache = PrepCache(ds, tmp_path)
    cache.run()
    artifact = _artifact(cache, "pca", "pca.npy")
    artifact.write_bytes(artifact.read_bytes()[:100])
    others = sorted(p.name for p in cache.stages_dir.iterdir() if p != artifact.parent)
    with pytest.raises(cache_mod.CacheMiss):
        PrepCache(ds, tmp_path).load_all(PrepParams())
    assert not artifact.parent.exists()
    assert sorted(p.name for p in cache.stages_dir.iterdir()) == others


def test_load_all_writes_nothing_to_the_cache(ds: Dataset, tmp_path: Path) -> None:
    cache = PrepCache(ds, tmp_path)
    cache.run()
    (cache.dir / "manifest.json").unlink()
    PrepCache(ds, tmp_path).load_all(PrepParams())
    assert not (cache.dir / "manifest.json").exists(), "the server only reads the cache"


def test_an_artifact_that_cannot_be_reached_is_not_deleted(ds: Dataset, tmp_path: Path) -> None:
    # A permission or I/O error says nothing about the file being damaged: keep the stage.
    PrepCache(ds, tmp_path).run()
    stage_dir = _artifact(PrepCache(ds, tmp_path), "pca", "pca.npy").parent
    denied = PermissionError(13, "Permission denied")
    with mock.patch.object(cache_mod, "_load_artifact", side_effect=denied), pytest.raises(PermissionError):
        PrepCache(ds, tmp_path).run()
    assert (stage_dir / "pca.npy").exists()


def test_a_stage_that_cannot_be_deleted_still_counts_as_taken_out(
    ds: Dataset, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    cache = PrepCache(ds, tmp_path)
    cache.run()
    artifact = _artifact(cache, "pca", "pca.npy")
    artifact.write_bytes(b"")
    capsys.readouterr()
    busy = OSError("device busy")
    with mock.patch.object(cache_mod.shutil, "rmtree", side_effect=busy), pytest.raises(cache_mod.CacheMiss):
        PrepCache(ds, tmp_path).load_all(PrepParams())
    assert not artifact.parent.exists(), "renamed out of the way even though deleting it failed"
    assert "device busy" in capsys.readouterr().err, "the leftover is reported"


def test_a_writer_that_loses_the_race_keeps_the_winners_copy(ds: Dataset, tmp_path: Path) -> None:
    cache = PrepCache(ds, tmp_path)
    cache.run()
    artifact = _artifact(cache, "pca", "pca.npy")
    final = artifact.parent
    winner = tmp_path / "winner"
    shutil.copytree(final, winner)
    artifact.unlink()  # this run has to recompute pca
    real_rename = Path.rename

    def racing_rename(self: Path, target: str | Path) -> Path:
        if Path(target) == final and self.name.startswith(".pca-") and not final.exists():
            real_rename(winner, final)  # another writer lands its complete copy first
        return real_rename(self, target)

    with mock.patch.object(Path, "rename", racing_rename):
        res = PrepCache(ds, tmp_path).run()
    assert {r.name: r.state for r in res.stages}["pca"] == "done"
    assert (final / "pca.npy").exists()
    assert not [p.name for p in cache.stages_dir.iterdir() if p.name.startswith(".")]
