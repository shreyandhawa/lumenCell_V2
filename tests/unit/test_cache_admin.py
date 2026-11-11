"""`lumencell cache list` and `lumencell cache prune`, on a fake cache under tmp_path.

Never run against the real cache: every test points LUMENCELL_CACHE_DIR and
LUMENCELL_DATA_DIR at tmp_path.
"""

import json
import os
import time
from datetime import timedelta
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

from lumencell.cli import app
from lumencell.prep.housekeeping import Removal, UnsafePath, apply_prune, list_cache, parse_age

runner = CliRunner()
DAY = 86400.0
FP_A, FP_B, FP_C, FP_D = "a" * 24, "b" * 24, "c" * 24, "d" * 24
MANIFEST_A = json.dumps({"dataset": {"name": "pancreas", "path": "/data/pancreas.h5ad"}, "stages": {}})
MANIFEST_B = json.dumps({"dataset": {"name": "bone_marrow", "path": "/data/bm.h5ad"}, "stages": {}})
MANIFEST_C = json.dumps({"pipeline": "1", "stages": {}})  # written before manifests named the dataset


def _write(path: Path, nbytes: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * nbytes)
    return path


def _age(path: Path, seconds: float, now: float) -> None:
    os.utime(path, (now - seconds, now - seconds), follow_symlinks=False)


def _age_tree(root: Path, seconds: float, now: float) -> None:
    for p in [*root.rglob("*"), root]:
        _age(p, seconds, now)


@pytest.fixture
def cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    root, state = tmp_path / "cache", tmp_path / "state"
    monkeypatch.setenv("LUMENCELL_CACHE_DIR", str(root))
    monkeypatch.setenv("LUMENCELL_DATA_DIR", str(state))
    now = time.time()
    v1 = root / "v1"

    a = v1 / FP_A  # last used 40 days ago
    _write(a / "stages" / "pca-0123456789abcdef" / "pca.npy", 1000)
    _write(a / "runs" / "r1" / "cells.lcb", 500)
    (a / "manifest.json").write_text(MANIFEST_A, encoding="utf-8")

    b = v1 / FP_B  # in use, with leftovers from cancelled jobs and one temp folder being written now
    _write(b / "stages" / "hvg-0123456789abcdef" / "hvg.npy", 200)
    (b / "manifest.json").write_text(MANIFEST_B, encoding="utf-8")
    old_stage_tmp = _write(b / "stages" / ".pca-k2j3h4x9" / "pca.npy", 300).parent
    old_manifest_tmp = _write(b / ".manifest-x1y2z3w4.tmp", 50)
    _write(b / "runs" / "r2" / "result.json", 10)
    old_run_tmp = _write(b / "runs" / "r2" / "result.json.tmp", 20)
    fresh_tmp = _write(b / "stages" / ".hvg-busy1234" / "hvg.npy", 40).parent

    c = v1 / FP_C  # last used 5 days ago, old-style manifest
    _write(c / "stages" / "normalize-0123456789abcdef" / "lognorm.npz", 100)
    (c / "manifest.json").write_text(MANIFEST_C, encoding="utf-8")

    finder = [_write(d / ".DS_Store", 6) for d in (b, b / "stages", b / "runs")]  # not ours, never removed
    not_a_dataset = _write(v1 / ("e" * 24) / "notes.txt", 8).parent  # a 24-hex name alone is not enough

    datasets = _write(root / "datasets" / "pancreas.h5ad", 70)  # benchmark downloads live here
    state_file = _write(state / f"{FP_A}.json", 5)

    for p in (old_stage_tmp, old_manifest_tmp, old_run_tmp, *finder):
        _age(p, 2 * 3600, now)
    _age_tree(a, 40 * DAY, now)
    _age_tree(c, 5 * DAY, now)
    return SimpleNamespace(
        root=root, v1=v1, a=a, b=b, c=c, now=now, datasets=datasets, state_file=state_file,
        old_temps=[old_stage_tmp, old_manifest_tmp, old_run_tmp], fresh_tmp=fresh_tmp,
        finder=finder, not_a_dataset=not_a_dataset,
    )  # fmt: skip


def _snapshot(base: Path) -> dict[str, tuple[int, float]]:
    return {str(p): (p.lstat().st_size, p.lstat().st_mtime) for p in base.rglob("*")}


def test_list_reports_each_dataset_with_size_and_last_use(cache: SimpleNamespace) -> None:
    entries = list_cache(cache.root)
    assert [e.fingerprint for e in entries] == [FP_B, FP_C, FP_A], "most recently used first"
    by = {e.fingerprint: e for e in entries}
    assert by[FP_A].name == "pancreas" and by[FP_A].source == "/data/pancreas.h5ad"
    assert by[FP_A].size == 1000 + 500 + len(MANIFEST_A)
    assert by[FP_B].size == 200 + len(MANIFEST_B) + 300 + 50 + 10 + 20 + 40 + 3 * 6
    assert by[FP_C].name is None
    assert abs(by[FP_A].last_used - (cache.now - 40 * DAY)) < 5

    res = runner.invoke(app, ["cache", "list"])
    assert res.exit_code == 0, res.output
    for text in ("pancreas", "bone_marrow", FP_C, "40 days ago", "5 days ago", str(cache.root)):
        assert text in res.output, text


def test_prune_older_than_removes_old_datasets_and_leftover_temps(cache: SimpleNamespace) -> None:
    res = runner.invoke(app, ["cache", "prune", "--older-than", "30d"])
    assert res.exit_code == 0, res.output
    assert not cache.a.exists()
    assert cache.b.exists() and cache.c.exists()
    assert not any(p.exists() for p in cache.old_temps), "leftovers from cancelled jobs are removed"
    assert cache.fresh_tmp.exists(), "a temp folder a job may be writing right now is kept"
    assert all(p.exists() for p in cache.finder), "only our own temporary names are leftovers"
    assert (cache.b / "runs" / "r2" / "result.json").exists()
    assert (cache.b / "stages" / "hvg-0123456789abcdef" / "hvg.npy").exists()
    assert cache.datasets.exists() and cache.state_file.exists()
    assert "removed" in res.output and "pancreas" in res.output


@pytest.mark.parametrize("flags", [["--all"], ["--older-than", "30d"]])
def test_prune_dry_run_removes_nothing(cache: SimpleNamespace, tmp_path: Path, flags: list[str]) -> None:
    before = _snapshot(tmp_path)
    res = runner.invoke(app, ["cache", "prune", *flags, "--dry-run"])
    assert res.exit_code == 0, res.output
    assert _snapshot(tmp_path) == before
    assert "would remove" in res.output and "pancreas" in res.output


def test_prune_all_removes_every_dataset_and_nothing_else(cache: SimpleNamespace) -> None:
    res = runner.invoke(app, ["cache", "prune", "--all"])
    assert res.exit_code == 0, res.output
    assert [p.name for p in cache.v1.iterdir()] == ["e" * 24], "only cache datasets are removed"
    assert cache.datasets.read_bytes() == b"x" * 70
    assert cache.state_file.read_bytes() == b"x" * 5


def test_prune_refuses_a_path_outside_the_cache_root(cache: SimpleNamespace, tmp_path: Path) -> None:
    outside = _write(tmp_path / "outside" / "keep.npy", 10).parent
    _write(outside / "manifest.json", 2)
    (cache.v1 / FP_D).symlink_to(outside, target_is_directory=True)
    res = runner.invoke(app, ["cache", "prune", "--all"])
    assert res.exit_code == 0, res.output
    assert "outside the cache" in res.output, "the refusal is reported"
    assert (outside / "keep.npy").exists() and (cache.v1 / FP_D).is_symlink()
    assert not cache.a.exists(), "one refused entry does not block the rest"


def test_removing_leftovers_keeps_a_datasets_last_use(cache: SimpleNamespace) -> None:
    leftover = _write(cache.c / ".manifest-a1b2c3d4.tmp", 5)  # directly in the dataset folder
    _age(leftover, 2 * 3600, cache.now)
    _age(cache.c, 5 * DAY, cache.now)
    assert runner.invoke(app, ["cache", "prune", "--older-than", "30d"]).exit_code == 0
    assert not leftover.exists()
    (c,) = [e for e in list_cache(cache.root) if e.fingerprint == FP_C]
    assert abs(c.last_used - (cache.now - 5 * DAY)) < 5, "cleaning up is not using it"


def test_prune_reports_what_it_could_not_remove(cache: SimpleNamespace) -> None:
    locked = cache.a / "runs"
    locked.chmod(0o500)  # its files cannot be deleted
    try:
        res = runner.invoke(app, ["cache", "prune", "--all"])
    finally:
        locked.chmod(0o700)
    assert res.exit_code == 1
    assert "could not remove" in res.output and "pancreas" in res.output
    assert not cache.b.exists() and not cache.c.exists(), "the others are still removed"
    assert "removed dataset bone_marrow" in res.output


def test_the_state_directory_is_never_removed(
    cache: SimpleNamespace, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("LUMENCELL_DATA_DIR", str(cache.root))  # state files next to v1/: fine
    assert runner.invoke(app, ["cache", "prune", "--all"]).exit_code == 0
    assert not cache.a.exists()
    d = cache.v1 / FP_D
    _write(d / "manifest.json", 2)
    inside = _write(d / "state" / f"{FP_D}.json", 5).parent
    monkeypatch.setenv("LUMENCELL_DATA_DIR", str(inside))  # state inside a dataset: refuse it
    res = runner.invoke(app, ["cache", "prune", "--all"])
    assert "saved state" in res.output
    assert (inside / f"{FP_D}.json").exists()


def test_apply_prune_only_removes_dataset_entries_under_the_cache(
    cache: SimpleNamespace, tmp_path: Path
) -> None:
    outside = _write(tmp_path / "outside" / "keep.npy", 10)
    for target in (outside, cache.datasets, cache.root / "v1", cache.state_file, cache.not_a_dataset):
        with pytest.raises(UnsafePath):
            apply_prune(cache.root, [Removal(target, 10, "test")])
        assert target.exists()


def test_prune_needs_older_than_or_all(cache: SimpleNamespace) -> None:
    assert runner.invoke(app, ["cache", "prune"]).exit_code != 0
    assert runner.invoke(app, ["cache", "prune", "--all", "--older-than", "1d"]).exit_code != 0
    assert runner.invoke(app, ["cache", "prune", "--older-than", "soon"]).exit_code != 0
    assert cache.a.exists()


@pytest.mark.parametrize(
    ("text", "age"),
    [("30d", timedelta(days=30)), ("12h", timedelta(hours=12)), ("2w", timedelta(weeks=2))],
)
def test_parse_age(text: str, age: timedelta) -> None:
    assert parse_age(text) == age


@pytest.mark.parametrize("text", ["30", "1y", "d", "0d", "-1d", "1.5d", ""])
def test_parse_age_refuses_anything_else(text: str) -> None:
    with pytest.raises(ValueError):
        parse_age(text)


def test_cache_commands_on_a_missing_cache_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LUMENCELL_CACHE_DIR", str(tmp_path / "nothing-here"))
    for args in (["cache", "list"], ["cache", "prune", "--all"]):
        res = runner.invoke(app, args)
        assert res.exit_code == 0, res.output
        assert "empty" in res.output
    assert not (tmp_path / "nothing-here").exists(), "looking does not create it"


def test_listing_survives_files_that_vanish_meanwhile(
    cache: SimpleNamespace, monkeypatch: pytest.MonkeyPatch
) -> None:
    from lumencell.prep import housekeeping

    real_walk = os.walk

    def walk_with_a_ghost(top: str, *args: object, **kw: object) -> object:
        for dirpath, dirs, files in real_walk(top, *args, **kw):  # type: ignore[call-overload]
            yield dirpath, dirs, [*files, "renamed-meanwhile.npy"]

    monkeypatch.setattr(housekeeping.os, "walk", walk_with_a_ghost)
    assert {e.fingerprint for e in list_cache(cache.root)} == {FP_A, FP_B, FP_C}


def test_a_cache_from_before_the_stamp_counts_its_newest_file_as_last_use(cache: SimpleNamespace) -> None:
    # Prepared 40 days ago, a run saved yesterday, and no last-used stamp (written before it existed).
    e = cache.v1 / ("f" * 24)
    (e / "manifest.json").parent.mkdir(parents=True)
    (e / "manifest.json").write_text(MANIFEST_C, encoding="utf-8")
    _write(e / "stages" / "pca-0123456789abcdef" / "pca.npy", 10)
    _age_tree(e, 40 * DAY, cache.now)
    run = _write(e / "runs" / "r9" / "result.json", 10)
    _age(run, DAY, cache.now)
    _age(e, 40 * DAY, cache.now)
    (entry,) = [x for x in list_cache(cache.root) if x.fingerprint == "f" * 24]
    assert abs(entry.last_used - (cache.now - DAY)) < 5
    assert runner.invoke(app, ["cache", "prune", "--older-than", "30d"]).exit_code == 0
    assert run.exists(), "a dataset used yesterday is kept"
