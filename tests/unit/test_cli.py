from pathlib import Path

import anndata as ad
from typer.testing import CliRunner

from lumencell import __version__
from lumencell.cli import app

runner = CliRunner()


def test_help_lists_commands() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "version" in result.output


def test_no_args_shows_help() -> None:
    result = runner.invoke(app, [])
    assert "Usage" in result.output


def test_version_prints_package_version() -> None:
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.output.strip() == __version__


def test_synth_writes_an_h5ad(tmp_path: Path) -> None:
    out = tmp_path / "data" / "s.h5ad"
    result = runner.invoke(
        app, ["synth", str(out), "--n-cells", "200", "--n-genes", "50", "--topology", "binary2"]
    )
    assert result.exit_code == 0, result.output
    adata = ad.read_h5ad(out)
    assert adata.shape == (200, 50)
    assert adata.uns["lumencell_truth"]["topology"] == "binary2"


def test_synth_rejects_unknown_topology(tmp_path: Path) -> None:
    result = runner.invoke(app, ["synth", str(tmp_path / "x.h5ad"), "--topology", "spiral"])
    assert result.exit_code != 0


def test_launch_reports_a_bad_time_limit_setting(tmp_path: Path, monkeypatch: object) -> None:
    import pytest

    mp = monkeypatch
    assert isinstance(mp, pytest.MonkeyPatch)
    mp.setenv("LUMENCELL_JOB_TIME_LIMITS", "prepare=soon")
    data = tmp_path / "d.h5ad"
    data.write_bytes(b"")
    result = runner.invoke(app, ["launch", str(data), "--no-browser"])
    assert result.exit_code != 0
    assert "LUMENCELL_JOB_TIME_LIMITS" in result.output
    assert isinstance(result.exception, SystemExit), "a message, not a traceback"
