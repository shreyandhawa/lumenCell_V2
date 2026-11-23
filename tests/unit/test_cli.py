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


def test_synth_can_add_clones(tmp_path: Path) -> None:
    base = ["--n-cells", "300", "--n-genes", "50"]
    plain, biased = tmp_path / "p.h5ad", tmp_path / "b.h5ad"
    assert runner.invoke(app, ["synth", str(plain), *base]).exit_code == 0
    result = runner.invoke(app, ["synth", str(biased), *base, "--fate-bias", "0.7"])
    assert result.exit_code == 0, result.output
    assert {"clone_id", "donor"} <= set(ad.read_h5ad(biased).obs.columns)
    assert not {"clone_id", "donor"} & set(ad.read_h5ad(plain).obs.columns)
    assert runner.invoke(app, ["synth", str(plain), "--fate-bias", "1.5"]).exit_code == 2


def test_synth_fate_bias_reaches_add_clones(tmp_path: Path, monkeypatch: object) -> None:
    import pytest

    import lumencell.bench.synth_clones as sc

    assert isinstance(monkeypatch, pytest.MonkeyPatch)
    seen: list[tuple[float, int]] = []
    real = sc.add_clones

    def spy(adata, layer, seed):  # type: ignore[no-untyped-def]
        seen.append((layer.fate_bias, seed))
        return real(adata, layer, seed)

    monkeypatch.setattr(sc, "add_clones", spy)
    base = ["--n-cells", "300", "--n-genes", "50"]
    out = tmp_path / "z.h5ad"
    zero = runner.invoke(app, ["synth", str(out), *base, "--fate-bias", "0"])
    assert zero.exit_code == 0 and {"clone_id", "donor"} <= set(ad.read_h5ad(out).obs.columns)
    assert runner.invoke(app, ["synth", str(out), *base, "--fate-bias", "0.7", "--seed", "3"]).exit_code == 0
    assert seen == [(0.0, 0), (0.7, 3)]
    assert runner.invoke(app, ["synth", str(out), "--fate-bias", "-0.1"]).exit_code == 2
