import json
import subprocess
import zipfile
from pathlib import Path

import pytest

WEB_DIST = Path(__file__).parent.parent / "src" / "lumencell" / "web_dist"


@pytest.mark.skipif(not (WEB_DIST / "index.html").exists(), reason="run `make web-dist` first")
def test_wheel_contains_the_web_app(tmp_path: Path) -> None:
    # Plain `uv build` builds the wheel from the sdist, where git-ignored files drop out.
    subprocess.run(["uv", "build", "-o", str(tmp_path)], check=True, capture_output=True)
    (wheel,) = tmp_path.glob("lumencell-*.whl")
    assert "lumencell/web_dist/index.html" in zipfile.ZipFile(wheel).namelist()


def test_installed_wheel_fills_the_reference_rows(tmp_path: Path) -> None:
    # Installed without dependencies: the reader imports only the standard library, so this checks
    # the rows the server reads from the wheel, not the repository's copy of the file.
    subprocess.run(["uv", "build", "-o", str(tmp_path)], check=True, capture_output=True)
    (wheel,) = tmp_path.glob("lumencell-*.whl")
    venv = tmp_path / "venv"
    subprocess.run(["uv", "venv", "-q", str(venv)], check=True, capture_output=True)
    python = venv / "bin" / "python"
    subprocess.run(
        ["uv", "pip", "install", "-q", "--python", str(python), "--no-deps", str(wheel)],
        check=True,
        capture_output=True,
    )
    code = (
        "import json, lumencell.bench.reference_data as data;"
        "print(data.__file__);"
        "print(json.dumps(data.rows_for('paga_dpt')))"
    )
    path, rows = subprocess.run(
        [str(python), "-c", code], cwd=tmp_path, check=True, capture_output=True, text=True
    ).stdout.split("\n", 1)
    assert Path(path).is_relative_to(venv), "the test must read the installed copy, not the repository"
    rows_json = json.loads(rows)
    assert any(r["dataset"].startswith("Bone marrow") and r["passed"] is False for r in rows_json)


def test_installed_wheel_loads_the_policy_files(tmp_path: Path) -> None:
    # The synthetic gates and baseline ship in the package; a wheel without them breaks `bench synth`
    # and the research repository that imports lumencell.bench.gates.
    subprocess.run(["uv", "build", "-o", str(tmp_path)], check=True, capture_output=True)
    (wheel,) = tmp_path.glob("lumencell-*.whl")
    names = zipfile.ZipFile(wheel).namelist()
    assert "lumencell/bench/policy/thresholds.toml" in names
    assert "lumencell/bench/policy/baseline.json" in names
    venv = tmp_path / "venv"
    subprocess.run(["uv", "venv", "-q", str(venv)], check=True, capture_output=True)
    python = venv / "bin" / "python"
    subprocess.run(
        ["uv", "pip", "install", "-q", "--python", str(python), "--no-deps", str(wheel)],
        check=True,
        capture_output=True,
    )
    code = (
        "import lumencell.bench.gates as g;"
        "t = g.load_thresholds(); b = g.load_baseline();"
        "print(g.THRESHOLDS); print(bool(t['synthetic']['gates'])); print(bool(b))"
    )
    path, gates, baseline = subprocess.run(
        [str(python), "-c", code], cwd=tmp_path, check=True, capture_output=True, text=True
    ).stdout.split("\n")[:3]
    assert Path(path).is_relative_to(venv), "the test must read the installed copy, not the repository"
    assert gates == baseline == "True"
