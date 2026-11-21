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
