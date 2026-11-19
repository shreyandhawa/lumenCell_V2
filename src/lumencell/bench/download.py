"""Verified downloads and derived datasets for the real-data benchmark.

A dataset is one or more ``Source`` files, each checked against a recorded size
and SHA-256, plus an optional build step that turns them into one derived
``.h5ad`` in the cache. Nothing here is called implicitly: ``lumencell bench
fetch`` shows the sources and asks first. Sources are never modified, and every
file is written to a temporary name and renamed into place. Each download goes into
its own fresh, empty directory first. Redirects must stay on https.
"""

import hashlib
import json
import os
import re
import shutil
import ssl
import tempfile
import urllib.request
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

import certifi
import platformdirs

USER_AGENT = "lumencell-bench"
SAFE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,200}")
_SHA256 = re.compile(r"[0-9a-f]{64}")


class FetchError(RuntimeError):
    pass


def datasets_dir() -> Path:
    override = os.environ.get("LUMENCELL_DATASETS_DIR")
    return Path(override) if override else Path(platformdirs.user_cache_dir("lumencell")) / "datasets"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass(frozen=True)
class Source:
    url: str
    filename: str
    size: int
    sha256: str | None

    def __post_init__(self) -> None:
        if not self.url.startswith("https://"):
            raise ValueError(f"source url must start with https://: {self.url!r}")
        if not SAFE_NAME.fullmatch(self.filename):
            raise ValueError(f"unsafe source filename {self.filename!r}")
        if self.size <= 0:
            raise ValueError(f"{self.filename}: size must be positive")
        if self.sha256 is not None and not _SHA256.fullmatch(self.sha256):
            raise ValueError(f"{self.filename}: sha256 must be 64 lowercase hex characters")


# Receives every source's local path keyed by Source.filename, the temporary
# path to write the .h5ad to, and the declared params.
Builder = Callable[[Mapping[str, Path], Path, Mapping[str, Any]], None]


@dataclass(frozen=True)
class Build:
    builder: Builder
    version: int
    params: Mapping[str, Any]  # JSON-serialisable


class Fetchable(Protocol):
    @property
    def id(self) -> str: ...
    @property
    def sources(self) -> tuple[Source, ...]: ...
    @property
    def build(self) -> Build | None: ...


def source_path(src: Source) -> Path:
    return datasets_dir() / src.filename


def is_downloaded(src: Source) -> bool:
    p = source_path(src)
    return p.exists() and p.stat().st_size == src.size


class _HttpsOnlyRedirect(urllib.request.HTTPRedirectHandler):
    """Follow redirects (hosts may change, e.g. figshare to S3) but never to a non-https target."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[no-untyped-def]
        if not newurl.lower().startswith("https://"):
            raise FetchError(f"refusing a redirect to a non-https address: {newurl!r}")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _tls_context() -> ssl.SSLContext:
    # certifi's CA bundle: some Python builds (e.g. python.org on macOS) ship without CA certificates.
    return ssl.create_default_context(cafile=certifi.where())


def _open(request: urllib.request.Request) -> Any:
    opener = urllib.request.build_opener(
        urllib.request.HTTPSHandler(context=_tls_context()), _HttpsOnlyRedirect()
    )
    return opener.open(request, timeout=60)


def fetch_source(src: Source, progress: Callable[[int, int], None] | None = None) -> Path:
    """Download one file (the caller has already asked the user). Returns its path.

    The bytes go into a fresh, empty directory made for this download, are checked there,
    and only then moved into place; the directory is always removed.
    """
    target = source_path(src)
    if is_downloaded(src):
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(src.url, headers={"User-Agent": USER_AGENT})
    workdir = Path(tempfile.mkdtemp(prefix=".download-", dir=target.parent))
    try:
        tmp = workdir / src.filename
        h = hashlib.sha256()
        done = 0
        with tmp.open("wb") as out, _open(request) as res:
            while chunk := res.read(1 << 20):
                out.write(chunk)
                h.update(chunk)
                done += len(chunk)
                if progress:
                    progress(done, src.size)
        if done != src.size:
            raise FetchError(f"{src.filename}: expected {src.size:,} bytes, got {done:,}")
        digest = h.hexdigest()
        if src.sha256 and digest != src.sha256:
            raise FetchError(f"{src.filename}: checksum mismatch ({digest} != {src.sha256})")
        os.replace(tmp, target)
        return target
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def build_key(item: Fetchable) -> str:
    """Twelve hex characters naming a derived file by everything that determines it."""
    b = item.build
    if b is None:
        raise ValueError(f"{item.id} has no build step")
    payload = {
        "id": item.id,
        "version": b.version,
        "params": b.params,
        "sources": [s.sha256 for s in item.sources],
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:12]


def local_path(item: Fetchable) -> Path:
    if item.build is None:
        return source_path(item.sources[0])
    return datasets_dir() / f"{item.id}-{build_key(item)}.h5ad"


def is_cached(item: Fetchable) -> bool:
    return all(is_downloaded(s) for s in item.sources) and local_path(item).exists()


def build_dataset(item: Fetchable) -> Path:
    """Run the build step on verified sources and rename the result into the cache."""
    b = item.build
    if b is None:
        raise ValueError(f"{item.id} has no build step")
    paths: dict[str, Path] = {}
    for s in item.sources:
        if not is_downloaded(s):
            raise FetchError(f"{item.id}: {s.filename} is not downloaded")
        if s.sha256 is None:
            raise FetchError(f"{item.id}: {s.filename} has no recorded checksum")
        p = source_path(s)
        if sha256_of(p) != s.sha256:
            raise FetchError(f"{item.id}: {s.filename} does not match its recorded checksum")
        paths[s.filename] = p
    target = local_path(item)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{item.id}-", suffix=".tmp.h5ad", dir=datasets_dir())
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        b.builder(paths, tmp, b.params)
        if tmp.stat().st_size == 0:
            raise FetchError(f"{item.id}: the build wrote an empty file")
        os.replace(tmp, target)
        return target
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def fetch(item: Fetchable, progress: Callable[[int, int], None] | None = None) -> Path:
    """Download what is missing, then build when there is a build step and no derived file."""
    for s in item.sources:
        if not is_downloaded(s):
            fetch_source(s, progress)
    if item.build is not None and not local_path(item).exists():
        build_dataset(item)
    return local_path(item)
