"""Inspect the Daniel 2022 matrix, barcode, feature and contig formats without reading a receptor value
(phase 3 Task 5a).

Run from the repo root: ``uv run python -I benchmarks/inspection/daniel2022_formats.py --dir <datasets dir>
[--capture-map <json: tar member -> barcode tag>]``. It reads the files it is given and writes nothing; the
tar is read in memory and never extracted.

Freeze rule: of a contig member only the header line and the ``READ_COLUMNS`` are read (an allow-list, so
a column the script has never seen is withheld too). Every other column is listed by name only; none of its
values is read, printed or used.

A merged barcode is ``<tag><16-nt core>-<n>``: its tag is everything before the last 16-nt core. A contig
barcode joins a tag when ``tag + barcode`` is a merged barcode of either matrix. With ``--capture-map`` each
mapped member gets one line ``join <member> <own tag> <hits>/<distinct> <share> <best other tag>
<hits>/<distinct> <share>``; the best other tag has the most hits (ties go to the first in sorted order).
"""

import argparse
import gzip
import io
import json
import re
import sys
import tarfile
import zlib
from collections import Counter
from pathlib import Path

import pandas as pd

from lumencell.bench.builders.common import BuildError, read_members
from lumencell.bench.builders.daniel2022 import SOURCES

MATRICES = {"main": "GSE188666_scRNA_LCMV", "d5": "GSE188666_scRNA_D5_LCMV"}
TAR = "GSE188666_RAW.tar"
FLAGS = ("chain", "productive", "high_confidence", "full_length", "is_cell")
READ_COLUMNS = ("barcode", *FLAGS)
N_VALUES = 100_000
MAX_LEVELS = 50
MERGED = re.compile(r"(.*)([ACGT]{16}-\d+)")
INTEGER = re.compile(r"[+-]?\d+")

_KNOWN = {s.filename for s in SOURCES}
assert TAR in _KNOWN and all(
    f"{stem}{end}" in _KNOWN
    for stem in MATRICES.values()
    for end in (".mtx.gz", "_barcodes.tsv.gz", "_features.tsv.gz")
)


def _lines(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n") for line in fh]


def matrix_head(path: Path) -> tuple[str, str, int, bool]:
    """Banner, size line, and whether the first N_VALUES stored values are integer text (and how many)."""
    with gzip.open(path, "rt") as fh:
        banner = fh.readline().rstrip("\n")
        size = ""
        for line in fh:
            if not line.startswith("%"):
                size = line.strip()
                break
        n, integer = 0, True
        for line in fh:
            if n == N_VALUES:
                break
            parts = line.split()
            integer &= len(parts) == 3 and bool(INTEGER.fullmatch(parts[2]))
            n += 1
    return banner, size, n, integer


def orientation(size: str, n_features: int, n_barcodes: int) -> str:
    dims = [int(x) for x in size.split()[:2]] if len(size.split()) == 3 else []
    if dims == [n_features, n_barcodes]:
        return "genes_x_cells"
    if dims == [n_barcodes, n_features]:
        return "cells_x_genes"
    return "neither"


def tag_of(barcode: str) -> str | None:
    m = MERGED.fullmatch(barcode)
    return m.group(1) if m else None


def inspect_matrices(d: Path) -> tuple[set[str], list[str]]:
    """Print each matrix's head and vocabularies; return the merged barcodes of both and the sorted tags."""
    merged: set[str] = set()
    tags: set[str] = set()
    features: dict[str, list[str]] = {}
    per_file: dict[str, set[str]] = {}
    for key, stem in MATRICES.items():
        mtx = d / f"{stem}.mtx.gz"
        barcodes = _lines(d / f"{stem}_barcodes.tsv.gz")
        feats = _lines(d / f"{stem}_features.tsv.gz")
        features[key] = feats
        banner, size, n, integer = matrix_head(mtx)
        print(f"== matrix {key} {mtx.name}")
        print(f"banner\t{banner}")
        print(f"size line\t{size}")
        print(f"barcodes\t{len(barcodes)}\tfirst\t" + "\t".join(barcodes[:3]))
        print(f"features\t{len(feats)}\tfirst\t" + "\t".join(feats[:3]))
        fields = sorted({len(f.split("\t")) for f in feats})
        print("feature fields per line\t" + ",".join(str(x) for x in fields))
        print(f"orientation\t{orientation(size, len(feats), len(barcodes))}")
        print(f"first {n} stored values integer\t{integer}")
        counts = Counter(tag_of(b) for b in barcodes)
        print("tag vocabulary (tag, count):")
        for t in sorted(counts, key=lambda x: (x is None, x or "")):
            print(f"  {t if t is not None else '<no 16-nt core>'}\t{counts[t]}")
        suffixes = Counter(b.rsplit("-", 1)[1] for b in barcodes if tag_of(b) is not None)
        print("barcode suffixes\t" + "\t".join(f"-{s}:{c}" for s, c in sorted(suffixes.items())))
        per_file[key] = set(barcodes)
        merged.update(barcodes)
        tags.update(t for t in counts if t is not None)
    print(f"D5 features equal main features\t{features['d5'] == features['main']}")
    print(f"D5 barcodes also among main barcodes\t{len(per_file['d5'] & per_file['main'])}")
    return merged, sorted(tags)


def _text(data: bytes, position: int) -> str:
    """Decode a member; bad gzip or UTF-8 is an error that names its position, never its bytes."""
    try:
        return (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode()
    except (OSError, EOFError, zlib.error, ValueError):
        raise BuildError(f"contig member {position} is not valid gzip or UTF-8 text") from None


def inspect_contigs(d: Path) -> dict[str, list[str]]:
    """Print each member's header and the pooled flag vocabularies; return each member's distinct barcodes."""
    tar = d / TAR
    # tarfile only lists the member names here; read_members reads the bytes, in memory.
    with tarfile.open(tar, "r:") as tf:
        names = [m.name for m in tf if m.isreg()]
    data, _ = read_members(tar, names)
    # Validate every member before printing anything from one: a header line that is not a header may be a
    # data row, so nothing of a member (its name or its first line) is shown until it is known to be a header.
    texts: dict[str, str] = {}
    for i, name in enumerate(names, 1):
        texts[name] = _text(data[name], i)
        if "barcode" not in texts[name].split("\n", 1)[0].rstrip("\r").split(","):
            raise BuildError(f"contig member {i} of {len(names)} has no barcode column in its header")
    print(f"== contigs {tar.name}")
    print(f"members\t{len(names)}")
    pooled: dict[str, Counter[str]] = {c: Counter() for c in FLAGS}
    absent: Counter[str] = Counter()
    withheld: set[str] = set()
    barcodes: dict[str, list[str]] = {}
    for name in names:
        text = texts[name]
        header = text.split("\n", 1)[0].rstrip("\r")
        columns = header.split(",")
        print(f"member\t{name}")
        print(f"header\t{header}")
        withheld.update(c for c in columns if c not in READ_COLUMNS)
        use = [c for c in READ_COLUMNS if c in columns]
        df = pd.read_csv(io.StringIO(text), usecols=use, dtype=str, keep_default_na=False)
        for c in FLAGS:
            if c in df:
                pooled[c].update(df[c])
            else:
                absent[c] += 1
        barcodes[name] = sorted(set(df["barcode"]))
    print("withheld columns (by name, never read)\t" + ",".join(sorted(withheld)))
    for c in FLAGS:
        print(f"distinct values {c} (pooled over {len(names) - absent[c]} members):")
        if len(pooled[c]) > MAX_LEVELS:
            print(f"  {len(pooled[c])} distinct values, not listed")
        else:
            for v, n in sorted(pooled[c].items()):
                print(f"  {v}\t{n}")
    if names:
        first = pd.read_csv(io.StringIO(texts[names[0]]), usecols=["barcode"], dtype=str, nrows=50)
        print("first member first barcodes\t" + "\t".join(list(dict.fromkeys(first["barcode"]))[:3]))
    return barcodes


def read_map(path: Path) -> dict[str, str]:
    try:
        m = json.loads(path.read_text())
    except (OSError, ValueError):
        sys.exit("error: the capture map is not a readable JSON file")
    if not isinstance(m, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in m.items()):
        sys.exit("error: the capture map must be a JSON object of tar member -> barcode tag strings")
    return m


def _share(hits: int, n: int) -> str:
    return f"{hits / n:.2f}" if n else "n/a"


def join(
    capture_map: dict[str, str], members: dict[str, list[str]], merged: set[str], tags: list[str]
) -> None:
    unknown = [k for k in capture_map if k not in members]
    if unknown:
        sys.exit(f"error: {len(unknown)} capture-map member(s) are not in the tar")
    print(f"tar members not in the capture map\t{sum(m not in capture_map for m in members)}")
    print("join (member, own tag, hits/distinct, share, best other tag, hits/distinct, share):")
    for member in sorted(capture_map):
        bcs = members[member]
        n = len(bcs)

        def hits(tag: str, bcs: list[str] = bcs) -> int:
            return sum(tag + b in merged for b in bcs)

        own = capture_map[member]
        h_own = hits(own)
        others = [(hits(t), t) for t in tags if t != own]
        h_best, best = min(others, key=lambda x: (-x[0], x[1])) if others else (0, "-")
        mine, other = f"{h_own}/{n}\t{_share(h_own, n)}", f"{h_best}/{n}\t{_share(h_best, n)}"
        print(f"join\t{member}\t{own}\t{mine}\t{best}\t{other}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--capture-map", type=Path)
    a = ap.parse_args()
    capture_map = read_map(a.capture_map) if a.capture_map is not None else None
    merged, tags = inspect_matrices(a.dir)
    try:
        members = inspect_contigs(a.dir)
    except BuildError as e:
        sys.exit(f"error: {e}")
    if capture_map is not None:
        join(capture_map, members, merged, tags)
    return 0


if __name__ == "__main__":
    sys.exit(main())
