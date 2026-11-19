"""Inspect the two Daniel 2022 metadata files without reading any receptor value (phase 3 Task 3).

Run from the repo root: ``uv run python -I benchmarks/inspection/daniel2022_metadata.py --dir <datasets dir>
[--labels COL --samples COL] [--compare]``. It reads the files it is given and writes nothing.

A column is withheld when its name matches ``RECEPTOR_COLUMN`` or contains, case-insensitively, one of
``WITHHELD_PARTS``. A withheld column is listed by name and dtype and dropped before anything else is
computed, so none of its values is printed, summarised or used.

``--compare`` (with ``--labels`` and ``--samples``) reads only the index and those two columns of the main
and the D5 file, and prints how many cells they share and the main-by-D5 label cross-tab of the shared cells.
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

from lumencell.bench.builders.daniel2022 import SOURCES
from lumencell.clonal.isolation import RECEPTOR_COLUMN

WITHHELD_PARTS = ("clon", "cdr3", "tcr", "tra", "trb", "chain", "vdj", "_gene")
CORE = re.compile(r"[ACGT]{16}")
MAX_LEVELS = 200


def withheld(name: str) -> bool:
    return bool(RECEPTOR_COLUMN.match(name)) or any(p in name.lower() for p in WITHHELD_PARTS)


def tag(index_value: str) -> str:
    m = CORE.search(index_value)
    return f"{index_value[: m.start()]}<core>{index_value[m.end() :]}" if m else "<no 16-nt core>"


def require_columns(columns: pd.Index, path: Path, *names: str | None) -> None:
    for name in names:
        if name not in columns:
            sys.exit(f"error: {name!r} is not a column of {path.name}")


def compare(main_path: Path, d5_path: Path, labels: str, samples: str) -> None:
    """Index overlap and label cross-tab of two files, reading only the index and two columns of each."""
    frames = []
    for path in (main_path, d5_path):
        header = pd.read_csv(path, sep="\t", nrows=0)
        require_columns(header.columns, path, labels, samples)
        frames.append(pd.read_csv(path, sep="\t", usecols=[labels, samples]))
    main, d5 = frames
    shared = main.index.intersection(d5.index)
    print(f"== comparison {main_path.name} vs {d5_path.name}")
    print(f"cells in both\t{len(shared)}")
    print(f"only in main\t{len(main.index.difference(d5.index))}")
    print(f"only in D5\t{len(d5.index.difference(main.index))}")
    print(f"main {labels} (rows) x D5 {labels} (columns), shared cells:")
    table = pd.crosstab(main.loc[shared, labels], d5.loc[shared, labels], dropna=False)
    print(table.to_csv(sep="\t"), end="")


def inspect(path: Path, labels: str | None, samples: str | None) -> None:
    # The GEO files are Seurat exports: the header has no field for the index, so pandas takes the
    # first column of each row as the index.
    df = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"== {path.name}")
    print(f"rows\t{len(df)}")
    print("first index values\t" + "\t".join(str(v) for v in df.index[:3]))
    print("index tag vocabulary (tag, count):")
    for t, n in pd.Series([tag(str(v)) for v in df.index]).value_counts().sort_index().items():
        print(f"  {t}\t{n}")
    print("columns (name, dtype, withheld):")
    for c in df.columns:
        print(f"  {c}\t{df[c].dtype}\t{'withheld' if withheld(str(c)) else '-'}")
    df = df[[c for c in df.columns if not withheld(str(c))]]
    for c in df.columns:
        col = df[c]
        if pd.api.types.is_numeric_dtype(col):
            print(f"numeric {c}: min {col.min()}  median {col.median()}  max {col.max()}")
        n_levels = col.nunique(dropna=False)
        if n_levels <= MAX_LEVELS:
            print(f"value counts {c} ({n_levels} distinct):")
            for v, n in col.value_counts(dropna=False).sort_index().items():
                print(f"  {v}\t{n}")
        elif not pd.api.types.is_numeric_dtype(col):
            print(f"{c}: {n_levels} distinct values, not listed")
    if labels is not None and samples is not None:
        require_columns(df.columns, path, labels, samples)
        print(f"cross-tab {samples} (rows) x {labels} (columns):")
        print(pd.crosstab(df[samples], df[labels], dropna=False).to_csv(sep="\t"), end="")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--labels")
    ap.add_argument("--samples")
    ap.add_argument("--compare", action="store_true")
    a = ap.parse_args()
    if (a.labels is None) != (a.samples is None):
        ap.error("--labels and --samples go together")
    for name in (a.labels, a.samples):
        if name is not None and withheld(name):
            ap.error(f"{name!r} is a withheld column and cannot be cross-tabulated")
    if a.compare and a.labels is None:
        ap.error("--compare needs --labels and --samples")
    files = [a.dir / src.filename for src in SOURCES if "_metadata." in src.filename]
    if a.compare:
        main_path, d5_path = sorted(files, key=lambda p: "_D5_" in p.name)
        compare(main_path, d5_path, a.labels, a.samples)
        return 0
    for path in files:
        inspect(path, a.labels, a.samples)
    return 0


if __name__ == "__main__":
    sys.exit(main())
