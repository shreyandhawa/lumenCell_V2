"""Inspect the two Daniel 2022 metadata files without reading any receptor value (phase 3 Task 3).

Run from the repo root: ``uv run python -I benchmarks/inspection/daniel2022_metadata.py --dir <datasets dir>
[--labels COL --samples COL]``. It reads the files it is given and writes nothing.

A column is withheld when its name matches ``RECEPTOR_COLUMN`` or contains, case-insensitively, one of
``WITHHELD_PARTS``. A withheld column is listed by name and dtype and dropped before anything else is
computed, so none of its values is printed, summarised or used.
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
        if labels not in df.columns or samples not in df.columns:
            print(f"cross-tab skipped: {labels!r} or {samples!r} is not a column of {path.name}")
            return
        print(f"cross-tab {samples} (rows) x {labels} (columns):")
        print(pd.crosstab(df[samples], df[labels], dropna=False).to_csv(sep="\t"), end="")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--labels")
    ap.add_argument("--samples")
    a = ap.parse_args()
    if (a.labels is None) != (a.samples is None):
        ap.error("--labels and --samples go together")
    for name in (a.labels, a.samples):
        if name is not None and withheld(name):
            ap.error(f"{name!r} is a withheld column and cannot be cross-tabulated")
    for src in SOURCES:
        if "_metadata." in src.filename:
            inspect(a.dir / src.filename, a.labels, a.samples)
    return 0


if __name__ == "__main__":
    sys.exit(main())
