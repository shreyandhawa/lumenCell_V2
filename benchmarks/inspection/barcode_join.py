"""Check that a clonotype table's cell barcodes join to an expression matrix's cell IDs.

Prints headers, dtypes, row count, the barcode column, three example barcodes and the share of
distinct table barcodes found among the cell IDs. Mapping rule: a table barcode ending in `-1` becomes
`<16-nt core><suffix>`; a barcode with any other ending (including a bare 16-nt one) does not join, is
counted in the denominator and is reported. Exits 1 when no column holds barcodes.
Prints no receptor value: only column names, dtypes, counts and barcodes.
Run as: uv run python -I benchmarks/inspection/barcode_join.py ...
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

BARCODE = re.compile(r"[ACGT]{16}(-\d+)?")
SPLIT = re.compile(r"[;,]")


def barcodes_in(values: pd.Series) -> list[str]:
    """Every barcode in a column, splitting `;`/`,` lists; [] unless over half the cells are barcodes."""
    cells = [str(v).strip() for v in values.dropna()]
    parts = [[p.strip() for p in SPLIT.split(c) if p.strip()] for c in cells]
    nonempty = [p for p in parts if p]
    if not nonempty or sum(all(BARCODE.fullmatch(x) for x in p) for p in nonempty) <= 0.5 * len(nonempty):
        return []
    return [x for p in nonempty for x in p if BARCODE.fullmatch(x)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--table", type=Path, required=True)
    ap.add_argument("--suffix", required=True)
    args = ap.parse_args()

    cells = {line.strip().strip('"').strip() for line in args.cells.read_text().splitlines() if line.strip()}
    table = pd.read_csv(args.table, dtype=str, keep_default_na=False, na_values=[""])
    print(f"rows {len(table)}, columns {len(table.columns)}")
    for col in table.columns:
        print(f"column {col!r} dtype {table[col].dtype}")

    best_col, best = None, []
    for col in table.columns:
        found = barcodes_in(table[col])
        if len(set(found)) > len(set(best)):
            best_col, best = col, found
    if best_col is None:
        print("error: no column holds cell barcodes", file=sys.stderr)
        sys.exit(1)
    distinct = sorted(set(best))
    joinable = [b for b in distinct if b.endswith("-1")]
    print(f"barcode column: {best_col}")
    print(f"distinct barcodes {len(distinct)}")
    print(f"not ending in -1 (do not join): {len(distinct) - len(joinable)}")
    print("examples " + ", ".join(distinct[:3]))
    mapped = {b[:16] + args.suffix for b in joinable}
    hit = len(mapped & cells)
    print(
        f"mapped to 16-nt core + {args.suffix!r}: {hit} of {len(distinct)} distinct barcodes "
        f"found among {len(cells)} cell IDs"
    )
    print(f"share {hit / len(distinct):.2f}")


if __name__ == "__main__":
    main()
