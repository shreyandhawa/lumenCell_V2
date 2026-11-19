import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "benchmarks" / "inspection" / "barcode_join.py"
A, B, C, D = "AAACCTGAGCGTCTAT", "AAAGATGAGAGTAAGG", "CCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGGG"


def run(tmp_path: Path, cells: list[str], table: str, suffix: str) -> str:
    cell_file = tmp_path / "cells.txt"
    cell_file.write_text("\n".join(cells) + "\n")
    table_file = tmp_path / "table.csv"
    table_file.write_text(table)
    out = subprocess.run(
        [
            sys.executable,
            "-I",
            str(SCRIPT),
            "--cells",
            str(cell_file),
            "--table",
            str(table_file),
            "--suffix",
            suffix,
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    return out.stdout


def test_join_share_by_suffix(tmp_path):
    table = f"barcode,n\n{A}-1,1\n{B}-1,2\n{C}-1,3\n"
    out = run(tmp_path, [f"{A}-2", f"{B}-2", f"{A}-1"], table, "-2")
    assert "share 0.67" in out


def test_barcode_lists_in_one_column_are_split(tmp_path):
    table = f'id,barcodes,n\n1,"{A}-1;{B}-1",2\n2,"{C}-1;{D}-1",2\n'
    out = run(tmp_path, [f"{A}-3", f"{B}-3"], table, "-3")
    assert "barcode column: barcodes" in out
    assert "distinct barcodes 4" in out
    assert "share 0.50" in out


def test_no_receptor_values_are_printed(tmp_path):
    table = f"barcode,cdr3,v_gene\n{A}-1,CASSPLANTEDF,TRBV-PLANTED\n{B}-1,CASSPLANTEDF,TRBV-PLANTED\n"
    out = run(tmp_path, [f"{A}-1"], table, "-1")
    assert "CASSPLANTEDF" not in out
    assert "TRBV-PLANTED" not in out
    assert "barcode column: barcode" in out
