import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "benchmarks" / "inspection" / "barcode_join.py"
A, B, C, D = "AAACCTGAGCGTCTAT", "AAAGATGAGAGTAAGG", "CCCCCCCCCCCCCCCC", "GGGGGGGGGGGGGGGG"


def run_raw(tmp_path: Path, cell_text: str, table: str, suffix: str) -> subprocess.CompletedProcess:
    cell_file = tmp_path / "cells.txt"
    cell_file.write_bytes(cell_text.encode())
    table_file = tmp_path / "table.csv"
    table_file.write_text(table)
    return subprocess.run(
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
    )


def run(tmp_path: Path, cells: list[str], table: str, suffix: str) -> str:
    out = run_raw(tmp_path, "\n".join(cells) + "\n", table, suffix)
    assert out.returncode == 0, out.stderr
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


def test_comma_separated_lists_are_split(tmp_path):
    table = f'id,barcodes\n1,"{A}-1,{B}-1"\n2,"{C}-1,{D}-1"\n'
    out = run(tmp_path, [f"{A}-3", f"{B}-3"], table, "-3")
    assert "barcode column: barcodes" in out
    assert "distinct barcodes 4" in out
    assert "share 0.50" in out


def test_padded_cell_values_are_trimmed(tmp_path):
    table = f'barcode\n"  {A}-1 "\n"{B}-1  "\n'
    out = run(tmp_path, [f"{A}-2", f"{B}-2"], table, "-2")
    assert "barcode column: barcode" in out
    assert "share 1.00" in out


def test_cells_file_with_crlf_trailing_spaces_and_quotes(tmp_path):
    table = f"barcode\n{A}-1\n{B}-1\n{C}-1\n"
    cells = f'"{A}-2"  \r\n" {B}-2 " \r\n  {C}-2\r\n\r\n'
    out = run_raw(tmp_path, cells, table, "-2")
    assert out.returncode == 0, out.stderr
    assert "share 1.00" in out.stdout


def test_column_at_exactly_half_barcodes_is_not_a_barcode_column(tmp_path):
    table = f"mixed\n{A}-1\nnot-a-barcode\n"
    out = run_raw(tmp_path, f"{A}-2\n", table, "-2")
    assert out.returncode == 1
    assert "no column" in out.stderr


def test_sixty_percent_column_is_chosen_and_twenty_percent_is_not(tmp_path):
    table = f"sixty,twenty\n{A}-1,{D}-1\n{B}-1,x\n{C}-1,x\nx,x\nx,y\n"
    out = run(tmp_path, [f"{A}-2", f"{B}-2", f"{C}-2"], table, "-2")
    assert "barcode column: sixty" in out
    assert "distinct barcodes 3" in out


def test_twenty_percent_column_alone_is_rejected(tmp_path):
    table = f"twenty\n{A}-1\nx\ny\nz\nw\n"
    out = run_raw(tmp_path, f"{A}-2\n", table, "-2")
    assert out.returncode == 1


def test_the_column_with_more_barcodes_wins(tmp_path):
    table = f"few,many\n{A}-1,{B}-1\n{A}-1,{C}-1\n{A}-1,{D}-1\n"
    out = run(tmp_path, [f"{B}-2"], table, "-2")
    assert "barcode column: many" in out
    assert "distinct barcodes 3" in out


def test_bare_16nt_barcodes_do_not_join_but_are_counted(tmp_path):
    table = f"barcode\n{A}\n{B}-1\n{C}-2\n"
    out = run(tmp_path, [f"{A}-3", f"{B}-3", f"{C}-3"], table, "-3")
    assert "distinct barcodes 3" in out
    assert "not ending in -1 (do not join): 2" in out
    assert "share 0.33" in out


def test_duplicate_barcodes_across_chains_count_once(tmp_path):
    table = f"chain,barcodes\nTRA,{A}-1;{B}-1\nTRB,{A}-1;{B}-1\nTRB,{A}-1\n"
    out = run(tmp_path, [f"{A}-2"], table, "-2")
    assert "distinct barcodes 2" in out
    assert "share 0.50" in out


def test_no_barcode_column_exits_non_zero_with_a_message(tmp_path):
    out = run_raw(tmp_path, f"{A}-2\n", "a,b\n1,2\n3,4\n", "-2")
    assert out.returncode == 1
    assert "no column holds cell barcodes" in out.stderr


def test_spaces_around_list_separators_are_trimmed(tmp_path):
    table = f'barcodes\n"{A}-1 ; {B}-1"\n"{C}-1 , {D}-1"\n'
    out = run(tmp_path, [f"{A}-2", f"{B}-2", f"{C}-2", f"{D}-2"], table, "-2")
    assert "distinct barcodes 4" in out
    assert "share 1.00" in out
