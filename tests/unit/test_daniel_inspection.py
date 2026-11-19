"""The Daniel 2022 inspection scripts list receptor columns by name and never print their values."""

import gzip
import re
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = str(Path(__file__).resolve().parents[2] / "benchmarks" / "inspection" / "daniel2022_metadata.py")

HEADER = ["orig.ident", "nCount_RNA", "celltype", "clonotype_id", "cdr3s_aa", "TRB_v_gene", "tissue"]


def _write(path: Path, sample: str, cores: list[str]) -> None:
    # The GEO layout: a Seurat export whose header has no field for the index.
    lines = ["\t".join(HEADER)]
    for i, core in enumerate(cores):
        label = "Texprog_fixture" if i % 2 == 0 else "Texterm_fixture"
        row = [f"{sample}_gex_{core}-1", sample, str(1000 + i), label]
        lines.append("\t".join([*row, "PLANTEDCLONE", "CASSPLANTEDF", "TRBV-PLANTED", "Spleen"]))
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


@pytest.fixture
def fixture_dir(tmp_path: Path) -> Path:
    cores = ["AAACCTGAGGCTAGAC", "AAACCTGAGTACGTTC", "AAACGGGAGCTGAAAT"]
    _write(tmp_path / "GSE188666_scRNA_LCMV_metadata.tsv.gz", "Spleen_Chronic_D8", cores * 2)
    _write(tmp_path / "GSE188666_scRNA_D5_LCMV_metadata.txt.gz", "Spleen_Chronic_D5", cores * 2)
    return tmp_path


def run(d: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-I", SCRIPT, "--dir", str(d), "--labels", "celltype", "--samples", "orig.ident"],
        capture_output=True,
        text=True,
        check=True,
    )


def test_the_inspection_never_prints_receptor_values(fixture_dir: Path) -> None:
    out = run(fixture_dir).stdout
    assert "PLANTED" not in out
    assert all(c in out for c in ("clonotype_id", "cdr3s_aa", "TRB_v_gene"))  # names are listed
    assert "Texprog_fixture" in out and "orig.ident" in out  # labels and samples are listed


def test_the_inspection_writes_nothing(fixture_dir: Path) -> None:
    before = sorted(fixture_dir.iterdir())
    run(fixture_dir)
    assert sorted(fixture_dir.iterdir()) == before


def test_the_index_tag_vocabulary_is_counted(fixture_dir: Path) -> None:
    out = run(fixture_dir).stdout
    assert "Spleen_Chronic_D8_gex_<core>-1\t6" in out


def test_a_withheld_column_is_refused_as_labels_or_samples(fixture_dir: Path) -> None:
    for flag in ("--labels", "--samples"):
        args = ["--labels", "celltype", "--samples", "orig.ident"]
        args[args.index(flag) + 1] = "clonotype_id"
        r = subprocess.run(
            [sys.executable, "-I", SCRIPT, "--dir", str(fixture_dir), *args], capture_output=True, text=True
        )
        assert r.returncode != 0 and "withheld" in r.stderr and "PLANTED" not in r.stdout + r.stderr


def _run_any(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, "-I", SCRIPT, *args], capture_output=True, text=True)


@pytest.mark.parametrize(
    "column",
    ["Clonotype", "ir_foo", "airr:x", "raw_clonotype_id", "my_cdr3", "Xclon", "foo_gene", "TRB_V", "CHAIN_x"],
)
def test_every_withholding_rule_hides_its_values(tmp_path: Path, column: str) -> None:
    for name in ("GSE188666_scRNA_LCMV_metadata.tsv.gz", "GSE188666_scRNA_D5_LCMV_metadata.txt.gz"):
        with gzip.open(tmp_path / name, "wt") as fh:
            fh.write(f"{column}\nAAACCTGAGGCTAGAC-1\tSECRETVALUE\nAAACCTGAGTACGTTC-1\tSECRETVALUE\n")
    r = _run_any(["--dir", str(tmp_path)])
    assert r.returncode == 0
    assert "SECRETVALUE" not in r.stdout + r.stderr
    assert re.search(rf"^  {re.escape(column)}\t\w+\twithheld$", r.stdout, re.M)
    assert "value counts" not in r.stdout and "numeric" not in r.stdout


@pytest.mark.parametrize("flag", ["--labels", "--samples"])
def test_a_missing_labels_or_samples_column_is_an_error(fixture_dir: Path, flag: str) -> None:
    args = ["--labels", "celltype", "--samples", "orig.ident"]
    args[args.index(flag) + 1] = "no_such_column"
    r = _run_any(["--dir", str(fixture_dir), *args])
    assert r.returncode != 0
    assert "no_such_column" in r.stderr and "not a column" in r.stderr


def test_the_output_shape_is_pinned(fixture_dir: Path) -> None:
    out = run(fixture_dir).stdout
    first = out.split("== ", 2)[1].splitlines()
    assert first[0] == "GSE188666_scRNA_LCMV_metadata.tsv.gz"
    assert first[1] == "rows\t6"
    assert first[2] == "first index values\t" + "\t".join(
        f"Spleen_Chronic_D8_gex_{c}-1" for c in ("AAACCTGAGGCTAGAC", "AAACCTGAGTACGTTC", "AAACGGGAGCTGAAAT")
    )
    order = [
        "index tag vocabulary (tag, count):",
        "columns (name, dtype, withheld):",
        "numeric nCount_RNA: min 1000  median 1002.5  max 1005",
        "value counts nCount_RNA (6 distinct):",
        "value counts celltype (2 distinct):",
        "  Texprog_fixture\t3",
        "  Texterm_fixture\t3",
        "value counts tissue (1 distinct):",
        "cross-tab orig.ident (rows) x celltype (columns):",
        "orig.ident\tTexprog_fixture\tTexterm_fixture",
    ]
    # orig.ident is the sample column; nCount_RNA is shifted by the index, so it holds the sample name.
    positions = [out.index(s) for s in order if s in out]
    assert positions == sorted(positions)
    assert (
        "cross-tab orig.ident (rows) x celltype (columns):\norig.ident\tTexprog_fixture\tTexterm_fixture\n"
        in out
    )
    assert "Spleen_Chronic_D8\t3\t3" in out


def test_the_cross_file_comparison_reads_only_the_index_and_two_columns(tmp_path: Path) -> None:
    def write(name: str, rows: list[tuple[str, str, str]]) -> None:
        with gzip.open(tmp_path / name, "wt") as fh:
            fh.write("orig.ident\tident\tcdr3s_aa\n")
            fh.writelines(f"{i}\t{s}\t{lab}\tSECRETVALUE\n" for i, s, lab in rows)

    write(
        "GSE188666_scRNA_LCMV_metadata.tsv.gz",
        [("c1-1", "S", "A"), ("c2-1", "S", "A"), ("c3-1", "S", "B"), ("c4-1", "S", "B")],
    )
    write(
        "GSE188666_scRNA_D5_LCMV_metadata.txt.gz",
        [("c1-1", "S", "X"), ("c2-1", "S", "Y"), ("c3-1", "S", "Y"), ("c9-1", "T", "Y")],
    )
    r = _run_any(["--dir", str(tmp_path), "--compare", "--labels", "ident", "--samples", "orig.ident"])
    assert r.returncode == 0, r.stderr
    assert "SECRETVALUE" not in r.stdout + r.stderr
    assert "cells in both\t3" in r.stdout
    assert "only in main\t1" in r.stdout and "only in D5\t1" in r.stdout
    assert (
        "main ident (rows) x D5 ident (columns), shared cells:\nident\tX\tY\nA\t1\t1\nB\t0\t1\n" in r.stdout
    )
    assert "cdr3s_aa" not in r.stdout


def test_the_comparison_refuses_a_withheld_column(fixture_dir: Path) -> None:
    r = _run_any(
        ["--dir", str(fixture_dir), "--compare", "--labels", "clonotype_id", "--samples", "orig.ident"]
    )
    assert r.returncode != 0 and "withheld" in r.stderr
