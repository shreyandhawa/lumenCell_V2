"""The Daniel 2022 inspection scripts list receptor columns by name and never print their values."""

import gzip
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
