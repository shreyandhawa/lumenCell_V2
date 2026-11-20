"""The Daniel 2022 inspection scripts list receptor columns by name and never print their values."""

import gzip
import io
import json
import re
import subprocess
import sys
import tarfile
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


# ---- The file-format script (phase 3 Task 5a) ----

FORMATS = str(Path(__file__).resolve().parents[2] / "benchmarks" / "inspection" / "daniel2022_formats.py")
# The real 10x header of every Daniel contig member.
CONTIG_HEADER = [
    "barcode",
    "is_cell",
    "contig_id",
    "high_confidence",
    "length",
    "chain",
    "v_gene",
    "d_gene",
    "j_gene",
    "c_gene",
    "full_length",
    "productive",
    "cdr3",
    "cdr3_nt",
    "reads",
    "umis",
    "raw_clonotype_id",
    "raw_consensus_id",
]
READ_COLUMNS = {"barcode", "is_cell", "chain", "productive", "high_confidence", "full_length"}
WITHHELD_CONTIG = [c for c in CONTIG_HEADER if c not in READ_COLUMNS]
C = [
    "AAACCTGAGGCTAGAC",
    "AAACCTGAGTACGTTC",
    "AAACGGGAGCTGAAAT",
    "AAAGATGAGAGTAAGG",
    "AACTCAGTCAGCTCGG",
    "AACTCCCAGTGGTCCC",
    "ACACCGGAGCTAGTCT",
    "ACGATGTAGGCATGGT",
]
MEMBER_X = "GSM1_Spleen_X_tcr_filtered_contig_annotations.csv.gz"
MEMBER_Y = "GSM2_Spleen_Y_tcr_filtered_contig_annotations.csv.gz"


def _gz(path: Path, lines: list[str]) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def _contigs(barcodes: list[str], planted: dict[str, str], header: list[str]) -> bytes:
    rows = [",".join(header)]
    for i, bc in enumerate(barcodes):
        for chain in ("TRA", "TRB"):
            row = {c: "0" for c in header}
            row.update(
                barcode=bc,
                is_cell="True",
                contig_id=f"{bc}_contig_{chain}",
                high_confidence="True" if i else "False",
                chain=chain,
                full_length="True",
                productive="None" if i == 2 else "True",
            )
            row.update(planted)
            rows.append(",".join(row[c] for c in header))
    return gzip.compress(("\n".join(rows) + "\n").encode())


def _tar(path: Path, members: dict[str, bytes]) -> None:
    with tarfile.open(path, "w") as tf:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))


PLANTED = {
    "cdr3": "CASSPLANTEDF",
    "cdr3_nt": "PLANTEDNT",
    "v_gene": "TRBV-PLANTED",
    "raw_clonotype_id": "clonotypePLANTED",
}


def _formats_dir(
    d: Path,
    planted: dict[str, str] = PLANTED,
    d5_cells_x_genes: bool = False,
    header: list[str] = CONTIG_HEADER,
    main_size: str | None = None,
    d5_genes: list[str] | None = None,
    extra_members: dict[str, bytes] | None = None,
) -> Path:
    genes = ["Xkr4", "Gm1992", "Cd8a"]
    main = [f"Spleen_X_gex_{c}-1" for c in C[:3]] + [f"Spleen_Y_gex_{c}-1" for c in C[3:6]]
    d5 = [f"Spleen_D5_gex_{c}-1" for c in C[6:8]]
    for stem, cells in (("GSE188666_scRNA_LCMV", main), ("GSE188666_scRNA_D5_LCMV", d5)):
        _gz(d / f"{stem}_barcodes.tsv.gz", cells)
        _gz(d / f"{stem}_features.tsv.gz", d5_genes if d5_genes and "D5" in stem else genes)
        size = f"{len(cells)} {len(genes)} 2" if d5_cells_x_genes and "D5" in stem else f"3 {len(cells)} 2"
        if main_size is not None and "D5" not in stem:
            size = main_size
        _gz(
            d / f"{stem}.mtx.gz",
            ["%%MatrixMarket matrix coordinate integer general", "%comment", size, "1 1 2", "3 2 7"],
        )
    # Member Y holds two Y cells and one X cell, so its own share is 2/3 and its best other share 1/3.
    _tar(
        d / "GSE188666_RAW.tar",
        {
            MEMBER_X: _contigs([f"{c}-1" for c in C[:3]], planted, header),
            MEMBER_Y: _contigs([f"{c}-1" for c in (C[3], C[4], C[0])], planted, header),
            **(extra_members or {}),
        },
    )
    return d


def run_formats(d: Path, *extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-I", FORMATS, "--dir", str(d), *extra], capture_output=True, text=True
    )


def _map(d: Path, mapping: dict[str, str]) -> str:
    p = d / "capture_map.json"
    p.write_text(json.dumps(mapping))
    return str(p)


def test_the_format_script_never_prints_receptor_values(tmp_path: Path) -> None:
    d = _formats_dir(tmp_path)
    m = _map(d, {MEMBER_X: "Spleen_X_gex_", MEMBER_Y: "Spleen_Y_gex_"})
    r = run_formats(d, "--capture-map", m)
    assert r.returncode == 0, r.stderr
    assert "PLANTED" not in r.stdout + r.stderr
    assert "distinct values chain (pooled over 2 members):\n  TRA\t6\n  TRB\t6\n" in r.stdout


@pytest.mark.parametrize("column", [*WITHHELD_CONTIG, "a_new_column"])
def test_every_withheld_contig_column_hides_its_values(tmp_path: Path, column: str) -> None:
    # A column the script has never heard of is withheld too: the read columns are an allow-list.
    header = CONTIG_HEADER if column in CONTIG_HEADER else [*CONTIG_HEADER, column]
    d = _formats_dir(tmp_path, {column: "SECRETVALUE"}, header=header)
    m = _map(d, {MEMBER_X: "Spleen_X_gex_", MEMBER_Y: "Spleen_Y_gex_"})
    r = run_formats(d, "--capture-map", m)
    assert r.returncode == 0, r.stderr
    assert "SECRETVALUE" not in r.stdout + r.stderr
    assert re.search(rf"^withheld columns \(by name, never read\)\t.*\b{re.escape(column)}\b", r.stdout, re.M)
    assert f"distinct values {column} " not in r.stdout


def test_the_join_share_is_reported_per_member(tmp_path: Path) -> None:
    d = _formats_dir(tmp_path)
    m = _map(d, {MEMBER_X: "Spleen_X_gex_", MEMBER_Y: "Spleen_Y_gex_"})
    out = run_formats(d, "--capture-map", m).stdout
    # A tie for the best other tag goes to the first tag in sorted order.
    assert f"join\t{MEMBER_X}\tSpleen_X_gex_\t3/3\t1.00\tSpleen_D5_gex_\t0/3\t0.00\n" in out
    assert f"join\t{MEMBER_Y}\tSpleen_Y_gex_\t2/3\t0.67\tSpleen_X_gex_\t1/3\t0.33\n" in out


def test_the_matrix_heads_and_vocabularies_are_printed(tmp_path: Path) -> None:
    out = run_formats(_formats_dir(tmp_path)).stdout
    assert "== matrix main GSE188666_scRNA_LCMV.mtx.gz\n" in out
    assert "banner\t%%MatrixMarket matrix coordinate integer general\n" in out
    assert "size line\t3 6 2\n" in out
    assert "barcodes\t6\tfirst\t" + "\t".join(f"Spleen_X_gex_{c}-1" for c in C[:3]) + "\n" in out
    assert "features\t3\tfirst\tXkr4\tGm1992\tCd8a\n" in out
    assert "feature fields per line\t1\n" in out
    assert "orientation\tgenes_x_cells\n" in out
    assert "first 2 stored values integer\tTrue\n" in out
    assert "tag vocabulary (tag, count):\n  Spleen_X_gex_\t3\n  Spleen_Y_gex_\t3\n" in out
    assert "D5 features equal main features\tTrue\n" in out
    assert "D5 barcodes also among main barcodes\t0\n" in out
    assert "first member first barcodes\t" + "\t".join(f"{c}-1" for c in C[:3]) + "\n" in out
    assert f"member\t{MEMBER_X}\nheader\t{','.join(CONTIG_HEADER)}\n" in out
    assert "distinct values productive (pooled over 2 members):\n  None\t4\n  True\t8\n" in out


def test_a_swapped_size_line_reads_as_cells_by_genes(tmp_path: Path) -> None:
    out = run_formats(_formats_dir(tmp_path, d5_cells_x_genes=True)).stdout
    d5 = out.split("== matrix d5 ", 1)[1]
    assert "orientation\tcells_x_genes\n" in d5


def test_a_non_integer_value_is_reported(tmp_path: Path) -> None:
    d = _formats_dir(tmp_path)
    _gz(
        d / "GSE188666_scRNA_LCMV.mtx.gz",
        ["%%MatrixMarket matrix coordinate real general", "3 6 2", "1 1 2", "3 2 0.5"],
    )
    assert "first 2 stored values integer\tFalse\n" in run_formats(d).stdout


def test_the_format_script_writes_nothing(tmp_path: Path) -> None:
    d = _formats_dir(tmp_path)
    m = _map(d, {MEMBER_X: "Spleen_X_gex_", MEMBER_Y: "Spleen_Y_gex_"})
    before = {p: p.read_bytes() for p in sorted(d.iterdir())}
    run_formats(d, "--capture-map", m)
    assert {p: p.read_bytes() for p in sorted(d.iterdir())} == before


def test_a_map_member_not_in_the_tar_is_refused_without_echoing_it(tmp_path: Path) -> None:
    d = _formats_dir(tmp_path)
    m = _map(d, {MEMBER_X: "Spleen_X_gex_", "GSM9_SECRETNAME.csv.gz": "Spleen_Y_gex_"})
    r = run_formats(d, "--capture-map", m)
    assert r.returncode != 0
    assert "1 capture-map member(s) are not in the tar" in r.stderr
    assert "SECRETNAME" not in r.stdout + r.stderr


def test_a_malformed_capture_map_is_refused(tmp_path: Path) -> None:
    d = _formats_dir(tmp_path)
    p = d / "bad.json"
    for text in ("[1, 2]", '{"a": 3}', "not json"):
        p.write_text(text)
        r = run_formats(d, "--capture-map", str(p))
        assert r.returncode != 0 and "capture map" in r.stderr


def _raw_member(text: bytes, gz: bool = True) -> bytes:
    return gzip.compress(text) if gz else text


def test_a_headerless_member_is_refused_before_any_of_it_is_printed(tmp_path: Path) -> None:
    row = ",".join(["AAACCTGAGGCTAGAC-1", "True", "c1", "True", "500", "TRB", "TRBV-PLANTED", "None"])
    row += ",None,TRBJ-PLANTED,C,True,True,CASSPLANTEDF,PLANTEDNT,10,2,clonotypePLANTED,consensus1\n"
    name = "GSM3_SECRETNAME_tcr_filtered_contig_annotations.csv.gz"
    d = _formats_dir(tmp_path, extra_members={name: _raw_member(row.encode())})
    r = run_formats(d)
    assert r.returncode != 0
    assert "PLANTED" not in r.stdout + r.stderr and "SECRETNAME" not in r.stdout + r.stderr
    assert "no barcode column" in r.stderr and "member 3 of 3" in r.stderr


@pytest.mark.parametrize(
    "data",
    [
        b"\x1f\x8b not really gzip, PLANTEDBYTES",
        b"\xff\xfe\xfd,PLANTEDBYTES\n",
        gzip.compress(b"\xff\xfe,PLANTEDBYTES\n"),
    ],
)
def test_a_member_that_is_not_gzip_or_utf8_gives_the_script_error(tmp_path: Path, data: bytes) -> None:
    name = "GSM3_SECRETNAME_tcr_filtered_contig_annotations.csv.gz"
    r = run_formats(_formats_dir(tmp_path, extra_members={name: data}))
    assert r.returncode != 0
    assert "not valid gzip or UTF-8" in r.stderr and "Traceback" not in r.stderr
    assert "PLANTED" not in r.stdout + r.stderr and "SECRETNAME" not in r.stdout + r.stderr


def test_d5_features_in_another_order_or_with_one_changed_name_are_not_equal(tmp_path: Path) -> None:
    for genes in (["Gm1992", "Xkr4", "Cd8a"], ["Xkr4", "Gm1992", "Cd8b"]):
        d = tmp_path / "-".join(genes)
        d.mkdir()
        out = run_formats(_formats_dir(d, d5_genes=genes)).stdout
        assert "D5 features equal main features\tFalse\n" in out


def test_the_join_runs_over_the_union_of_both_barcode_files(tmp_path: Path) -> None:
    # Member Z's cells exist only in the D5 file, under a D5-only tag: joining against main alone gives 0.
    z = "GSM3_Spleen_D5_tcr_filtered_contig_annotations.csv.gz"
    zdata = _contigs([f"{c}-1" for c in C[6:8]], PLANTED, CONTIG_HEADER)
    d = _formats_dir(tmp_path, extra_members={z: zdata})
    m = _map(d, {MEMBER_X: "Spleen_X_gex_", MEMBER_Y: "Spleen_Y_gex_", z: "Spleen_D5_gex_"})
    out = run_formats(d, "--capture-map", m).stdout
    assert f"join\t{z}\tSpleen_D5_gex_\t2/2\t1.00\t" in out


def test_a_member_with_no_barcodes_has_no_share(tmp_path: Path) -> None:
    z = "GSM3_Spleen_Z_tcr_filtered_contig_annotations.csv.gz"
    empty = _raw_member((",".join(CONTIG_HEADER) + "\n").encode())
    d = _formats_dir(tmp_path, extra_members={z: empty})
    m = _map(d, {MEMBER_X: "Spleen_X_gex_", MEMBER_Y: "Spleen_Y_gex_", z: "Spleen_X_gex_"})
    out = run_formats(d, "--capture-map", m).stdout
    assert f"join\t{z}\tSpleen_X_gex_\t0/0\tn/a\tSpleen_D5_gex_\t0/0\tn/a\n" in out


def test_a_map_covering_one_of_two_members_reports_one_unmapped(tmp_path: Path) -> None:
    d = _formats_dir(tmp_path)
    out = run_formats(d, "--capture-map", _map(d, {MEMBER_X: "Spleen_X_gex_"})).stdout
    assert "tar members not in the capture map\t1\n" in out


@pytest.mark.parametrize("size", ["3 6", "3 5 2", "4 6 2", "6 4 2"])
def test_an_unmatched_size_line_is_neither(tmp_path: Path, size: str) -> None:
    out = run_formats(_formats_dir(tmp_path, main_size=size)).stdout
    main = out.split("== matrix main ", 1)[1].split("== matrix d5 ", 1)[0]
    assert "orientation\tneither\n" in main


def _big_matrix(d: Path, last: str) -> None:
    lines = ["%%MatrixMarket matrix coordinate integer general", "3 6 100001", *["1 1 2"] * 100_000, last]
    _gz(d / "GSE188666_scRNA_LCMV.mtx.gz", lines)


def test_the_integer_check_stops_at_exactly_100000_values(tmp_path: Path) -> None:
    d = _formats_dir(tmp_path)
    _big_matrix(d, "1 1 0.5")  # the 100,001st value is not looked at
    assert "first 100000 stored values integer\tTrue\n" in run_formats(d).stdout


def test_the_integer_check_sees_the_100000th_value(tmp_path: Path) -> None:
    d = _formats_dir(tmp_path)
    lines = ["%%MatrixMarket matrix coordinate integer general", "3 6 100000", *["1 1 2"] * 99_999, "1 1 0.5"]
    _gz(d / "GSE188666_scRNA_LCMV.mtx.gz", lines)
    assert "first 100000 stored values integer\tFalse\n" in run_formats(d).stdout


@pytest.mark.real_data
@pytest.mark.slow
def test_the_capture_map_joins_on_the_real_files(tmp_path: Path) -> None:
    from lumencell.bench.builders.daniel2022 import CAPTURES, SOURCES
    from lumencell.bench.download import datasets_dir

    d = datasets_dir()
    if not all((d / s.filename).is_file() for s in SOURCES):
        pytest.skip("run: lumencell bench fetch daniel2022")
    m = tmp_path / "map.json"
    m.write_text(json.dumps({c.member: c.tag for c in CAPTURES}))
    r = run_formats(d, "--capture-map", str(m))
    assert r.returncode == 0, r.stderr
    rows = [line.split("\t") for line in r.stdout.splitlines() if line.startswith("join\t")]
    assert len(rows) == len(CAPTURES) == 23
    for _, _member, _tag, own, _share, _other, best, _ in rows:
        own_hits, n = map(int, own.split("/"))
        best_hits, _n = map(int, best.split("/"))
        assert 2 * own_hits >= n and own_hits > best_hits
