"""A tiny GEO layout of Daniel 2022 in the formats Task 5a recorded (note sections 2 and 9).

Matrices are MatrixMarket coordinate integer, genes x cells; features are one gene-symbol column;
barcodes are ``<tag><16-nt core>-1``; the metadata is a Seurat export (one header field fewer than each
row, the first field of a row being the cell id); contig members are gzipped 10x
``filtered_contig_annotations`` tables with bare ``<core>-1`` barcodes. Every receptor value is made up.
"""

import gzip
import io
import tarfile
from pathlib import Path

import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp

from lumencell.bench.builders import daniel2022
from lumencell.bench.builders.daniel2022 import Capture

GSM_A, GSM_B, GSM_D = "GSM_A", "GSM_B", "GSM_D"
TAG_A, TAG_B, TAG_D, TAG_X = "A_gex_", "B_gex_", "D_gex_", "X_gex_"
MEMBER_END = "_tcr_filtered_contig_annotations.csv.gz"

CORE = "AAACCTGAGAAGGCCT"  # in GSM_A, GSM_B and (D5) GSM_D
THREE_TRB = "AAACCTGAGCCCAAAA"  # GSM_A: three distinct qualifying TRBs
BAD_TRA = "AAACCTGAGCCCCCCC"  # GSM_A: its only TRA is productive="False"
NO_CONTIG = "AAACGGGCATTGTGCA"  # GSM_A: a cell with no contig
NOT_IN_META = "AAACGGGTCACAGGCC"  # GSM_A: in the matrix, not in the metadata
MISSED_CORE = "TTTCCTGAGAAGGCCT"  # GSM_A: qualifying contigs, no cell
TRA_MISS = "ACACCTGAGAAGGCCT"  # GSM_A: only a qualifying TRA, no cell (a join miss)
NOT_CELL = "AAACGGGGTAGCTAGT"  # GSM_A: a matrix cell whose only contigs have is_cell=False
NOT_CELL_MISS = "AAACGGGGTTTTTTTT"  # GSM_A: is_cell=False contigs, no cell (not a join miss)
D_MISS = "GTGCCTGAGAAGGCCT"  # GSM_D: qualifying contigs, no cell (a join miss)
B2 = "CCCCCTGAGAAGGCCT"  # GSM_B: a second paired cell
UNQUALIFIED_B = "GGGCCTGAGAAGGCCT"  # GSM_B: no cell, and no qualifying contig (not a join miss)

GENES = ["Xkr4", "Cd8a", "Trbc2", "Trav1", "Tcf7", "Havcr2"]
STRAY = "../stray.csv.gz"

PARAMS = {"captures": [GSM_A, GSM_B], "namespace": {GSM_A: "m1", GSM_B: "m2"}}

CAPTURES = (
    Capture(GSM_A, f"{GSM_A}_A{MEMBER_END}", TAG_A, "main", "LCMV-c13", 21, "spleen", "CD8+ gp33+ T cells"),
    Capture(
        GSM_B, f"{GSM_B}_B{MEMBER_END}", TAG_B, "main", "LCMV-Armstrong", 8, "lung", "CD8+ gp33- T cells"
    ),
    Capture(GSM_D, f"{GSM_D}_D{MEMBER_END}", TAG_D, "d5", "LCMV-c13", 5, "liver", "CD8+ gp33+ T cells"),
)

# (tag, core, label); the main matrix's 8 cells and the D5 matrix's 3, in file order.
MAIN_CELLS = [
    (TAG_A, CORE, "Tpex"),
    (TAG_B, CORE, "Tex"),
    (TAG_A, THREE_TRB, "Tex"),
    (TAG_A, BAD_TRA, "Teff"),
    (TAG_A, NO_CONTIG, ""),  # a blank label
    (TAG_A, NOT_IN_META, "Tpex"),
    (TAG_B, B2, "Teff"),
    (TAG_X, CORE, "Tex"),  # a tag of no capture
    (TAG_A, NOT_CELL, "Tex"),
]
D5_CELLS = [
    (TAG_D, CORE, "D5_early"),
    (TAG_A, CORE, "D5_repeat"),  # the D5 matrix repeats main cells, with its own labels
    (TAG_B, CORE, "D5_repeat"),
]

_HEADER = (
    "barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,"
    "full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id"
)


def _contig(
    core: str,
    chain: str,
    nt: str,
    clonotype: str,
    productive: str = "True",
    full: str = "True",
    is_cell: str = "True",
) -> str:
    v = "TRAV1" if chain == "TRA" else "TRBV1"
    c = "TRAC" if chain == "TRA" else "TRBC2"
    return (
        f"{core}-1,{is_cell},{core}-1_contig_1,True,500,{chain},{v},None,None,{c},"
        f"{full},{productive},CAS,{nt},100,5,{clonotype},{clonotype}_consensus_1"
    )


CONTIGS = {
    GSM_A: [
        _contig(CORE, "TRA", "TGTGCTAGA", "clonotype1"),
        _contig(CORE, "TRB", "TGTGCCAGC", "clonotype1"),
        _contig(THREE_TRB, "TRA", "TGTGCAACA", "clonotype2"),
        _contig(THREE_TRB, "TRB", "TGTAAATTT", "clonotype2"),
        _contig(THREE_TRB, "TRB", "TGTCCCTTT", "clonotype2"),
        _contig(THREE_TRB, "TRB", "TGTGGGTTT", "clonotype2"),
        _contig(BAD_TRA, "TRA", "TGTTTTAAA", "clonotype3", productive="False"),
        _contig(BAD_TRA, "TRB", "TGTACAGGG", "clonotype3"),
        _contig(MISSED_CORE, "TRA", "TGTACCAAA", "clonotype4"),
        _contig(MISSED_CORE, "TRB", "TGTCCAGGG", "clonotype4"),
        _contig(TRA_MISS, "TRA", "TGTACCTTT", "clonotype5"),
        _contig(NOT_CELL, "TRA", "TGTGCTCCC", "clonotype6", is_cell="False"),
        _contig(NOT_CELL, "TRB", "TGTGCCCCC", "clonotype6", is_cell="False"),
        _contig(NOT_CELL_MISS, "TRA", "TGTGCTGGG", "clonotype7", is_cell="False"),
        _contig(NOT_CELL_MISS, "TRB", "TGTGCCGGG", "clonotype7", is_cell="False"),
    ],
    GSM_B: [
        _contig(CORE, "TRA", "TGTGGAAAA", "clonotype1"),
        _contig(CORE, "TRB", "TGTCGAGGG", "clonotype1"),
        _contig(B2, "TRA", "TGTAGAAAA", "clonotype2"),
        _contig(B2, "TRB", "TGTAGCGGG", "clonotype2"),
        _contig(UNQUALIFIED_B, "TRB", "TGTAAACCC", "clonotype3", full="False"),
        _contig(UNQUALIFIED_B, "Multi", "TGTAAACCC", "clonotype3"),
        # The core of a GSM_A cell, not of any GSM_B cell: a GSM_B join miss.
        _contig(THREE_TRB, "TRA", "TGTCAGAAA", "clonotype4"),
        _contig(THREE_TRB, "TRB", "TGTCAGGGG", "clonotype4"),
    ],
    GSM_D: [
        _contig(CORE, "TRA", "TGTGCTAGA", "clonotype1"),
        _contig(CORE, "TRB", "TGTGCCAGC", "clonotype1"),
        _contig(D_MISS, "TRA", "TGTGATAAA", "clonotype2"),
        _contig(D_MISS, "TRB", "TGTGATGGG", "clonotype2"),
    ],
}


def patch_captures(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(daniel2022, "CAPTURES", CAPTURES)
    monkeypatch.setattr(daniel2022, "LABEL_COLUMN", "ident")


def gz(path: Path, text: str) -> None:
    path.write_bytes(gzip.compress(text.encode()))


def _matrix(path: Path, genes_by_cells: np.ndarray) -> None:
    buf = io.BytesIO()
    scipy.io.mmwrite(buf, sp.coo_matrix(genes_by_cells.astype(np.int64)))
    path.write_bytes(gzip.compress(buf.getvalue()))


def _metadata(path: Path, cells: list[tuple[str, str, str]]) -> None:
    # A Seurat export: 4 header fields, 5 fields per row. tcr_cdr3s_nt is planted and never read.
    lines = ["orig.ident\tnCount_RNA\ttcr_cdr3s_nt\tident"]
    lines += [f"{tag}{core}-1\t{tag[:-5]}\t100\tTGTPLANTED\t{label}" for tag, core, label in cells]
    gz(path, "\n".join(lines) + "\n")


def _member(data: bytes, name: str) -> tuple[tarfile.TarInfo, io.BytesIO]:
    info = tarfile.TarInfo(name)
    info.size = len(data)
    return info, io.BytesIO(data)


def write_layout(
    d: Path,
    *,
    d5_genes: list[str] | None = None,
    features: list[str] | None = None,
    contigs: dict[str, list[str]] | None = None,
) -> tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]]:
    """Write the nine source files into ``d``. Returns the paths by file name and, per (GSM, core) of a
    capture cell, its expected counts in gene order."""
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    main = rng.integers(0, 9, size=(len(GENES), len(MAIN_CELLS)))
    d5 = rng.integers(20, 29, size=(len(GENES), len(D5_CELLS)))  # never equal to a main count
    by_tag = {c.tag: c.gsm for c in CAPTURES}
    counts: dict[tuple[str, str], np.ndarray] = {}
    for j, (tag, core, _) in enumerate(MAIN_CELLS):
        if tag in by_tag and by_tag[tag] != GSM_D:
            counts[(by_tag[tag], core)] = main[:, j].astype(np.float32)
    for j, (tag, core, _) in enumerate(D5_CELLS):
        if tag == TAG_D:
            counts[(GSM_D, core)] = d5[:, j].astype(np.float32)

    p = "GSE188666_scRNA_"
    feature_lines = features if features is not None else GENES
    for stem, cells, m, genes in (
        ("LCMV", MAIN_CELLS, main, feature_lines),
        ("D5_LCMV", D5_CELLS, d5, d5_genes if d5_genes is not None else feature_lines),
    ):
        _matrix(d / f"{p}{stem}.mtx.gz", m)
        gz(d / f"{p}{stem}_barcodes.tsv.gz", "".join(f"{t}{c}-1\n" for t, c, _ in cells))
        gz(d / f"{p}{stem}_features.tsv.gz", "".join(f"{g}\n" for g in genes))
    _metadata(d / f"{p}LCMV_metadata.tsv.gz", [c for c in MAIN_CELLS if c[1] != NOT_IN_META])
    _metadata(d / f"{p}D5_LCMV_metadata.txt.gz", D5_CELLS)

    with tarfile.open(d / "GSE188666_RAW.tar", "w") as tf:
        for c in CAPTURES:
            text = "\n".join([_HEADER, *(contigs or CONTIGS)[c.gsm]]) + "\n"
            tf.addfile(*_member(gzip.compress(text.encode()), c.member))
        tf.addfile(*_member(gzip.compress(f"{_HEADER}\n".encode()), STRAY))

    paths = {s.filename: d / s.filename for s in daniel2022.SOURCES}
    assert all(q.exists() for q in paths.values())
    return paths, counts
