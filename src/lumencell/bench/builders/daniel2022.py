"""Daniel et al. 2022 (chronic LCMV, paired scRNA and TCR): the GEO files and their formats."""

from dataclasses import dataclass
from typing import Literal

from lumencell.bench.download import Source

GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188666/suppl/"

_FILES: tuple[tuple[str, int], ...] = (
    ("GSE188666_scRNA_LCMV.mtx.gz", 602240886),
    ("GSE188666_scRNA_LCMV_barcodes.tsv.gz", 487511),
    ("GSE188666_scRNA_LCMV_features.tsv.gz", 113716),
    ("GSE188666_scRNA_LCMV_metadata.tsv.gz", 4846117),
    ("GSE188666_scRNA_D5_LCMV.mtx.gz", 544512106),
    ("GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz", 389725),
    ("GSE188666_scRNA_D5_LCMV_features.tsv.gz", 113719),
    ("GSE188666_scRNA_D5_LCMV_metadata.txt.gz", 3976768),
    ("GSE188666_RAW.tar", 12113920),
)

# SHA-256 as recorded on the first download (2026-11-19). The two mtx sizes are 574 and 519 MiB.
_SHA256 = {
    "GSE188666_scRNA_LCMV.mtx.gz": ("be4a6e2825d9cdffbb4bebb094e43be00ed388faa11a5bac95c2418bdd576ed6"),
    "GSE188666_scRNA_LCMV_barcodes.tsv.gz": (
        "c75bc8d63c793b4c74215d46d40d243118b5fb8f4e1aff248af0451ee15196fe"
    ),
    "GSE188666_scRNA_LCMV_features.tsv.gz": (
        "dee1b30b23a9ceae7117d71b894ec7e4d5441fc934fbcb242493ef57b2e06e2d"
    ),
    "GSE188666_scRNA_LCMV_metadata.tsv.gz": (
        "0594a0463b0d667307736fa45453d9aae9ed761ce105d53b3a24c30dffa46393"
    ),
    "GSE188666_scRNA_D5_LCMV.mtx.gz": ("14bc8bf47f68d11b27dbabfc8ec6613385f58a83578a0a81b7187817602657b3"),
    "GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz": (
        "3ec9ca5d7650e5ab428121161c7533d1da42d72ca1187dd9941e8de0207151c8"
    ),
    "GSE188666_scRNA_D5_LCMV_features.tsv.gz": (
        "e3367683bbf603ef45124593bcc9a2868ffe393e422df1a5a035c24e7e76a576"
    ),
    "GSE188666_scRNA_D5_LCMV_metadata.txt.gz": (
        "9f8943dea624e8f56c894421eaad3b723b9e2698a3f8d148122fa516bab7c4ec"
    ),
    "GSE188666_RAW.tar": ("db4932008b1b9e603c2f076aa7cf60302f38d8b66108188f3827a121a676c9a9"),
}

SOURCES: tuple[Source, ...] = tuple(Source(GEO_SUPPL + n, n, size, _SHA256[n]) for n, size in _FILES)

LICENCE = (
    "Data: GEO GSE188666 (SuperSeries GSE188670); "
    "UNVERIFIED: https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html is behind a CAPTCHA, "
    "not bypassed (read attempt 2026-11-19). "
    "Paper: Daniel et al. 2022, Nat Immunol 23:1614-1627, doi:10.1038/s41590-022-01337-5; "
    "UNVERIFIED: https://pmc.ncbi.nlm.nih.gov/articles/PMC11225711/ is behind a CAPTCHA, "
    "not bypassed (read attempt 2026-11-19)."
)


@dataclass(frozen=True)
class Capture:
    """One TCR library (a tar member) and the expression capture whose cells it holds.

    `gsm` is the TCR GSM; `member` its tar member. `tag` is the cell-id prefix of that capture in the
    matrix barcodes: a contig barcode `b` is the cell `tag + b`. `arm`, `day`, `tissue` and `sort_gate`
    describe those cells, as the GSM records write them (`sort_gate` is `cell type:`).
    """

    gsm: str
    member: str
    tag: str
    matrix: Literal["main", "d5"]
    arm: str
    day: int
    tissue: str
    sort_gate: str


_MEMBER_END = "_tcr_filtered_contig_annotations.csv.gz"


def _c(
    gsm: str, name: str, matrix: Literal["main", "d5"], arm: str, day: int, tissue: str, gate: str
) -> Capture:
    # For every member the join puts its barcodes under the tag of its own file name (see below).
    return Capture(gsm, f"{gsm}_{name}{_MEMBER_END}", f"{name}_gex_", matrix, arm, day, tissue, gate)


# From the Task 3 note (research/2026-11-19-daniel2022-metadata.md) and the Task 5a join check there.
# Every member's barcodes sit under its own tag (share 0.58-0.96; every other tag at most 0.02).
# GSM6381667 and GSM6381668: the titles say TetPos_rep1 and TetNeg_rep2, the file names TetNeg_rep2 and
# TetPos_rep1. The join follows the file names (0.75 and 0.84 against 0.01 by title), so tag and labels
# follow the file names: the labels describe the cells the barcodes are.
# GSM5687829-32: paired by title, as the note does (0.85-0.96; the sequencer pairing gets at most 0.02).
CAPTURES: tuple[Capture, ...] = (
    _c(
        "GSM5687816",
        "Liver_Chronic_D21_CD8_TetNeg",
        "main",
        "LCMV-c13",
        21,
        "liver",
        "CD8+ gp33- T cells",
    ),
    _c(
        "GSM5687817",
        "Liver_Chronic_D21_CD8_TetPos",
        "main",
        "LCMV-c13",
        21,
        "liver",
        "CD8+ gp33+ T cells",
    ),
    _c(
        "GSM5687818",
        "Lung_Chronic_D21_CD8_TetNeg",
        "main",
        "LCMV-c13",
        21,
        "lung",
        "CD8+ gp33- T cells",
    ),
    _c(
        "GSM5687819",
        "Lung_Chronic_D21_CD8_TetPos",
        "main",
        "LCMV-c13",
        21,
        "lung",
        "CD8+ gp33+ T cells",
    ),
    _c(
        "GSM5687820",
        "Spleen_Acute_D21_CD8_TetNeg",
        "main",
        "LCMV-Armstrong",
        21,
        "spleen",
        "CD8+ gp33- T cells",
    ),
    _c(
        "GSM5687821",
        "Spleen_Acute_D21_CD8_TetPos",
        "main",
        "LCMV-Armstrong",
        21,
        "spleen",
        "CD8+ gp33+ T cells",
    ),
    _c(
        "GSM5687822",
        "Spleen_Acute_D8_CD8_TetNeg",
        "main",
        "LCMV-Armstrong",
        8,
        "spleen",
        "CD8+ gp33- T cells",
    ),
    _c(
        "GSM5687823",
        "Spleen_Acute_D8_CD8_TetPos",
        "main",
        "LCMV-Armstrong",
        8,
        "spleen",
        "CD8+ gp33+ T cells",
    ),
    _c(
        "GSM5687824",
        "Spleen_Chronic_D21_CD8_CX3CR1Neg_Slamf6Pos",
        "main",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ CX3CR1- SLAMF6+ T cells",
    ),
    _c(
        "GSM5687825",
        "Spleen_Chronic_D21_CD8_CX3CR1Neg",
        "main",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ CX3CR1- T cells",
    ),
    _c(
        "GSM5687826",
        "Spleen_Chronic_D21_CD8_CX3CR1Pos",
        "main",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ CX3CR1+ T cells",
    ),
    _c(
        "GSM5687827",
        "Spleen_Chronic_D21_CD8_PD1Pos",
        "main",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ PD1+ T cells",
    ),
    _c(
        "GSM5687828",
        "Spleen_Chronic_D21_CD8_TetMid",
        "main",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ gp33 mid T cells",
    ),
    _c(
        "GSM5687829",
        "Spleen_Chronic_D21_CD8_TetNeg_rep1",
        "main",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ gp33- T cells",
    ),
    _c(
        "GSM5687830",
        "Spleen_Chronic_D21_CD8_TetNeg",
        "main",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ gp33- T cells",
    ),
    _c(
        "GSM5687831",
        "Spleen_Chronic_D21_CD8_TetPos_rep1",
        "main",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ gp33+ T cells",
    ),
    _c(
        "GSM5687832",
        "Spleen_Chronic_D21_CD8_TetPos",
        "main",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ gp33+ T cells",
    ),
    _c(
        "GSM5687833",
        "Spleen_Chronic_D8_CD8_TetNeg",
        "main",
        "LCMV-c13",
        8,
        "spleen",
        "CD8+ gp33- T cells",
    ),
    _c(
        "GSM5687834",
        "Spleen_Chronic_D8_CD8_TetPos",
        "main",
        "LCMV-c13",
        8,
        "spleen",
        "CD8+ gp33+ T cells",
    ),
    _c(
        "GSM6381666",
        "Spleen_Chronic_D5_CD8_TetNeg_rep1",
        "d5",
        "LCMV-c13",
        5,
        "spleen",
        "CD8+ gp33- T cells",
    ),
    _c(
        "GSM6381667",
        "Spleen_Chronic_D5_CD8_TetNeg_rep2",
        "d5",
        "LCMV-c13",
        5,
        "spleen",
        "CD8+ gp33- T cells",
    ),
    _c(
        "GSM6381668",
        "Spleen_Chronic_D5_CD8_TetPos_rep1",
        "d5",
        "LCMV-c13",
        5,
        "spleen",
        "CD8+ gp33+ T cells",
    ),
    _c(
        "GSM6381669",
        "Spleen_Chronic_D5_CD8_TetPos_rep2",
        "d5",
        "LCMV-c13",
        5,
        "spleen",
        "CD8+ gp33+ T cells",
    ),
)

LABEL_COLUMN: str | None = "ident"
# The real contig headers already use the CONTIG_COLUMNS spellings.
COLUMN_RENAMES: dict[str, str] = {}
# Size lines 31053 x 96750 and 31053 x 77614: rows are the 31,053 features, columns the barcodes.
MATRIX_ORIENTATION: dict[str, Literal["genes_x_cells", "cells_x_genes"]] = {
    "main": "genes_x_cells",
    "d5": "genes_x_cells",
}
