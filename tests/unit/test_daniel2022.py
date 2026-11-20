from lumencell.bench.builders.common import CONTIG_COLUMNS
from lumencell.bench.builders.daniel2022 import (
    CAPTURES,
    COLUMN_RENAMES,
    GEO_SUPPL,
    LABEL_COLUMN,
    LICENCE,
    MATRIX_ORIENTATION,
    SOURCES,
)

EXPECTED = (
    "GSE188666_scRNA_LCMV.mtx.gz",
    "GSE188666_scRNA_LCMV_barcodes.tsv.gz",
    "GSE188666_scRNA_LCMV_features.tsv.gz",
    "GSE188666_scRNA_LCMV_metadata.tsv.gz",
    "GSE188666_scRNA_D5_LCMV.mtx.gz",
    "GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz",
    "GSE188666_scRNA_D5_LCMV_features.tsv.gz",
    "GSE188666_scRNA_D5_LCMV_metadata.txt.gz",
    "GSE188666_RAW.tar",
)


def test_sources_are_the_nine_geo_files() -> None:
    assert tuple(s.filename for s in SOURCES) == EXPECTED
    assert all(s.url == GEO_SUPPL + s.filename and s.sha256 and len(s.sha256) == 64 for s in SOURCES)
    assert 1.0e9 < sum(s.size for s in SOURCES) < 1.2e9 and not any("GSE188670" in s.url for s in SOURCES)


def test_the_licence_is_recorded() -> None:
    assert "GSE188666" in LICENCE and "10.1038/s41590-022-01337-5" in LICENCE


# Sizes and SHA-256 as first downloaded on 2026-11-19 and recomputed from the cache by the reviewer.
EXPECTED_PINS = {
    "GSE188666_scRNA_LCMV.mtx.gz": (
        602240886,
        "be4a6e2825d9cdffbb4bebb094e43be00ed388faa11a5bac95c2418bdd576ed6",
    ),
    "GSE188666_scRNA_LCMV_barcodes.tsv.gz": (
        487511,
        "c75bc8d63c793b4c74215d46d40d243118b5fb8f4e1aff248af0451ee15196fe",
    ),
    "GSE188666_scRNA_LCMV_features.tsv.gz": (
        113716,
        "dee1b30b23a9ceae7117d71b894ec7e4d5441fc934fbcb242493ef57b2e06e2d",
    ),
    "GSE188666_scRNA_LCMV_metadata.tsv.gz": (
        4846117,
        "0594a0463b0d667307736fa45453d9aae9ed761ce105d53b3a24c30dffa46393",
    ),
    "GSE188666_scRNA_D5_LCMV.mtx.gz": (
        544512106,
        "14bc8bf47f68d11b27dbabfc8ec6613385f58a83578a0a81b7187817602657b3",
    ),
    "GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz": (
        389725,
        "3ec9ca5d7650e5ab428121161c7533d1da42d72ca1187dd9941e8de0207151c8",
    ),
    "GSE188666_scRNA_D5_LCMV_features.tsv.gz": (
        113719,
        "e3367683bbf603ef45124593bcc9a2868ffe393e422df1a5a035c24e7e76a576",
    ),
    "GSE188666_scRNA_D5_LCMV_metadata.txt.gz": (
        3976768,
        "9f8943dea624e8f56c894421eaad3b723b9e2698a3f8d148122fa516bab7c4ec",
    ),
    "GSE188666_RAW.tar": (12113920, "db4932008b1b9e603c2f076aa7cf60302f38d8b66108188f3827a121a676c9a9"),
}


def test_every_source_keeps_its_recorded_size_and_checksum() -> None:
    assert {s.filename: (s.size, s.sha256) for s in SOURCES} == EXPECTED_PINS


def test_the_files_come_from_geo() -> None:
    assert GEO_SUPPL == "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188666/suppl/"


def test_the_licence_says_what_could_not_be_read() -> None:
    assert LICENCE.count("UNVERIFIED") == 2
    assert "not bypassed" in LICENCE


def test_captures_match_the_inspection() -> None:
    # 23, not 24: the number the Task 3 note records (one TCR member per capture).
    assert len(CAPTURES) == 23
    assert len({c.gsm for c in CAPTURES}) == len({c.member for c in CAPTURES}) == 23
    assert len({c.tag for c in CAPTURES}) == 23  # one expression capture per TCR member
    assert all(c.member.endswith("_tcr_filtered_contig_annotations.csv.gz") for c in CAPTURES)
    assert all(c.member.startswith(c.gsm + "_") for c in CAPTURES)
    assert all(c.tag.endswith("_gex_") for c in CAPTURES)
    assert {c.matrix for c in CAPTURES} == {"main", "d5"}
    assert all((c.matrix == "d5") == (c.day == 5) for c in CAPTURES)
    assert set(MATRIX_ORIENTATION) == {"main", "d5"}
    assert set(MATRIX_ORIENTATION.values()) <= {"genes_x_cells", "cells_x_genes"}
    assert LABEL_COLUMN == "ident"
    assert set(COLUMN_RENAMES.values()) <= set(CONTIG_COLUMNS)


def test_the_two_record_conflicts_follow_the_join() -> None:
    # Settled by the barcode shares of Task 5a (the note's Join check), not by the GSM titles.
    by_gsm = {c.gsm: c for c in CAPTURES}
    d5 = "Spleen_Chronic_D5_CD8_"
    # Titled TetPos_rep1 / TetNeg_rep2; their barcodes are the TetNeg_rep2 / TetPos_rep1 cells.
    assert (by_gsm["GSM6381667"].tag, by_gsm["GSM6381667"].sort_gate) == (
        f"{d5}TetNeg_rep2_gex_",
        "CD8+ gp33- T cells",
    )
    assert (by_gsm["GSM6381668"].tag, by_gsm["GSM6381668"].sort_gate) == (
        f"{d5}TetPos_rep1_gex_",
        "CD8+ gp33+ T cells",
    )
    # The D21 spleen replicates pair by title, not by sequencer.
    d21 = "Spleen_Chronic_D21_CD8_"
    assert [by_gsm[f"GSM56878{n}"].tag for n in (29, 30, 31, 32)] == [
        f"{d21}{s}_gex_" for s in ("TetNeg_rep1", "TetNeg", "TetPos_rep1", "TetPos")
    ]
