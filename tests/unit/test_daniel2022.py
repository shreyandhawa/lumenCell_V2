import dataclasses
import gzip
import io
import re
import tracemalloc
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.io
import scipy.sparse
import scipy.sparse as sp

from lumencell.bench.builders import daniel2022
from lumencell.bench.builders.common import CONTIG_COLUMNS, BuildError
from lumencell.bench.builders.daniel2022 import (
    BUILD_VERSION,
    CAPTURES,
    COLUMN_RENAMES,
    GEO_SUPPL,
    LABEL_COLUMN,
    LICENCE,
    MATRIX_ORIENTATION,
    SOURCES,
    UNSAFE_MEMBER,
    build,
)
from lumencell.clonal.isolation import receptor_gene_mask, scoring_only_columns
from lumencell.data.loader import load
from tests.unit.daniel_fixture import (
    B2,
    BAD_TRA,
    CONTIGS,
    CORE,
    GENES,
    MISSED_CORE,
    NO_CONTIG,
    NOT_CELL,
    NOT_CELL_MISS,
    PARAMS,
    THREE_TRB,
    gz,
    patch_captures,
    write_layout,
)
from tests.unit.daniel_fixture import (
    CAPTURES as FIX_CAPTURES,
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


# ---- The 23 captures, written out by hand from the Task 3 note (sections 1, 4 and 10) ----
# gsm -> (file-name stem, matrix, arm, day, tissue, sort gate). Member is `<gsm>_<stem>` + _END,
# and tag `<stem>_gex_`. GSM6381667/68 follow their file names (note section 10), not their titles.
_C13, _ARM = "LCMV-c13", "LCMV-Armstrong"
_NEG, _POS, _MID = "CD8+ gp33- T cells", "CD8+ gp33+ T cells", "CD8+ gp33 mid T cells"
PINNED = {
    "GSM5687816": ("Liver_Chronic_D21_CD8_TetNeg", "main", _C13, 21, "liver", _NEG),
    "GSM5687817": ("Liver_Chronic_D21_CD8_TetPos", "main", _C13, 21, "liver", _POS),
    "GSM5687818": ("Lung_Chronic_D21_CD8_TetNeg", "main", _C13, 21, "lung", _NEG),
    "GSM5687819": ("Lung_Chronic_D21_CD8_TetPos", "main", _C13, 21, "lung", _POS),
    "GSM5687820": ("Spleen_Acute_D21_CD8_TetNeg", "main", _ARM, 21, "spleen", _NEG),
    "GSM5687821": ("Spleen_Acute_D21_CD8_TetPos", "main", _ARM, 21, "spleen", _POS),
    "GSM5687822": ("Spleen_Acute_D8_CD8_TetNeg", "main", _ARM, 8, "spleen", _NEG),
    "GSM5687823": ("Spleen_Acute_D8_CD8_TetPos", "main", _ARM, 8, "spleen", _POS),
    "GSM5687824": (
        "Spleen_Chronic_D21_CD8_CX3CR1Neg_Slamf6Pos",
        "main",
        _C13,
        21,
        "spleen",
        "CD8+ CX3CR1- SLAMF6+ T cells",
    ),
    "GSM5687825": ("Spleen_Chronic_D21_CD8_CX3CR1Neg", "main", _C13, 21, "spleen", "CD8+ CX3CR1- T cells"),
    "GSM5687826": ("Spleen_Chronic_D21_CD8_CX3CR1Pos", "main", _C13, 21, "spleen", "CD8+ CX3CR1+ T cells"),
    "GSM5687827": ("Spleen_Chronic_D21_CD8_PD1Pos", "main", _C13, 21, "spleen", "CD8+ PD1+ T cells"),
    "GSM5687828": ("Spleen_Chronic_D21_CD8_TetMid", "main", _C13, 21, "spleen", _MID),
    "GSM5687829": ("Spleen_Chronic_D21_CD8_TetNeg_rep1", "main", _C13, 21, "spleen", _NEG),
    "GSM5687830": ("Spleen_Chronic_D21_CD8_TetNeg", "main", _C13, 21, "spleen", _NEG),
    "GSM5687831": ("Spleen_Chronic_D21_CD8_TetPos_rep1", "main", _C13, 21, "spleen", _POS),
    "GSM5687832": ("Spleen_Chronic_D21_CD8_TetPos", "main", _C13, 21, "spleen", _POS),
    "GSM5687833": ("Spleen_Chronic_D8_CD8_TetNeg", "main", _C13, 8, "spleen", _NEG),
    "GSM5687834": ("Spleen_Chronic_D8_CD8_TetPos", "main", _C13, 8, "spleen", _POS),
    "GSM6381666": ("Spleen_Chronic_D5_CD8_TetNeg_rep1", "d5", _C13, 5, "spleen", _NEG),
    "GSM6381667": ("Spleen_Chronic_D5_CD8_TetNeg_rep2", "d5", _C13, 5, "spleen", _NEG),
    "GSM6381668": ("Spleen_Chronic_D5_CD8_TetPos_rep1", "d5", _C13, 5, "spleen", _POS),
    "GSM6381669": ("Spleen_Chronic_D5_CD8_TetPos_rep2", "d5", _C13, 5, "spleen", _POS),
}
_END = "_tcr_filtered_contig_annotations.csv.gz"


def test_all_captures_are_pinned_literally() -> None:
    assert len(PINNED) == len(CAPTURES) == 23
    got = {c.gsm: (c.member, c.tag, c.matrix, c.arm, c.day, c.tissue, c.sort_gate) for c in CAPTURES}
    want = {
        g: (f"{g}_{stem}{_END}", f"{stem}_gex_", matrix, arm, day, tissue, gate)
        for g, (stem, matrix, arm, day, tissue, gate) in PINNED.items()
    }
    assert got == want


_GATES = {
    "TetPos": "CD8+ gp33+ T cells",
    "TetNeg": "CD8+ gp33- T cells",
    "TetMid": "CD8+ gp33 mid T cells",
    "CX3CR1Neg": "CD8+ CX3CR1- T cells",
    "CX3CR1Neg_Slamf6Pos": "CD8+ CX3CR1- SLAMF6+ T cells",
    "CX3CR1Pos": "CD8+ CX3CR1+ T cells",
    "PD1Pos": "CD8+ PD1+ T cells",
}


def test_tissue_day_arm_and_gate_agree_with_the_member_name() -> None:
    for c in CAPTURES:
        m = re.fullmatch(rf"{c.gsm}_(\w+?)_(Acute|Chronic)_D(\d+)_CD8_(\w+?)(?:_rep\d)?{_END}", c.member)
        assert m, c.member
        tissue, arm, day, gate = m.groups()
        assert c.tissue == tissue.lower()
        assert c.day == int(day)
        assert c.arm == {"Acute": "LCMV-Armstrong", "Chronic": "LCMV-c13"}[arm]
        assert c.sort_gate == _GATES[gate]
        assert c.tag == c.member.removeprefix(c.gsm + "_").removesuffix(_END) + "_gex_"


def test_the_matrix_orientation_and_column_renames_are_pinned() -> None:
    assert MATRIX_ORIENTATION == {"main": "genes_x_cells", "d5": "genes_x_cells"}
    assert COLUMN_RENAMES == {}


# ---- Task 5b: the build, on the synthetic layout of tests/unit/daniel_fixture.py ----


@pytest.fixture
def layout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]]:
    patch_captures(monkeypatch)
    return write_layout(tmp_path / "src")


@pytest.fixture
def counts(
    layout: tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]],
) -> dict[tuple[str, str], np.ndarray]:
    return layout[1]


def _build(paths: dict[str, Path], out_dir: Path, params: dict[str, Any] = PARAMS) -> ad.AnnData:
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "daniel.h5ad"
    build(paths, out, params)
    return ad.read_h5ad(out)


@pytest.fixture
def built(layout: tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]], tmp_path: Path) -> ad.AnnData:
    return _build(layout[0], tmp_path / "out")


def test_duplicate_barcodes_across_captures_stay_apart(
    built: ad.AnnData, counts: dict[tuple[str, str], np.ndarray]
) -> None:
    a, b = "GSM_A_AAACCTGAGAAGGCCT", "GSM_B_AAACCTGAGAAGGCCT"
    assert built.obs.loc[a, "clone_id"] != built.obs.loc[b, "clone_id"]
    np.testing.assert_array_equal(built[a].X.toarray().ravel(), counts[("GSM_A", "AAACCTGAGAAGGCCT")])
    np.testing.assert_array_equal(built[b].X.toarray().ravel(), counts[("GSM_B", "AAACCTGAGAAGGCCT")])


def test_a_barcode_join_miss_is_counted_and_dropped(built: ad.AnnData) -> None:
    assert built.uns["lumencell_build"]["join_misses"] == {"GSM_A": 2, "GSM_B": 1}
    assert not any(MISSED_CORE in n for n in built.obs_names)


def test_clone_columns_are_scoring_only(built: ad.AnnData) -> None:
    so = scoring_only_columns(built.obs.columns, ("donor",))
    assert {"clone_id", "clone_id_trb", "raw_clonotype_id", "receptor_n_tra", "receptor_n_trb", "donor"} <= so
    assert not {"state", "day", "tissue", "capture", "sort_gate", "arm"} & so


def test_the_kept_cells_are_those_in_the_metadata_and_a_listed_capture(built: ad.AnnData) -> None:
    assert list(built.obs_names) == [
        f"GSM_A_{c}" for c in (CORE, THREE_TRB, BAD_TRA, NO_CONTIG, NOT_CELL)
    ] + [f"GSM_B_{c}" for c in (CORE, B2)]
    assert built.X.dtype == np.float32 and sp.issparse(built.X) and built.X.format == "csr"
    assert list(built.var_names) == GENES and "gene_symbols" not in built.var
    obs = built.obs
    assert obs.loc[f"GSM_A_{CORE}", "state"] == "Tpex" and obs.loc[f"GSM_B_{B2}", "state"] == "Teff"
    a = obs.loc[f"GSM_A_{CORE}"]
    assert (a.capture, a.arm, a.day, a.tissue, a.sort_gate, a.donor) == (
        "GSM_A",
        "LCMV-c13",
        21,
        "spleen",
        "CD8+ gp33+ T cells",
        "m1",
    )
    assert obs.loc[f"GSM_B_{CORE}", "day"] == 8 and obs.loc[f"GSM_B_{CORE}", "donor"] == "m2"
    u = built.uns["lumencell_build"]
    assert u["version"] == BUILD_VERSION == 1 and u["n_cells"] == 7
    assert u["contig_barcodes"] == {"GSM_A": 7, "GSM_B": 4}
    assert (u["labelled_cells"], u["labelled_cells_trb"]) == (3, 4)
    assert list(u["params"]["captures"]) == ["GSM_A", "GSM_B"]
    assert dict(u["params"]["namespace"]) == {"GSM_A": "m1", "GSM_B": "m2"}


def test_cells_without_contigs_are_unlabelled(built: ad.AnnData) -> None:
    row = built.obs.loc[f"GSM_A_{NO_CONTIG}"]
    assert pd.isna(row.clone_id) and pd.isna(row.clone_id_trb) and pd.isna(row.raw_clonotype_id)
    assert row.receptor_n_tra == 0 and row.receptor_n_trb == 0
    assert built.obs["receptor_n_tra"].dtype == np.int64 and built.obs["receptor_n_trb"].dtype == np.int64


def test_the_clone_call_runs_through_the_build(built: ad.AnnData) -> None:
    three, bad = built.obs.loc[f"GSM_A_{THREE_TRB}"], built.obs.loc[f"GSM_A_{BAD_TRA}"]
    assert pd.isna(three.clone_id) and pd.isna(bad.clone_id)
    assert pd.notna(bad.clone_id_trb) and bad.receptor_n_tra == 0
    assert three.receptor_n_trb == 3 and pd.isna(three.clone_id_trb)
    paired = built.obs.loc[f"GSM_A_{CORE}"]
    assert str(paired.clone_id).startswith("c") and paired.raw_clonotype_id == "GSM_A:clonotype1"


def test_only_listed_captures_are_kept_and_namespaced(
    layout: tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]], tmp_path: Path
) -> None:
    one = _build(layout[0], tmp_path / "one", {"captures": ["GSM_A"], "namespace": {"GSM_A": "m1"}})
    assert set(one.obs["capture"]) == {"GSM_A"} and set(one.obs["donor"]) == {"m1"}
    assert one.uns["lumencell_build"]["join_misses"] == {"GSM_A": 2}
    assert list(one.uns["lumencell_build"]["ignored_members"]) == [UNSAFE_MEMBER]


@pytest.mark.parametrize(
    "params",
    [
        {"captures": [], "namespace": {}},
        {"captures": ["GSM_A", "GSM_Z"], "namespace": {"GSM_A": "m1", "GSM_Z": "m9"}},
        {"captures": ["GSM_A", "GSM_B"], "namespace": {"GSM_A": "m1"}},
        {"captures": ["GSM_A"], "namespace": {"GSM_A": "m1", "GSM_B": "m2"}},
        {"captures": ["GSM_A", "GSM_A"], "namespace": {"GSM_A": "m1"}},
        {"captures": ["GSM_A"], "namespace": {"GSM_A": ""}},
        {"captures": "GSM_A", "namespace": {"GSM_A": "m1"}},
        {"captures": ["GSM_A"]},
    ],
    ids=[
        "empty",
        "unknown",
        "namespace-missing",
        "namespace-extra",
        "duplicate",
        "blank-donor",
        "string",
        "no-ns",
    ],
)
def test_bad_params_are_refused(
    layout: tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]], tmp_path: Path, params: dict[str, Any]
) -> None:
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    with pytest.raises(BuildError) as e:
        build(layout[0], out_dir / "daniel.h5ad", params)
    assert "GSM_Z" not in str(e.value) and list(out_dir.iterdir()) == []


def test_the_stray_member_is_ignored_and_nothing_is_written_beside_out(
    layout: tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]], tmp_path: Path
) -> None:
    before = sorted(p.name for p in (tmp_path / "src").iterdir())
    built = _build(layout[0], tmp_path / "out")
    ignored = list(built.uns["lumencell_build"]["ignored_members"])
    assert ignored == [UNSAFE_MEMBER] and not any("stray" in str(n) for n in ignored)
    assert [p.name for p in (tmp_path / "out").iterdir()] == ["daniel.h5ad"]
    assert sorted(p.name for p in (tmp_path / "src").iterdir()) == before
    assert not (tmp_path / "stray.csv.gz").exists()


def test_the_output_loads_as_counts_with_receptor_symbols(
    layout: tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]], tmp_path: Path
) -> None:
    _build(layout[0], tmp_path / "out")
    ds = load(tmp_path / "out" / "daniel.h5ad")
    assert ds.x_kind == "counts"
    assert receptor_gene_mask(ds.gene_symbols).sum() == 2


def test_a_d5_capture_reads_the_d5_matrix_and_labels(
    layout: tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]], tmp_path: Path
) -> None:
    params = {"captures": ["GSM_A", "GSM_D"], "namespace": {"GSM_A": "m1", "GSM_D": "m1"}}
    both = _build(layout[0], tmp_path / "both", params)
    counts = layout[1]
    a, d = f"GSM_A_{CORE}", f"GSM_D_{CORE}"
    # GSM_A reads the main matrix even though the D5 matrix repeats its cell (no dedupe, no D5 copy).
    np.testing.assert_array_equal(both[a].X.toarray().ravel(), counts[("GSM_A", CORE)])
    np.testing.assert_array_equal(both[d].X.toarray().ravel(), counts[("GSM_D", CORE)])
    assert both.obs.loc[a, "state"] == "Tpex" and both.obs.loc[d, "state"] == "D5_early"
    assert both.obs.loc[d, "day"] == 5 and both.n_obs == 6
    # Same receptors in two captures: the same id string, kept apart by capture and donor columns.
    assert both.obs.loc[a, "clone_id"] == both.obs.loc[d, "clone_id"]
    assert both.uns["lumencell_build"]["join_misses"] == {"GSM_A": 2, "GSM_D": 1}


def test_matrices_with_different_features_are_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    patch_captures(monkeypatch)
    paths, _ = write_layout(tmp_path / "src", d5_genes=[*GENES[:-1], "Pdcd1"])
    params = {"captures": ["GSM_A", "GSM_D"], "namespace": {"GSM_A": "m1", "GSM_D": "m2"}}
    with pytest.raises(BuildError, match="feature"):
        _build(paths, tmp_path / "out", params)
    assert list((tmp_path / "out").iterdir()) == []
    # One matrix alone needs no comparison.
    assert _build(paths, tmp_path / "d5", {"captures": ["GSM_D"], "namespace": {"GSM_D": "m2"}}).n_obs == 1


def test_features_with_an_id_and_a_symbol_column_keep_both(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    patch_captures(monkeypatch)
    ids = [f"ENSMUSG{i:011d}\t{g}\tGene Expression" for i, g in enumerate(GENES)]
    paths, _ = write_layout(tmp_path / "src", features=ids)
    out = _build(paths, tmp_path / "out")
    assert list(out.var_names) == [f"ENSMUSG{i:011d}" for i in range(len(GENES))]
    assert list(out.var["gene_symbols"]) == GENES


def test_a_failed_build_leaves_no_output(
    layout: tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]], tmp_path: Path
) -> None:
    paths = dict(layout[0])
    bad = tmp_path / "bad.mtx.gz"
    bad.write_bytes(gzip.compress(b"%%MatrixMarket matrix coordinate integer general\n6 8 1\n1 1\n"))
    paths["GSE188666_scRNA_LCMV.mtx.gz"] = bad
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    with pytest.raises(BuildError):
        build(paths, out_dir / "daniel.h5ad", PARAMS)
    assert list(out_dir.iterdir()) == []


def test_without_a_label_column_there_is_no_state(
    layout: tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(daniel2022, "LABEL_COLUMN", None)
    out = _build(layout[0], tmp_path / "out")
    assert "state" not in out.obs and out.n_obs == 7


# ---- Task 5b fix round 1 ----

Layout = tuple[dict[str, Path], dict[tuple[str, str], np.ndarray]]
PARAMS_AD = {"captures": ["GSM_A", "GSM_D"], "namespace": {"GSM_A": "m1", "GSM_D": "m2"}}
P = "GSE188666_scRNA_"
MAIN_MTX, MAIN_BARCODES, MAIN_FEATURES, MAIN_META = (
    f"{P}LCMV.mtx.gz",
    f"{P}LCMV_barcodes.tsv.gz",
    f"{P}LCMV_features.tsv.gz",
    f"{P}LCMV_metadata.tsv.gz",
)
D5_FEATURES = f"{P}D5_LCMV_features.tsv.gz"


def _edit(path: Path, fn: Any) -> None:
    gz(path, fn(gzip.decompress(path.read_bytes()).decode()))


def _mtx(path: Path, m: np.ndarray) -> None:
    buf = io.BytesIO()
    scipy.io.mmwrite(buf, scipy.sparse.coo_matrix(m))
    path.write_bytes(gzip.compress(buf.getvalue()))


def _refused(paths: dict[str, Path], tmp_path: Path, params: dict[str, Any] = PARAMS) -> str:
    out_dir = tmp_path / "refused"
    out_dir.mkdir()
    with pytest.raises(BuildError) as e:
        build(paths, out_dir / "daniel.h5ad", params)
    assert list(out_dir.iterdir()) == []
    msg = str(e.value)
    assert str(tmp_path) not in msg
    return msg


def test_a_miss_whose_only_qualifying_contig_is_a_tra_is_counted(built: ad.AnnData) -> None:
    # MISSED_CORE (TRA + TRB) and TRA_MISS (a lone TRA): both are misses in GSM_A.
    assert built.uns["lumencell_build"]["join_misses"]["GSM_A"] == 2


def test_misses_are_counted_under_their_own_capture_and_never_leak_across(
    layout: Layout, tmp_path: Path
) -> None:
    # GSM_B has contigs on THREE_TRB, the core of a GSM_A cell only: a GSM_B miss, not matched to GSM_A.
    b = _build(layout[0], tmp_path / "b")
    assert b.uns["lumencell_build"]["join_misses"] == {"GSM_A": 2, "GSM_B": 1}
    assert f"GSM_B_{THREE_TRB}" not in b.obs_names and f"GSM_A_{THREE_TRB}" in b.obs_names
    d = _build(layout[0], tmp_path / "d", PARAMS_AD)
    assert d.uns["lumencell_build"]["join_misses"] == {"GSM_A": 2, "GSM_D": 1}


def test_is_cell_false_contigs_give_no_clone_and_no_join_miss(built: ad.AnnData) -> None:
    row = built.obs.loc[f"GSM_A_{NOT_CELL}"]
    assert pd.isna(row.clone_id) and pd.isna(row.clone_id_trb) and pd.isna(row.raw_clonotype_id) is False
    assert row.receptor_n_tra == 0 and row.receptor_n_trb == 0
    # NOT_CELL_MISS has no cell and only is_cell=False contigs: counted as a barcode, not as a miss.
    assert built.uns["lumencell_build"]["join_misses"]["GSM_A"] == 2
    assert not any(NOT_CELL_MISS in n for n in built.obs_names)


def test_capture_is_the_gsm_prefix_of_every_obs_name(built: ad.AnnData) -> None:
    prefix = [n.rsplit("_", 1)[0] for n in built.obs_names]
    assert list(built.obs["capture"].astype(str)) == prefix


def test_arm_day_tissue_gate_and_donor_follow_each_capture(layout: Layout, tmp_path: Path) -> None:
    both = _build(
        layout[0],
        tmp_path / "all",
        {"captures": ["GSM_A", "GSM_B", "GSM_D"], "namespace": {"GSM_A": "m1", "GSM_B": "m2", "GSM_D": "m3"}},
    )
    by = {c.gsm: c for c in FIX_CAPTURES}
    donor = {"GSM_A": "m1", "GSM_B": "m2", "GSM_D": "m3"}
    assert len({(c.arm, c.tissue) for c in FIX_CAPTURES}) == 3
    for name, row in both.obs.iterrows():
        c = by[str(row["capture"])]
        got = (row["arm"], row["day"], row["tissue"], row["sort_gate"], row["donor"])
        assert got == (c.arm, c.day, c.tissue, c.sort_gate, donor[c.gsm]), name


def test_the_obs_day_is_int64(built: ad.AnnData) -> None:
    assert built.obs["day"].dtype == np.int64


def test_a_blank_label_becomes_nan(built: ad.AnnData) -> None:
    assert pd.isna(built.obs.loc[f"GSM_A_{NO_CONTIG}", "state"])
    assert built.obs["state"].isna().sum() == 1


# -- the P12 refusals: a listed capture that joins nothing --


def test_a_listed_capture_with_no_kept_cell_is_refused(
    layout: Layout, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A wrong tag: no GSM_A cell is in the matrix.
    wrong = tuple(
        c if c.gsm != "GSM_A" else dataclasses.replace(c, tag="Q_gex_") for c in daniel2022.CAPTURES
    )
    monkeypatch.setattr(daniel2022, "CAPTURES", wrong)
    msg = _refused(layout[0], tmp_path)
    assert "GSM_A" in msg and "m1" not in msg


def test_a_listed_capture_whose_metadata_keeps_nothing_is_refused(layout: Layout, tmp_path: Path) -> None:
    _edit(
        layout[0][MAIN_META],
        lambda t: "\n".join(r for r in t.splitlines() if "B_gex_" in r or "orig" in r) + "\n",
    )
    assert "GSM_A" in _refused(layout[0], tmp_path)


def test_a_capture_whose_qualifying_barcodes_all_miss_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    patch_captures(monkeypatch)
    changed = {**CONTIGS, "GSM_A": [r.replace("-1,True,", "-2,True,", 1) for r in CONTIGS["GSM_A"]]}
    paths, _ = write_layout(tmp_path / "src", contigs=changed)
    msg = _refused(paths, tmp_path)
    assert "GSM_A" in msg and "GSM_B" not in msg


# -- the refusals of malformed inputs --


def test_a_metadata_cell_absent_from_its_matrix_is_refused(layout: Layout, tmp_path: Path) -> None:
    _edit(
        layout[0][MAIN_META],
        lambda t: t + "A_gex_ACGTACGTACGTACGT-1\tA\t100\tTGTPLANTED\tTex\n",
    )
    msg = _refused(layout[0], tmp_path)
    assert "metadata" in msg and "GSM_A" in msg


@pytest.mark.parametrize("bad", [1.5, -1], ids=["non-integer", "negative"])
def test_non_integer_and_negative_counts_are_refused(layout: Layout, tmp_path: Path, bad: float) -> None:
    m = np.ones((len(GENES), 9))
    m[2, 3] = bad
    _mtx(layout[0][MAIN_MTX], m)
    assert "counts" in _refused(layout[0], tmp_path)


@pytest.mark.parametrize("shape", [(len(GENES) - 1, 9), (len(GENES), 8)], ids=["features", "barcodes"])
def test_a_matrix_shape_that_disagrees_with_barcodes_and_features_is_refused(
    layout: Layout, tmp_path: Path, shape: tuple[int, int]
) -> None:
    _mtx(layout[0][MAIN_MTX], np.ones(shape))
    assert " x " in _refused(layout[0], tmp_path)


def test_duplicate_matrix_barcodes_are_refused(layout: Layout, tmp_path: Path) -> None:
    _edit(layout[0][MAIN_BARCODES], lambda t: t.replace("B_gex_CCCCCTGAGAAGGCCT", "A_gex_AAACCTGAGAAGGCCT"))
    assert "barcode" in _refused(layout[0], tmp_path)


def test_duplicate_feature_names_are_refused(layout: Layout, tmp_path: Path) -> None:
    _edit(layout[0][MAIN_FEATURES], lambda t: t.replace("Cd8a", "Xkr4"))
    assert "repeats" in _refused(layout[0], tmp_path)


def test_non_unique_metadata_ids_are_refused(layout: Layout, tmp_path: Path) -> None:
    _edit(layout[0][MAIN_META], lambda t: t + t.splitlines()[1] + "\n")
    assert "cell id" in _refused(layout[0], tmp_path)


def test_metadata_without_an_id_column_is_refused(layout: Layout, tmp_path: Path) -> None:
    # Every row now has as many fields as the header: no field is left over to be the cell id.
    _edit(layout[0][MAIN_META], lambda t: t.replace("orig.ident", "id\torig.ident", 1))
    assert "cell id" in _refused(layout[0], tmp_path)


def test_features_mixing_one_and_two_fields_are_refused(layout: Layout, tmp_path: Path) -> None:
    _edit(layout[0][MAIN_FEATURES], lambda t: t.replace("Cd8a", "ENSM1\tCd8a", 1))
    assert "mixes" in _refused(layout[0], tmp_path)


def test_d5_features_in_a_different_order_are_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    patch_captures(monkeypatch)
    paths, _ = write_layout(tmp_path / "src", d5_genes=[*GENES[1:], GENES[0]])
    assert "feature" in _refused(paths, tmp_path, PARAMS_AD)


def test_d5_features_with_the_same_ids_but_other_symbols_are_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    patch_captures(monkeypatch)
    ids = [f"ENSMUSG{i:011d}\t{g}" for i, g in enumerate(GENES)]
    paths, _ = write_layout(tmp_path / "src", features=ids)
    _edit(paths[D5_FEATURES], lambda t: t.replace("\tCd8a", "\tCd8b"))
    assert "feature" in _refused(paths, tmp_path, PARAMS_AD)


@pytest.mark.parametrize(
    "params",
    [
        {"captures": ["GSM_A"], "namespace": ["m1"]},
        {"captures": ["GSM_A"], "namespace": {"GSM_A": 5}},
        {"captures": ["GSM_A"], "namespace": {"GSM_A": "  "}},
    ],
    ids=["list", "non-str-donor", "whitespace-donor"],
)
def test_malformed_namespaces_are_refused(layout: Layout, tmp_path: Path, params: dict[str, Any]) -> None:
    _refused(layout[0], tmp_path, params)


def test_errors_never_carry_donor_values_or_paths(layout: Layout, tmp_path: Path) -> None:
    msg = _refused(
        layout[0], tmp_path, {"captures": ["GSM_A", "GSM_B"], "namespace": {"GSM_A": "m1", "GSM_X": "m2"}}
    )
    assert "m1" not in msg and "m2" not in msg
    bad = tmp_path / "bad.mtx.gz"
    bad.write_bytes(gzip.compress(b"%%MatrixMarket matrix coordinate integer general\n6 9 1\n1 1\n"))
    paths = {**layout[0], MAIN_MTX: bad}
    (tmp_path / "refused").rmdir()
    msg = _refused(paths, tmp_path)
    assert "m1" not in msg and "m2" not in msg and "bad.mtx" not in msg


# -- the atomic write --


def test_a_failing_write_leaves_no_output_and_no_temp_file(
    layout: Layout, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def partial(self: ad.AnnData, filename: Path, *a: Any, **k: Any) -> None:
        Path(filename).write_bytes(b"partial")
        raise OSError("disk full")

    monkeypatch.setattr(ad.AnnData, "write_h5ad", partial)
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    with pytest.raises(OSError, match="disk full"):
        build(layout[0], out_dir / "daniel.h5ad", PARAMS)
    assert list(out_dir.iterdir()) == []
    # An earlier output survives a failed rebuild unchanged.
    (out_dir / "daniel.h5ad").write_bytes(b"earlier")
    with pytest.raises(OSError, match="disk full"):
        build(layout[0], out_dir / "daniel.h5ad", PARAMS)
    assert [p.name for p in out_dir.iterdir()] == ["daniel.h5ad"]
    assert (out_dir / "daniel.h5ad").read_bytes() == b"earlier"


# -- memory of the matrix read (P13) --


def test_reading_an_integer_matrix_makes_no_float_temporaries(tmp_path: Path) -> None:
    rng = np.random.default_rng(1)
    n_genes, n_cells = 400, 5000
    m = sp.random(
        n_genes,
        n_cells,
        density=0.3,
        format="coo",
        random_state=rng,
        data_rvs=lambda k: rng.integers(1, 9, k),
    )
    buf = io.BytesIO()
    scipy.io.mmwrite(buf, m)
    mtx = tmp_path / "GSE188666_scRNA_LCMV.mtx.gz"
    mtx.write_bytes(gzip.compress(buf.getvalue()))
    gz(
        tmp_path / "GSE188666_scRNA_LCMV_barcodes.tsv.gz",
        "".join(f"A_gex_{i:016d}-1\n" for i in range(n_cells)),
    )
    gz(tmp_path / "GSE188666_scRNA_LCMV_features.tsv.gz", "".join(f"g{i}\n" for i in range(n_genes)))
    paths = {p.name: p for p in tmp_path.iterdir()}
    nnz = m.nnz
    tracemalloc.start()
    try:
        daniel2022._read_matrix(paths, "main")
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    # Measured 28.7 bytes a value now (mmread's 16 plus the CSR's 12); the old float64 check took it to 37.7.
    assert peak < nnz * 33
