"""Daniel et al. 2022 (chronic LCMV, paired scRNA and TCR): the GEO files, their formats and the build.

The build joins each listed TCR capture's contigs to the cells of its own expression capture (a contig
barcode ``b`` is the matrix cell ``tag + b``), calls clones with ``call_clones`` and writes one ``.h5ad``.
"""

import gzip
import os
import re
import tempfile
import zlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from anndata import AnnData

from lumencell.bench.builders.common import CONTIG_COLUMNS, BuildError, call_clones, read_members
from lumencell.bench.download import SAFE_NAME, Source

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


# ---- The build (Task 5b) ----

BUILD_VERSION = 1
TAR = "GSE188666_RAW.tar"
_STEM: dict[str, str] = {"main": "GSE188666_scRNA_LCMV", "d5": "GSE188666_scRNA_D5_LCMV"}
_METADATA: dict[str, str] = {
    "main": "GSE188666_scRNA_LCMV_metadata.tsv.gz",
    "d5": "GSE188666_scRNA_D5_LCMV_metadata.txt.gz",
}
# A tar member name that is not a plain file name is listed as this, never as itself.
UNSAFE_MEMBER = "<unsafe member name>"
# After its tag, a cell id is the 16-nt core and "-1" (note section 9). obs names keep the core only.
_CORE_AND_SUFFIX = r"([ACGT]{16})-1"
_RECEPTOR_IDS = ("clone_id", "clone_id_trb", "raw_clonotype_id")
_RECEPTOR_COUNTS = ("receptor_n_tra", "receptor_n_trb")
_READ_ERRORS = (
    pd.errors.ParserError,
    pd.errors.EmptyDataError,
    ValueError,
    OSError,
    EOFError,
    zlib.error,
)


@dataclass(frozen=True)
class _Matrix:
    x: sp.csr_matrix  # cells x genes, as stored (integer counts)
    barcodes: list[str]
    features: list[str]
    symbols: list[str] | None


def _check_params(params: Mapping[str, Any]) -> tuple[list[Capture], dict[str, str]]:
    """The listed captures and their namespace. Errors give counts, never the values."""
    if set(params) != {"captures", "namespace"}:
        raise BuildError("params must have exactly the keys captures and namespace")
    caps, ns = params["captures"], params["namespace"]
    if isinstance(caps, str) or not isinstance(caps, Sequence) or not all(isinstance(g, str) for g in caps):
        raise BuildError("params captures must be a list of GSM strings")
    if not caps:
        raise BuildError("params captures is empty")
    if len(set(caps)) != len(caps):
        raise BuildError(f"params captures repeats {len(caps) - len(set(caps))} GSM(s)")
    by_gsm = {c.gsm: c for c in CAPTURES}
    unknown = sum(g not in by_gsm for g in caps)
    if unknown:
        raise BuildError(f"{unknown} listed capture(s) are not in CAPTURES")
    if not isinstance(ns, Mapping):
        raise BuildError("params namespace must map each listed GSM to a string")
    missing, extra = sum(g not in ns for g in caps), sum(k not in set(caps) for k in ns)
    if missing or extra:
        raise BuildError(
            f"the namespace must name exactly the listed captures: {missing} missing, {extra} not listed"
        )
    if not all(isinstance(ns[g], str) and ns[g].strip() for g in caps):
        raise BuildError("every namespace value must be a non-empty string")
    return [by_gsm[g] for g in caps], {g: ns[g] for g in caps}


def _lines(path: Path, name: str) -> list[str]:
    try:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            return [line.rstrip("\r\n") for line in fh]
    except _READ_ERRORS as e:
        raise BuildError(f"{name} could not be read ({type(e).__name__})") from e


def _features(lines: list[str], name: str) -> tuple[list[str], list[str] | None]:
    """Feature names, and the symbols when they have their own column (the second field)."""
    if not lines:
        raise BuildError(f"{name} lists no features")
    fields = [line.split("\t") for line in lines]
    widths = {len(f) for f in fields}
    if widths == {1}:
        names, symbols = lines, None
    elif min(widths) >= 2:
        names, symbols = [f[0] for f in fields], [f[1] for f in fields]
    else:
        raise BuildError(f"{name} mixes lines of one field with lines of several")
    if any(not n for n in names):
        raise BuildError(f"{name} has {sum(not n for n in names)} blank feature name(s)")
    if len(set(names)) != len(names):
        raise BuildError(f"{name} repeats {len(names) - len(set(names))} feature name(s)")
    return names, symbols


def _read_matrix(paths: Mapping[str, Path], which: Literal["main", "d5"]) -> _Matrix:
    stem = _STEM[which]
    barcodes = _lines(paths[f"{stem}_barcodes.tsv.gz"], f"{stem}_barcodes.tsv.gz")
    if not barcodes or len(set(barcodes)) != len(barcodes):
        raise BuildError(f"{stem}_barcodes.tsv.gz is empty or repeats a barcode")
    features, symbols = _features(
        _lines(paths[f"{stem}_features.tsv.gz"], f"{stem}_features.tsv.gz"), f"{stem}_features.tsv.gz"
    )
    name = f"{stem}.mtx.gz"
    try:
        with gzip.open(paths[name], "rb") as fh:
            m = scipy.io.mmread(fh, spmatrix=False)
    except _READ_ERRORS as e:
        raise BuildError(f"{name} could not be read as a MatrixMarket file ({type(e).__name__})") from e
    if not sp.issparse(m):
        raise BuildError(f"{name} is not a coordinate matrix")
    genes_x_cells = MATRIX_ORIENTATION[which] == "genes_x_cells"
    want = (len(features), len(barcodes)) if genes_x_cells else (len(barcodes), len(features))
    if tuple(m.shape) != want:
        raise BuildError(
            f"{name} is {m.shape[0]} x {m.shape[1]}; its barcodes and features make {want[0]} x {want[1]}"
        )
    x = sp.csr_matrix(m.T if genes_x_cells else m)
    del m  # the coordinate arrays are as large as the CSR; drop them before the value check
    if not _are_counts(x.data):
        raise BuildError(f"{name} does not hold non-negative integer counts")
    return _Matrix(x, barcodes, features, symbols)


def _are_counts(d: np.ndarray, step: int = 1 << 20) -> bool:
    """Whether every value is a finite, non-negative whole number, without full-size temporaries."""
    if d.size == 0:
        return True
    if d.dtype.kind in "iu":  # an integer dtype is whole and finite already
        return bool(d.dtype.kind == "u" or d.min() >= 0)
    for i in range(0, d.size, step):
        c = d[i : i + step]
        if not np.isfinite(c).all() or c.min() < 0 or (c != np.rint(c)).any():
            return False
    return True


def _read_metadata(path: Path, name: str) -> pd.Series:
    """The label of each cell id (all None when LABEL_COLUMN is None), reading the index and one column.

    The files are Seurat exports: the header has one field fewer than each row, and the first field of a
    row is the cell id, the same string as the matrix barcode (note sections 2 and 9).
    """
    try:
        header = pd.read_csv(path, sep="\t", nrows=0, compression="gzip")
        if len(header.columns) == 0:
            raise BuildError(f"{name} has no columns")
        if LABEL_COLUMN is not None and LABEL_COLUMN not in header.columns:
            raise BuildError(f"{name} has no column LABEL_COLUMN")
        col = LABEL_COLUMN if LABEL_COLUMN is not None else str(header.columns[0])
        df = pd.read_csv(path, sep="\t", usecols=[col], dtype=str, keep_default_na=False, compression="gzip")
    except _READ_ERRORS as e:
        raise BuildError(f"{name} could not be read ({type(e).__name__})") from e
    if isinstance(df.index, pd.RangeIndex) or not df.index.is_unique:
        raise BuildError(f"{name} does not have one cell id per row in a column of its own")
    labels = df[col].where(df[col] != "") if LABEL_COLUMN is not None else pd.Series(None, index=df.index)
    return pd.Series(labels.to_numpy(dtype=object), index=df.index.astype(str))


def _contigs(members: Mapping[str, bytes], listed: Sequence[Capture]) -> pd.DataFrame:
    """Every contig row of the listed captures, with ``capture`` set to the GSM."""
    frames = []
    for c in listed:
        data = members[c.member]
        try:
            text = gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data
            t = pd.read_csv(BytesIO(text), dtype=str, keep_default_na=False)
        except _READ_ERRORS as e:
            raise BuildError(f"{c.gsm}: the contig table could not be read ({type(e).__name__})") from e
        t = t.rename(columns=COLUMN_RENAMES)
        absent = [k for k in CONTIG_COLUMNS if k not in t.columns]
        if absent:
            raise BuildError(f"{c.gsm}: the contig table lacks {', '.join(absent)}")
        keep = [*CONTIG_COLUMNS, *(["is_cell"] if "is_cell" in t.columns else [])]
        frames.append(t[keep].assign(capture=c.gsm))
    contigs = pd.concat(frames, ignore_index=True)
    if "is_cell" in contigs:
        # Only a member without the column has NaN here (blanks read as ""): no is_cell filter for it.
        contigs["is_cell"] = contigs["is_cell"].fillna("True")
    return contigs


def build(paths: Mapping[str, Path], out: Path, params: Mapping[str, Any]) -> None:
    """Build the listed captures into ``out``, written by atomic rename (a ``Builder``).

    Kept cells are those of a listed capture (matrix id ``tag + core + "-1"``) that the metadata of the
    capture's own matrix lists. Each capture is read from its own matrix only, so the D5 matrix's copies
    of main-matrix cells are never read and nothing is deduplicated. obs names are ``<gsm>_<core>``.
    """
    listed, namespace = _check_params(params)
    used: list[Literal["main", "d5"]] = [w for w in ("main", "d5") if any(c.matrix == w for c in listed)]
    matrices = {w: _read_matrix(paths, w) for w in used}
    if len(matrices) == 2:
        a, b = matrices["main"], matrices["d5"]
        if a.features != b.features or a.symbols != b.symbols:
            raise BuildError("the main and D5 matrices list different features")
    metadata = {w: _read_metadata(paths[_METADATA[w]], _METADATA[w]) for w in used}
    members, others = read_members(paths[TAR], [c.member for c in listed])
    known = {c.member for c in CAPTURES}
    ignored = [n if SAFE_NAME.fullmatch(n) else UNSAFE_MEMBER for n in others if n not in known]
    calls = call_clones(_contigs(members, listed))

    obs_parts: list[pd.DataFrame] = []
    x_parts: list[sp.csr_matrix] = []
    kept: set[tuple[str, str]] = set()
    for c in listed:
        m, meta = matrices[c.matrix], metadata[c.matrix]
        cell = re.compile(re.escape(c.tag) + _CORE_AND_SUFFIX)
        rows: list[int] = []
        cores: list[str] = []
        for i, bc in enumerate(m.barcodes):
            hit = cell.fullmatch(bc)
            if hit is not None and bc in meta.index:
                rows.append(i)
                cores.append(hit.group(1))
        if not rows:
            raise BuildError(
                f"{c.gsm}: no cell of this capture is in both the {c.matrix} matrix and its metadata"
            )
        absent = sum(cell.fullmatch(bc) is not None for bc in meta.index) - len(rows)
        if absent:
            raise BuildError(f"{c.gsm}: {absent} metadata cell(s) are not in the {c.matrix} matrix")
        contig_keys = [f"{core}-1" for core in cores]  # the contig barcode b of the cell tag + b
        kept.update((c.gsm, k) for k in contig_keys)
        cl = calls.reindex(pd.MultiIndex.from_arrays([[c.gsm] * len(rows), contig_keys]))
        part = pd.DataFrame(
            {
                "capture": c.gsm,
                "arm": c.arm,
                "day": np.int64(c.day),
                "tissue": c.tissue,
                "sort_gate": c.sort_gate,
            },
            index=pd.Index([f"{c.gsm}_{core}" for core in cores]),
        )
        if LABEL_COLUMN is not None:
            part["state"] = meta.loc[[m.barcodes[i] for i in rows]].to_numpy(dtype=object)
        part["donor"] = namespace[c.gsm]
        for col in _RECEPTOR_IDS:
            part[col] = cl[col].to_numpy(dtype=object)
        for col in _RECEPTOR_COUNTS:
            part[col] = cl[col].fillna(0).to_numpy(dtype=np.int64)
        obs_parts.append(part)
        x_parts.append(m.x[rows])
    obs = pd.concat(obs_parts)
    for col in ("capture", "arm", "tissue", "sort_gate", "donor", *(["state"] if "state" in obs else [])):
        obs[col] = pd.Categorical(obs[col])
    obs["day"] = obs["day"].astype(np.int64)
    for col in _RECEPTOR_IDS:
        obs[col] = obs[col].astype(object)

    gsm_of = calls.index.get_level_values("capture")
    qualifying = (calls["receptor_n_tra"] + calls["receptor_n_trb"]).to_numpy() > 0
    misses = {c.gsm: 0 for c in listed}
    n_qualifying = dict(misses)
    for key in calls.index[qualifying]:
        n_qualifying[key[0]] += 1
        if key not in kept:
            misses[key[0]] += 1
    for c in listed:
        if n_qualifying[c.gsm] and misses[c.gsm] == n_qualifying[c.gsm]:
            raise BuildError(
                f"{c.gsm}: none of its {n_qualifying[c.gsm]} qualifying contig barcodes joins a kept cell"
            )

    first = matrices[used[0]]
    var = pd.DataFrame(index=pd.Index(first.features))
    if first.symbols is not None:
        var["gene_symbols"] = first.symbols
    adata = AnnData(X=sp.vstack(x_parts, format="csr").astype(np.float32), obs=obs, var=var)
    adata.uns["lumencell_build"] = {
        "version": BUILD_VERSION,
        "params": {"captures": [c.gsm for c in listed], "namespace": namespace},
        "n_cells": int(adata.n_obs),
        "contig_barcodes": {c.gsm: int((gsm_of == c.gsm).sum()) for c in listed},
        "join_misses": misses,
        "ignored_members": ignored,
        "labelled_cells": int(obs["clone_id"].notna().sum()),
        "labelled_cells_trb": int(obs["clone_id_trb"].notna().sum()),
    }
    _write_atomic(adata, out)


def _write_atomic(adata: AnnData, out: Path) -> None:
    fd, name = tempfile.mkstemp(prefix=".daniel2022-", suffix=".tmp.h5ad", dir=out.parent)
    os.close(fd)
    tmp = Path(name)
    try:
        adata.write_h5ad(tmp)
        os.replace(tmp, out)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
