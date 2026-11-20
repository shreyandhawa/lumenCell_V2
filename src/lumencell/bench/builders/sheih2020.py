"""Sheih et al. 2020 (CAR-T cells in four patients, scRNA and TCR): the GEO files and the build.

The expression matrix is a genes x cells csv with cell ids ``<16 nt>-<aggr id>``. The TCR tables
(inside ``GSE125881_RAW.tar``) hold one row per clonotype x chain with ``;``-separated barcode
lists (``<16 nt>-1``) and no contig flags, so the clones are called from the per-cell chains the
tables list, with the declared rule (paired, CDR3 nucleotides, no filters, at most two per chain)
and the same hashing as the 10x builders (``call_clones``).
"""

import csv
import gzip
import re
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from lumencell.bench.builders.common import MAX_PER_CHAIN, BuildError, call_clones, listable, read_members
from lumencell.bench.download import Source

BUILD_VERSION = 2
CHUNK_ROWS = 500
JOINED_SHARE_MIN = 0.5

GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125881/suppl/"
MATRIX = "GSE125881_raw.expMatrix.csv.gz"
TAR = "GSE125881_RAW.tar"

# Sizes from HEAD requests (2026-11-20, equal to Task 9's); SHA-256 as recorded on download.
SOURCES: tuple[Source, ...] = (
    Source(
        GEO_SUPPL + MATRIX,
        MATRIX,
        114_573_707,
        "243e6e2dd09158e08dcaecfbf35b23122ac3096c1039355b7590474e76d8c02c",
    ),
    Source(
        GEO_SUPPL + TAR,
        TAR,
        2_048_000,
        "769fa5b71b12b7efe433f51164d4370e0af00c20bebc647029cc4fbd9a943487",
    ),
)

LICENCE = (
    "Article: CC BY 4.0 (Sheih et al. 2020, Nat Commun 11:219, doi:10.1038/s41467-019-13880-1, "
    "PMC6954177). Data: public GEO accession GSE125881; UNVERIFIED: the GEO terms page is behind a "
    "CAPTCHA, not bypassed (read attempt 2026-11-19), so cite the paper and the accession."
)


@dataclass(frozen=True)
class Sample:
    gsm: str  # the expression (GE) GSM; the VDJ GSM is this plus 16
    member: str  # name of the clonotype table inside the tar
    patient: str
    timepoint: str  # GEO's "time point" as written
    day: int | None  # the sampling day in the label; None for the infusion product (not a sampling day)
    aggr_id: int  # the suffix of this sample's cell ids in the matrix


def _s(gsm: int, patient: str, timepoint: str, aggr: int) -> Sample:
    day = None if timepoint == "IP" else int(timepoint[1:])
    member = f"GSM{gsm + 16}_{patient}_{timepoint}-unpaired-clonotypes.csv.gz"
    return Sample(f"GSM{gsm}", member, patient, timepoint, day, aggr)


# Task 9's GSM table (docs/superpowers/research/2026-11-19-sheih-verification.md).
REAL_SAMPLES: tuple[Sample, ...] = (
    _s(3583861, "CLL-1", "IP", 1),
    _s(3583862, "CLL-1", "d21", 2),
    _s(3583863, "CLL-1", "d38", 3),
    _s(3583864, "CLL-1", "d112", 4),
    _s(3583865, "CLL-2", "IP", 8),
    _s(3583866, "CLL-2", "d12", 9),
    _s(3583867, "CLL-2", "d29", 10),
    _s(3583868, "CLL-2", "d83", 11),
    _s(3583869, "NHL-6", "IP", 16),
    _s(3583870, "NHL-6", "d12", 5),
    _s(3583871, "NHL-6", "d29", 6),
    _s(3583872, "NHL-6", "d102", 7),
    _s(3583873, "NHL-7", "IP", 12),
    _s(3583874, "NHL-7", "d12", 13),
    _s(3583875, "NHL-7", "d28", 14),
    _s(3583876, "NHL-7", "d89", 15),
)
SAMPLES: tuple[Sample, ...] = REAL_SAMPLES

_CELL_ID = re.compile(r"([ACGT]{16})-([0-9]{1,3})")
_TABLE_BARCODE = re.compile(r"([ACGT]{16})-1")
_NT = re.compile(r"[ACGTN]+")
_REQUIRED = ("chain", "clonotype_ids", "barcodes", "cdr3_nts")
_DECLARED = {"chains": "paired", "cdr3": "nt", "filters": [], "max_per_chain": MAX_PER_CHAIN}


def _check_params(params: Mapping[str, Any]) -> None:
    """The builder implements one rule; a declaration it cannot honour is refused, never approximated."""
    ok = (
        set(params) == set(_DECLARED)
        and params["chains"] == "paired"
        and params["cdr3"] == "nt"
        and list(params["filters"]) == []
        and params["max_per_chain"] == MAX_PER_CHAIN
    )
    if not ok:
        raise BuildError(
            "params must equal the declared clone definition: paired, nt, no filters, 2 per chain"
        )


def _cell_ids(matrix: Path) -> list[str]:
    with gzip.open(matrix, "rt", newline="") as fh:
        try:
            header = next(csv.reader(fh))
        except StopIteration:
            raise BuildError("the expression matrix is empty") from None
    if not header or header[0] != "" or len(header) < 2:
        raise BuildError("the expression matrix header is not an empty field followed by cell ids")
    cells = header[1:]
    if not all(_CELL_ID.fullmatch(c) for c in cells):
        raise BuildError("the expression matrix has cell ids that are not <16 nt>-<number>")
    if len(set(cells)) != len(cells):
        raise BuildError("the expression matrix has duplicate cell ids")
    return cells


def _table_calls(
    members: Mapping[str, bytes], cells: set[str]
) -> tuple[pd.DataFrame, dict[str, int], dict[str, int], dict[str, int]]:
    """Per-cell chains of every sample's table, joined to the matrix's cell ids.

    Returns the contig-style rows of the joined barcodes and, per sample GSM, the number of distinct
    barcodes, how many of them match no cell of the sample, and the rows that are not usable (a chain
    that is not TRA or TRB, or a CDR3 field that is not one bare ACGTN string; any such row refuses the
    build, so the counts are all 0 in a built file).
    """
    frames: list[pd.DataFrame] = []
    barcodes: dict[str, int] = {}
    misses: dict[str, int] = {}
    unusable: dict[str, int] = {}
    for s in SAMPLES:
        try:
            t = pd.read_csv(BytesIO(members[s.member]), compression="gzip", dtype=str, keep_default_na=False)
        except (pd.errors.ParserError, pd.errors.EmptyDataError, OSError, EOFError, UnicodeDecodeError) as e:
            raise BuildError(f"{s.gsm}: the clonotype table could not be read ({type(e).__name__})") from e
        if any(c not in t.columns for c in _REQUIRED):
            raise BuildError(f"{s.gsm}: the clonotype table lacks one of {', '.join(_REQUIRED)}")
        usable = t["chain"].isin(["TRA", "TRB"]) & t["cdr3_nts"].str.fullmatch(_NT.pattern)
        unusable[s.gsm] = int((~usable).sum())
        if unusable[s.gsm]:
            # A dropped row could make a cell with three TRAs look like it has two, and label it.
            raise BuildError(
                f"{s.gsm}: {unusable[s.gsm]} clonotype row(s) are not TRA/TRB with one bare ACGTN CDR3"
            )
        t = t[usable]
        ex = pd.DataFrame(
            {
                "tb": t["barcodes"].str.split(";"),
                "chain": t["chain"],
                "cdr3_nt": t["cdr3_nts"],
            }
        ).explode("tb")
        ex["tb"] = ex["tb"].str.strip()
        ex = ex[ex["tb"] != ""]
        ex["barcode"] = [
            f"{m.group(1)}-{s.aggr_id}" if (m := _TABLE_BARCODE.fullmatch(b)) else "" for b in ex["tb"]
        ]
        distinct = ex["tb"].nunique()
        if distinct == 0:
            raise BuildError(f"{s.gsm}: the clonotype table has no usable barcodes")
        joined = ex["barcode"].isin(cells)
        barcodes[s.gsm] = int(distinct)
        misses[s.gsm] = int(distinct - ex.loc[joined, "tb"].nunique())
        j = ex[joined]
        frames.append(
            pd.DataFrame(
                {
                    "capture": s.gsm,
                    "barcode": j["barcode"].to_numpy(),
                    "chain": j["chain"].to_numpy(),
                    "productive": "True",
                    "high_confidence": "True",
                    "full_length": "True",
                    "cdr3_nt": j["cdr3_nt"].to_numpy(),
                    "raw_clonotype_id": None,  # per-chain ids in an unpaired table: not a cell's clonotype
                }
            )
        )
    return pd.concat(frames, ignore_index=True), barcodes, misses, unusable


def _read_counts(matrix: Path, n_cells: int, chunk_rows: int) -> tuple[sp.csr_matrix, list[str]]:
    """The genes x cells csv as cells x genes CSR float32, read ``chunk_rows`` genes at a time."""
    blocks: list[sp.csr_matrix] = []
    genes: list[str] = []
    try:
        with gzip.open(matrix, "rt", newline="") as fh:
            fh.readline()  # the header, already read and checked
            reader = pd.read_csv(
                fh,
                header=None,
                index_col=0,
                chunksize=chunk_rows,
                converters={0: str},
                keep_default_na=False,
            )
            for chunk in reader:
                if chunk.shape[1] != n_cells:
                    raise BuildError("an expression matrix row does not have one value per cell")
                block = sp.csr_matrix(chunk.to_numpy(dtype=np.float32))
                data = block.data
                if data.size and (
                    not np.isfinite(data).all() or data.min() < 0 or (data != np.rint(data)).any()
                ):
                    raise BuildError("the expression matrix does not hold non-negative integer counts")
                blocks.append(block)
                genes.extend(map(str, chunk.index))
    except (pd.errors.ParserError, ValueError, OSError, EOFError) as e:
        if isinstance(e, BuildError):
            raise
        raise BuildError(f"the expression matrix could not be read ({type(e).__name__})") from e
    if not blocks:
        raise BuildError("the expression matrix has no gene rows")
    if len(set(genes)) != len(genes):
        raise BuildError("the expression matrix has duplicate gene names")
    genes_by_cells = sp.vstack(blocks, format="csr")
    return genes_by_cells.T.tocsr().astype(np.float32), genes


def build(
    paths: Mapping[str, Path], out: Path, params: Mapping[str, Any], *, chunk_rows: int = CHUNK_ROWS
) -> None:
    _check_params(params)
    if chunk_rows < 1:
        raise BuildError("chunk_rows must be positive")
    matrix = paths[MATRIX]
    cells = _cell_ids(matrix)
    by_aggr = {s.aggr_id: s for s in SAMPLES}
    suffix = [int(_CELL_ID.fullmatch(c).group(2)) for c in cells]  # type: ignore[union-attr]
    if any(a not in by_aggr for a in suffix):
        raise BuildError("the expression matrix has cells whose sample suffix is not a known sample")
    members, ignored = read_members(paths[TAR], {s.member for s in SAMPLES})

    contigs, barcodes, misses, unusable = _table_calls(members, set(cells))
    share = {g: 1 - misses[g] / barcodes[g] for g in barcodes}
    low = [g for g, v in share.items() if v < JOINED_SHARE_MIN]
    if low:
        raise BuildError(
            f"joined share below {JOINED_SHARE_MIN}: " + ", ".join(f"{g} {share[g]:.2f}" for g in low)
        )

    x, genes = _read_counts(matrix, len(cells), chunk_rows)
    sample = [by_aggr[a] for a in suffix]
    obs = pd.DataFrame(
        {
            "patient": pd.Categorical([s.patient for s in sample]),
            "timepoint": pd.Categorical([s.timepoint for s in sample]),
            # NaN for the infusion product: it is not a sampling day, and 0 would correlate by construction
            "day": np.array([np.nan if s.day is None else s.day for s in sample], dtype=np.float64),
        },
        index=pd.Index(cells, name=None),
    )
    called = call_clones(contigs).drop(columns="raw_clonotype_id")
    keys = pd.MultiIndex.from_arrays([[s.gsm for s in sample], cells])
    calls = called.reindex(keys)
    for col in ("clone_id", "clone_id_trb"):
        obs[col] = calls[col].astype(object).to_numpy()
    for col in ("receptor_n_tra", "receptor_n_trb"):
        obs[col] = calls[col].fillna(0).astype(np.int64).to_numpy()

    var = pd.DataFrame(index=pd.Index(genes))
    adata = AnnData(X=x, obs=obs, var=var)
    adata.uns["lumencell_build"] = {
        "version": BUILD_VERSION,
        "params": {k: params[k] for k in _DECLARED},
        "n_cells": len(cells),
        "contig_barcodes": barcodes,
        "join_misses": misses,
        "joined_share": share,
        "cdr3_unusable_rows": unusable,
        "ignored_members": listable(ignored),
        "labelled_cells": int(obs["clone_id"].notna().sum()),
        "labelled_cells_trb": int(obs["clone_id_trb"].notna().sum()),
        "cells_per_sample": {g: int(n) for g, n in Counter(s.gsm for s in sample).items()},
    }
    adata.write_h5ad(out)
