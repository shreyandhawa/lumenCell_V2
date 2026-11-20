"""Shared pieces of the dataset builders: strict contig flags, in-memory tar
members, and the one clone call every builder uses.

Pure functions; nothing here touches the network or writes a file.
"""

import hashlib
import tarfile
from collections.abc import Collection, Iterable
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd

from lumencell.bench.download import SAFE_NAME


class BuildError(RuntimeError):
    """A source could not be turned into a dataset."""


CONTIG_COLUMNS = (
    "barcode",
    "chain",
    "productive",
    "high_confidence",
    "full_length",
    "cdr3_nt",
    "raw_clonotype_id",
)
MAX_PER_CHAIN = 2
# A tar member name that is not a plain file name is listed as this, never as itself.
UNSAFE_MEMBER = "<unsafe member name>"
_NT = r"[ACGTN]+"


def parse_flag(values: pd.Series) -> npt.NDArray[np.bool_]:
    """True only for bool True and the strings true/True/TRUE (stripped)."""

    def one(v: object) -> bool:
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        return isinstance(v, str) and v.strip() in ("True", "true", "TRUE")

    return np.fromiter((one(v) for v in values), dtype=np.bool_, count=len(values))


def read_members(tar_path: Path, wanted: Collection[str]) -> tuple[dict[str, bytes], list[str]]:
    """Bytes of the wanted regular-file members, and the sorted names of all others.

    Member names are only compared to `wanted`, never used as paths; nothing is extracted.
    Wanted names are never in the second list. For a repeated name the last member wins.
    """
    want = set(wanted)
    found: dict[str, bytes] = {}
    other: list[str] = []
    try:
        with tarfile.open(tar_path, "r:") as tf:
            for m in tf:
                if m.name not in want:
                    other.append(m.name)
                elif m.isreg():
                    f = tf.extractfile(m)
                    if f is None:
                        raise BuildError("tar member could not be read")
                    found[m.name] = f.read()
                else:
                    found.pop(m.name, None)  # last one wins: a later non-file hides an earlier file
    except tarfile.TarError as e:
        raise BuildError(f"{tar_path} is not a readable tar file ({type(e).__name__})") from e
    missing = sorted(want - found.keys())
    if missing:
        raise BuildError(f"{len(missing)} wanted tar member(s) missing or not regular files")
    return found, sorted(other)


def listable(names: Iterable[str]) -> list[str]:
    """Member names fit to store or print: a plain file name as itself, any other as ``UNSAFE_MEMBER``."""
    return [n if SAFE_NAME.fullmatch(n) else UNSAFE_MEMBER for n in names]


def _hash(prefix: str, parts: set[str]) -> str:
    return prefix + hashlib.sha256(";".join(sorted(parts)).encode()).hexdigest()[:16]


def call_clones(contigs: pd.DataFrame) -> pd.DataFrame:
    """One row per (capture, barcode): paired and TRB-only clone ids and chain counts.

    A contig qualifies when productive, high_confidence, full_length (and is_cell, if
    present) all parse True, chain is TRA or TRB, and cdr3_nt is fully ACGTN.
    `clone_id` needs 1-2 distinct TRA and 1-2 distinct TRB; `clone_id_trb` needs 1-2 TRB.
    No namespace is in the id strings: the caller keys clones by (capture, id).
    """
    if len(contigs) == 0:
        raise BuildError("no contigs")
    need = ["capture", *CONTIG_COLUMNS]
    absent = [c for c in need if c not in contigs.columns]
    if absent:
        raise BuildError(f"contig table lacks columns: {', '.join(absent)}")
    d = contigs.reset_index(drop=True)
    ok = np.ones(len(d), dtype=bool)
    for col in ("productive", "high_confidence", "full_length", *(["is_cell"] if "is_cell" in d else [])):
        ok &= parse_flag(d[col])
    ok &= d["chain"].isin(["TRA", "TRB"]).to_numpy()
    ok &= d["cdr3_nt"].astype("string").str.fullmatch(_NT).fillna(False).to_numpy(dtype=bool)
    q = d[ok]
    seqs: dict[tuple[str, str], dict[str, set[str]]] = {}
    for cap, bc, chain, nt in zip(q["capture"], q["barcode"], q["chain"], q["cdr3_nt"], strict=True):
        seqs.setdefault((cap, bc), {"TRA": set(), "TRB": set()})[chain].add(nt)

    raw: dict[tuple[str, str], object] = {}
    for cap, bc, rc in zip(d["capture"], d["barcode"], d["raw_clonotype_id"], strict=True):
        key = (cap, bc)
        raw.setdefault(key, np.nan)
        if pd.notna(rc) and str(rc).strip() not in ("", "None") and pd.isna(raw[key]):
            raw[key] = f"{cap}:{rc}"

    index = pd.MultiIndex.from_tuples(sorted(raw), names=["capture", "barcode"])
    rows = []
    for key in index:
        s = seqs.get(key, {"TRA": set(), "TRB": set()})
        n_a, n_b = len(s["TRA"]), len(s["TRB"])
        trb = {f"TRB:{x}" for x in s["TRB"]}
        tra = {f"TRA:{x}" for x in s["TRA"]}
        rows.append(
            {
                "clone_id": _hash("c", tra | trb)
                if 1 <= n_a <= MAX_PER_CHAIN and 1 <= n_b <= MAX_PER_CHAIN
                else np.nan,
                "clone_id_trb": _hash("b", trb) if 1 <= n_b <= MAX_PER_CHAIN else np.nan,
                "raw_clonotype_id": raw[key],
                "receptor_n_tra": n_a,
                "receptor_n_trb": n_b,
            }
        )
    out = pd.DataFrame(rows, index=index)
    for c in ("clone_id", "clone_id_trb", "raw_clonotype_id"):
        out[c] = out[c].astype(object)
    out["receptor_n_tra"] = out["receptor_n_tra"].astype(int)
    out["receptor_n_trb"] = out["receptor_n_trb"].astype(int)
    return out
