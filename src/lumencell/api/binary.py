"""LCB1: the binary container for per-cell arrays.

Every endpoint that returns per-cell or per-gene data (embeddings, obs columns,
expression, trajectory cells, heatmaps) uses this one format, so the browser
decodes everything with one function and gets zero-copy TypedArray views.

Layout (all integers little-endian)::

    offset  size          content
    0       4             magic b"LCB1"
    4       4             uint32 format version (1)
    8       4             uint32 manifest byte length (M)
    12      4             uint32 reserved, 0
    16      M             UTF-8 JSON manifest
    16+M    pad           zero bytes up to the next multiple of 8 (= data start)
    ...                   array blocks, each starting on an 8-byte boundary

The manifest is ``{"arrays": [...], "meta": {...}}``. Each array entry is
``{"name", "dtype", "shape", "offset", "byteLength"}`` where ``offset`` is
relative to the data start. Arrays are C-ordered (row-major), so a
``(lineages, cells)`` array gives each lineage a contiguous slice.

The 8-byte alignment matters: ``new Float32Array(buffer, offset, n)`` throws
unless ``offset`` is a multiple of the element size.
"""

import json
import struct
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

MAGIC = b"LCB1"
VERSION = 1
HEADER = struct.Struct("<4sIII")
ALIGN = 8
MEDIA_TYPE = "application/vnd.lumencell.lcb1"

# Explicit little-endian dtypes: the wire format must not depend on host byte order.
DTYPES: dict[str, np.dtype[Any]] = {
    "u8": np.dtype("<u1"),
    "i8": np.dtype("<i1"),
    "u16": np.dtype("<u2"),
    "i16": np.dtype("<i2"),
    "u32": np.dtype("<u4"),
    "i32": np.dtype("<i4"),
    "f32": np.dtype("<f4"),
}


class LcbError(ValueError):
    """A buffer is not a valid LCB1 container, or an array cannot be encoded."""


@dataclass(frozen=True)
class Decoded:
    arrays: dict[str, npt.NDArray[Any]]
    meta: dict[str, Any]
    version: int


def _align(n: int) -> int:
    return (n + ALIGN - 1) // ALIGN * ALIGN


def encode(arrays: Mapping[str, tuple[str, npt.ArrayLike]], meta: Mapping[str, Any] | None = None) -> bytes:
    """Encode named arrays as an LCB1 buffer.

    ``arrays`` maps a name to ``(dtype, values)``. The dtype is explicit rather
    than inferred so callers decide the wire type (e.g. category codes as u8 when
    there are few categories). Values are cast with ``astype``; a cast that
    would silently wrap integers raises instead.
    """
    entries: list[dict[str, Any]] = []
    blocks: list[bytes] = []
    cursor = 0
    for name, (dtype, values) in arrays.items():
        if dtype not in DTYPES:
            raise LcbError(f"array {name!r}: unsupported dtype {dtype!r}; use one of {sorted(DTYPES)}")
        target = DTYPES[dtype]
        src = np.asarray(values)
        if src.dtype.kind in "iub" and target.kind in "iu" and src.size:
            info = np.iinfo(target)
            lo, hi = int(src.min()), int(src.max())
            if lo < info.min or hi > info.max:
                raise LcbError(f"array {name!r}: values [{lo}, {hi}] do not fit in {dtype}")
        if src.dtype.kind == "f" and target.kind in "iu":
            raise LcbError(f"array {name!r}: refusing to cast floating values to {dtype}")
        data = np.ascontiguousarray(src.astype(target, copy=False)).tobytes()
        pad = _align(cursor) - cursor
        if pad:
            blocks.append(b"\x00" * pad)
            cursor += pad
        entries.append(
            {
                "name": name,
                "dtype": dtype,
                "shape": [int(d) for d in src.shape],
                "offset": cursor,
                "byteLength": len(data),
            }
        )
        blocks.append(data)
        cursor += len(data)

    names = [e["name"] for e in entries]
    if len(set(names)) != len(names):
        raise LcbError("array names must be unique")

    manifest = json.dumps(
        {"arrays": entries, "meta": dict(meta or {})},
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    head = HEADER.pack(MAGIC, VERSION, len(manifest), 0) + manifest
    head += b"\x00" * (_align(len(head)) - len(head))
    return head + b"".join(blocks)


def decode(buf: bytes | bytearray | memoryview) -> Decoded:
    """Decode an LCB1 buffer. The returned arrays are read-only views into ``buf``."""
    mv = memoryview(buf).cast("B")
    if len(mv) < HEADER.size:
        raise LcbError("buffer is shorter than the LCB1 header")
    magic, version, manifest_len, _reserved = HEADER.unpack_from(mv, 0)
    if magic != MAGIC:
        raise LcbError(f"bad magic {bytes(magic)!r}")
    if version != VERSION:
        raise LcbError(f"unsupported LCB1 version {version}")
    manifest_end = HEADER.size + manifest_len
    if len(mv) < manifest_end:
        raise LcbError("buffer is truncated inside the manifest")
    try:
        manifest = json.loads(bytes(mv[HEADER.size : manifest_end]).decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise LcbError(f"manifest is not valid UTF-8 JSON: {exc}") from exc
    data_start = _align(manifest_end)

    arrays: dict[str, npt.NDArray[Any]] = {}
    for entry in manifest.get("arrays", []):
        name, dtype = entry["name"], entry["dtype"]
        if dtype not in DTYPES:
            raise LcbError(f"array {name!r}: unknown dtype {dtype!r}")
        np_dtype = DTYPES[dtype]
        shape = tuple(int(d) for d in entry["shape"])
        count = int(np.prod(shape, dtype=np.int64)) if shape else 1
        start = data_start + int(entry["offset"])
        if start % np_dtype.itemsize:
            raise LcbError(f"array {name!r} is not aligned to its element size")
        if int(entry["byteLength"]) != count * np_dtype.itemsize:
            raise LcbError(f"array {name!r}: byteLength does not match shape and dtype")
        if start + count * np_dtype.itemsize > len(mv):
            raise LcbError(f"array {name!r} runs past the end of the buffer")
        arr = np.frombuffer(mv, dtype=np_dtype, count=count, offset=start).reshape(shape)
        arrays[name] = arr
    return Decoded(arrays=arrays, meta=dict(manifest.get("meta", {})), version=int(version))
