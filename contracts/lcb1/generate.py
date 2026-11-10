"""(Re)generate the LCB1 cross-language fixture.

    uv run python contracts/lcb1/generate.py

Writes ``fixture.lcb`` (the encoded buffer) and ``fixture.json`` (the values both
decoders must produce). The Python and TypeScript test suites both read these
files, so a change to either codec that alters the bytes fails both suites until
the fixture is regenerated on purpose.
"""

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from lumencell.api import binary

HERE = Path(__file__).resolve().parent

F32_MAX = float(np.finfo(np.float32).max)


def fixture_arrays() -> dict[str, tuple[str, Any]]:
    # Ordering is deliberate: odd-length u8/i8/i16 blocks sit in front of wider
    # types, so the encoder has to insert alignment padding between blocks.
    return {
        "codes_u8": ("u8", np.array([0, 1, 255], dtype=np.uint8)),
        "pseudotime": ("f32", np.array([0.0, 1.5, -2.25, np.nan, np.inf, -np.inf, F32_MAX])),
        "small_i8": ("i8", np.array([-128, 0, 127])),
        "primary_lineage": ("i16", np.array([-2, -1, 0, 3, 32767])),
        "wide_u16": ("u16", np.array([0, 65535])),
        "segment": ("i32", np.array([-2147483648, 0, 2147483647])),
        "counts_u32": ("u32", np.array([0, 4294967295], dtype=np.uint64)),
        "matrix": ("f32", np.arange(6, dtype=np.float64).reshape(2, 3) + 0.5),
        "membership": ("u8", np.array([[1, 0, 1], [0, 1, 1]], dtype=bool)),
        "empty": ("f32", np.zeros((0,), dtype=np.float32)),
    }


FIXTURE_META: dict[str, Any] = {"n_obs": 7, "label": "Café β-cells ✓", "nested": {"ok": True}}


def _json_value(v: float | int) -> float | int | str:
    if isinstance(v, float):
        if math.isnan(v):
            return "NaN"
        if math.isinf(v):
            return "Infinity" if v > 0 else "-Infinity"
    return v


def expected(arrays: dict[str, tuple[str, Any]], meta: dict[str, Any]) -> dict[str, Any]:
    out: dict[str, Any] = {"meta": meta, "arrays": {}}
    for name, (dtype, values) in arrays.items():
        cast = np.asarray(values).astype(binary.DTYPES[dtype])
        out["arrays"][name] = {
            "dtype": dtype,
            "shape": list(cast.shape),
            "values": [_json_value(v) for v in cast.ravel().tolist()],
        }
    return out


def main() -> None:
    arrays = fixture_arrays()
    (HERE / "fixture.lcb").write_bytes(binary.encode(arrays, FIXTURE_META))
    (HERE / "fixture.json").write_text(
        json.dumps(expected(arrays, FIXTURE_META), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"wrote {HERE / 'fixture.lcb'} and {HERE / 'fixture.json'}")


if __name__ == "__main__":
    main()
