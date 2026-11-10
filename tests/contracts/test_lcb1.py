import importlib.util
import json
import math
import struct
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra import numpy as hnp

from lumencell.api import binary

FIXTURE_DIR = Path(__file__).resolve().parents[2] / "contracts" / "lcb1"


def _load_generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("lcb1_generate", FIXTURE_DIR / "generate.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _from_json(v: Any) -> float:
    return {"NaN": math.nan, "Infinity": math.inf, "-Infinity": -math.inf}.get(v, v)


# ---- shared fixture -----------------------------------------------------------


def test_fixture_decodes_to_expected_values() -> None:
    decoded = binary.decode((FIXTURE_DIR / "fixture.lcb").read_bytes())
    want = json.loads((FIXTURE_DIR / "fixture.json").read_text(encoding="utf-8"))

    assert decoded.meta == want["meta"]
    assert list(decoded.arrays) == list(want["arrays"])
    for name, spec in want["arrays"].items():
        arr = decoded.arrays[name]
        assert arr.dtype == binary.DTYPES[spec["dtype"]], name
        assert list(arr.shape) == spec["shape"], name
        expected = np.array([_from_json(v) for v in spec["values"]], dtype=np.float64)
        np.testing.assert_array_equal(arr.ravel().astype(np.float64), expected, err_msg=name)


def test_encoder_reproduces_fixture_bytes() -> None:
    gen = _load_generator()
    rebuilt = binary.encode(gen.fixture_arrays(), gen.FIXTURE_META)
    assert rebuilt == (FIXTURE_DIR / "fixture.lcb").read_bytes(), (
        "LCB1 encoding changed; regenerate with `uv run python contracts/lcb1/generate.py` "
        "if the change is intended"
    )


def test_fixture_blocks_are_8_byte_aligned() -> None:
    buf = (FIXTURE_DIR / "fixture.lcb").read_bytes()
    _, _, manifest_len, _ = binary.HEADER.unpack_from(buf, 0)
    data_start = binary.HEADER.size + manifest_len
    data_start += (-data_start) % 8
    manifest = json.loads(buf[binary.HEADER.size : binary.HEADER.size + manifest_len])
    assert data_start % 8 == 0
    for entry in manifest["arrays"]:
        assert (data_start + entry["offset"]) % 8 == 0, entry["name"]


# ---- encoder behaviour ----------------------------------------------------------


def test_header_layout() -> None:
    buf = binary.encode({"x": ("f32", [1.0])})
    magic, version, manifest_len, reserved = struct.unpack_from("<4sIII", buf, 0)
    assert (magic, version, reserved) == (b"LCB1", 1, 0)
    assert manifest_len > 0


def test_output_is_little_endian_for_big_endian_input() -> None:
    big = np.array([1.0, 2.0], dtype=">f4")
    decoded = binary.decode(binary.encode({"x": ("f32", big)}))
    assert decoded.arrays["x"].dtype.byteorder in ("<", "=")
    np.testing.assert_array_equal(decoded.arrays["x"], [1.0, 2.0])


def test_non_contiguous_input() -> None:
    m = np.arange(12, dtype=np.float32).reshape(3, 4)[:, ::2]
    decoded = binary.decode(binary.encode({"m": ("f32", m)}))
    np.testing.assert_array_equal(decoded.arrays["m"], m)


def test_row_major_rows_are_contiguous() -> None:
    m = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    decoded = binary.decode(binary.encode({"m": ("f32", m)}))
    np.testing.assert_array_equal(decoded.arrays["m"][1], [4, 5, 6])


@pytest.mark.parametrize(
    ("dtype", "values"),
    [("u8", [256]), ("u8", [-1]), ("i8", [128]), ("i16", [40000]), ("u16", [-1]), ("i32", [2**31])],
)
def test_integer_overflow_is_rejected(dtype: str, values: list[int]) -> None:
    with pytest.raises(binary.LcbError, match="do not fit"):
        binary.encode({"x": (dtype, np.array(values, dtype=np.int64))})


def test_float_to_int_cast_is_rejected() -> None:
    with pytest.raises(binary.LcbError, match="floating"):
        binary.encode({"x": ("i32", np.array([1.5]))})


def test_unsupported_dtype_is_rejected() -> None:
    with pytest.raises(binary.LcbError, match="unsupported dtype"):
        binary.encode({"x": ("f64", [1.0])})


def test_nan_in_meta_is_rejected() -> None:
    with pytest.raises(ValueError):
        binary.encode({}, {"bad": float("nan")})


# ---- decoder validation ---------------------------------------------------------


def test_bad_magic() -> None:
    buf = bytearray(binary.encode({"x": ("f32", [1.0])}))
    buf[:4] = b"NOPE"
    with pytest.raises(binary.LcbError, match="magic"):
        binary.decode(bytes(buf))


def test_truncated_buffer() -> None:
    buf = binary.encode({"x": ("f32", np.zeros(16))})
    with pytest.raises(binary.LcbError, match="past the end"):
        binary.decode(buf[:-4])


def test_short_buffer() -> None:
    with pytest.raises(binary.LcbError, match="shorter"):
        binary.decode(b"LCB1")


def test_unknown_version() -> None:
    buf = bytearray(binary.encode({}))
    struct.pack_into("<I", buf, 4, 99)
    with pytest.raises(binary.LcbError, match="version"):
        binary.decode(bytes(buf))


# ---- property: anything we can encode round-trips exactly -----------------------

_dtype_names = st.sampled_from(sorted(binary.DTYPES))


@st.composite
def _named_arrays(draw: st.DrawFn) -> dict[str, tuple[str, np.ndarray]]:
    n = draw(st.integers(min_value=0, max_value=5))
    out: dict[str, tuple[str, np.ndarray]] = {}
    for i in range(n):
        name = draw(_dtype_names)
        arr = draw(
            hnp.arrays(
                dtype=binary.DTYPES[name],
                shape=hnp.array_shapes(min_dims=1, max_dims=3, min_side=0, max_side=7),
            )
        )
        out[f"a{i}"] = (name, arr)
    return out


@given(_named_arrays())
def test_round_trip_property(arrays: dict[str, tuple[str, np.ndarray]]) -> None:
    decoded = binary.decode(binary.encode(arrays, {"k": 1}))
    assert decoded.meta == {"k": 1}
    for name, (_dtype, arr) in arrays.items():
        got = decoded.arrays[name]
        assert got.shape == arr.shape
        np.testing.assert_array_equal(got, arr)
