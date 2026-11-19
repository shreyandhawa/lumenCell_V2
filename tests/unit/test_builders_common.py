import gzip
import hashlib
import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lumencell.bench.builders.common import (
    CONTIG_COLUMNS,
    BuildError,
    call_clones,
    parse_flag,
    read_members,
)


def contig(capture: str, barcode: str, chain: str, nt: str, **flags: object) -> dict[str, object]:
    row: dict[str, object] = {
        "capture": capture,
        "barcode": barcode,
        "chain": chain,
        "productive": "True",
        "high_confidence": "True",
        "full_length": "True",
        "cdr3_nt": nt,
        "raw_clonotype_id": "clonotype1",
    }
    row.update(flags)
    return row


def frame(*rows: dict[str, object]) -> pd.DataFrame:
    return pd.DataFrame(list(rows), columns=["capture", *CONTIG_COLUMNS])


def h(prefix: str, parts: list[str]) -> str:
    return prefix + hashlib.sha256(";".join(sorted(parts)).encode()).hexdigest()[:16]


def test_flag_strings_are_parsed_strictly() -> None:  # Review Focus 1
    s = pd.Series([True, "True", "true", " TRUE ", False, "False", "None", "", np.nan, "1"], dtype=object)
    assert parse_flag(s).tolist() == [True] * 4 + [False] * 6


def test_a_cell_with_three_trbs_is_unlabelled() -> None:
    out = call_clones(
        frame(
            contig("A", "b1", "TRA", "AAA"), *(contig("A", "b1", "TRB", nt) for nt in ("CCC", "GGG", "TTT"))
        )
    )
    row = out.loc[("A", "b1")]
    assert pd.isna(row.clone_id) and pd.isna(row.clone_id_trb) and row.receptor_n_trb == 3


def test_two_of_each_chain_are_allowed_in_any_order() -> None:
    a = [("TRA", "AAA"), ("TRA", "ACA"), ("TRB", "CCC"), ("TRB", "CGC")]
    out = call_clones(
        frame(
            *(contig("A", "b1", c, n) for c, n in a),
            *(contig("A", "b2", c, n) for c, n in reversed(a)),
            contig("A", "b3", "TRA", "TTT"),
            contig("A", "b3", "TRB", "CCC"),
            contig("A", "b3", "TRB", "CGC"),
        )
    )
    assert out.loc[("A", "b1")].clone_id == out.loc[("A", "b2")].clone_id
    assert out.loc[("A", "b1")].clone_id == h("c", [f"{c}:{n}" for c, n in a])
    assert out.loc[("A", "b3")].clone_id != out.loc[("A", "b1")].clone_id
    assert out.loc[("A", "b3")].clone_id_trb == out.loc[("A", "b1")].clone_id_trb
    assert out.loc[("A", "b1")].clone_id_trb == h("b", ["TRB:CCC", "TRB:CGC"])
    assert out.loc[("A", "b1")].receptor_n_tra == 2


@pytest.mark.parametrize(
    "override",
    [
        {"productive": "False"},
        {"high_confidence": "False"},
        {"full_length": "False"},
        {"cdr3_nt": "None"},
        {"chain": "IGH"},
        {"cdr3_nt": "AAXA"},
        {"is_cell": "False"},
    ],
)
def test_non_qualifying_contigs_are_ignored(override: dict[str, object]) -> None:
    tra = {**contig("A", "b1", "TRA", "AAA"), **override}
    trb = contig("A", "b1", "TRB", "CCC")
    df = frame(tra, trb)
    if "is_cell" in override:
        df["is_cell"] = ["False", "True"]
    row = call_clones(df).loc[("A", "b1")]
    assert pd.isna(row.clone_id) and row.clone_id_trb == h("b", ["TRB:CCC"])
    assert row.receptor_n_tra == 0 and row.receptor_n_trb == 1


def test_duplicate_contig_rows_count_once() -> None:
    row = call_clones(frame(contig("A", "b1", "TRB", "CCC"), contig("A", "b1", "TRB", "CCC"))).loc[
        ("A", "b1")
    ]
    assert row.receptor_n_trb == 1 and pd.notna(row.clone_id_trb)


def test_strings_carry_no_namespace() -> None:
    rows = [contig(c, "b1", ch, nt) for c in ("A", "B") for ch, nt in (("TRA", "AAA"), ("TRB", "CCC"))]
    out = call_clones(frame(*rows))
    assert out.loc[("A", "b1")].clone_id == out.loc[("B", "b1")].clone_id
    assert out.loc[("A", "b1")].raw_clonotype_id == "A:clonotype1"
    assert out.loc[("B", "b1")].raw_clonotype_id == "B:clonotype1"


def test_missing_raw_clonotype_is_nan_and_unlabelled_barcodes_are_kept() -> None:
    out = call_clones(
        frame(contig("A", "b1", "TRB", "CCC", raw_clonotype_id="None"), contig("A", "b2", "IGH", "CCC"))
    )
    assert pd.isna(out.loc[("A", "b1")].raw_clonotype_id)
    assert list(out.index) == [("A", "b1"), ("A", "b2")]
    assert out.loc[("A", "b2")].receptor_n_trb == 0 and pd.isna(out.loc[("A", "b2")].clone_id_trb)


def test_missing_columns_raise() -> None:
    with pytest.raises(BuildError, match="cdr3_nt"):
        call_clones(frame(contig("A", "b1", "TRB", "CCC")).drop(columns="cdr3_nt"))


def test_tar_members_are_read_in_memory(tmp_path: Path) -> None:
    payload = gzip.compress(b"a,b\n1,2\n")
    tar_path = tmp_path / "t.tar"
    with tarfile.open(tar_path, "w") as tf:
        info = tarfile.TarInfo("x.csv.gz")
        info.size = len(payload)
        tf.addfile(info, io.BytesIO(payload))
        stray = tarfile.TarInfo("../stray.csv.gz")
        stray.size = 1
        tf.addfile(stray, io.BytesIO(b"z"))
        d = tarfile.TarInfo("d")
        d.type = tarfile.DIRTYPE
        tf.addfile(d)
    before = sorted(p.name for p in tmp_path.iterdir())
    got, ignored = read_members(tar_path, ["x.csv.gz"])
    assert got == {"x.csv.gz": payload}
    assert ignored == ["../stray.csv.gz", "d"]
    assert sorted(p.name for p in tmp_path.iterdir()) == before
    with pytest.raises(BuildError):
        read_members(tar_path, ["missing"])
    with pytest.raises(BuildError):
        read_members(tar_path, ["d"])


def test_zero_contig_rows_raise_a_build_error() -> None:
    with pytest.raises(BuildError, match="no contigs"):
        call_clones(frame())


@pytest.mark.parametrize("blank", ["", "  ", "\t", "None", np.nan, None])
def test_blank_raw_clonotype_ids_are_missing(blank: object) -> None:
    out = call_clones(frame(contig("A", "b1", "TRB", "CCC", raw_clonotype_id=blank)))
    assert pd.isna(out.loc[("A", "b1")].raw_clonotype_id)


def test_raw_clonotype_comes_from_the_first_non_missing_contig() -> None:
    def raw(*ids: object) -> object:
        rows = [contig("A", "b1", "TRB", "CCC", raw_clonotype_id=i) for i in ids]
        return call_clones(frame(*rows)).loc[("A", "b1")].raw_clonotype_id

    assert raw("c1", "c2") == "A:c1"
    assert raw("None", "c2") == "A:c2"
    assert raw("", np.nan, "c2") == "A:c2"


def test_three_distinct_tras_with_one_trb_leave_only_the_trb_clone() -> None:
    rows = [contig("A", "b1", "TRA", nt) for nt in ("AAA", "ACA", "AGA")] + [contig("A", "b1", "TRB", "CCC")]
    row = call_clones(frame(*rows)).loc[("A", "b1")]
    assert pd.isna(row.clone_id) and row.clone_id_trb == h("b", ["TRB:CCC"]) and row.receptor_n_tra == 3


def test_the_same_tra_on_two_contigs_counts_once() -> None:
    rows = [contig("A", "b1", "TRA", "AAA"), contig("A", "b1", "TRA", "AAA"), contig("A", "b1", "TRB", "CCC")]
    row = call_clones(frame(*rows)).loc[("A", "b1")]
    assert row.receptor_n_tra == 1 and row.clone_id == h("c", ["TRA:AAA", "TRB:CCC"])


@pytest.mark.parametrize("bad", ["yes", "Yes", "T", "TRUE1", "tRue", "1.0"])
def test_flag_parser_rejects_near_misses(bad: str) -> None:
    assert parse_flag(pd.Series([bad], dtype=object)).tolist() == [False]


def test_lowercase_cdr3_is_ignored_and_ana_qualifies() -> None:
    rows = [contig("A", "b1", "TRB", "ccc"), contig("A", "b2", "TRB", "ANA")]
    out = call_clones(frame(*rows))
    assert out.loc[("A", "b1")].receptor_n_trb == 0
    assert out.loc[("A", "b2")].clone_id_trb == h("b", ["TRB:ANA"])


def test_output_index_is_sorted_for_any_input_order() -> None:
    rows = [contig(c, b, "TRB", "CCC") for c, b in (("B", "z"), ("A", "y"), ("B", "a"), ("A", "b"))]
    out = call_clones(frame(*rows))
    assert list(out.index) == [("A", "b"), ("A", "y"), ("B", "a"), ("B", "z")]


def test_is_cell_nan_none_and_numpy_bools_go_through_call_clones() -> None:
    rows = [contig("A", f"b{i}", "TRB", "CCC") for i in range(4)]
    df = frame(*rows)
    df["is_cell"] = pd.Series([np.True_, np.False_, np.nan, "None"], dtype=object)
    out = call_clones(df)
    assert out["receptor_n_trb"].tolist() == [1, 0, 0, 0]


def test_error_messages_never_carry_member_names(tmp_path: Path) -> None:
    tar_path = tmp_path / "t.tar"
    with tarfile.open(tar_path, "w") as tf:
        info = tarfile.TarInfo("SECRET_stray.csv")
        info.size = 1
        tf.addfile(info, io.BytesIO(b"z"))
    with pytest.raises(BuildError) as e:
        read_members(tar_path, ["SECRET_wanted.csv"])
    assert "SECRET" not in str(e.value)


def test_a_non_tar_file_raises_a_build_error_naming_the_path(tmp_path: Path) -> None:
    p = tmp_path / "notatar.tar"
    p.write_bytes(b"this is not a tar" * 100)
    with pytest.raises(BuildError, match=r"notatar\.tar"):
        read_members(p, ["x"])


def test_repeated_wanted_names_follow_last_one_wins(tmp_path: Path) -> None:
    def make(name: str, *members: tuple[bytes | None, bool]) -> Path:
        path = tmp_path / name
        with tarfile.open(path, "w") as tf:
            for data, is_file in members:
                info = tarfile.TarInfo("w")
                if is_file:
                    assert data is not None
                    info.size = len(data)
                    tf.addfile(info, io.BytesIO(data))
                else:
                    info.type = tarfile.DIRTYPE
                    tf.addfile(info)
        return path

    got, other = read_members(make("a.tar", (b"1", True), (b"2", True)), ["w"])
    assert got == {"w": b"2"} and other == []
    got, other = read_members(make("b.tar", (None, False), (b"2", True)), ["w"])
    assert got == {"w": b"2"} and other == []
    with pytest.raises(BuildError):
        read_members(make("c.tar", (b"1", True), (None, False)), ["w"])
