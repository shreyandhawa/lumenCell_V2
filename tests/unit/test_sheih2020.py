import gzip
import io
import tarfile
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from lumencell.bench import download as dl
from lumencell.bench.builders import sheih2020
from lumencell.bench.builders.common import BuildError
from lumencell.bench.datasets import CASE_STUDIES, CaseStudyDataset
from lumencell.bench.download import Build, Source, build_dataset, local_path
from lumencell.bench.gates import load_thresholds
from lumencell.clonal.isolation import receptor_gene_mask, scoring_only_columns
from lumencell.clonal.labels import encode_clones
from lumencell.data.loader import load

GENES = ["G1", "G2", "TRBC2", "TRAV1", "G5"]
# Cell ids: <16 nt>-<aggr id>. c1 and c4 share a barcode core, one in each sample.
CORE = "AAAACCCCGGGGTTTT"
CELLS = {
    "c1": f"{CORE}-1",
    "c2": "AAAACCCCGGGGTTTA-1",
    "c3": "AAAACCCCGGGGTTAT-1",
    "c4": f"{CORE}-2",
    "c5": "AAAACCCCGGGGTTAA-2",
    "c6": "AAAACCCCGGGGTATT-2",
}
MISS = "CCCCCCCCCCCCCCCC-1"
COUNTS = np.array(
    [[0, 1, 0, 2, 0, 3], [4, 0, 5, 0, 6, 0], [1, 1, 0, 0, 2, 2], [0, 0, 7, 7, 0, 1], [9, 8, 0, 0, 0, 0]]
)  # genes x cells
SAMPLES = (
    sheih2020.Sample("GSM_A", "GSM_VA_P1_IP.csv.gz", "P1", "IP", 0, 1),
    sheih2020.Sample("GSM_B", "GSM_VB_P2_d12.csv.gz", "P2", "d12", 12, 2),
)
HEADER = ["chain", "clonotype_ids", "consensus_ids", "barcodes", "contig_ids", "cdr3s", "cdr3_nts"]


def _row(chain: str, nt: str, barcodes: list[str]) -> dict[str, str]:
    return {
        "chain": chain,
        "clonotype_ids": "clonotype1",
        "consensus_ids": "c1",
        "barcodes": ";".join(barcodes),
        "contig_ids": "x",
        "cdr3s": "CASS",
        "cdr3_nts": nt,
    }


def _bc(*cells: str) -> list[str]:
    return [CELLS[c].rsplit("-", 1)[0] + "-1" for c in cells]


TRA_X, TRB_Y, TRA_X2, TRB_Y2, TRB_Z = "TGTGCAGCA", "TGCGCCAGC", "TGTGCTTCC", "TGCAGTGCT", "TGCCCCCCC"
TABLE_A = [
    _row("TRA", TRA_X, [*_bc("c1"), MISS]),
    _row("TRB", TRB_Y, [*_bc("c1"), MISS]),
    _row("TRA", TRA_X2, _bc("c2")),
    _row("TRB", TRB_Y2, _bc("c2")),
    _row("TRB", TRB_Z, _bc("c3")),
]
TABLE_B = [
    _row("TRA", TRA_X, _bc("c4")),  # the same receptor as c1, in another patient
    _row("TRB", TRB_Y, _bc("c4")),
    _row("TRA", TRA_X2, _bc("c5")),
    _row("TRB", "TGAAAAAAA", _bc("c5")),
    _row("TRB", "TGAAAAAAC", _bc("c5")),
    _row("TRB", "TGAAAAAAG", _bc("c5")),
]


def _gz(frame: pd.DataFrame) -> bytes:
    return gzip.compress(frame.to_csv(index=False).encode())


def _write_matrix(path: Path, cells: dict[str, str] = CELLS) -> None:
    head = '""' + "".join(f',"{c}"' for c in cells.values())
    lines = [head] + [f'"{g}",' + ",".join(str(v) for v in row) for g, row in zip(GENES, COUNTS, strict=True)]
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))


def _write_tar(path: Path, tables: dict[str, list[dict[str, str]]], stray: bool = True) -> None:
    with tarfile.open(path, "w") as tf:
        members = {name: _gz(pd.DataFrame(rows, columns=HEADER)) for name, rows in tables.items()}
        if stray:
            members["../stray.csv.gz"] = b"not a table"
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))


@pytest.fixture(autouse=True)
def _two_samples(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sheih2020, "SAMPLES", SAMPLES)


@pytest.fixture
def sources(tmp_path: Path) -> dict[str, Path]:
    d = tmp_path / "src"
    d.mkdir()
    _write_matrix(d / "m.csv.gz")
    _write_tar(d / "t.tar", {SAMPLES[0].member: TABLE_A, SAMPLES[1].member: TABLE_B})
    return {sheih2020.MATRIX: d / "m.csv.gz", sheih2020.TAR: d / "t.tar"}


PARAMS: dict[str, Any] = {"chains": "paired", "cdr3": "nt", "filters": [], "max_per_chain": 2}


def _build(sources: dict[str, Path], out: Path, params: dict[str, Any] = PARAMS, **kw: Any) -> ad.AnnData:
    sheih2020.build(sources, out, params, **kw)
    return ad.read_h5ad(out)


@pytest.fixture
def built(sources: dict[str, Path], tmp_path: Path) -> ad.AnnData:
    return _build(sources, tmp_path / "out.h5ad", chunk_rows=2)


def test_the_matrix_is_read_in_chunks_without_changing_values(built: ad.AnnData) -> None:
    assert built.shape == (6, 5) and list(built.var_names) == GENES
    np.testing.assert_array_equal(built.X.toarray(), COUNTS.T)
    assert built.X.dtype == np.float32 and built.X.format == "csr"
    assert list(built.obs_names) == list(CELLS.values())


def test_chunk_size_does_not_change_the_result(sources: dict[str, Path], tmp_path: Path) -> None:
    a = _build(sources, tmp_path / "a.h5ad", chunk_rows=1)
    b = _build(sources, tmp_path / "b.h5ad", chunk_rows=500)
    np.testing.assert_array_equal(a.X.toarray(), b.X.toarray())
    pd.testing.assert_frame_equal(a.obs, b.obs)


def test_cells_map_to_patient_and_time_point_by_suffix(built: ad.AnnData) -> None:
    assert built.obs["patient"].astype(str).tolist() == ["P1"] * 3 + ["P2"] * 3
    assert built.obs["timepoint"].astype(str).tolist() == ["IP"] * 3 + ["d12"] * 3
    assert built.obs["day"].tolist() == [0, 0, 0, 12, 12, 12]
    assert built.obs["sample"].astype(str).tolist() == ["GSM_A"] * 3 + ["GSM_B"] * 3


def test_the_same_receptor_in_two_patients_is_two_clones(built: ad.AnnData) -> None:
    o = built.obs
    assert o.loc[CELLS["c1"], "clone_id"] == o.loc[CELLS["c4"], "clone_id"]
    assert o.loc[CELLS["c1"], "patient"] != o.loc[CELLS["c4"], "patient"]
    codes = encode_clones(
        o["clone_id"].to_numpy(dtype=object),
        o["patient"].astype(str).to_numpy(dtype=object),
    )
    assert codes.clone[0] >= 0 and codes.clone[3] >= 0 and codes.clone[0] != codes.clone[3]
    assert o.loc[CELLS["c1"], "clone_id"] != o.loc[CELLS["c2"], "clone_id"]


def test_the_declared_chains_are_used(built: ad.AnnData) -> None:
    o = built.obs
    # c3 has a TRB only: no paired clone, but a TRB-only one. c5 has three TRBs: unlabelled.
    assert pd.isna(o.loc[CELLS["c3"], "clone_id"]) and pd.notna(o.loc[CELLS["c3"], "clone_id_trb"])
    assert (o.loc[CELLS["c3"], "receptor_n_tra"], o.loc[CELLS["c3"], "receptor_n_trb"]) == (0, 1)
    assert pd.isna(o.loc[CELLS["c5"], "clone_id"]) and pd.isna(o.loc[CELLS["c5"], "clone_id_trb"])
    assert o.loc[CELLS["c5"], "receptor_n_trb"] == 3


def test_a_cell_without_a_table_row_is_unlabelled_with_zero_counts(built: ad.AnnData) -> None:
    r = built.obs.loc[CELLS["c6"]]
    assert pd.isna(r["clone_id"]) and (r["receptor_n_tra"], r["receptor_n_trb"]) == (0, 0)


def test_the_clone_id_format_is_the_common_call(built: ad.AnnData) -> None:
    from lumencell.bench.builders.common import call_clones

    contigs = pd.DataFrame(
        {
            "capture": "GSM_A",
            "barcode": [CELLS["c1"]] * 2,
            "chain": ["TRA", "TRB"],
            "productive": "True",
            "high_confidence": "True",
            "full_length": "True",
            "cdr3_nt": [TRA_X, TRB_Y],
            "raw_clonotype_id": np.nan,
        }
    )
    expected = call_clones(contigs).loc[("GSM_A", CELLS["c1"]), "clone_id"]
    assert built.obs.loc[CELLS["c1"], "clone_id"] == expected


def test_a_join_miss_is_counted(built: ad.AnnData) -> None:
    b = built.uns["lumencell_build"]
    assert b["join_misses"] == {"GSM_A": 1, "GSM_B": 0}
    assert b["contig_barcodes"] == {"GSM_A": 4, "GSM_B": 2}
    assert b["n_cells"] == 6 and b["labelled_cells"] == 3 and b["labelled_cells_trb"] == 4
    assert b["version"] == sheih2020.BUILD_VERSION
    p = b["params"]  # an empty list reads back as an empty array
    assert (p["chains"], p["cdr3"], p["max_per_chain"], len(p["filters"])) == ("paired", "nt", 2, 0)
    assert b["joined_share"]["GSM_A"] == pytest.approx(0.75) and b["joined_share"]["GSM_B"] == 1.0
    assert b["ignored_members"] == ["../stray.csv.gz"]
    assert not any(MISS.split("-")[0] in n for n in built.obs_names)


def test_a_multi_value_cdr3_field_is_not_used_and_is_counted(
    sources: dict[str, Path], tmp_path: Path
) -> None:
    rows = [*TABLE_B, _row("TRB", "TGAAAAAAA;TGAAAAAAC", _bc("c6")), _row("TRA", "", _bc("c6"))]
    _write_tar(sources[sheih2020.TAR], {SAMPLES[0].member: TABLE_A, SAMPLES[1].member: rows})
    b = _build(sources, tmp_path / "o.h5ad")
    assert b.uns["lumencell_build"]["cdr3_unusable_rows"] == {"GSM_A": 0, "GSM_B": 2}
    assert b.obs.loc[CELLS["c6"], "receptor_n_trb"] == 0 and pd.isna(b.obs.loc[CELLS["c6"], "clone_id"])


def test_a_sample_below_the_joined_share_floor_refuses_to_build(
    sources: dict[str, Path], tmp_path: Path
) -> None:
    rows = [_row("TRA", TRA_X, [MISS, "GGGGGGGGGGGGGGGG-1", *_bc("c4")]), _row("TRB", TRB_Y, [MISS])]
    _write_tar(sources[sheih2020.TAR], {SAMPLES[0].member: TABLE_A, SAMPLES[1].member: rows})
    out = tmp_path / "o.h5ad"
    with pytest.raises(BuildError, match="joined share"):
        sheih2020.build(sources, out, PARAMS)
    assert not out.exists()


def test_a_cell_with_an_unknown_sample_suffix_is_refused(sources: dict[str, Path], tmp_path: Path) -> None:
    _write_matrix(sources[sheih2020.MATRIX], {**CELLS, "c6": "AAAACCCCGGGGTATT-9"})
    with pytest.raises(BuildError):
        sheih2020.build(sources, tmp_path / "o.h5ad", PARAMS)


@pytest.mark.parametrize(
    "bad",
    [
        {**PARAMS, "chains": "trb"},
        {**PARAMS, "cdr3": "aa"},
        {**PARAMS, "filters": ["productive"]},
        {**PARAMS, "max_per_chain": 3},
        {"chains": "paired"},
    ],
)
def test_params_the_builder_cannot_honour_are_refused(
    sources: dict[str, Path], tmp_path: Path, bad: dict[str, Any]
) -> None:
    with pytest.raises(BuildError):
        sheih2020.build(sources, tmp_path / "o.h5ad", bad)


def test_hostile_values_never_reach_an_error_message(sources: dict[str, Path], tmp_path: Path) -> None:
    _write_matrix(sources[sheih2020.MATRIX], {**CELLS, "c6": "../../etc/passwd;rm -rf-2"})
    with pytest.raises(BuildError) as e:
        sheih2020.build(sources, tmp_path / "o.h5ad", PARAMS)
    assert "passwd" not in str(e.value) and "rm -rf" not in str(e.value)


def test_clone_columns_are_scoring_only(built: ad.AnnData) -> None:
    so = scoring_only_columns(built.obs.columns, ("patient",))
    assert {"clone_id", "clone_id_trb", "receptor_n_tra", "receptor_n_trb", "patient"} <= so
    assert not {"timepoint", "day", "sample"} & so
    assert "raw_clonotype_id" not in built.obs.columns


def test_the_output_loads_as_counts_with_receptor_symbols(built: ad.AnnData, tmp_path: Path) -> None:
    built.write_h5ad(tmp_path / "again.h5ad")
    ds = load(tmp_path / "again.h5ad")
    assert ds.x_kind == "counts" and receptor_gene_mask(ds.gene_symbols).sum() == 2


def test_nothing_is_written_beside_the_output(sources: dict[str, Path], tmp_path: Path) -> None:
    out_dir = tmp_path / "outdir"
    out_dir.mkdir()
    _build(sources, out_dir / "o.h5ad")
    assert [p.name for p in out_dir.iterdir()] == ["o.h5ad"]


def test_params_mirror_the_declaration() -> None:
    cs = load_thresholds()["case_study"]["sheih2020"]
    assert CASE_STUDIES["sheih2020"].build.params == {
        k: cs[k] for k in ("chains", "cdr3", "filters", "max_per_chain")
    }


def test_case_study_sources_are_complete() -> None:
    cs = CASE_STUDIES["sheih2020"]
    assert isinstance(cs, CaseStudyDataset) and cs.licence.strip() and cs.build.version >= 1
    assert {s.filename for s in cs.sources} == {sheih2020.MATRIX, sheih2020.TAR}
    for s in cs.sources:
        assert s.sha256 and len(s.sha256) == 64 and s.url.startswith("https://") and s.size > 0


def test_the_samples_follow_the_note() -> None:
    real = sheih2020.REAL_SAMPLES  # the real table, not the two-sample patch
    assert len(real) == 16 and sorted(s.aggr_id for s in real) == list(range(1, 17))
    assert {s.patient for s in real} == {"CLL-1", "CLL-2", "NHL-6", "NHL-7"}
    assert {s.timepoint for s in real} == {
        "IP", "d12", "d21", "d28", "d29", "d38", "d83", "d89", "d102", "d112",
    }  # fmt: skip
    assert all((s.timepoint == "IP") == (s.day == 0) for s in real)
    assert len({s.gsm for s in real}) == 16 and len({s.member for s in real}) == 16


def test_build_dataset_writes_the_derived_file_by_atomic_rename(
    sources: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache = tmp_path / "cache"
    cache.mkdir()
    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(cache))
    srcs = []
    for name, p in sources.items():
        (cache / name).write_bytes(p.read_bytes())
        srcs.append(Source("https://example.org/" + name, name, p.stat().st_size, dl.sha256_of(p)))
    item = CaseStudyDataset(
        id="tiny", title="t", citation="c", sources=tuple(srcs), licence="l",
        build=Build(sheih2020.build, 1, PARAMS),
    )  # fmt: skip
    out = build_dataset(item)
    assert out == local_path(item) and out.exists()
    assert sorted(p.name for p in cache.iterdir() if p.suffix == ".h5ad") == [out.name]
    assert not [p for p in cache.iterdir() if p.name.endswith(".tmp.h5ad")]


def test_bench_fetch_finds_a_case_study_by_its_id_but_all_skips_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from typer.testing import CliRunner

    from lumencell.cli import app

    monkeypatch.setenv("LUMENCELL_DATASETS_DIR", str(tmp_path))
    out = CliRunner().invoke(app, ["bench", "fetch", "sheih2020"], input="n\n")
    assert out.exit_code == 0 and "skipped" in out.output, out.output
    assert all(s.url in out.output for s in CASE_STUDIES["sheih2020"].sources)
    assert not list(tmp_path.iterdir())
    everything = CliRunner().invoke(app, ["bench", "fetch", "all"], input="n\n" * 4)
    assert "sheih2020" not in everything.output
    assert CliRunner().invoke(app, ["bench", "fetch", "nope"]).exit_code == 2
