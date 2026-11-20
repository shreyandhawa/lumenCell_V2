"""The case-study declaration (phase 3 Task 10), written before any analysis.

The facts it is checked against come from Task 9's note
(docs/superpowers/research/2026-11-19-sheih-verification.md), never from clone values.
"""

from datetime import date
from typing import Any

import pytest

from lumencell.api.models import ClonalSettings
from lumencell.bench.builders.common import CONTIG_COLUMNS
from lumencell.bench.gates import THRESHOLDS, load_thresholds
from lumencell.engine.roots import ClusterRoot, parse_root
from lumencell.trust.clonal import settings_columns

# Task 9's note: the 13 columns of every *-unpaired-clonotypes.csv.gz table.
SHEIH_TABLE_COLUMNS = (
    "chain",
    "clonotype_ids",
    "consensus_ids",
    "barcodes",
    "contig_ids",
    "cdr3s",
    "cdr3_nts",
    "v_genes",
    "d_genes",
    "j_genes",
    "c_genes",
    "frequency",
    "proportion",
)
# Task 9's note, GSM table: GEO's "time point" values as written, per patient.
SHEIH_TIMEPOINTS = {
    "CLL-1": ("IP", "d21", "d38", "d112"),
    "CLL-2": ("IP", "d12", "d29", "d83"),
    "NHL-6": ("IP", "d12", "d29", "d102"),
    "NHL-7": ("IP", "d12", "d28", "d89"),
}
CONTIG_FLAGS = ("productive", "high_confidence", "full_length")


@pytest.fixture
def cs() -> dict[str, Any]:
    section: dict[str, Any] = load_thresholds()["case_study"]["sheih2020"]
    return section


def test_case_study_is_declared() -> None:
    cs = load_thresholds()["case_study"]["sheih2020"]
    s = ClonalSettings(
        clone_key=cs["clone_key"],
        donor_key=cs["donor_key"],
        single_donor=cs["single_donor"],
        strata_keys=cs["strata_keys"],
    )
    assert (
        s.donor_key == "patient" and s.strata_keys == ["timepoint"] and "timepoint" not in settings_columns(s)
    )
    root = parse_root(cs["root"])
    assert isinstance(root, ClusterRoot) and root.key == "timepoint"
    assert cs["methods"] == ["paga_dpt", "palantir", "scfates_ppt"] and cs["max_per_chain"] == 2
    assert (
        cs["chains"] in ("paired", "trb")
        and cs["cdr3"] in ("nt", "aa")
        and date.fromisoformat(cs["declared"])
    )
    assert cs["unit_key"] == "patient" == cs["donor_key"] and cs["pooled_diagnostic"] is True


def test_the_section_holds_exactly_the_declared_keys(cs: dict[str, Any]) -> None:
    assert set(cs) == {
        "declared",
        "accession",
        "paper",
        "methods",
        "root",
        "clone_key",
        "donor_key",
        "single_donor",
        "strata_keys",
        "chains",
        "cdr3",
        "filters",
        "max_per_chain",
        "unit_key",
        "pooled_diagnostic",
    }


def test_the_source_is_named(cs: dict[str, Any]) -> None:
    assert cs["declared"] == "2026-11-19"
    assert cs["accession"] == "GSE125881"
    assert cs["paper"] == "Sheih et al. 2020, Nat Commun 11:219, doi:10.1038/s41467-019-13880-1"


def test_the_clone_and_namespace_columns(cs: dict[str, Any]) -> None:
    assert cs["clone_key"] == "clone_id" and cs["donor_key"] == "patient" and cs["single_donor"] is False
    assert cs["max_per_chain"] == 2


def test_the_root_is_the_infusion_product_as_geo_writes_it(cs: dict[str, Any]) -> None:
    root = parse_root(cs["root"])
    assert isinstance(root, ClusterRoot)
    assert root == ClusterRoot("timepoint", "IP")
    # The root is a stratum, so methods may see it, and every patient has that time point.
    assert root.key in cs["strata_keys"]
    assert all(root.value in tps for tps in SHEIH_TIMEPOINTS.values())


def test_filters_are_the_contig_flags_the_tables_carry(cs: dict[str, Any]) -> None:
    assert set(CONTIG_FLAGS) <= set(CONTIG_COLUMNS)  # the flags call_clones knows
    carried = [f for f in CONTIG_FLAGS if f in SHEIH_TABLE_COLUMNS]
    assert cs["filters"] == carried == []


def test_chains_follow_the_documented_table_format(cs: dict[str, Any]) -> None:
    # Every row names its chain (TRA or TRB) and lists its cells' barcodes, so both chains
    # can be grouped by cell barcode: the plan's rule gives paired.
    groupable = {"chain", "barcodes"} <= set(SHEIH_TABLE_COLUMNS)
    assert cs["chains"] == ("paired" if groupable else "trb") == "paired"


def test_cdr3_is_nucleotide_when_the_tables_carry_it(cs: dict[str, Any]) -> None:
    assert cs["cdr3"] == ("nt" if "cdr3_nts" in SHEIH_TABLE_COLUMNS else "aa") == "nt"


def test_the_declaration_has_its_dated_changelog_entry(cs: dict[str, Any]) -> None:
    """The changelog entry carries the section's date, the assumed root and the per-patient cost."""
    text = THRESHOLDS.read_text()
    head = f"# {cs['declared']}  Case study declared: [case_study.sheih2020]"
    assert text.count(head) == 1
    entry = text[text.index(head) :]
    assert "Assumption: IP as the start is assumed" in entry
    assert "12 runs with stability instead of 3" in entry
