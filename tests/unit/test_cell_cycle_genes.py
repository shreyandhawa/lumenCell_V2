"""The shipped Tirosh 2016 lists, the mouse ortholog table and the symbol matching."""

from __future__ import annotations

import re
from importlib.resources import files

import pytest

from lumencell.prep.cell_cycle_genes import (
    GENE_SETS,
    load_gene_set,
    match_genes,
    resolve_species,
)

DATA = files("lumencell.prep") / "cell_cycle_data"


def _text(name: str) -> str:
    return (DATA / name).read_text(encoding="utf-8")


def _rows(name: str) -> list[list[str]]:
    lines = [ln for ln in _text(name).splitlines() if ln and not ln.startswith("#")]
    return [ln.split("\t") for ln in lines[1:]]


def test_tirosh_lists_have_43_s_and_54_g2m_genes():
    h = load_gene_set("tirosh2016", "human")
    assert (len(h.s), len(h.g2m)) == (43, 54)
    flat = [sym for gene in h.s + h.g2m for sym in gene]
    assert len(flat) == len(set(flat))
    assert h.unmapped == ()


def test_no_symbol_is_repeated_in_the_shipped_tsv():
    rows = _rows("tirosh2016.tsv")
    assert len(rows) == 97
    published = [r[1] for r in rows]
    assert len(published) == len(set(published))
    current = [r[2] for r in rows if r[2]]
    assert len(current) == len(set(current))
    assert not set(current) & set(published)
    orth = _rows("human_mouse_orthologs.tsv")
    assert len({r[0] for r in orth}) == len(orth)
    assert len({r[1] for r in orth}) == len(orth)


def test_known_members():
    h = load_gene_set("tirosh2016", "human")
    assert ("MCM5",) in h.s and ("PCNA",) in h.s and ("MLF1IP", "CENPU") in h.s
    assert ("TOP2A",) in h.g2m and ("HN1", "JPT1") in h.g2m
    assert ("FAM64A", "PIMREG") in h.g2m
    m = load_gene_set("tirosh2016", "mouse")
    assert ("Mcm5",) in m.s and ("Top2a",) in m.g2m


def test_every_human_gene_is_mapped_or_listed_as_unmapped():
    m = load_gene_set("tirosh2016", "mouse")
    assert len(m.s) + len(m.g2m) + len(m.unmapped) == 97


def test_headers_name_their_sources():
    for name, needle in (
        ("tirosh2016.tsv", "10.1126/science.aad0501"),
        ("human_mouse_orthologs.tsv", "HOM_MouseHumanSequence.rpt"),
    ):
        head = "\n".join(ln for ln in _text(name).splitlines() if ln.startswith("#"))
        assert needle in head
        assert re.search(r"\d{4}-\d{2}-\d{2}", head)
        assert re.search(r"\b[0-9a-f]{64}\b", head)
        assert any(ln.startswith("# verified:") for ln in head.splitlines())


def test_resolve_species_from_symbol_case():
    h = load_gene_set("tirosh2016", "human")
    upper = [g[0] for g in h.s + h.g2m]
    assert resolve_species(upper, "auto") == "human"
    assert resolve_species([s.capitalize() for s in upper], "auto") == "mouse"
    assert resolve_species(["G00000", "G00001"], "auto") is None
    assert resolve_species(["PCNA", "Mcm5"], "auto") is None
    assert resolve_species(["G00000"], "mouse") == "mouse"


def test_matching_is_case_insensitive_and_counts_genes_once():
    h = load_gene_set("tirosh2016", "human")
    assert match_genes(["pcna", "PCNA", "CENPU", "x"], h.s) == [0, 2]


def test_unknown_gene_set_raises():
    assert GENE_SETS == ("tirosh2016",)
    with pytest.raises(ValueError):
        load_gene_set("nope", "human")


def test_renamed_genes_map_to_their_mouse_orthologs():
    m = load_gene_set("tirosh2016", "mouse")
    assert ("Cenpu",) in m.s
    assert ("Pimreg",) in m.g2m and ("Jpt1",) in m.g2m
    assert m.unmapped == ()
    assert len(m.s) == 43 and len(m.g2m) == 54


def test_an_unknown_species_raises():
    with pytest.raises(ValueError, match="species"):
        load_gene_set("tirosh2016", "rat")  # type: ignore[arg-type]


def test_the_first_var_index_of_any_accepted_symbol_wins():
    h = load_gene_set("tirosh2016", "human")
    assert match_genes(["CENPU", "MLF1IP"], h.s) == [0]
    assert match_genes(["MLF1IP", "x", "CENPU"], h.s) == [0]


def test_a_tsv_without_its_header_or_with_a_short_row_is_refused(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    import lumencell.prep.cell_cycle_genes as ccg

    class _Dir:
        def __init__(self, text: str) -> None:
            self.text = text

        def __truediv__(self, _: str) -> _Dir:
            return self

        def read_text(self, encoding: str) -> str:
            return self.text

    monkeypatch.setattr(ccg, "files", lambda _pkg: _Dir("# c\nS\tMCM5\t\n"))
    with pytest.raises(ValueError, match="header"):
        ccg._rows("tirosh2016.tsv")
    monkeypatch.setattr(ccg, "files", lambda _pkg: _Dir("phase\tsymbol\tcurrent_symbol\nS\tMCM5\n"))
    with pytest.raises(ValueError, match="line 2"):
        ccg._rows("tirosh2016.tsv")
