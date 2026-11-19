"""The Tirosh 2016 cell-cycle gene lists and their mouse orthologs, read from shipped TSVs."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from typing import Literal

GENE_SETS = ("tirosh2016",)
Species = Literal["human", "mouse"]
_SOURCE = "Tirosh et al. 2016, Science 352:189, doi:10.1126/science.aad0501, Table S5"
_GeneList = tuple[tuple[str, ...], ...]


@dataclass(frozen=True)
class CellCycleGenes:
    """Each gene is a tuple of the symbols it may appear under."""

    gene_set: str
    species: Species
    s: _GeneList
    g2m: _GeneList
    source: str
    unmapped: tuple[str, ...]


def _rows(name: str) -> list[list[str]]:
    text = (files("lumencell.prep") / "cell_cycle_data" / name).read_text(encoding="utf-8")
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    return [ln.split("\t") for ln in lines[1:]]


@cache
def _load(gene_set: str, species: Species) -> CellCycleGenes:
    s: list[tuple[str, ...]] = []
    g2m: list[tuple[str, ...]] = []
    unmapped: list[str] = []
    mouse = {h: m for h, m, _ in _rows("human_mouse_orthologs.tsv")} if species == "mouse" else {}
    for phase, symbol, current in _rows("tirosh2016.tsv"):
        if species == "human":
            gene = (symbol, current) if current else (symbol,)
        else:
            human = current or symbol
            if human not in mouse:
                unmapped.append(human)
                continue
            gene = (mouse[human],)
        (s if phase == "S" else g2m).append(gene)
    return CellCycleGenes(gene_set, species, tuple(s), tuple(g2m), _SOURCE, tuple(unmapped))


def load_gene_set(gene_set: str, species: Species) -> CellCycleGenes:
    if gene_set not in GENE_SETS:
        raise ValueError(f"unknown cell-cycle gene set {gene_set!r}; expected one of {GENE_SETS}")
    if species not in ("human", "mouse"):
        raise ValueError(f"unknown species {species!r}")
    return _load(gene_set, species)


def _lowered(genes: CellCycleGenes) -> set[str]:
    return {sym.lower() for gene in genes.s + genes.g2m for sym in gene}


def resolve_species(
    symbols: Sequence[str],
    requested: Literal["auto", "human", "mouse"],
    gene_set: str = "tirosh2016",
) -> Species | None:
    """An explicit species is returned as given. For `auto`, upper-case matches vote human,
    any-lower-case matches vote mouse, and a tie (including no match) gives None."""
    if requested != "auto":
        return requested
    known = _lowered(load_gene_set(gene_set, "human")) | _lowered(load_gene_set(gene_set, "mouse"))
    upper = mixed = 0
    for sym in symbols:
        if sym.lower() in known:
            if sym == sym.upper():
                upper += 1
            else:
                mixed += 1
    if upper == mixed:
        return None
    return "human" if upper > mixed else "mouse"


def match_genes(symbols: Sequence[str], genes: _GeneList) -> list[int]:
    """For each gene the var index of the first symbol that matches any accepted symbol,
    case-insensitively. Unmatched genes are left out and each var index appears once."""
    first: dict[str, int] = {}
    for i, sym in enumerate(symbols):
        first.setdefault(sym.lower(), i)
    out: list[int] = []
    used: set[int] = set()
    for gene in genes:
        hits = [first[s.lower()] for s in gene if s.lower() in first]
        if hits and min(hits) not in used:
            used.add(min(hits))
            out.append(min(hits))
    return out
