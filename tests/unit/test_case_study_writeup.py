"""The Sheih 2020 write-up (phase 3 Task 13): its numbers match the committed results, and its text
claims no direction from clones, cites its sources, states the failed verification and the cap, and
leads with the per-patient runs.
"""

import json
import re

from lumencell.bench.case_study import render_markdown
from lumencell.bench.gates import REPO
from lumencell.trust.constants import VERDICT_VERIFIED_ON

DOC, START, END = (
    REPO / "docs/case-studies/sheih-2020.md",
    "<!-- case-study:start -->",
    "<!-- case-study:end -->",
)
DIRECTION = (
    "give rise",
    "gives rise",
    "gave rise",
    "differentiat",
    "precursor",
    "derived from",
    "become",
    "became",
    "turn into",
    "turns into",
    "->",
    "→",
)


def test_generated_section_matches_the_committed_results() -> None:  # Review Focus 5
    latest = sorted((REPO / "benchmarks/case_studies/results").glob("sheih2020-*.json"))[-1]
    section = DOC.read_text(encoding="utf-8").split(START, 1)[1].split(END, 1)[0]
    assert section.strip() == render_markdown(json.loads(latest.read_text(encoding="utf-8"))).strip()


CLONAL_WORDS = ("clon", "tcr", "receptor", "shar", "coupl")


def test_no_sentence_about_clones_claims_direction() -> None:
    head, rest = DOC.read_text(encoding="utf-8").split(START, 1)
    for sentence in re.split(r"(?<=[.!?])\s+", head + rest.split(END, 1)[1]):
        if any(w in sentence.lower() for w in CLONAL_WORDS):
            assert not any(p in sentence.lower() for p in DIRECTION), sentence


def test_the_paper_the_data_and_the_methods_are_cited() -> None:
    text = DOC.read_text(encoding="utf-8")
    assert "10.1038/s41467-019-13880-1" in text and "GSE125881" in text  # Sheih et al. 2020 and its data
    for name, doi in (
        ("STARTRAC", "10.1038/s41586-018-0694-x"),
        ("LRT", "10.1371/journal.pcbi.1011300"),
        ("CoSpar", "10.1038/s41587-022-01209-1"),
    ):  # spec §10, checked
        assert name in text and doi in text, name


def test_the_failed_verification_and_the_cap_are_stated() -> None:  # (amended 2026-11-19, phase 2 X-1)
    text = DOC.read_text(encoding="utf-8")
    assert "failed its verification" in text and VERDICT_VERIFIED_ON in text  # from lumencell.trust.constants
    latest = sorted((REPO / "benchmarks/case_studies/results").glob("sheih2020-*.json"))[-1]
    report = json.loads(latest.read_text(encoding="utf-8"))
    levels = {m["tree_support"]["level"] for u in report["units"].values() for m in u["methods"].values()}
    assert "supported" not in levels  # the cap makes it impossible; the text must not claim it


def test_per_patient_results_lead_and_pooled_is_only_a_diagnostic() -> None:
    text = DOC.read_text(encoding="utf-8")
    assert text.index("Per patient (primary)") < text.index("Pooled diagnostic")
    head, rest = text.split(START, 1)
    for sentence in re.split(r"(?<=[.!?])\s+", head + rest.split(END, 1)[1]):
        if "pooled" in sentence.lower():
            assert "diagnostic" in sentence.lower(), sentence


# Fix round 1: guards for the rest of the plan's rules for the text.


def _parts() -> tuple[str, str, str]:
    """The text before the generated section, the generated section, and the text after it."""
    head, rest = DOC.read_text(encoding="utf-8").split(START, 1)
    generated, tail = rest.split(END, 1)
    return head, generated, tail


def _sentences(text: str) -> list[str]:
    return re.split(r"(?<=[.!?])\s+", text)


def _section(text: str, heading: str) -> str:
    """One `## ` section's body, up to the next `## ` heading."""
    body = text.split(f"\n## {heading}\n", 1)[1]
    return re.split(r"^## ", body, maxsplit=1, flags=re.M)[0]


NUMBER = re.compile(r"(?<![\w.#-])-?\d[\d,]*(?:\.\d+)?(?!\w)")
NOT_IN_TABLES = {
    "0.01",  # the declared C5 bar, [clonal] in benchmarks/thresholds.toml
    "1.0",  # coverage, each method's "coverage" in the results file
    "4.3",  # a section number: spec §4.3
}


def test_every_decimal_in_the_reading_appears_in_the_generated_section() -> None:
    _, generated, tail = _parts()
    in_tables = {n.lstrip("-") for n in NUMBER.findall(generated)}
    reading = _section(tail, "Reading")
    decimals = {n.lstrip("-") for n in NUMBER.findall(reading) if "." in n}
    assert not decimals - in_tables - NOT_IN_TABLES, sorted(decimals - in_tables - NOT_IN_TABLES)


DIRECTION_MORE = (
    "expanded from",
    "originat",
    "progress",
    "lineage of origin",
    "descend",
    "arise",
    "upstream",
    "downstream",
    "precede",
    "transition",
    "convert",
    "lead to",
    "leads to",
    "ancestor",
)


def test_no_sentence_about_clones_uses_a_broader_direction_word() -> None:
    head, _, tail = _parts()
    for sentence in _sentences(head + tail):
        low = sentence.lower()
        if any(w in low for w in CLONAL_WORDS):
            assert not any(p in low for p in DIRECTION_MORE), sentence


SUPPORT = re.compile(r"\bsupport(?:s|ed|ing)?\b")
NEGATION = re.compile(r"\b(?:not|never|no|none|cannot)\b")


def test_provisional_is_never_called_support() -> None:
    head, _, tail = _parts()
    for sentence in _sentences(head + tail):
        low = sentence.lower()
        if "provisional" not in low:
            continue
        for match in SUPPORT.finditer(low):
            clause = re.split(r"[,;:]", low[: match.start()])[-1]  # from the last clause break to the word
            assert NEGATION.search(clause), sentence


def test_the_prose_itself_states_the_failed_verification() -> None:  # not only the generated note
    head, _, tail = _parts()
    assert "failed its verification" in head + tail and VERDICT_VERIFIED_ON in head + tail


def test_the_plans_five_limits_are_stated() -> None:
    limits = _section(_parts()[2], "Limits")
    for item in (
        "quarter of the cells",
        "no batch correction",
        "Unpaired chains",
        "Only four patients",
        "infusion product is manufactured",
    ):
        assert item in limits, item


SECTIONS = (
    "Question",
    "Data",
    "Declaration",
    "What the paper reports",
    "Results",
    "Reading",
    "Limits",
    "Reproduce",
)


def test_the_plans_sections_unverified_marks_and_the_root_statement() -> None:
    text = DOC.read_text(encoding="utf-8")
    headings = re.findall(r"^## (.+)$", text, flags=re.M)
    assert all(s in headings for s in SECTIONS), headings
    head, _, tail = _parts()
    assert "(unverified" in head + tail
    assert re.search(r"pseudotime direction\W+comes from the declared root", (head + tail).lower())


def test_the_reading_comes_after_the_results() -> None:
    text = DOC.read_text(encoding="utf-8")
    assert text.index("\n## Results\n") < text.index(START) < text.index(END) < text.index("\n## Reading\n")
