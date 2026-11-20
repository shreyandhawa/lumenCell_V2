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
