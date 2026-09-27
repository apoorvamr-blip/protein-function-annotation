"""Contradiction-detection tests."""

from __future__ import annotations

from protein_annotator.integration.contradiction import competing_candidates, detect_contradictions

from .scenarios import (
    SETTINGS,
    conflicting_summary,
    empty_summary,
    enzyme_without_domain_summary,
    strong_summary,
)


def kinds(contradictions):
    return {c.kind for c in contradictions}


def test_consistent_evidence_has_no_major_contradictions():
    summary, interpro = strong_summary()
    found = detect_contradictions(summary, SETTINGS, interpro)
    assert not [c for c in found if c.severity == "major"]


def test_conflicting_homologs_flagged_with_competitors():
    summary = conflicting_summary()
    assert competing_candidates(summary, SETTINGS)
    found = detect_contradictions(summary, SETTINGS, None)
    conflict = next(c for c in found if c.kind == "conflicting_homolog_annotations")
    assert conflict.severity == "major"
    assert "Conflicting functional evidence detected" in conflict.description
    text = " ".join(conflict.competing_interpretations)
    assert "Hemoglobin subunit alpha" in text and "DNA gyrase subunit B" in text


def test_enzyme_without_catalytic_domain_flagged():
    summary, interpro = enzyme_without_domain_summary()
    found = detect_contradictions(summary, SETTINGS, interpro)
    assert {"enzyme_without_catalytic_domain", "homology_domain_mismatch"} <= kinds(found)
    assert found[0].severity == "major"  # sorted most severe first


def test_no_evidence_no_contradictions():
    assert detect_contradictions(empty_summary(), SETTINGS, None) == []
