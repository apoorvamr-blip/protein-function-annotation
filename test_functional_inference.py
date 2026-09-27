"""Functional-inference tests: confidence categories and conservative wording."""

from __future__ import annotations

from protein_annotator.integration.contradiction import detect_contradictions
from protein_annotator.integration.evidence import build_evidence_summary
from protein_annotator.integration.functional_inference import assign_confidence, infer_function
from protein_annotator.models.annotation import ConfidenceCategory, EvidenceType

from .conftest import make_interpro
from .scenarios import (
    GLOBIN_DOMAINS,
    SETTINGS,
    blast,
    conflicting_summary,
    empty_summary,
    enzyme_without_domain_summary,
    globin_hits_and_entries,
    strong_summary,
    weak_summary,
)


def annotate(summary, interpro=None):
    contradictions = detect_contradictions(summary, SETTINGS, interpro)
    return infer_function(summary, contradictions, SETTINGS, interpro=interpro)


def test_strong_evidence():
    summary, interpro = strong_summary()
    annotation = annotate(summary, interpro)
    assert annotation.confidence_category is ConfidenceCategory.STRONG
    assert annotation.predicted_function == "Hemoglobin subunit alpha"
    assert all(line.startswith("✓") for line in annotation.confidence_rationale)
    types = {item.evidence_type for item in annotation.evidence}
    assert EvidenceType.HOMOLOGY in types and EvidenceType.DOMAIN in types
    assert annotation.enzyme_assessment.status == "no_evidence"
    assert not annotation.conflicting_evidence
    assert "definitely" not in annotation.interpretation.lower()


def test_strong_requires_domain_support():
    hits, entries = globin_hits_and_entries()
    summary = build_evidence_summary(blast(hits), None, entries, SETTINGS, uniprot_available=True)
    category, rationale = assign_confidence(summary, [])
    assert category is ConfidenceCategory.MODERATE
    assert any("Not strong" in r and "domains" in r for r in rationale)
    assert annotate(summary).predicted_function.startswith("Probable ")


def test_weak_evidence_is_limited_and_putative():
    annotation = annotate(weak_summary())
    assert annotation.confidence_category is ConfidenceCategory.LIMITED
    assert annotation.predicted_function.startswith("Putative ")
    assert any(r.startswith("✗") for r in annotation.confidence_rationale)


def test_conflicting_evidence_is_not_resolved_blindly():
    annotation = annotate(conflicting_summary())
    assert annotation.conflicting_evidence
    assert annotation.confidence_category is ConfidenceCategory.LIMITED
    assert "(conflicting evidence)" in annotation.predicted_function
    assert annotation.alternative_functions
    assert any("competing" in alt for alt in annotation.alternative_functions)


def test_no_evidence_is_insufficient():
    annotation = annotate(empty_summary())
    assert annotation.confidence_category is ConfidenceCategory.INSUFFICIENT
    assert annotation.predicted_function.startswith("Uncharacterized protein")
    assert annotation.limitations


def test_domain_only_evidence_is_limited():
    interpro = make_interpro(GLOBIN_DOMAINS)
    summary = build_evidence_summary(None, interpro, {}, SETTINGS, uniprot_available=False)
    annotation = annotate(summary, interpro)
    assert annotation.confidence_category is ConfidenceCategory.LIMITED
    assert "Globin" in annotation.predicted_function
    assert any(item.evidence_type is EvidenceType.DOMAIN for item in annotation.evidence)


def test_enzyme_contradiction_downgrades_confidence():
    summary, interpro = enzyme_without_domain_summary()
    annotation = annotate(summary, interpro)
    assert annotation.confidence_category is ConfidenceCategory.LIMITED
    assert annotation.enzyme_assessment.status == "homology_only"
    assert any(c.kind == "enzyme_without_catalytic_domain" for c in annotation.contradictions)


def test_categories_are_not_probabilities():
    annotation = annotate(strong_summary()[0], strong_summary()[1])
    dumped = annotation.model_dump_json()
    assert "%" not in annotation.confidence_category.value
    assert "probability" not in dumped.lower()
