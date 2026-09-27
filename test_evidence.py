"""Evidence-integration tests: agreement, weak, conflicting and absent evidence."""

from __future__ import annotations

import pytest

from protein_annotator.integration.evidence import (
    build_candidates,
    description_tokens,
    hit_weight,
    is_uninformative,
    jaccard,
    normalize_description,
    select_significant_hits,
    suggests_enzyme,
)

from .conftest import make_entry, make_hit
from .scenarios import (
    SETTINGS,
    conflicting_summary,
    empty_summary,
    enzyme_without_domain_summary,
    strong_summary,
    weak_summary,
)


def test_normalize_description_uses_recname():
    text = "RecName: Full=Hemoglobin subunit alpha; AltName: Full=Alpha-globin [Homo sapiens]"
    assert normalize_description(text) == "Hemoglobin subunit alpha"
    assert normalize_description("hypothetical protein [Escherichia coli]").lower().startswith("hypothetical")


@pytest.mark.parametrize("label", ["hypothetical protein", "uncharacterized protein", "Unnamed protein product"])
def test_uninformative_labels(label):
    assert is_uninformative(label)


def test_informative_label():
    assert not is_uninformative("Hemoglobin subunit alpha")


def test_enzyme_keyword_detection():
    assert suggests_enzyme(description_tokens("Serine protease trypsin"))
    assert not suggests_enzyme(description_tokens("Hemoglobin subunit alpha"))


def test_jaccard():
    assert jaccard({"a", "b"}, {"a", "b"}) == 1.0
    assert jaccard({"a"}, {"b"}) == 0.0
    assert jaccard(set(), set()) == 0.0


def test_significant_hits_filter_evalue_and_coverage():
    hits = [
        make_hit("A", "x", evalue=1e-10, coverage=0.9),
        make_hit("B", "x", evalue=1e-2, coverage=0.9),   # not significant
        make_hit("C", "x", evalue=1e-10, coverage=0.2),  # low coverage
    ]
    adequate, low = select_significant_hits(hits, SETTINGS)
    assert [h.accession for h in adequate] == ["A"]
    assert [h.accession for h in low] == ["C"]


def test_hit_weight_rewards_curation_and_experiment():
    hit = make_hit("A", "x", identity=80, coverage=1.0, evalue=1e-60)
    plain = hit_weight(hit)
    reviewed = hit_weight(hit, make_entry("A", "x", reviewed=True))
    experimental = hit_weight(hit, make_entry("A", "x", reviewed=True, experimental=True))
    assert plain == pytest.approx(0.8)
    assert reviewed == pytest.approx(1.0)
    assert experimental == pytest.approx(1.2)


def test_clustering_groups_similar_descriptions():
    hits = [
        make_hit("A", "Hemoglobin subunit alpha [x]", bitscore=300),
        make_hit("B", "Hemoglobin subunit alpha-1 [y]", bitscore=290),
        make_hit("C", "hypothetical protein [z]", bitscore=280),
        make_hit("D", "Myoglobin [w]", bitscore=100),
    ]
    candidates, uninformative = build_candidates(hits, {}, SETTINGS)
    assert uninformative == 1
    assert candidates[0].hit_count == 2
    assert candidates[0].hit_accessions == ["A", "B"]
    assert any(c.label == "Myoglobin" for c in candidates[1:])


def test_strong_agreement_summary():
    summary, _ = strong_summary()
    top = summary.candidates[0]
    assert len(summary.candidates) == 1
    assert top.label == "Hemoglobin subunit alpha"
    assert top.hit_count == 3
    # every entry carries experimental (IDA) GO annotations, which count as experimental support
    assert top.reviewed_count == 3 and top.experimental_count == 3
    assert top.domain_support is True
    assert summary.enzyme.status == "no_evidence"
    go = {g.id: g for g in summary.supported_go_terms}
    assert go["GO:0005344"].experimental and go["GO:0005344"].homolog_count == 3
    assert "GO:0020037" not in go or go["GO:0020037"].sources  # InterPro2GO only if present


def test_weak_evidence_summary():
    summary = weak_summary()
    assert summary.significant_hit_count == 1
    assert summary.candidates[0].domain_support is None
    assert summary.interpro_available is False


def test_conflicting_summary_has_two_candidates():
    summary = conflicting_summary()
    labels = {c.label for c in summary.candidates}
    assert labels == {"Hemoglobin subunit alpha", "DNA gyrase subunit B"}


def test_no_evidence_summary():
    summary = empty_summary()
    assert summary.candidates == []
    assert summary.total_hits == 0
    assert summary.enzyme.status == "no_evidence"


def test_enzyme_homology_without_catalytic_domain():
    summary, _ = enzyme_without_domain_summary()
    assert summary.enzyme.status == "homology_only"
    assert summary.enzyme.ec_numbers == ["3.4.21.4"]
    assert summary.candidates[0].domain_support is False
