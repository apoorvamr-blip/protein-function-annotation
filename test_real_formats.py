"""Regression tests against genuine service output (see tests/fixtures/README.md for provenance).

These guard against the parsers drifting from the real formats, which hand-built
fixtures alone cannot detect.
"""

from __future__ import annotations

import json

import pytest
from Bio.Blast import NCBIXML

from protein_annotator.analysis.blast import parse_blast_xml
from protein_annotator.analysis.uniprot import mapping_source_database, parse_uniprot_entry
from protein_annotator.config import InferenceSettings
from protein_annotator.integration.contradiction import detect_contradictions
from protein_annotator.integration.evidence import build_evidence_summary
from protein_annotator.integration.functional_inference import infer_function
from protein_annotator.models.annotation import ConfidenceCategory


@pytest.fixture
def real_blast(fixtures_dir):
    return fixtures_dir / "real_ncbi_blastp_nr.xml"


def test_blast_parser_agrees_with_biopython(real_blast):
    ours = parse_blast_xml(real_blast.read_text())
    with open(real_blast) as handle:
        record = next(NCBIXML.parse(handle))
    assert ours.database == "nr"
    assert ours.query_length == record.query_length == 239
    assert len(ours.hits) == len(record.alignments) == 10
    for hit, alignment in zip(ours.hits, record.alignments, strict=True):
        best = min(alignment.hsps, key=lambda h: (h.expect, -h.bits))
        assert hit.evalue == pytest.approx(best.expect)
        assert hit.bitscore == pytest.approx(best.bits)
        assert hit.identity == pytest.approx(100 * best.identities / best.align_length, abs=0.01)
        assert hit.accession == alignment.accession.split(".")[0]


def test_real_nr_hits_are_mapped_via_refseq(real_blast):
    hits = parse_blast_xml(real_blast.read_text()).hits
    assert all(h.uniprot_accession is None for h in hits)  # nr hits need ID mapping
    refseq = [h for h in hits if h.accession.startswith("WP_")]
    assert refseq and all(mapping_source_database(h.sequence_id) == "RefSeq_Protein" for h in refseq)


def test_real_blast_integration_is_conservative(real_blast):
    settings = InferenceSettings()
    blast = parse_blast_xml(real_blast.read_text())
    summary = build_evidence_summary(blast, None, {}, settings, uniprot_available=False)
    annotation = infer_function(summary, detect_contradictions(summary, settings), settings)
    assert summary.candidates[0].hit_count == 10
    assert "sigma factor SigE" in annotation.predicted_function
    # No curated homologs and no domain analysis -> cannot be "strong evidence".
    assert annotation.confidence_category is ConfidenceCategory.MODERATE
    assert annotation.predicted_function.startswith("Probable ")


def test_real_uniprot_entry(fixtures_dir):
    data = json.loads((fixtures_dir / "real_uniprot_P30542.json").read_text())
    entry = parse_uniprot_entry(data)
    assert entry.accession == "P30542"
    assert entry.entry_name == "AA1R_HUMAN"
    assert entry.reviewed and entry.annotation_status == "Reviewed / curated"
    assert entry.protein_name == "Adenosine receptor A1"
    assert entry.gene == "ADORA1"
    assert entry.organism == "Homo sapiens"
    assert entry.function and entry.function[0].startswith("Receptor for adenosine")
    assert entry.ec_numbers == []  # a receptor, not an enzyme
    assert "Cell membrane" in entry.subcellular_locations
    assert "PF00001" in entry.pfam_ids and entry.interpro_ids
    go = {g.id: g for g in entry.go_terms}
    assert len(go) > 50
    assert go["GO:0030673"].category == "cellular_component"
    assert go["GO:0030673"].name == "axolemma"
    assert go["GO:0030673"].evidence == "IEA"
    assert entry.has_experimental_function  # experimental GO evidence codes present
