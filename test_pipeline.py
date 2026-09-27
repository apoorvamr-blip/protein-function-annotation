"""End-to-end pipeline tests with injected fake service clients (no network)."""

from __future__ import annotations

import json

import pytest

from protein_annotator.analysis.blast import parse_blast_xml
from protein_annotator.analysis.interpro import parse_interproscan_json
from protein_annotator.analysis.uniprot import parse_uniprot_entry
from protein_annotator.config import AnnotatorConfig
from protein_annotator.io.fasta import read_fasta
from protein_annotator.models.annotation import ConfidenceCategory, StageState
from protein_annotator.pipeline.runner import AnnotationPipeline
from protein_annotator.pipeline.stages import PipelineOptions

from .conftest import make_entry


class FakeBlast:
    def __init__(self, result):
        self.result = result
        self.calls = 0

    def search(self, protein, **_):
        self.calls += 1
        return self.result


class FakeInterPro:
    def __init__(self, result=None, error=None):
        self.result, self.error = result, error

    def annotate(self, protein, **_):
        if self.error:
            raise self.error
        return self.result


class FakeUniProt:
    def __init__(self, entries):
        self.entries = entries
        self.requested: list[str] = []

    def fetch_entries(self, accessions):
        self.requested.extend(accessions)
        found = {a: self.entries[a] for a in accessions if a in self.entries}
        return found, [f"{a}: missing" for a in accessions if a not in self.entries]

    def map_to_uniprot(self, identifiers):
        return {}


@pytest.fixture
def clients(fixtures_dir):
    blast = parse_blast_xml((fixtures_dir / "blast_example.xml").read_text())
    interpro = parse_interproscan_json((fixtures_dir / "interpro_example.json").read_text())
    hba = parse_uniprot_entry(json.loads((fixtures_dir / "uniprot_P69905.json").read_text()))
    entries = {
        "P69905": hba,
        "P69907": make_entry("P69907", "Hemoglobin subunit alpha", interpro=["IPR000971"], pfam=["PF00042"]),
        "P01942": make_entry("P01942", "Hemoglobin subunit alpha", interpro=["IPR000971"], pfam=["PF00042"]),
        "P02144": make_entry("P02144", "Myoglobin", interpro=["IPR000971"], pfam=["PF00042"]),
    }
    return FakeBlast(blast), FakeInterPro(interpro), FakeUniProt(entries)


@pytest.fixture
def protein(fixtures_dir):
    return read_fasta(fixtures_dir / "example_protein.fasta")[0]


def config(**kw):
    return AnnotatorConfig(email="user@example.org", use_cache=False, **kw)


def test_full_pipeline_with_fixtures(clients, protein):
    blast, interpro, uniprot = clients
    pipeline = AnnotationPipeline(config(), blast_client=blast, interpro_client=interpro,
                                  uniprot_client=uniprot)
    result = pipeline.run(protein)
    states = {s.name: s.state for s in result.stages}
    # The RefSeq hit (WP_000000001) cannot be mapped to UniProtKB by the fake client -> warning.
    assert states.pop("UniProt") is StageState.WARNING
    assert all(state is StageState.SUCCESS for state in states.values()), result.stages
    annotation = result.annotation
    assert annotation is not None
    assert "hemoglobin" in annotation.predicted_function.lower()
    assert annotation.confidence_category is ConfidenceCategory.STRONG
    assert annotation.enzyme_assessment.status == "no_evidence"
    assert result.properties is not None and result.properties.length == 142
    assert "P69905" in uniprot.requested


def test_failed_optional_stage_does_not_stop_pipeline(clients, protein):
    blast, _, uniprot = clients
    pipeline = AnnotationPipeline(
        config(), blast_client=blast, interpro_client=FakeInterPro(error=RuntimeError("EBI down")),
        uniprot_client=uniprot,
    )
    result = pipeline.run(protein)
    interpro = next(s for s in result.stages if s.name == "InterPro")
    assert interpro.state is StageState.FAILED and "EBI down" in interpro.message
    assert result.annotation is not None
    assert result.annotation.confidence_category is not ConfidenceCategory.STRONG
    assert any("InterPro" in limit or "domain" in limit.lower() for limit in result.annotation.limitations)


def test_skips_and_missing_email(clients, protein):
    blast, interpro, uniprot = clients
    pipeline = AnnotationPipeline(
        AnnotatorConfig(email=None, use_cache=False), PipelineOptions(skip_blast=True),
        blast_client=blast, interpro_client=interpro, uniprot_client=uniprot,
    )
    result = pipeline.run(protein)
    states = {s.name: s.state for s in result.stages}
    assert states["BLASTP"] is StageState.SKIPPED
    assert states["InterPro"] is StageState.SKIPPED  # no e-mail configured
    assert states["UniProt"] is StageState.SKIPPED
    assert blast.calls == 0
    assert result.annotation.confidence_category is ConfidenceCategory.INSUFFICIENT
