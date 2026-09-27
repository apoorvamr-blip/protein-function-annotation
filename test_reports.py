"""Report-generation tests (JSON, CSV, Markdown)."""

from __future__ import annotations

import csv
import json

from protein_annotator.io.reports import render_markdown, write_reports
from protein_annotator.pipeline.runner import AnnotationPipeline

from .test_pipeline import clients, config, protein  # noqa: F401  (pytest fixtures)

EXPECTED = {
    "annotation.json", "annotation.md", "blast_hits.csv", "domains.csv",
    "sequence_statistics.json", "complete_results.json",
}


def run(clients, protein):  # noqa: F811
    blast, interpro, uniprot = clients
    return AnnotationPipeline(config(), blast_client=blast, interpro_client=interpro,
                              uniprot_client=uniprot).run(protein)


def test_write_reports(tmp_path, clients, protein):  # noqa: F811
    result = run(clients, protein)
    written = write_reports(result, tmp_path / "out")
    names = {p.name for p in (tmp_path / "out").iterdir()}
    assert EXPECTED <= names
    assert all(path.exists() for path in written.values())

    annotation = json.loads((tmp_path / "out" / "annotation.json").read_text())
    for key in ("predicted_function", "supporting_evidence", "alternative_functions",
                "confidence_category", "limitations"):
        assert key in annotation

    with open(tmp_path / "out" / "blast_hits.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 5 and rows[0]["accession"] == "P69905"

    with open(tmp_path / "out" / "domains.csv", newline="") as handle:
        assert any(r["accession"] == "PF00042" for r in csv.DictReader(handle))

    stats = json.loads((tmp_path / "out" / "sequence_statistics.json").read_text())
    assert stats["properties"]["length"] == 142
    assert stats["properties"]["molecular_weight"] > 15000

    complete = json.loads((tmp_path / "out" / "complete_results.json").read_text())
    assert complete["sequence"]["sha256"] == protein.sha256


def test_markdown_sections(clients, protein):  # noqa: F811
    text = render_markdown(run(clients, protein))
    for heading in ("# Protein Functional Annotation Report", "## Sequence", "## Homology (BLASTP)",
                    "## Domains (InterProScan)", "## UniProt evidence", "## GO annotations",
                    "## Functional interpretation", "### Supporting evidence", "### Limitations"):
        assert heading in text
    assert "Reviewed / curated" in text
