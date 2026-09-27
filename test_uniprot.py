"""UniProt REST parsing and client tests using saved/mock responses."""

from __future__ import annotations

import json

import pytest

from protein_annotator.analysis.uniprot import (
    UniProtClient,
    is_uniprot_accession,
    mapping_source_database,
    parse_uniprot_entry,
)
from protein_annotator.config import UniProtSettings
from protein_annotator.utils.http import HttpError

from .conftest import FakeResponse

BASE = "https://uniprot.example"


@pytest.fixture
def entry_json(fixtures_dir):
    return (fixtures_dir / "uniprot_P69905.json").read_text()


def _client(fake_http):
    return UniProtClient(fake_http, UniProtSettings(), base_url=BASE, sleep=lambda _s: None)


def test_parse_reviewed_entry(entry_json):
    entry = parse_uniprot_entry(json.loads(entry_json))
    assert entry.accession == "P69905"
    assert entry.reviewed is True
    assert entry.annotation_status == "Reviewed / curated"
    assert entry.protein_name == "Hemoglobin subunit alpha"
    assert entry.gene == "HBA1"
    assert entry.organism == "Homo sapiens"
    assert entry.has_experimental_function
    assert entry.ec_numbers == []
    assert "PF00042" in entry.pfam_ids and "IPR000971" in entry.interpro_ids
    go = {g.id: g for g in entry.go_terms}
    assert go["GO:0005344"].category == "molecular_function"
    assert go["GO:0005344"].evidence == "IDA"
    assert "Cytoplasm" in entry.subcellular_locations


def test_unreviewed_status(entry_json):
    data = json.loads(entry_json)
    data["entryType"] = "UniProtKB unreviewed (TrEMBL)"
    entry = parse_uniprot_entry(data)
    assert entry.reviewed is False
    assert entry.annotation_status == "Unreviewed / computational"


def test_inactive_entry_returns_none():
    assert parse_uniprot_entry({"primaryAccession": "P00000", "entryType": "Inactive"}) is None


def test_accession_helpers():
    assert is_uniprot_accession("P69905")
    assert is_uniprot_accession("A0A023GPI8")
    assert not is_uniprot_accession("WP_000000001")
    assert mapping_source_database("WP_000000001.1") == "RefSeq_Protein"
    assert mapping_source_database("AAA12345.1") == "EMBL-GenBank-DDBJ_CDS"


def test_fetch_entry(fake_http, entry_json):
    fake_http.route("GET", "/uniprotkb/P69905.json", FakeResponse(entry_json))
    entry = _client(fake_http).fetch_entry("P69905")
    assert entry is not None and entry.accession == "P69905"


def test_fetch_entry_not_found(fake_http):
    def not_found(**_):
        raise HttpError("HTTP 404", status_code=404)

    fake_http.route("GET", "/uniprotkb/", not_found)
    entries, failures = _client(fake_http).fetch_entries(["Q00000"])
    assert entries == {}
    assert failures and "not found" in failures[0]


def test_id_mapping(fake_http):
    fake_http.route("POST", "/idmapping/run", FakeResponse(payload={"jobId": "map1"}))
    fake_http.route("GET", "/idmapping/status/map1", FakeResponse(payload={"results": []}))
    fake_http.route("GET", "/idmapping/results/map1", FakeResponse(payload={
        "results": [
            {"from": "WP_000000001.1", "to": {"primaryAccession": "A0A000AAA1"}},
            {"from": "WP_000000001.1", "to": {"primaryAccession": "A0A000AAA2"}},
        ]
    }))
    mapping = _client(fake_http).map_to_uniprot(["WP_000000001.1"])
    assert mapping == {"WP_000000001.1": "A0A000AAA1"}
    assert fake_http.calls[0][2]["data"]["from"] == "RefSeq_Protein"
