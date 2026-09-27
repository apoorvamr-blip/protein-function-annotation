"""InterProScan (EMBL-EBI) parsing and client tests using saved/mock responses."""

from __future__ import annotations

import itertools

import pytest

from protein_annotator.analysis.domains import (
    functional_domains,
    integrated_entries,
    interpro_go_terms,
    query_domain_ids,
)
from protein_annotator.analysis.interpro import (
    InterProClient,
    InterProError,
    InterProTimeoutError,
    parse_interproscan_json,
)
from protein_annotator.config import InterProSettings
from protein_annotator.models.sequence import ProteinSequence

from .conftest import FakeResponse

BASE = "https://ebi.example/iprscan5"
PROTEIN = ProteinSequence(identifier="q", sequence="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF")


@pytest.fixture
def interpro_json(fixtures_dir):
    return (fixtures_dir / "interpro_example.json").read_text()


def test_parse_fixture(interpro_json):
    result = parse_interproscan_json(interpro_json, job_id="iprscan5-R1")
    assert result.job_id == "iprscan5-R1"
    pfam = [d for d in result.domains if d.database == "PFAM"]
    assert len(pfam) == 1
    assert pfam[0].accession == "PF00042"
    assert pfam[0].interpro_accession == "IPR000971"
    assert (pfam[0].start, pfam[0].end) == (27, 138)
    assert pfam[0].go_terms[0].id == "GO:0020037"
    assert pfam[0].go_terms[0].category == "molecular_function"


def test_domain_helpers(interpro_json):
    result = parse_interproscan_json(interpro_json)
    assert query_domain_ids(result) == {"IPR000971", "PF00042"}  # InterPro + Pfam IDs
    assert "IPR000971" in integrated_entries(result)
    # disorder predictions (MobiDB-lite) are not functional domains
    assert all(d.database != "MOBIDB_LITE" for d in functional_domains(result))
    assert [g.id for g in interpro_go_terms(result)] == ["GO:0020037"]


# Structure copied from the JSON example in the official InterProScan 5 documentation
# (interproscan-docs, OutputFormats: "JavaScript Object Notation (JSON)").
DOCUMENTED_EXAMPLE = {
    "interproscan-version": "5.26-65.0",
    "results": [{
        "sequence": "MSKIGKSIRLERIIDRKTRKTVIVPMDHG",
        "md5": "88d47cc807fe8e977130b0cc93e0bd61",
        "matches": [{
            "signature": {
                "accession": "PIRSF038992", "name": "Aldolase_Ia", "description": None, "type": None,
                "signatureLibraryRelease": {"library": "PIRSF", "version": "3.01"},
                "models": {"PIRSF038992": {"accession": "PIRSF038992", "name": "Aldolase_Ia",
                                           "description": None, "key": "PIRSF038992"}},
                "entry": {
                    "accession": "IPR002915", "name": "DeoC/FbaB/lacD_aldolase",
                    "description": "DeoC/FbaB/ lacD aldolase", "type": "FAMILY",
                    "goXRefs": [{"identifier": "GO:0016829", "name": "lyase activity",
                                 "databaseName": "GO", "category": "MOLECULAR_FUNCTION"}],
                    "pathwayXRefs": [{"identifier": "R-HSA-71336",
                                      "name": "Pentose phosphate pathway (hexose monophosphate shunt)",
                                      "databaseName": "Reactome"}],
                },
            },
            "locations": [{"start": 1, "end": 265, "hmmStart": 2, "hmmEnd": 262, "hmmBounds": "INCOMPLETE",
                           "evalue": 3.3e-94, "score": 302.6, "envelopeStart": 1, "envelopeEnd": 265}],
            "evalue": 3.0e-94,
            "score": 302.7,
        }],
    }],
}


def test_parse_documented_json_format():
    result = parse_interproscan_json(DOCUMENTED_EXAMPLE)
    assert result.interproscan_version == "5.26-65.0"
    (domain,) = result.domains
    assert (domain.database, domain.accession) == ("PIRSF", "PIRSF038992")
    assert (domain.interpro_accession, domain.interpro_type) == ("IPR002915", "FAMILY")
    assert (domain.start, domain.end, domain.evalue) == (1, 265, 3.3e-94)
    assert [(g.id, g.category) for g in domain.go_terms] == [("GO:0016829", "molecular_function")]
    assert domain.pathways == ["Reactome:R-HSA-71336 Pentose phosphate pathway (hexose monophosphate shunt)"]


def test_parse_empty_results():
    result = parse_interproscan_json({"interproscan-version": "5", "results": [{"matches": []}]})
    assert result.domains == []
    assert functional_domains(result) == []


def _client(fake_http, *, email="user@example.org", timeout=100.0):
    clock = itertools.count(0, 15)
    return InterProClient(
        fake_http, InterProSettings(timeout=timeout, poll_interval=15), base_url=BASE, email=email,
        sleep=lambda _s: None, clock=lambda: float(next(clock)),
    )


def test_requires_email(fake_http):
    with pytest.raises(InterProError, match="e-mail"):
        _client(fake_http, email=None).annotate(PROTEIN)
    assert fake_http.calls == []


def test_submit_poll_fetch(fake_http, interpro_json):
    statuses = iter(["RUNNING", "FINISHED"])
    fake_http.route("POST", "/run", FakeResponse("iprscan5-R20260924-000001-0001-12345678-p1m"))
    fake_http.route("GET", "/status/", lambda **_: FakeResponse(next(statuses)))
    fake_http.route("GET", "/result/", FakeResponse(interpro_json))

    result = _client(fake_http).annotate(PROTEIN)
    assert result.job_id.startswith("iprscan5-")
    assert any(d.accession == "PF00042" for d in result.domains)
    post = fake_http.calls[0][2]["data"]
    assert post["email"] == "user@example.org"
    assert post["goterms"] == "true" and post["pathways"] == "true"
    assert "stype" not in post  # the service defaults to protein, as in EBI's official client


def test_applications_sent_comma_separated(fake_http, interpro_json):
    fake_http.route("POST", "/run", FakeResponse("iprscan5-job"))
    fake_http.route("GET", "/status/", FakeResponse("FINISHED"))
    fake_http.route("GET", "/result/", FakeResponse(interpro_json))
    clock = itertools.count(0, 15)
    client = InterProClient(
        fake_http, InterProSettings(applications=("PfamA", "CDD")), base_url=BASE,
        email="user@example.org", sleep=lambda _s: None, clock=lambda: float(next(clock)),
    )
    client.annotate(PROTEIN)
    assert fake_http.calls[0][2]["data"]["appl"] == "PfamA,CDD"


def test_failed_job(fake_http):
    fake_http.route("POST", "/run", FakeResponse("iprscan5-job"))
    fake_http.route("GET", "/status/", FakeResponse("ERROR"))
    with pytest.raises(InterProError, match="ERROR"):
        _client(fake_http).annotate(PROTEIN)


def test_timeout(fake_http):
    fake_http.route("POST", "/run", FakeResponse("iprscan5-job"))
    fake_http.route("GET", "/status/", FakeResponse("RUNNING"))
    with pytest.raises(InterProTimeoutError):
        _client(fake_http, timeout=30).annotate(PROTEIN)


def test_rejects_unexpected_job_id(fake_http):
    fake_http.route("POST", "/run", FakeResponse("<html>Service error</html>"))
    with pytest.raises(InterProError):
        _client(fake_http).annotate(PROTEIN)
