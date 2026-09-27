"""BLASTP parsing and client tests (no network: NCBI is replaced by FakeHttp)."""

from __future__ import annotations

import itertools

import pytest

from protein_annotator.analysis.blast import (
    BlastClient,
    BlastError,
    BlastParseError,
    BlastTimeoutError,
    extract_organism,
    merged_coverage,
    parse_blast_xml,
    parse_hit_identifiers,
    parse_put_response,
    parse_search_info,
)
from protein_annotator.config import BlastSettings
from protein_annotator.models.sequence import ProteinSequence
from protein_annotator.utils.caching import ResultCache

from .conftest import FakeResponse

BLAST_URL = "https://blast.example/Blast.cgi"
PROTEIN = ProteinSequence(identifier="q", sequence="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF")


@pytest.fixture
def blast_xml(fixtures_dir):
    return (fixtures_dir / "blast_example.xml").read_text()


def test_parse_fixture_hits(blast_xml):
    result = parse_blast_xml(blast_xml)
    assert result.database == "swissprot"
    assert result.query_length == 142
    assert [h.accession for h in result.hits] == ["P69905", "P69907", "P01942", "P02144", "WP_000000001"]
    top = result.hits[0]
    assert top.rank == 1
    assert top.identity == 100.0
    assert top.query_coverage == 1.0
    assert top.evalue == pytest.approx(1e-99)
    assert top.organism == "Homo sapiens"
    assert top.uniprot_accession == "P69905"
    assert top.description.startswith("RecName: Full=Hemoglobin subunit alpha")
    assert (top.query_start, top.query_end) == (1, 142)


def test_multi_hsp_coverage_and_refseq_ids(blast_xml):
    hit = parse_blast_xml(blast_xml).hits[-1]
    assert hit.num_hsps == 2
    assert hit.uniprot_accession is None
    assert hit.sequence_id == "WP_000000001.1"
    assert 0.0 < hit.query_coverage < 0.5


def test_invalid_xml_raises():
    with pytest.raises(BlastParseError):
        parse_blast_xml("<not-closed")
    with pytest.raises(BlastParseError):
        parse_blast_xml("<Other/>")


def test_merged_coverage_merges_overlaps():
    assert merged_coverage([(1, 50), (40, 100)], 100) == pytest.approx(1.0)
    assert merged_coverage([(1, 10), (21, 30)], 100) == pytest.approx(0.2)
    assert merged_coverage([], 100) == 0.0


def test_identifier_parsing():
    acc, versioned, uniprot = parse_hit_identifiers("sp|P69905.2|HBA_HUMAN", "P69905")
    assert (acc, versioned, uniprot) == ("P69905", "P69905.2", "P69905")
    acc, versioned, uniprot = parse_hit_identifiers("ref|WP_000000001.1|", "WP_000000001")
    assert acc == "WP_000000001" and uniprot is None


def test_extract_organism():
    assert extract_organism("Myoglobin [Homo sapiens]") == "Homo sapiens"
    assert extract_organism("no organism here") is None


def test_put_and_search_info_parsing():
    rid, rtoe = parse_put_response("<!--QBlastInfoBegin\n    RID = ABC123XYZ\n    RTOE = 27\nQBlastInfoEnd-->")
    assert (rid, rtoe) == ("ABC123XYZ", 27)
    with pytest.raises(BlastError):
        parse_put_response("<html>error</html>")
    assert parse_search_info("Status=WAITING") == ("WAITING", None)
    assert parse_search_info("Status=READY\nThereAreHits=yes") == ("READY", True)


def _client(fake_http, *, timeout=600.0, cache=None):
    clock = itertools.count(0, 60)
    return BlastClient(
        fake_http, BlastSettings(timeout=timeout, poll_interval=60), base_url=BLAST_URL,
        tool="test-tool", email="user@example.org", cache=cache,
        sleep=lambda _s: None, clock=lambda: float(next(clock)),
    )


def test_client_submit_poll_fetch(fake_http, blast_xml):
    statuses = iter(["Status=WAITING", "Status=READY\nThereAreHits=yes"])
    fake_http.route("POST", BLAST_URL, FakeResponse("RID = RID42\nRTOE = 5\n"))

    def get_handler(url, params=None, **_):
        if params.get("FORMAT_OBJECT") == "SearchInfo":
            return FakeResponse(next(statuses))
        return FakeResponse(blast_xml)

    fake_http.route("GET", BLAST_URL, get_handler)
    result = _client(fake_http).search(PROTEIN, database="swissprot", evalue=1e-5, max_hits=10)

    assert result.rid == "RID42"
    assert len(result.hits) == 5
    method, _, kwargs = fake_http.calls[0]
    assert method == "POST"
    assert kwargs["data"]["PROGRAM"] == "blastp"  # never BLASTX for protein input
    assert kwargs["data"]["CMD"] == "Put"
    assert kwargs["data"]["FORMAT_TYPE"] == "XML"
    # NCBI usage guidelines: identify the client with the lower-case "email" and "tool" parameters.
    assert kwargs["data"]["email"] == "user@example.org"
    assert kwargs["data"]["tool"] == "test-tool"
    status_params = next(c[2]["params"] for c in fake_http.calls if c[0] == "GET")
    assert status_params["tool"] == "test-tool" and status_params["email"] == "user@example.org"


def test_client_reports_failed_job(fake_http):
    fake_http.route("POST", BLAST_URL, FakeResponse("RID = RIDX\nRTOE = 1\n"))
    fake_http.route("GET", BLAST_URL, FakeResponse("Status=FAILED"))
    with pytest.raises(BlastError, match="failed"):
        _client(fake_http).search(PROTEIN)


def test_client_times_out(fake_http):
    fake_http.route("POST", BLAST_URL, FakeResponse("RID = RIDT\nRTOE = 1\n"))
    fake_http.route("GET", BLAST_URL, FakeResponse("Status=WAITING"))
    with pytest.raises(BlastTimeoutError):
        _client(fake_http, timeout=180.0).search(PROTEIN)


def test_client_uses_cache(fake_http, blast_xml, tmp_path):
    cache = ResultCache(tmp_path)
    fake_http.route("POST", BLAST_URL, FakeResponse("RID = RIDC\nRTOE = 1\n"))

    def get_handler(url, params=None, **_):
        if params.get("FORMAT_OBJECT") == "SearchInfo":
            return FakeResponse("Status=READY\nThereAreHits=yes")
        return FakeResponse(blast_xml)

    fake_http.route("GET", BLAST_URL, get_handler)
    first = _client(fake_http, cache=cache).search(PROTEIN)
    n_calls = len(fake_http.calls)
    second = _client(fake_http, cache=cache).search(PROTEIN)
    assert not first.from_cache and second.from_cache
    assert len(fake_http.calls) == n_calls  # no new requests
    assert [h.accession for h in second.hits] == [h.accession for h in first.hits]
