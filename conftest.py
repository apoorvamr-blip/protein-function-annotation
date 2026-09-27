"""Shared test helpers. No test performs real network access."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from protein_annotator.models.blast import BlastHit
from protein_annotator.models.domain import DomainHit, GOTerm, InterProResult
from protein_annotator.models.uniprot import UniProtEntry

FIXTURES = Path(__file__).parent / "fixtures"


class FakeResponse:
    def __init__(self, text: str = "", status_code: int = 200, payload: Any = None):
        self.text = text if payload is None else json.dumps(payload)
        self.status_code = status_code
        self._payload = payload
        self.headers: dict[str, str] = {}

    def json(self) -> Any:
        return self._payload if self._payload is not None else json.loads(self.text)


class FakeHttp:
    """Routes requests to handler callables by (method, url-substring)."""

    def __init__(self) -> None:
        self.routes: list[tuple[str, str, Callable[..., FakeResponse]]] = []
        self.calls: list[tuple[str, str, dict[str, Any]]] = []

    def route(self, method: str, fragment: str, handler: Callable[..., FakeResponse] | FakeResponse) -> None:
        func = handler if callable(handler) else (lambda **_: handler)
        self.routes.append((method, fragment, func))

    def _dispatch(self, method: str, url: str, **kwargs: Any) -> FakeResponse:
        self.calls.append((method, url, kwargs))
        for route_method, fragment, handler in self.routes:
            if route_method == method and fragment in url:
                return handler(url=url, **kwargs)
        raise AssertionError(f"Unexpected {method} {url}")

    def get(self, url: str, **kwargs: Any) -> FakeResponse:
        return self._dispatch("GET", url, **kwargs)

    def post(self, url: str, **kwargs: Any) -> FakeResponse:
        return self._dispatch("POST", url, **kwargs)

    def get_json(self, url: str, **kwargs: Any) -> Any:
        return self.get(url, **kwargs).json()


def make_hit(accession: str, description: str, *, identity: float = 90.0, coverage: float = 1.0,
             evalue: float = 1e-50, bitscore: float = 200.0, rank: int = 1,
             uniprot: str | None = None, organism: str | None = None) -> BlastHit:
    return BlastHit(
        rank=rank, hit_id=f"sp|{accession}.1|X", accession=accession, sequence_id=f"{accession}.1",
        description=description, organism=organism, identity=identity, alignment_length=100,
        query_coverage=coverage, evalue=evalue, bitscore=bitscore, query_start=1, query_end=100,
        subject_start=1, subject_end=100, uniprot_accession=uniprot if uniprot is not None else accession,
    )


def make_entry(accession: str, name: str, *, reviewed: bool = True, experimental: bool = False,
               ec: list[str] | None = None, interpro: list[str] | None = None,
               pfam: list[str] | None = None, go: list[tuple[str, str, str, str]] | None = None,
               function: str | None = None, organism: str = "Test organism") -> UniProtEntry:
    return UniProtEntry(
        accession=accession, reviewed=reviewed, protein_name=name, organism=organism,
        function=[function] if function else [], function_experimental=experimental,
        ec_numbers=ec or [], interpro_ids=interpro or [], pfam_ids=pfam or [],
        go_terms=[GOTerm(id=i, name=n, category=c, source=f"UniProtKB:{accession}", evidence=e)
                  for i, n, c, e in (go or [])],
    )


def make_interpro(domains: list[tuple[str, str, str, str, str]]) -> InterProResult:
    """domains: (database, signature, interpro_accession, name, type)."""
    return InterProResult(
        job_id="job", domains=[
            DomainHit(database=db, accession=sig, name=name, start=1, end=100,
                      interpro_accession=ipr or None, interpro_name=name if ipr else None,
                      interpro_type=typ or None)
            for db, sig, ipr, name, typ in domains
        ],
    )


@pytest.fixture
def fixtures_dir() -> Path:
    return FIXTURES


@pytest.fixture
def fake_http() -> FakeHttp:
    return FakeHttp()
