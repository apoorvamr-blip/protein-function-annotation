"""Result-cache tests."""

from __future__ import annotations

from protein_annotator.utils.caching import ResultCache, sequence_hash, stable_hash


def test_sequence_hash_normalizes():
    assert sequence_hash("mkt ay\n") == sequence_hash("MKTAY")
    assert len(sequence_hash("MKTAY")) == 64


def test_stable_hash_is_order_independent():
    assert stable_hash({"a": 1, "b": 2}) == stable_hash({"b": 2, "a": 1})
    assert stable_hash({"a": 1}) != stable_hash({"a": 2})


def test_roundtrip_and_disabled(tmp_path):
    cache = ResultCache(tmp_path)
    key = cache.key("blast", sequence_hash="x", database="swissprot")
    assert cache.get("blast", key) is None
    cache.set("blast", key, "<BlastOutput/>")
    assert cache.get("blast", key) == "<BlastOutput/>"
    assert not list(tmp_path.rglob("*.tmp"))

    disabled = ResultCache(tmp_path, enabled=False)
    assert disabled.get("blast", key) is None


def test_ttl_expiry(tmp_path):
    now = [1000.0]
    cache = ResultCache(tmp_path, ttl_seconds=10, clock=lambda: now[0])
    cache.set("uniprot", "k", "payload")
    assert cache.get("uniprot", "k") == "payload"
    now[0] += 11
    assert cache.get("uniprot", "k") is None


def test_corrupt_entry_is_discarded(tmp_path):
    cache = ResultCache(tmp_path)
    cache.set("interpro", "k", "{}")
    path = next(tmp_path.rglob("k.json"))
    path.write_text("not json")
    assert cache.get("interpro", "k") is None
    assert not path.exists()
