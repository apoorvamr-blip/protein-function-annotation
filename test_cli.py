"""CLI tests. Only offline modes are exercised (all remote services skipped)."""

from __future__ import annotations

import json

import pytest

from protein_annotator.cli import EXIT_INPUT_ERROR, EXIT_OK, main


@pytest.fixture(autouse=True)
def isolated_env(monkeypatch, tmp_path):
    for name in ("ANNOTATOR_EMAIL", "ANNOTATOR_CACHE_DIR"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.chdir(tmp_path)


def test_help(capsys):
    with pytest.raises(SystemExit) as info:
        main(["--help"])
    assert info.value.code == 0
    out = capsys.readouterr().out
    for option in ("--input", "--output", "--email", "--blast-database", "--evalue", "--max-hits",
                   "--timeout", "--skip-blast", "--skip-interpro", "--skip-uniprot", "--json", "--verbose"):
        assert option in out


def test_nucleotide_input_rejected(tmp_path, capsys):
    fasta = tmp_path / "dna.fasta"
    fasta.write_text(">dna\nATGGCTAGCTAGCTAGCTAGGCTAGCTAGCATCGATCGATCGA\n")
    assert main(["--input", str(fasta), "--no-cache"]) == EXIT_INPUT_ERROR
    assert "nucleotide" in capsys.readouterr().err


def test_missing_file(tmp_path):
    assert main(["--input", str(tmp_path / "missing.fasta"), "--no-cache"]) == EXIT_INPUT_ERROR


def test_offline_run_writes_reports(tmp_path, fixtures_dir, capsys):
    out = tmp_path / "results"
    code = main(["--input", str(fixtures_dir / "example_protein.fasta"), "--output", str(out),
                 "--skip-blast", "--skip-interpro", "--skip-uniprot", "--no-cache", "--json"])
    assert code == EXIT_OK
    assert (out / "annotation.md").exists() and (out / "complete_results.json").exists()
    payload = json.loads(capsys.readouterr().out)
    assert payload["confidence_category"] == "insufficient evidence"


def test_multi_record_input_uses_subdirectories(tmp_path):
    fasta = tmp_path / "two.faa"
    fasta.write_text(">p1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ\n>p2\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF\n")
    out = tmp_path / "res"
    assert main(["--input", str(fasta), "--output", str(out), "--skip-blast", "--skip-interpro",
                 "--skip-uniprot", "--no-cache"]) == EXIT_OK
    assert (out / "p1" / "annotation.json").exists() and (out / "p2" / "annotation.json").exists()
