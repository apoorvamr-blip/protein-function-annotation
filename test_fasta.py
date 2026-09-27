import pytest

from protein_annotator.io.fasta import (
    FastaValidationError,
    NucleotideSequenceError,
    parse_fasta_text,
    read_fasta,
    safe_identifier,
    validate_protein_sequence,
)


def test_valid_single_record():
    records = parse_fasta_text(">prot1 some description\nMKTAYIAKQRQISFVKSHFSRQ\n")
    assert len(records) == 1
    assert records[0].identifier == "prot1"
    assert records[0].description == "some description"
    assert records[0].sequence == "MKTAYIAKQRQISFVKSHFSRQ"


def test_multiline_and_whitespace():
    records = parse_fasta_text(">p\nMKTAY IAKQ\n  RQISF\tVKSHF\n\nSRQ\n")
    assert records[0].sequence == "MKTAYIAKQRQISFVKSHFSRQ"


def test_lowercase_is_uppercased_and_stop_removed():
    records = parse_fasta_text(">p\nmktayiakqrqisfvkshfsrq*\n")
    assert records[0].sequence == "MKTAYIAKQRQISFVKSHFSRQ"


def test_multiple_records():
    records = parse_fasta_text(">a\nMKTAYIAKQR\n>b\nGGHHKKLLMM\n")
    assert [r.identifier for r in records] == ["a", "b"]


def test_invalid_characters_rejected():
    with pytest.raises(FastaValidationError, match="invalid amino-acid"):
        parse_fasta_text(">p\nMKT1AY#IAK\n")


def test_empty_sequence_rejected():
    with pytest.raises(FastaValidationError, match="empty sequence"):
        parse_fasta_text(">p\n\n>q\nMKT\n")


def test_missing_header_rejected():
    with pytest.raises(FastaValidationError, match="before a FASTA header"):
        parse_fasta_text("MKTAYIAKQR\n")


def test_no_records():
    with pytest.raises(FastaValidationError, match="No FASTA records"):
        parse_fasta_text("\n\n")


def test_duplicate_identifiers():
    with pytest.raises(FastaValidationError, match="Duplicate"):
        parse_fasta_text(">a\nMKT\n>a\nMKT\n")


def test_nucleotide_detected():
    with pytest.raises(NucleotideSequenceError) as info:
        parse_fasta_text(">dna\nATGGCTAGCTAGGATCCATGCATGCAAGCTT\n")
    assert "expects an amino-acid protein sequence" in str(info.value)
    assert "primarily nucleotide characters" in str(info.value)


def test_ambiguous_residues_allowed():
    assert validate_protein_sequence("MKXBZ") == "MKXBZ"


def test_read_fasta_extensions(tmp_path):
    for ext in (".fa", ".fasta", ".faa"):
        path = tmp_path / f"x{ext}"
        path.write_text(">p\nMKTAYIAKQR\n")
        assert read_fasta(path)[0].sequence == "MKTAYIAKQR"
    bad = tmp_path / "x.txt"
    bad.write_text(">p\nMKT\n")
    with pytest.raises(FastaValidationError, match="Unsupported file extension"):
        read_fasta(bad)


def test_read_fasta_missing_file(tmp_path):
    with pytest.raises(FastaValidationError, match="not found"):
        read_fasta(tmp_path / "missing.fasta")


def test_example_file(fixtures_dir):
    records = read_fasta(fixtures_dir / "example_protein.fasta")
    assert records[0].length == 142


def test_sequence_hash_is_normalised():
    a = parse_fasta_text(">a\nmktay\n")[0]
    b = parse_fasta_text(">b\nMKT AY\n")[0]
    assert a.sha256 == b.sha256


def test_safe_identifier():
    assert safe_identifier("sp|P69905|HBA_HUMAN") == "sp_P69905_HBA_HUMAN"
