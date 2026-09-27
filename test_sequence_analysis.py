import pytest

from protein_annotator.analysis.sequence_analysis import compute_properties
from protein_annotator.io.fasta import read_fasta
from protein_annotator.models.sequence import ProteinSequence


@pytest.fixture
def example(fixtures_dir):
    return read_fasta(fixtures_dir / "example_protein.fasta")[0]


def test_length_and_composition(example):
    props = compute_properties(example)
    assert props.length == 142
    assert sum(props.composition.values()) == 142
    assert props.composition["A"] == example.sequence.count("A")
    assert abs(sum(props.composition_percent.values()) - 100) < 0.5


def test_molecular_weight(example):
    # Average mass of the 142-residue chain is ~15.26 kDa.
    props = compute_properties(example)
    assert 15000 < props.molecular_weight < 15500


def test_isoelectric_point(example):
    props = compute_properties(example)
    assert 7.5 < props.isoelectric_point < 9.5


def test_simple_known_values():
    props = compute_properties(ProteinSequence(identifier="g", sequence="GGGG"))
    # Glycine residue 57.05 Da ×4 + water 18.02 Da
    assert props.molecular_weight == pytest.approx(246.22, abs=0.1)
    assert props.positive_residues == 0 and props.negative_residues == 0


def test_charged_counts():
    props = compute_properties(ProteinSequence(identifier="c", sequence="KRDEHAAA"))
    assert props.positive_residues == 2
    assert props.negative_residues == 2
    assert props.histidine_residues == 1


def test_non_standard_residues_noted():
    props = compute_properties(ProteinSequence(identifier="x", sequence="MKTAYXXAKQ"))
    assert props.non_standard_residues == {"X": 2}
    assert props.analysed_length == 8
    assert any("Non-standard" in n for n in props.notes)
