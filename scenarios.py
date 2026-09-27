"""Reusable evidence scenarios for integration/inference tests (hand-constructed, not real results)."""

from __future__ import annotations

from protein_annotator.config import InferenceSettings
from protein_annotator.integration.evidence import build_evidence_summary
from protein_annotator.models.blast import BlastResult

from .conftest import make_entry, make_hit, make_interpro

SETTINGS = InferenceSettings()

GLOBIN_DOMAINS = [
    ("PFAM", "PF00042", "IPR000971", "Globin", "DOMAIN"),
    ("PROSITE_PROFILES", "PS01033", "IPR000971", "Globin", "DOMAIN"),
]
GLOBIN_GO = [
    ("GO:0005344", "oxygen carrier activity", "molecular_function", "IDA"),
    ("GO:0015671", "oxygen transport", "biological_process", "IDA"),
]


def blast(hits):
    return BlastResult(database="swissprot", query_length=142, hits=hits)


def globin_hits_and_entries(n=3, identity=95.0):
    accessions = ["P69905", "P69907", "P01942", "P01946"][:n]
    hits = [make_hit(a, f"RecName: Full=Hemoglobin subunit alpha [Organism {i}]", identity=identity - i,
                     bitscore=300 - 10 * i, rank=i + 1) for i, a in enumerate(accessions)]
    entries = {a: make_entry(a, "Hemoglobin subunit alpha", experimental=(i == 0),
                             interpro=["IPR000971", "IPR002338"], pfam=["PF00042"], go=GLOBIN_GO)
               for i, a in enumerate(accessions)}
    return hits, entries


def strong_summary():
    hits, entries = globin_hits_and_entries()
    interpro = make_interpro(GLOBIN_DOMAINS)
    return build_evidence_summary(blast(hits), interpro, entries, SETTINGS, uniprot_available=True), interpro


def weak_summary():
    hits = [make_hit("Q9XYZ1", "Hemoglobin-like protein [Some organism]", identity=32.0, coverage=0.6,
                     evalue=1e-8, bitscore=60)]
    return build_evidence_summary(blast(hits), None, {}, SETTINGS, uniprot_available=False)


def conflicting_summary():
    hits = [
        make_hit("P11111", "Hemoglobin subunit alpha [A]", identity=45, coverage=0.9, evalue=1e-30, bitscore=120),
        make_hit("P22222", "Hemoglobin subunit alpha [B]", identity=44, coverage=0.9, evalue=1e-29, bitscore=118),
        make_hit("P33333", "DNA gyrase subunit B [C]", identity=46, coverage=0.9, evalue=1e-31, bitscore=125),
        make_hit("P44444", "DNA gyrase subunit B [D]", identity=45, coverage=0.9, evalue=1e-30, bitscore=121),
    ]
    return build_evidence_summary(blast(hits), None, {}, SETTINGS, uniprot_available=False)


def empty_summary():
    return build_evidence_summary(None, None, {}, SETTINGS, uniprot_available=False)


def enzyme_without_domain_summary():
    accessions = ["P00760", "P00761", "P00762"]
    hits = [make_hit(a, "Trypsin [Organism]", identity=60 - i, bitscore=200 - i, rank=i + 1)
            for i, a in enumerate(accessions)]
    entries = {a: make_entry(a, "Serine protease trypsin", ec=["3.4.21.4"], interpro=["IPR001254"],
                             pfam=["PF00089"]) for a in accessions}
    interpro = make_interpro(GLOBIN_DOMAINS)  # query carries an unrelated domain only
    summary = build_evidence_summary(blast(hits), interpro, entries, SETTINGS, uniprot_available=True)
    return summary, interpro
