# protein-function-annotator

Evidence-integrating functional annotation of protein sequences using **NCBI BLASTP**,
**EMBL-EBI InterProScan** and the **UniProt REST API**.

Given a protein FASTA sequence, the pipeline searches for homologs, identifies conserved domains,
retrieves curated annotations of the homologs, and then integrates these independent lines of
evidence into a conservative, explicitly justified functional annotation — including alternative
interpretations, contradictions, a documented confidence category and a list of limitations.

It answers questions such as:

* What is this protein likely to do, and how strong is the evidence?
* Which proteins are homologous to it, and are they experimentally characterized?
* Which conserved domains and families does it contain?
* Which GO terms are supported, and by what evidence?
* Is there evidence for enzymatic activity (EC numbers, catalytic domains, active sites)?
* Do the sources disagree, and what are the competing interpretations?

No local BLAST, InterProScan, HMMER or Pfam installation is needed: the computationally expensive
searches run on the public services. Results are never fabricated — if a service is unavailable the
report says so and the confidence category is reduced accordingly.

---

## Installation

Requires Python 3.11 or newer.

```bash
git clone https://github.com/<username>/protein-function-annotator.git
cd protein-function-annotator
python -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -e .                     # or: pip install -e ".[dev]" for tests and linting
```

## Configuration

No API keys are required. InterProScan at EMBL-EBI **requires a contact e-mail address**; the same
address is sent to NCBI with BLAST requests, as NCBI recommends.

```bash
cp .env.example .env
# edit .env and set ANNOTATOR_EMAIL=you@your-institution.org
```

or pass `--email you@your-institution.org` on the command line. Without an e-mail, the InterPro
stage is skipped (the pipeline still runs, with reduced confidence).

| Variable | Default | Purpose |
|---|---|---|
| `ANNOTATOR_EMAIL` | – | Contact e-mail for EBI (required for InterProScan) and NCBI |
| `BLAST_DATABASE` | `swissprot` | NCBI protein database (`swissprot`, `refseq_protein`, `nr`, `pdb`, …) |
| `BLAST_EVALUE` | `1e-5` | BLASTP E-value cut-off |
| `BLAST_MAX_HITS` | `25` | Maximum hits retrieved |
| `BLAST_TIMEOUT` / `INTERPRO_TIMEOUT` | `1800` | Maximum seconds to wait for a remote job |
| `BLAST_POLL_INTERVAL` | `60` | Seconds between status checks (NCBI minimum 60) |
| `INTERPRO_POLL_INTERVAL` | `15` | Seconds between InterProScan status checks |
| `INTERPRO_APPLICATIONS` | all | Comma-separated member databases (valid names: `iprscan5/parameterdetails/appl`, see docs/api.md) |
| `UNIPROT_TIMEOUT` / `UNIPROT_MAX_ENTRIES` | `60` / `15` | UniProt request timeout / entries fetched |
| `HTTP_TIMEOUT` / `HTTP_MAX_RETRIES` | `60` / `5` | Per-request timeout / retry budget |
| `ANNOTATOR_CACHE_DIR` | `~/.cache/protein-annotator` | Result cache location |
| `ANNOTATOR_USE_CACHE` / `ANNOTATOR_CACHE_TTL_DAYS` | `true` / `30` | Cache switch / expiry |

`.env` is listed in `.gitignore`; never commit it.

## Usage

```bash
protein-annotator --input examples/example_protein.fasta --output results/ --email you@example.org
```

Equivalent without installing the entry point:

```bash
python scripts/run_annotation.py --input examples/example_protein.fasta
```

Options:

```
--input, -i        Protein FASTA file (.fa, .fasta, .faa)                    [required]
--output, -o       Output directory (default: results/)
--email            Contact e-mail (overrides ANNOTATOR_EMAIL)
--blast-database   NCBI database (default: swissprot)
--evalue           BLASTP E-value cut-off (default: 1e-5)
--max-hits         Maximum BLASTP hits (default: 25)
--timeout          Maximum seconds per remote job
--skip-blast / --skip-interpro / --skip-uniprot
--json             Also print the annotation JSON to stdout
--no-cache         Do not read or write the result cache
--cache-dir DIR    Cache directory
--env-file FILE    Load configuration from a specific .env file
--verbose, -v      Debug logging
--version
```

Example progress output:

```
Protein Function Annotator
==========================

✓ FASTA validated (1 record(s))
✓ Sequence validation completed (142 residues, SHA-256 14725a105989…)
✓ Sequence properties completed (MW 15.3 kDa, pI 8.72)
⏳ Running BLASTP...
✓ BLASTP completed (25 hit(s) against swissprot)
⏳ Running InterProScan...
✓ InterPro completed (9 signature match(es), 4 InterPro entr(ies))
⏳ Retrieving UniProt annotations...
✓ UniProt completed (15 entr(ies) retrieved (15 reviewed))
⏳ Integrating evidence...
✓ Evidence integration completed (2 functional candidate(s))
⏳ Checking for contradictions...
✓ Contradiction analysis completed (0 contradiction(s), 0 major)
⏳ Generating functional interpretation...
✓ Functional inference completed (strong evidence)

Likely function: Hemoglobin subunit alpha
Confidence: strong evidence

Results written to:
results/
```

(Counts above are illustrative; actual numbers depend on current database contents. `--verbose` additionally shows job IDs and polling status.)

Remote jobs usually take a few minutes: BLASTP against `swissprot` typically 1–5 min (longer
against `nr`), InterProScan 1–5 min. Re-running the same sequence with the same parameters uses the
cache and completes in seconds. Multi-record FASTA files are processed sequentially, with one
sub-directory per record.

Exit codes: `0` success, `1` pipeline error, `2` invalid input/configuration, `130` interrupted.

### Library use

```python
from protein_annotator.config import AnnotatorConfig
from protein_annotator.io.fasta import read_fasta
from protein_annotator.io.reports import write_reports
from protein_annotator.pipeline import AnnotationPipeline

config = AnnotatorConfig.from_env()
with AnnotationPipeline(config) as pipeline:
    result = pipeline.run(read_fasta("examples/example_protein.fasta")[0])
print(result.annotation.predicted_function, result.annotation.confidence_category.value)
write_reports(result, "results/")
```

## Outputs

| File | Content |
|---|---|
| `annotation.json` | `predicted_function`, `interpretation`, `supporting_evidence`, typed `evidence` items, `alternative_functions`, `confidence_category`, `confidence_rationale`, `contradictions`, `enzyme_assessment`, `supported_go_terms`, `limitations` |
| `annotation.md` | Human-readable scientific report (sequence, homology, domains, UniProt evidence, GO, interpretation, limitations) |
| `blast_hits.csv` | Rank, accession, UniProt accession, description, organism, identity, positives, alignment length, query coverage, E-value, bit score, query/subject coordinates |
| `domains.csv` | Member-database signature, InterPro entry and type, coordinates, scores, GO terms, sites |
| `go_terms.csv` | Supported GO terms with sources, evidence codes and homolog counts |
| `sequence_statistics.json` | Sequence and physicochemical properties |
| `complete_results.json` | Everything above plus raw parsed service results, stage statuses and run metadata |

## How it works

```
FASTA → validation → sequence properties → BLASTP → InterProScan → UniProt enrichment
      → evidence integration → contradiction analysis → functional inference → reports
```

1. **Validation** — robust FASTA parsing; lower-case and whitespace normalised; invalid characters
   and empty sequences rejected; nucleotide sequences detected and rejected with the message
   *"This pipeline expects an amino-acid protein sequence. The supplied sequence appears to contain
   primarily nucleotide characters."* BLASTX is never used on protein input.
2. **Sequence properties** — Biopython ProtParam (MW, pI, GRAVY, instability index, composition…).
3. **BLASTP** — NCBI BLAST URL API: submit, poll (≤ once per minute), fetch XML, parse.
4. **InterProScan** — EMBL-EBI REST: submit, poll, fetch JSON; InterPro entries, Pfam and other
   member signatures, coordinates, GO terms, pathways, sites.
5. **UniProt** — curated annotations for significant hits (non-UniProt accessions mapped via ID
   mapping); reviewed vs. unreviewed status made explicit.
6. **Evidence integration** — homologs clustered by functional description; each cluster weighted by
   identity, coverage, E-value, curation and experimental support; domains, GO terms and EC numbers
   checked against each candidate.
7. **Contradictions** — conflicting homolog annotations, missing catalytic domains, domain/homology
   mismatches, inconsistent EC numbers, top-hit vs. consensus disagreement.
8. **Inference** — confidence category from explicit criteria; conservative wording.

Details: [docs/architecture.md](docs/architecture.md), [docs/methodology.md](docs/methodology.md),
[docs/api.md](docs/api.md).

## Interpreting results

| Category | Meaning | Wording |
|---|---|---|
| **strong evidence** | ≥3 agreeing homologs, ≥40 % identity, ≥70 % coverage, characterized homologs, supporting domains, no major contradictions | *"The strongest computational evidence supports …"* |
| **moderate evidence** | ≥2 agreeing homologs (or one close full-length homolog), ≥30 % identity, ≥50 % coverage, no contradicting domains | *"Probable …"* |
| **limited evidence** | A candidate exists but weaker criteria fail, conflicting evidence, or domain-only annotation | *"Putative …"* — treat as a hypothesis |
| **insufficient evidence** | No informative homologs and no functional domains | *"Uncharacterized protein (function not confidently assigned)"* |

These are heuristic evidence tiers, **not probabilities**; the report lists which criterion passed
or failed. Evidence statements are tagged as *direct*, *homology*, *domain* or *computational*.
Read the `contradictions` and `limitations` sections before relying on an annotation.

## Example data

`examples/example_protein.fasta` is **human hemoglobin subunit alpha** (UniProtKB **P69905**,
gene *HBA1*, 142 residues), a thoroughly characterized protein chosen so that results can be
checked against expert knowledge. Its FASTA header deliberately does not name the protein.

Expected general outcome (not hard-coded anywhere in the pipeline): a globin / hemoglobin alpha
chain annotation, globin domain signatures (InterPro IPR000971, Pfam PF00042), GO terms related to
oxygen carrier activity / oxygen transport / heme binding, and **no** evidence of enzymatic activity.
Because the query is itself in Swiss-Prot, the top hit should be the same protein at 100 % identity;
the report then adds a conditionally worded *direct evidence* statement.

`tests/fixtures/` contains small hand-constructed fixtures (illustrative values, not real results)
plus unmodified genuine NCBI BLAST and UniProt responses used to check the parsers against the
real formats; provenance is documented in `tests/fixtures/README.md`.

## Development

```bash
pip install -e ".[dev]"
pytest              # all offline: services mocked with saved and genuine responses
ruff check src tests scripts
make example        # run the example against the live services
```

## Limitations

* Remote services impose queueing and rate limits; this tool is designed for individual proteins or
  small batches, not proteome-scale annotation (use local BLAST+/InterProScan for that).
* Annotation transfer by homology can propagate database errors and cannot detect loss of
  function, neofunctionalisation or moonlighting.
* Confidence thresholds are documented heuristics, not a calibrated statistical model.
* Service APIs may change; endpoints are configurable in `config.py` (`Endpoints`).

## Data sources and citation

Please cite the underlying resources when publishing results: BLAST (Altschul et al., 1997;
Camacho et al., 2009),
InterProScan (Jones et al., 2014) and InterPro (Paysan-Lafosse et al., 2023), UniProt (The UniProt
Consortium, 2023), Biopython (Cock et al., 2009). Use of the services is subject to the terms of
NCBI, EMBL-EBI and UniProt.

## License

MIT — see [LICENSE](LICENSE).
