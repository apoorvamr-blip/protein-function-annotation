import io
import os
from pathlib import Path

import streamlit as st

from protein_annotator.config import AnnotatorConfig
from protein_annotator.io.fasta import FastaValidationError, parse_fasta_text
from protein_annotator.io.reports import render_markdown
from protein_annotator.models.annotation import StageState
from protein_annotator.pipeline.runner import AnnotationPipeline
from protein_annotator.pipeline.stages import PipelineOptions


st.set_page_config(
    page_title="Protein Function Annotator",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 Protein Function Annotator")
st.caption(
    "Evidence-integrating functional annotation using NCBI BLASTP, "
    "EMBL-EBI InterProScan, and UniProt."
)

with st.sidebar:
    st.header("Settings")
    email_default = os.getenv("ANNOTATOR_EMAIL", "")
    email = st.text_input(
        "Contact e-mail",
        value=email_default,
        help="Used as contact information for NCBI/EMBL-EBI requests. "
             "Set ANNOTATOR_EMAIL in Streamlit Secrets to avoid entering it each run.",
    )

    blast_database = st.selectbox(
        "BLAST database",
        ["swissprot", "nr", "refseq_protein", "pdb"],
        index=0,
    )
    max_hits = st.slider("Maximum BLAST hits", 5, 50, 25, 5)
    evalue = st.number_input(
        "BLAST E-value cutoff",
        min_value=1e-100,
        max_value=1.0,
        value=1e-5,
        format="%.1e",
    )

    st.subheader("Pipeline stages")
    skip_blast = st.checkbox("Skip BLASTP", value=False)
    skip_interpro = st.checkbox("Skip InterProScan", value=False)
    skip_uniprot = st.checkbox("Skip UniProt enrichment", value=False)

    st.divider()
    st.caption(
        "Remote services can take several minutes and may enforce rate limits. "
        "For repeated sequences, the built-in cache can reduce repeated requests."
    )

uploaded = st.file_uploader(
    "Upload a protein FASTA file",
    type=["fa", "fasta", "faa", "txt"],
    help="You can also paste a FASTA sequence below.",
)

pasted = st.text_area(
    "Or paste FASTA",
    height=180,
    placeholder=">my_protein\nMKTIIALSYIFCLVFADYKDDDDK...",
)

run = st.button("🚀 Annotate protein", type="primary", use_container_width=True)

if run:
    if uploaded is not None:
        fasta_text = uploaded.getvalue().decode("utf-8-sig")
    elif pasted.strip():
        fasta_text = pasted
    else:
        st.error("Please upload a FASTA file or paste a FASTA sequence.")
        st.stop()

    try:
        records = parse_fasta_text(fasta_text)
    except (UnicodeDecodeError, FastaValidationError) as exc:
        st.error(f"Invalid FASTA: {exc}")
        st.stop()

    if not email.strip() and not skip_interpro:
        st.warning(
            "No contact e-mail is configured, so InterProScan will be skipped. "
            "Enter an e-mail or set ANNOTATOR_EMAIL in Streamlit Secrets."
        )

    config = AnnotatorConfig.from_env()
    config = config.with_overrides(
        email=email.strip() or config.email,
        blast=config.blast.__class__(
            database=blast_database,
            evalue=evalue,
            max_hits=max_hits,
            timeout=config.blast.timeout,
            poll_interval=config.blast.poll_interval,
        ),
    )
    config.validate()

    options = PipelineOptions(
        skip_blast=skip_blast,
        skip_interpro=skip_interpro,
        skip_uniprot=skip_uniprot,
    )

    progress_box = st.empty()
    status_box = st.status("Running annotation pipeline…", expanded=True)

    class StreamlitObserver:
        def stage_started(self, name: str) -> None:
            progress_box.info(f"Running **{name}**…")

        def stage_finished(self, status) -> None:
            icon = {
                StageState.SUCCESS: "✅",
                StageState.WARNING: "⚠️",
                StageState.FAILED: "❌",
                StageState.SKIPPED: "⏭️",
            }.get(status.state, "•")
            status_box.write(f"{icon} **{status.name}** — {status.message}")

        def progress(self, message: str) -> None:
            progress_box.caption(message)

    results = []
    try:
        with AnnotationPipeline(
            config,
            options,
            observer=StreamlitObserver(),
        ) as pipeline:
            for record in records:
                status_box.write(f"Processing **{record.identifier}** ({record.length} aa)…")
                results.append(pipeline.run(record))
        status_box.update(label="Annotation complete", state="complete", expanded=False)
        progress_box.empty()
    except Exception as exc:
        status_box.update(label="Pipeline failed", state="error", expanded=True)
        st.exception(exc)
        st.stop()

    for result in results:
        st.divider()
        st.subheader(f"Result: {result.sequence.identifier}")

        ann = result.annotation
        if ann:
            c1, c2, c3 = st.columns(3)
            c1.metric("Length", f"{result.sequence.length} aa")
            c2.metric("Confidence", ann.confidence_category.value)
            c3.metric(
                "BLAST hits",
                len(result.blast.hits) if result.blast else 0,
            )

            st.markdown("### Predicted function")
            st.success(ann.predicted_function)
            st.write(ann.interpretation)

            if ann.supporting_evidence:
                with st.expander("Supporting evidence", expanded=True):
                    for item in ann.supporting_evidence:
                        st.write(f"• {item}")

            if ann.confidence_rationale:
                with st.expander("Why this confidence category"):
                    for item in ann.confidence_rationale:
                        st.write(f"• {item}")

            if ann.contradictions:
                with st.expander("Contradictions"):
                    for contradiction in ann.contradictions:
                        st.warning(
                            f"**{contradiction.severity.upper()} — "
                            f"{contradiction.kind}**: {contradiction.description}"
                        )
            else:
                st.info("No contradictions detected.")

            if ann.limitations:
                with st.expander("Limitations"):
                    for item in ann.limitations:
                        st.write(f"• {item}")

        if result.properties:
            with st.expander("Sequence properties"):
                p = result.properties
                st.json(p.model_dump(mode="json"))

        if result.blast and result.blast.hits:
            with st.expander("BLASTP hits"):
                rows = [
                    {
                        "Rank": h.rank,
                        "Accession": h.accession,
                        "Description": h.description,
                        "Organism": h.organism,
                        "Identity": h.identity,
                        "Coverage": h.query_coverage,
                        "E-value": h.evalue,
                        "Bitscore": h.bitscore,
                    }
                    for h in result.blast.hits
                ]
                st.dataframe(rows, use_container_width=True, hide_index=True)

        if result.interpro and result.interpro.domains:
            with st.expander("InterPro domains"):
                rows = [
                    {
                        "Database": d.database,
                        "Accession": d.accession,
                        "Name": d.name or d.description,
                        "Start": d.start,
                        "End": d.end,
                        "InterPro": d.interpro_accession,
                    }
                    for d in result.interpro.domains
                ]
                st.dataframe(rows, use_container_width=True, hide_index=True)

        if result.uniprot_entries:
            with st.expander("UniProt evidence"):
                rows = [
                    {
                        "Accession": e.accession,
                        "Protein": e.protein_name,
                        "Organism": e.organism,
                        "Reviewed": e.reviewed,
                        "EC": ", ".join(e.ec_numbers),
                    }
                    for e in result.uniprot_entries
                ]
                st.dataframe(rows, use_container_width=True, hide_index=True)

        report_md = render_markdown(result)
        st.download_button(
            "⬇️ Download Markdown report",
            data=report_md,
            file_name=f"{result.sequence.identifier}_annotation.md",
            mime="text/markdown",
        )

        complete_json = result.model_dump_json(indent=2)
        st.download_button(
            "⬇️ Download complete JSON",
            data=complete_json,
            file_name=f"{result.sequence.identifier}_results.json",
            mime="application/json",
        )
