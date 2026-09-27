"""Command-line interface: ``protein-annotator``."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

from rich.console import Console

from protein_annotator import __version__
from protein_annotator.config import NCBI_MIN_POLL_INTERVAL, AnnotatorConfig, ConfigurationError
from protein_annotator.io.fasta import FastaValidationError, read_fasta, safe_identifier
from protein_annotator.io.reports import write_reports
from protein_annotator.models.annotation import StageState, StageStatus
from protein_annotator.pipeline.runner import AnnotationPipeline, PipelineError
from protein_annotator.pipeline.stages import PipelineOptions
from protein_annotator.utils.logging import configure_logging

EXIT_OK = 0
EXIT_PIPELINE_ERROR = 1
EXIT_INPUT_ERROR = 2
EXIT_INTERRUPTED = 130

RUNNING_LABELS = {
    "BLASTP": "Running BLASTP",
    "InterPro": "Running InterProScan",
    "UniProt": "Retrieving UniProt annotations",
    "Evidence integration": "Integrating evidence",
    "Contradiction analysis": "Checking for contradictions",
    "Functional inference": "Generating functional interpretation",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="protein-annotator",
        description=(
            "Evidence-integrating functional annotation of a protein sequence using NCBI BLASTP, "
            "EMBL-EBI InterProScan and UniProt."
        ),
    )
    parser.add_argument("-i", "--input", required=True, type=Path,
                        help="Protein FASTA file (.fa, .fasta, .faa)")
    parser.add_argument("-o", "--output", type=Path, default=Path("results"),
                        help="Output directory (default: results/)")
    parser.add_argument("--email", help="Contact e-mail for NCBI/EBI (overrides ANNOTATOR_EMAIL)")
    parser.add_argument("--blast-database",
                        help="NCBI protein database (default: swissprot; e.g. nr, refseq_protein, pdb)")
    parser.add_argument("--evalue", type=float, help="BLASTP E-value cutoff (default: 1e-5)")
    parser.add_argument("--max-hits", type=int, help="Maximum BLASTP hits to retrieve (default: 25)")
    parser.add_argument("--timeout", type=float,
                        help="Maximum seconds to wait for each remote job (BLAST, InterProScan)")
    parser.add_argument("--skip-blast", action="store_true", help="Do not run BLASTP")
    parser.add_argument("--skip-interpro", action="store_true", help="Do not run InterProScan")
    parser.add_argument("--skip-uniprot", action="store_true", help="Do not query UniProt")
    parser.add_argument("--json", action="store_true",
                        help="Also print the final annotation JSON to stdout")
    parser.add_argument("--no-cache", action="store_true", help="Disable the result cache")
    parser.add_argument("--cache-dir", type=Path, help="Cache directory (default: ~/.cache/protein-annotator)")
    parser.add_argument("--env-file", type=Path, help="Path to a .env file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def build_config(args: argparse.Namespace) -> AnnotatorConfig:
    config = AnnotatorConfig.from_env(args.env_file)
    blast = config.blast
    if args.blast_database:
        blast = replace(blast, database=args.blast_database)
    if args.evalue is not None:
        blast = replace(blast, evalue=args.evalue)
    if args.max_hits is not None:
        blast = replace(blast, max_hits=args.max_hits)
    interpro = config.interpro
    if args.timeout is not None:
        blast = replace(blast, timeout=args.timeout, poll_interval=max(blast.poll_interval, NCBI_MIN_POLL_INTERVAL))
        interpro = replace(interpro, timeout=args.timeout)
    inference = config.inference
    if args.evalue is not None:
        inference = replace(inference, significance_evalue=min(args.evalue, inference.significance_evalue))
    config = config.with_overrides(
        email=args.email or config.email,
        blast=blast,
        interpro=interpro,
        inference=inference,
        use_cache=config.use_cache and not args.no_cache,
        cache_dir=args.cache_dir or config.cache_dir,
    )
    config.validate()
    return config


class ConsoleObserver:
    """Render stage progress with Rich."""

    def __init__(self, console: Console, verbose: bool = False) -> None:
        self.console = console
        self.verbose = verbose

    def stage_started(self, name: str) -> None:
        label = RUNNING_LABELS.get(name)
        if label:
            self.console.print(f"⏳ {label}...")

    def stage_finished(self, status: StageStatus) -> None:
        detail = f" [dim]({status.message})[/dim]" if status.message else ""
        if status.state is StageState.SUCCESS:
            self.console.print(f"[green]✓[/green] {status.name} completed{detail}")
        elif status.state is StageState.WARNING:
            self.console.print(f"[yellow]⚠️[/yellow]  {status.name} completed with warnings{detail}")
        elif status.state is StageState.SKIPPED:
            self.console.print(f"[dim]– {status.name} skipped{detail}[/dim]")
        else:
            self.console.print(f"[red]✗[/red] {status.name} unavailable{detail}")

    def progress(self, message: str) -> None:
        if self.verbose:
            self.console.print(f"  [dim]{message}[/dim]")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    console = Console(stderr=True, highlight=False)
    configure_logging(args.verbose, console=console)

    try:
        config = build_config(args)
    except ConfigurationError as exc:
        console.print(f"[red]Configuration error:[/red] {exc}")
        return EXIT_INPUT_ERROR

    console.print("[bold]Protein Function Annotator[/bold]\n==========================\n")
    try:
        records = read_fasta(args.input)
    except FastaValidationError as exc:
        console.print(f"[red]✗ Invalid input:[/red] {exc}")
        return EXIT_INPUT_ERROR
    console.print(f"[green]✓[/green] FASTA validated ({len(records)} record(s))")
    if not config.email and not args.skip_interpro:
        console.print("[yellow]⚠️  No contact e-mail configured; InterProScan (EMBL-EBI) will be "
                      "skipped. Set ANNOTATOR_EMAIL or use --email.[/yellow]")

    options = PipelineOptions(
        skip_blast=args.skip_blast, skip_interpro=args.skip_interpro, skip_uniprot=args.skip_uniprot
    )
    annotations = []
    try:
        with AnnotationPipeline(config, options, observer=ConsoleObserver(console, args.verbose)) as pipeline:
            for record in records:
                if len(records) > 1:
                    console.print(f"\n[bold]» {record.identifier}[/bold]")
                result = pipeline.run(record)
                out_dir = args.output if len(records) == 1 else args.output / safe_identifier(record.identifier)
                write_reports(result, out_dir)
                if result.annotation:
                    annotations.append({"sequence_id": record.identifier,
                                        **result.annotation.model_dump(mode="json")})
                    console.print(
                        f"\n[bold]Likely function:[/bold] {result.annotation.predicted_function}\n"
                        f"[bold]Confidence:[/bold] {result.annotation.confidence_category.value}"
                    )
    except PipelineError as exc:
        console.print(f"[red]✗ Pipeline error:[/red] {exc}")
        return EXIT_PIPELINE_ERROR
    except KeyboardInterrupt:
        console.print("\n[red]Interrupted.[/red]")
        return EXIT_INTERRUPTED

    console.print(f"\nResults written to:\n{args.output}/")
    if args.json:
        payload = annotations[0] if len(annotations) == 1 else annotations
        sys.stdout.write(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    return EXIT_OK


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
