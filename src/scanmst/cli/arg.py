"""Parse command line arguments."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Any

from scanmst import __version__


@dataclass
class DefaultOptions:
    """Cli default options."""

    input: str
    ref: str
    gtf: str
    output: str
    blat_two_bit: str
    output_sequence_choices: tuple[str, str, str] = ("consensus", "reference", "both")
    output_sequence_choice: str = "consensus"
    blat_closed: bool = True
    blat_sleep: bool = True
    blat_port: int = 88888
    aligner: tuple[str, str] = ("blat", "")
    support_reads: int = 1
    splice_bin: int = 5
    mapq: int = 20
    noncanonical: bool = False
    bound: bool = True
    graph: bool = False
    log: str = "warning"
    species: str = "human"
    species_choices: tuple[str, str] = ("human", "mouse")
    thread: int = 1
    min_soft_seg_len: int = 200
    max_allowed_nm: int = 100
    max_allowed_micro_insertion: int = 50
    min_required_insertion_length: int = 100
    ident_cutoff: float = 0.90
    prune_threshold: int = 10  # for merging conditions
    soft_len: int = 5
    mismatch: int = 3
    alignment_fraction: float = 0.8
    long_indel_length: int = 10
    substitutions_fraction: float = 0.05
    indel_fraction: float = 0.001
    circular_rna: str = "remove"
    circular_rna_choices: tuple[str, ...] = ("remove", "keep", "extract")
    # junctions within one annotated exon filter
    exon_filter: bool = True
    rt_switching_filter_len: int = 10
    ignore_circle: bool = False
    rescue_sr: bool = False
    refine: bool = False


COLOR = "bold magenta"

BANNER = {
    "   _____                 __  ___ _____ ______": COLOR,
    "  / ___/_________ _____ /  |/  // ___//_  __/": COLOR,
    "  \\__ \\/ ___/ __ `/ __ \\/ /|_/ / \\__ \\ / /": COLOR,
    " ___/ / /__/ /_/ / / / / /  / / ___/ // /": COLOR,
    "/____/\\___/\\__,_/_/ /_/_/  /_/ /____//_/": COLOR,
}


def print_banner() -> None:
    """Print banner."""
    from rich.console import Console

    console = Console()
    for line, color in BANNER.items():
        console.print(line, style=color)


class RichArgParser(argparse.ArgumentParser):
    """RichArgParser."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """RichArgParser."""
        from rich.console import Console

        self.console = Console()
        super().__init__(*args, **kwargs)

    @staticmethod
    def _color_message(message: str, color: str = "green") -> str:
        """Color message."""
        import re

        pattern = re.compile(r"(?<!\w)(?P<arg>(?:--[\w-]+|-h))(?!\w)")
        return pattern.sub(lambda m: f"[bold {color}]{m.group('arg')}[/]", message)

    def _print_message(self, message: str | None, _file: Any = None) -> None:
        if message:
            self.console.print(self._color_message(message))


class RichHelpFormatter(argparse.HelpFormatter):
    """RichHelpFormatter."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """RichHelpFormatter."""
        super().__init__(*args, max_help_position=42, **kwargs)  # type: ignore


def parse_args() -> argparse.ArgumentParser:
    """Parse command line arguments."""
    parser = RichArgParser(
        description="[red]scanmst[/] :rocket: Multi-segment transcript (MST) identification using transcriptomic long reads data",
        formatter_class=RichHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    parser.add_argument(
        "--input",
        action="store",
        dest="input",
        help="Input alignment BAM file, which must contain both cs and SA tags.",
        required=True,
    )
    parser.add_argument(
        "--ref",
        action="store",
        dest="ref",
        help="Reference genome in FASTA format (with fai index)",
        required=True,
    )
    parser.add_argument(
        "--gtf",
        action="store",
        dest="gtf",
        help="Gene annotations in GTF format",
        required=True,
    )
    parser.add_argument(
        "--output",
        action="store",
        dest="output",
        help="Output file prefix",
        required=True,
    )
    parser.add_argument(
        "--output-seq",
        action="store",
        dest="output_sequence_choice",
        help="Output sequence type (default: %(default)s)",
        choices=DefaultOptions.output_sequence_choices,
        default=DefaultOptions.output_sequence_choice,
    )
    parser.add_argument(
        "--sr",
        action="store",
        dest="support_reads",
        type=int,
        help="The minimum number of supporting reads required for calling MST. (default: %(default)s)",
        default=DefaultOptions.support_reads,
    )
    parser.add_argument(
        "--splice-bin",
        action="store",
        dest="splice_bin",
        type=int,
        help="Bin size for searching canonical splice sites. (default: %(default)s)",
        default=DefaultOptions.splice_bin,
    )
    parser.add_argument(
        "--mapq",
        action="store",
        dest="mapq",
        type=int,
        help="Minimum MAPQ of reads required for calling MST. (default: %(default)s)",
        default=DefaultOptions.mapq,
    )

    parser.add_argument(
        "--log-level",
        action="store",
        dest="log",
        choices=["info", "debug", "trace", "warning"],  # "warning", "error", "critical"
        default=DefaultOptions.log,
        help="Set log level (default: %(default)s)",
    )
    parser.add_argument(
        "--thread",
        action="store",
        dest="thread",
        type=int,
        default=DefaultOptions.thread,
        help="Set the thread number (default: %(default)s)",
    )
    parser.add_argument(
        "--aligner",
        dest="aligner",
        type=str,
        choices=DefaultOptions.aligner,
        help="Aligner used for additional realignment to recover missing chimeric alignments. (default: %(default)s)",
        required=False,
    )
    parser.add_argument(
        "--blat-identity",
        action="store",
        dest="ident_cutoff",
        type=float,
        help="BLAT identity cutoff (default: %(default)s)",
        default=DefaultOptions.ident_cutoff,
    )
    parser.add_argument(
        "--blat-2bit",
        action="store",
        dest="blat_two_bit",
        help="Reference genome in 2bit format for BLAT aligner",
        required=False,
    )
    parser.add_argument(
        "--blat-nclosed",
        action="store_false",
        dest="blat_closed",
        default=DefaultOptions.blat_closed,
        help="Close BLAT server when the job is complete (default: %(default)s)",
    )
    parser.add_argument(
        "--blat-nsleep",
        action="store_false",
        dest="blat_sleep",
        default=DefaultOptions.blat_sleep,
        help="Whether to sleep randomly before starting BLAT server (default: %(default)s)",
    )

    parser.add_argument(
        "--blat-port",
        action="store",
        dest="blat_port",
        type=int,
        help="Port for BLAT server (default: %(default)s)",
        default=DefaultOptions.blat_port,
    )
    parser.add_argument(
        "--species",
        action="store",
        dest="species",
        help="Name of the species for the reference genome (default: %(default)s)",
        choices=DefaultOptions.species_choices,
        default=DefaultOptions.species,
    )
    parser.add_argument(
        "--circular-rna-filter",
        action="store",
        dest="circular_rna",
        help="The way of dealing with putative circular RNAs (default: %(default)s)",
        choices=DefaultOptions.circular_rna_choices,
        default=DefaultOptions.circular_rna,
    )
    parser.add_argument(
        "--off-exon-filter",
        action="store_false",
        dest="exon_filter",
        default=DefaultOptions.exon_filter,
        help="Turn on exon filter (default: %(default)s)",
    )

    parser.add_argument(
        "--rt-switching-filter",
        action="store",
        dest="rt_switching_filter_len",
        type=int,
        default=DefaultOptions.rt_switching_filter_len,
        help="Set the length threshold for RT switching filter. (default length: %(default)s)",
    )
    parser.add_argument(
        "--ncan",
        action="store_true",
        dest="noncanonical",
        default=DefaultOptions.noncanonical,
        help="Considering non-canonical splice sites (default: %(default)s)",
    )
    parser.add_argument(
        "--graph",
        action="store_true",
        dest="graph",
        default=DefaultOptions.graph,
        help="Whether to output transcript segment graph. (default: %(default)s)",
    )
    parser.add_argument(
        "--refine",
        action="store_true",
        dest="refine",
        default=DefaultOptions.refine,
        help="Whether to refine the transcript segment graph after construction. (default: %(default)s)",
    )
    parser.add_argument(
        "--prune-threshold",
        action="store",
        dest="prune_threshold",
        type=int,
        default=DefaultOptions.prune_threshold,
        help="Length threshold for pruning the transcript segment graph (default: %(default)s)",
    )

    parser.add_argument(
        "--max-allowed-nm",
        action="store",
        dest="max_allowed_nm",
        type=int,
        help="Maximum allowed edit distance (NM tag). (default: %(default)s)",
        default=DefaultOptions.max_allowed_nm,
    )

    parser.add_argument(
        "--max-allowed-ins",
        action="store",
        dest="max_allowed_ins",
        type=int,
        help="Maximum allowed micro-insertion length (default: %(default)s)",
        default=DefaultOptions.max_allowed_micro_insertion,
    )
    parser.add_argument(
        "--min-required-ins",
        action="store",
        dest="min_required_ins",
        type=int,
        help="Minimum required insertion length in read to infer chimeric alignment (default: %(default)s)",
        default=DefaultOptions.min_required_insertion_length,
    )
    parser.add_argument(
        "--min-soft-seg-len",
        action="store",
        dest="min_soft_seg_len",
        type=int,
        help="Minimum length of soft-clipped portion required to trigger BLAT alignment. (default: %(default)s)",
        default=DefaultOptions.min_soft_seg_len,
    )
    # Reads filter parameters
    parser.add_argument(
        "--long-indel-length",
        action="store",
        dest="long_indel_length",
        type=int,
        default=DefaultOptions.long_indel_length,
        help="Length cutoff for defining long indels in reads. (default: %(default)s)",
    )
    parser.add_argument(
        "--indel-fraction",
        action="store",
        dest="indel_fraction",
        type=float,
        default=DefaultOptions.indel_fraction,
        help="Maximum allowed fraction of long indels in the reads. (default: %(default)s)",
    )
    parser.add_argument(
        "--substitution-fraction",
        action="store",
        dest="substitutions_fraction",
        type=float,
        default=DefaultOptions.substitutions_fraction,
        help="Maximum allowed fraction of substitutions in the reads (default: %(default)s)",
    )
    # SR Rescuer parameters
    parser.add_argument(
        "--rescue-sr",
        action="store_true",
        dest="rescue_sr",
        default=DefaultOptions.rescue_sr,
        help="Whether to rescue SR for segment links (default: %(default)s)",
    )
    parser.add_argument(
        "--soft-len",
        action="store",
        dest="soft_len",
        type=int,
        help="Minimum length of soft-clipped portion to be rescued (default: %(default)s)",
        default=DefaultOptions.soft_len,
    )

    parser.add_argument(
        "--mismatch",
        action="store",
        dest="mismatch",
        type=int,
        help="Maximum number of mismatched bases allowed in a rescued segment (default: %(default)s)",
        default=DefaultOptions.mismatch,
    )
    parser.add_argument(
        "--alignment-fraction",
        action="store",
        dest="alignment_fraction",
        type=float,
        help="Minimum fraction of the sequence that must align in Smith-Waterman local alignment. (default: %(default)s)",
        default=DefaultOptions.alignment_fraction,
    )
    parser.add_argument(
        "--nbound",
        action="store_false",
        dest="bound",
        default=DefaultOptions.bound,
        help="Whether to add maximum increment limit using average reads depth when rescuing SR (default: %(default)s)",
    )

    parser.add_argument(
        "--ignore-circle",
        action="store_true",
        dest="ignore_circle",
        default=DefaultOptions.ignore_circle,
        help="Whether to export result when the transcript segment graph contains a circle (default: %(default)s)",
    )

    return parser
