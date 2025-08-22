"""Parse command line arguments."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Any

from scannls import __version__


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
    substitutions_num: int = 20
    substitutions_fraction: float = 0.2
    indel_fraction: float = 0.2
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
    "   _____                 _   ___________   ______  ________": COLOR,
    "  / ___/________ _____  / | / / ____/ /  /_  __/ / /_  __/": COLOR,
    "  \\__ \\/ ___/ __ `/ __ \\/  |/ / /   / /    / /  / __/ /": COLOR,
    " ___/ / /__/ /_/ / / / / /|  / /___/ /____/ /  / / / /": COLOR,
    "/____/\\___/\\__,_/_/ /_/_/ |_/\\____/_____/_/  /_/ /_/": COLOR,
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
        description="[red]scannls[/] :rocket: Non-co-linear transcripts (NCLT) identification using transcriptomic long reads data",
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
        help="input BAM file",
        required=True,
    )
    parser.add_argument(
        "--ref",
        action="store",
        dest="ref",
        help="reference genome in FASTA format (with fai index)",
        required=True,
    )
    parser.add_argument(
        "--gtf",
        action="store",
        dest="gtf",
        help="gene annotations in GTF format",
        required=True,
    )
    parser.add_argument(
        "--output",
        action="store",
        dest="output",
        help="output prefix",
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
        help="minimum number of support reads for reporting NCLT (default: %(default)s)",
        default=DefaultOptions.support_reads,
    )
    parser.add_argument(
        "--splice-bin",
        action="store",
        dest="splice_bin",
        type=int,
        help="splice site bin size (default: %(default)s)",
        default=DefaultOptions.splice_bin,
    )
    parser.add_argument(
        "--mapq",
        action="store",
        dest="mapq",
        type=int,
        help="minimum MAPQ of reads for calling NCLT (default: %(default)s)",
        default=DefaultOptions.mapq,
    )

    parser.add_argument(
        "--log-level",
        action="store",
        dest="log",
        choices=["info", "debug", "trace", "warning"],  # "warning", "error", "critical"
        default=DefaultOptions.log,
        help="set log level (default: %(default)s)",
    )
    parser.add_argument(
        "--thread",
        action="store",
        dest="thread",
        type=int,
        default=DefaultOptions.thread,
        help="set the thread number (default: %(default)s)",
    )
    parser.add_argument(
        "--aligner",
        dest="aligner",
        type=str,
        choices=DefaultOptions.aligner,
        help="aligner to use for mapping reads (default: %(default)s)",
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
        help="reference genome in 2bit format for blat aligner",
        required=False,
    )
    parser.add_argument(
        "--blat-nclosed",
        action="store_false",
        dest="blat_closed",
        default=DefaultOptions.blat_closed,
        help="close BLAT server when job has done (default: %(default)s)",
    )
    parser.add_argument(
        "--blat-nsleep",
        action="store_false",
        dest="blat_sleep",
        default=DefaultOptions.blat_sleep,
        help="if sleep randomly before starting BLAT server (default: %(default)s)",
    )

    parser.add_argument(
        "--blat-port",
        action="store",
        dest="blat_port",
        type=int,
        help="port for BLAT server (default: %(default)s)",
        default=DefaultOptions.blat_port,
    )
    parser.add_argument(
        "--species",
        action="store",
        dest="species",
        help="species name for reference genome (default: %(default)s)",
        choices=DefaultOptions.species_choices,
        default=DefaultOptions.species,
    )
    parser.add_argument(
        "--circular-rna-filter",
        action="store",
        dest="circular_rna",
        help="The way of dealing with circular RNAs (default: %(default)s)",
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
        help="Set RT switching filter (default length: %(default)s)",
    )
    parser.add_argument(
        "--ncan",
        action="store_true",
        dest="noncanonical",
        default=DefaultOptions.noncanonical,
        help="considering Non-canonical spliced sites  (default: %(default)s)",
    )
    parser.add_argument(
        "--graph",
        action="store_true",
        dest="graph",
        default=DefaultOptions.graph,
        help="if output graph (default: %(default)s)",
    )
    parser.add_argument(
        "--refine",
        action="store_true",
        dest="refine",
        default=DefaultOptions.refine,
        help="if refine the graph (default: %(default)s)",
    )
    parser.add_argument(
        "--nbound",
        action="store_false",
        dest="bound",
        default=DefaultOptions.bound,
        help="if add maximum increment limit using average reads depth when rescuing sr (default: %(default)s)",
    )
    parser.add_argument(
        "--max-allowed-nm",
        action="store",
        dest="max_allowed_nm",
        type=int,
        help="maximum allowed NM to keep AS tag (default: %(default)s)",
        default=DefaultOptions.max_allowed_nm,
    )

    parser.add_argument(
        "--max-allowed-ins",
        action="store",
        dest="max_allowed_ins",
        type=int,
        help="maximum allowed micro-insertion length (default: %(default)s)",
        default=DefaultOptions.max_allowed_micro_insertion,
    )

    parser.add_argument(
        "--min-required-ins",
        action="store",
        dest="min_required_ins",
        type=int,
        help="minimum required insertion length in read (default: %(default)s)",
        default=DefaultOptions.min_required_insertion_length,
    )

    # Reads filter parameters
    parser.add_argument(
        "--long-indel-length",
        action="store",
        dest="long_indel_length",
        type=int,
        default=DefaultOptions.long_indel_length,
        help="the length cutoff of defining long indel in the reads (default: %(default)s)",
    )
    parser.add_argument(
        "--substitution-num",
        action="store",
        dest="substitutions_num",
        type=int,
        default=DefaultOptions.substitutions_num,
        help="the allowed maximum substitution number in the reads (default: %(default)s)",
    )
    parser.add_argument(
        "--indel-fraction",
        action="store",
        dest="indel_fraction",
        type=float,
        default=DefaultOptions.indel_fraction,
        help="the allowed maximum long indel fraction in the reads (default: %(default)s)",
    )
    parser.add_argument(
        "--prune-threshold",
        action="store",
        dest="prune_threshold",
        type=int,
        default=DefaultOptions.prune_threshold,
        help="splice graph pruning length threshold (default: %(default)s)",
    )

    # SR Rescuer parameters
    parser.add_argument(
        "--soft-len",
        action="store",
        dest="soft_len",
        type=int,
        help="minimum softclipped segment length to be rescued (default: %(default)s)",
        default=DefaultOptions.soft_len,
    )
    parser.add_argument(
        "--mismatch",
        action="store",
        dest="mismatch",
        type=int,
        help="maximum allowed mismatch bases of rescued segment (default: %(default)s)",
        default=DefaultOptions.mismatch,
    )
    parser.add_argument(
        "--min-soft-seg-len",
        action="store",
        dest="min_soft_seg_len",
        type=int,
        help="minimum softclipped segment length to trigger BLAT alignment (default: %(default)s)",
        default=DefaultOptions.min_soft_seg_len,
    )
    parser.add_argument(
        "--alignment-fraction",
        action="store",
        dest="alignment_fraction",
        type=float,
        help="minimal fraction of aligned part for smith waterman local alignment (default: %(default)s)",
        default=DefaultOptions.alignment_fraction,
    )
    parser.add_argument(
        "--substitution-fraction",
        action="store",
        dest="substitutions_fraction",
        type=float,
        default=DefaultOptions.substitutions_fraction,
        help="the allowed maximum substitution fraction in the reads (default: %(default)s)",
    )
    parser.add_argument(
        "--ignore-circle",
        action="store_true",
        dest="ignore_circle",
        default=DefaultOptions.ignore_circle,
        help="if export result if the nlgraph has a circle  (default: %(default)s)",
    )
    parser.add_argument(
        "--rescue-sr",
        action="store_true",
        dest="rescue_sr",
        default=DefaultOptions.rescue_sr,
        help="if rescuing sr for edge  (default: %(default)s)",
    )

    return parser
