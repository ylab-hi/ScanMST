# !/usr/bin/env python
"""Parse command line arguments.

@Filename:    arg.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        3/25/22 9:01 AM
"""
import argparse
import textwrap
from typing import Any
from typing import Optional

from scannls import __version__


# https://github.com/plasma-umass/scalene/blob/master/scalene/scalene_parseargs.py#L11:7


class RichArgParser(argparse.ArgumentParser):
    """RichArgParser."""

    def __init__(self, *args: Any, **kwargs: Any):
        """RichArgParser."""
        from rich.console import Console

        self.console = Console()
        super().__init__(*args, **kwargs)

    def _print_message(self, message: Optional[str], file: Any = None) -> None:
        if message:
            self.console.print(message)


def parse_args() -> argparse.ArgumentParser:
    """Parse command line arguments."""
    parser = RichArgParser(
        description="[red]ScanNLS[/red]: Nonlinear splicing "
        "(NLS) events identification using transcriptomic"
        " long-reads data",
        epilog=textwrap.dedent(
            """Authors: Ting-You Wang and Yangyang Li, Hormel Institute,
            University of Minnesota, 2022"""
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    parser.add_argument(
        "--input",
        action="store",
        dest="input",
        help="Input BAM file",
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
        "--sr",
        action="store",
        dest="support_reads",
        type=int,
        help="minimum number of support reads for reporting NLS (default: %(default)s)",
        default=1,
    )
    parser.add_argument(
        "--splice-bin",
        action="store",
        dest="splice_bin",
        type=int,
        help="splice site bin size (default: %(default)s)",
        default=5,
    )
    parser.add_argument(
        "--mapq",
        action="store",
        dest="mapq",
        type=int,
        help="minimum MAPQ of reads for calling NLS (default: %(default)s)",
        default=15,
    )
    parser.add_argument(
        "--non-can",
        action="store_true",
        dest="noncanonical",
        default=False,
        help="Considering Non-canonical spliced sites",
    )
    parser.add_argument(
        "--log-level",
        action="store",
        dest="log",
        choices=["info", "debug", "trace", "warning", "error", "critical"],
        default="info",
        help="set log level (default: %(default)s)",
    )
    parser.add_argument(
        "--parallel",
        action="store",
        dest="parallel",
        type=int,
        default=1,
        help="set working mode in processor (default: %(default)s)",
    )
    parser.add_argument(
        "--2bit",
        action="store",
        dest="two_bit",
        help="reference genome in 2bit format",
        required=True,
    )
    parser.add_argument(
        "--nclosed",
        action="store_false",
        dest="closed",
        default=True,
        help="close BLAT server when job has done (default: %(default)s)",
    )
    parser.add_argument(
        "--nsleep",
        action="store_false",
        dest="nsleep",
        default=True,
        help="If sleep randomly before starting BLAT server (default: %(default)s)",
    )
    parser.add_argument(
        "--port",
        action="store",
        dest="port",
        type=int,
        help="port for BLAT server (default: %(default)s)",
        default=88888,
    )
    parser.add_argument(
        "--min-soft-seg-len",
        action="store",
        dest="min_soft_seg_len",
        type=int,
        help="minimum softclipped segment length to trigger BLAT alignment (default: %(default)s)",
        default=200,
    )
    parser.add_argument(
        "--max-allowed-nm",
        action="store",
        dest="max_allowed_nm",
        type=int,
        help="Maximum allowed NM to keep AS tag (default: %(default)s)",
        default=60,
    )
    parser.add_argument(
        "--identity",
        action="store",
        dest="ident_cutoff",
        type=float,
        help="blat_ident_pct_cutoff (default: %(default)s)",
        default=0.99,
    )
    parser.add_argument(
        "--tmp",
        action="store",
        dest="tmp_dir",
        type=str,
        help="BLAT temporary directory (default: %(default)s)",
        default="/tmp",
    )
    # Reads filter parameters
    parser.add_argument(
        "--long-indel-length",
        action="store",
        dest="long_indel_length",
        type=int,
        default=5,
        help="The length cutoff of defining long indels in the reads (default: %(default)s)",
    )
    parser.add_argument(
        "--substitution-num",
        action="store",
        dest="substitutions_num",
        type=int,
        default=5,
        help="The allowed maximum substitution number in the reads (default: %(default)s)",
    )
    parser.add_argument(
        "--substitution-fraction",
        action="store",
        dest="substitutions_fraction",
        type=float,
        default=0.2,
        help="The allowed maximum substitution fraction in the reads (default: %(default)s)",
    )
    parser.add_argument(
        "--indel-fraction",
        action="store",
        dest="indel_fraction",
        type=float,
        default=0.2,
        help="The allowed maximum long indel fraction in the reads (default: %(default)s)",
    )
    # SR Rescuer parameters
    parser.add_argument(
        "--soft-len",
        action="store",
        dest="soft_len",
        type=int,
        help="minimum softclipped segment length to be rescued (default: %(default)s)",
        default=5,
    )
    parser.add_argument(
        "--mismatch",
        action="store",
        dest="mismatch",
        type=int,
        help="maximum allowed mismatch bases of rescued segment (default: %(default)s)",
        default=3,
    )
    parser.add_argument(
        "--alignment-fraction",
        action="store",
        dest="alignment_fraction",
        type=float,
        help="minimal fraction of aligned part for smith-waterman local alignment (default: %(default)s)",
        default=0.8,
    )

    return parser
