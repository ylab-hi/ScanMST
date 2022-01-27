# !/usr/bin/env python
"""CLi for scannls_simulator.

@Filename:    cli.py
@license:     MIT Licence
@Time:        1/11/22 4:28 PM
"""
import argparse
import sys
import textwrap
import time
from dataclasses import dataclass
from typing import Union

from loguru import logger

from .main import multi_hop_generator
from .main import single_hop_generator
from scannls import __version__
from scannls import external_tool_checking


@dataclass
class SimulatorOptions:
    """Cli options for testing."""

    input: str
    ref: str
    gtf: str
    output: str
    num: int
    hops: int = 2
    choice: str = "s"
    log: str = "info"
    shift: int = 7
    max_len: int = 200
    nls_type: str = "TDUP"


def parse_args() -> argparse.ArgumentParser:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="ScanNLS: Nonlinear splicing (NLS) events identification using transcriptomic"
        " long-reads data",
        epilog=textwrap.dedent(
            """Authors: Ting-You Wang and Yangyang Li, Hormel Institute,
            University of Minnesota, 2021"""
        ),
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "-g",
        "--gtf",
        action="store",
        dest="gtf",
        help="gene annotations in GTF format",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--ref",
        action="store",
        dest="ref",
        help="reference genome in FASTA format (with fai index)",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        help="output prefix (default: %(default)s)",
        default="output",
    )
    parser.add_argument(
        "-c",
        "--choice",
        action="store",
        dest="choice",
        choices=["s", "m"],
        help="Choose single-hop NLS (s) or multi-hop NLS (m) (default: %(default)s)",
        default="s",
    )
    parser.add_argument(
        "-n",
        "--num",
        action="store",
        dest="num",
        type=int,
        help="The number of transcripts for NLS generation",
    )
    parser.add_argument(
        "--hops",
        action="store",
        dest="hops",
        type=int,
        help="The number of hops in multi-hop NLS generation (default: %(default)s)",
        default=2,
    )
    parser.add_argument(
        "-s",
        "--shift",
        action="store",
        dest="shift",
        type=int,
        help="The shift size (default: %(default)s)",
        default=7,
    )
    parser.add_argument(
        "-m",
        "--max_len",
        action="store",
        dest="max_len",
        type=int,
        help="maximum sequence length of MT (intergenic and intronic) (default: %(default)s)",
        default=200,
    )
    parser.add_argument(
        "-t",
        "--type",
        action="store",
        dest="nls_type",
        choices=["TDUP", "IDUP", "INV", "TRA"],
        default="TDUP",
        help="The NLS type in one-hop NLS generation",
    )
    parser.add_argument(
        "--log_level",
        action="store",
        dest="log",
        choices=["info", "debug", "trace", "warning", "error", "critical"],
        default="info",
        help="set log level (default: %(default)s)",
    )

    return parser


def cli(options: Union[argparse.Namespace, SimulatorOptions]) -> None:
    """Cli function."""
    # add logger
    logger.remove()
    logger.add(
        sys.stdout,
        level=options.log.upper(),
        enqueue=True,
        colorize=True,
        backtrace=True,
        diagnose=True,
    )

    # check external tools used
    external_tool_checking(logger=logger)

    logger.info(f"{options.gtf=}")
    start = time.time()
    if options.choice == "s":
        single_hop_generator(
            annotation_gtf=options.gtf,
            reference=options.ref,
            num_of_transcripts=options.num,
            sv_type=options.nls_type,
            output_prefix=options.output,
            shift=options.shift,
            max_len=options.max_len,
            logger=logger,
        )
    else:
        multi_hop_generator(
            annotation_gtf=options.gtf,
            reference=options.ref,
            num_of_transcripts=options.num,
            num_of_hops=options.hops,
            output_prefix=options.output,
            shift=options.shift,
            max_len=options.max_len,
            logger=logger,
        )

    logger.info("ScanNLS simulator running done")
    end = time.time()
    logger.info(f"ScanNLS simulator takes {end - start} seconds.")
