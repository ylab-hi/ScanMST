# !/usr/bin/env python
"""CLi for scannls.

@Filename:    cli.py
@license:     MIT Licence
@Time:        1/11/22 4:28 PM
"""
import argparse
import sys
import textwrap
import time
from typing import Union

from loguru import logger

from . import __version__
from . import Blat
from . import CliqueFinder
from . import FastaWriter
from . import GTFWriter
from . import Options
from . import SpliceGraph
from . import SRRescuer
from . import VCFWriter
from .core.main import scanbam_run
from .utils import external_tool_checking


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
        "-i",
        "--input",
        action="store",
        dest="input",
        help="Input BAM file",
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
        "-g",
        "--gtf",
        action="store",
        dest="gtf",
        help="gene annotations in GTF format",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        help="output BAM file",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--splice_bin",
        action="store",
        dest="splice_bin",
        type=int,
        help="minimal observation count for ITD (default: %(default)s)",
        default=5,
    )
    parser.add_argument(
        "-m",
        "--mapq",
        action="store",
        dest="mapq",
        type=int,
        help="minimal MAPQ in BAM for calling NLS (default: %(default)s)",
        default=15,
    )
    parser.add_argument(
        "-n",
        "--noncanonical",
        action="store_true",
        dest="noncanonical",
        default=False,
        help="Considering Non-canonical spliced sites",
    )
    parser.add_argument(
        "--log_level",
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
        "-c",
        "--closed",
        action="store_true",
        dest="closed",
        default=True,
        help="Considering Non-canonical spliced sites",
    )
    parser.add_argument(
        "-p",
        "--port",
        action="store",
        dest="port",
        type=int,
        help="port for BLAT server (default: %(default)s)",
        default=88888,
    )
    parser.add_argument(
        "--min_soft_seg_len",
        action="store",
        dest="min_soft_seg_len",
        type=int,
        help="minimum softclipped segment length to trigger BLAT alignment (default: %(default)s)",
        default=200,
    )
    parser.add_argument(
        "--max_allowed_nm",
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
    # SR Rescuer parameters
    parser.add_argument(
        "--soft_len",
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
        "-a",
        "--alignment_fraction",
        action="store",
        dest="alignment_fraction",
        type=float,
        help="minimal fraction of aligned part for smith-waterman local alignment (default: %(default)s)",
        default=0.8,
    )

    return parser


def cli(options: Union[argparse.Namespace, Options]) -> None:
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

    if options.parallel > 1:
        logger.info("scannls starts running in parallel mode")
    else:
        logger.info("scannls starts running in normal mode")

    logger.info(f"{options.input=}")
    start = time.time()
    blat = Blat(options.two_bit, logger, options.port, options.tmp_dir)
    blat.start_server()
    blat_info = blat.log_file_path, blat.is_start_server
    # CIGAR string refinement or add SV tag
    motif_required = not options.noncanonical

    intact_series_list, in_bam_io_object = scanbam_run(
        two_bit=options.two_bit,
        port=options.port,
        tmp_dir=options.tmp_dir,
        blat_info=blat_info,
        in_bam_path=options.input,
        mapq_cutoff=options.mapq,
        ref_genome=options.ref,
        gtf=options.gtf,
        splice_bin=options.splice_bin,
        blat=blat,
        logger=logger,
        motif_required=motif_required,
        parallel=options.parallel,
        max_allowed_nm=options.max_allowed_nm,
        min_soft_seg_len=options.min_soft_seg_len,
        blat_ident_pct_cutoff=options.ident_cutoff,
    )

    fasta_writer = FastaWriter("test.fasta", options.ref, logger)
    with fasta_writer.open() as _:
        fasta_writer.write_data(intact_series_list[0])

    gtf_writer = GTFWriter("test.gtf", logger)
    with gtf_writer.open() as _:
        gtf_writer.write_data(intact_series_list[0])

    vcf_writer = VCFWriter(
        "test.vcf", options.ref, in_bam_io_object, options.output, logger
    )
    with vcf_writer.open() as _:
        vcf_writer.write_header()
        vcf_writer.write_data(intact_series_list[0])

    logger.info(f"Total Series: {len(intact_series_list)}")
    rescuer = SRRescuer(
        in_bam_io_object,
        options.mapq,
        options.soft_len,
        options.mismatch,
        options.alignment_fraction,
        logger,
    )
    splice_graph = SpliceGraph(logger)
    clique_finder = CliqueFinder(intact_series_list, logger)
    # cliques is generator
    cliques = clique_finder.find_clique()
    for clique in cliques:
        for i in splice_graph(clique, rescuer):
            logger.debug(f"Series{i}")
    in_bam_io_object.close()

    if options.closed:
        blat.stop_server()

    logger.info("ScanNLS build running done")
    end = time.time()
    logger.info(f"ScanNLS build takes {end - start} seconds.")
