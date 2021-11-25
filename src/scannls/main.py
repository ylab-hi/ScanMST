#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import sys
import textwrap
import time
from pathlib import Path

from loguru import logger

from . import __version__
from .classes import Blat
from .classes import LengthAction
from .draft.main import scan_bam
from .draft.main2 import BamScanner
from .draft.main2 import scan_run
from scannls.utils import external_tool_checking


def parse_args():
    parser = argparse.ArgumentParser(
        description="ScanNLS: Nonlinear splicing (NLS) events identification using transcriptomic long-reads data",
        epilog=textwrap.dedent(
            """Authors: Ting-You Wang and Yangyang Li, Hormel Institute, University of Minnesota, 2021"""
        ),
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s {}".format(__version__)
    )
    sub_parsers = parser.add_subparsers(help="sub-command help", dest="sub_command")

    draft_parser = sub_parsers.add_parser(
        "draft",
        help="add additional tags to build.BAM",
        description="%(prog)s -i input_bam_file -o output_bam_file -r ref_genome_fasta -g gtf_file [opts]",
        epilog=textwrap.dedent(
            """Authors: Ting-You Wang and Yangyang Li, Hormel Institute, University of Minnesota, 2021"""
        ),
    )

    draft_parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="Input BAM file",
        required=True,
    )
    draft_parser.add_argument(
        "-r",
        "--ref",
        action="store",
        dest="ref",
        help="reference genome in FASTA format (with fai index)",
        required=True,
    )
    draft_parser.add_argument(
        "-g",
        "--gtf",
        action="store",
        dest="gtf",
        help="gene annotations in GTF format",
        required=True,
    )
    draft_parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        help="output BAM file",
        required=True,
    )
    draft_parser.add_argument(
        "-s",
        "--splice_bin",
        action="store",
        dest="splice_bin",
        type=int,
        help="minimal observation count for ITD (default: %(default)s)",
        default=5,
    )
    draft_parser.add_argument(
        "-m",
        "--mapq",
        action="store",
        dest="mapq",
        type=int,
        help="minimal MAPQ in BAM for calling NLS (default: %(default)s)",
        default=15,
    )
    draft_parser.add_argument(
        "-n",
        "--noncanonical",
        action="store_true",
        dest="noncanonical",
        default=False,
        help="Considering Non-canonical spliced sites",
    )
    draft_parser.add_argument(
        "--log_level",
        action="store",
        dest="log",
        choices=["info", "debug", "trace", "warning", "error", "critical"],
        default="info",
        help="set log level (default: %(default)s)",
    )
    draft_parser.add_argument(
        "--parallel",
        action="store",
        dest="parallel",
        type=int,
        default=1,
        help="set working mode in processor (default: %(default)s)",
    )
    draft_parser.add_argument(
        "--2bit",
        action="store",
        dest="two_bit",
        help="reference genome in 2bit format",
        required=True,
    )
    draft_parser.add_argument(
        "-p",
        "--port",
        action="store",
        dest="port",
        type=int,
        help="port for BLAT server (default: %(default)s)",
        default=88888,
    )
    draft_parser.add_argument(
        "--min_soft_seg_len",
        action="store",
        dest="min_soft_seg_len",
        type=int,
        help="minimum softclipped segement length to trigger BLAT alignment (default: %(default)s)",
        default=200,
    )
    draft_parser.add_argument(
        "--max_allowed_nm",
        action="store",
        dest="max_allowed_nm",
        type=int,
        help="Maximum allowed NM to keep AS tag (default: %(default)s)",
        default=100,
    )
    draft_parser.add_argument(
        "--identity",
        action="store",
        dest="ident_cutoff",
        type=float,
        help="blat_ident_pct_cutoff (default: %(default)s)",
        default=0.99,
    )
    draft_parser.add_argument(
        "--tmp",
        action="store",
        dest="tmp_dir",
        type=str,
        help="BLAT temporary directory (default: %(default)s)",
        default="/tmp",
    )

    call_parser = sub_parsers.add_parser(
        "call",
        help="call NLS events from build.BAM",
        description="%(prog)s -i input_bam_file_from_ScanNLS_build -o output_vcf_file_prefix [opts]",
        epilog=textwrap.dedent(
            """Authors: Ting-You Wang and Yangyang Li, Hormel Institute, University of Minnesota, 2021"""
        ),
    )

    call_parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="Input BAM file",
        required=True,
    )
    call_parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        help="output file prefix",
        required=True,
    )
    call_parser.add_argument(
        "-a",
        "--alignment_fraction",
        action="store",
        dest="alignment_fraction",
        type=float,
        help="minimal fraction of aligned part for smith-waterman local alignment (default: %(default)s)",
        default=0.8,
    )
    call_parser.add_argument(
        "-c",
        "--sr",
        action="store",
        dest="sr",
        type=int,
        help="minimal observation supporting reads for SV (default: %(default)s)",
        default=4,
    )
    call_parser.add_argument(
        "-d",
        "--depth",
        action="store",
        dest="depth",
        type=int,
        help="minimal depth to call SV (default: %(default)s)",
        default=10,
    )
    call_parser.add_argument(
        "-p",
        "--pso",
        action="store",
        dest="pso",
        type=float,
        help="minimal variant allele frequency (default: %(default)s)",
        default=0.1,
    )
    call_parser.add_argument(
        "-l",
        "--length",
        action=LengthAction,
        dest="length",
        type=int,
        help="minimal length (>=1) of SV to report (default: %(default)s)",
        default=1000,
    )
    call_parser.add_argument(
        "-m",
        "--mapq",
        action="store",
        dest="mapq",
        type=int,
        help="minimal MAPQ of read from BAM file to call NLS (default: %(default)s)",
        default=15,
    )
    call_parser.add_argument(
        "-n",
        action="store",
        dest="mismatch",
        type=int,
        help="maximum mismatch bases of pairwise local alignment (default: %(default)s)",
        default=3,
    )
    call_parser.add_argument(
        "-s",
        "--seed",
        action="store",
        dest="seed",
        type=int,
        help="maximum seed observation (reads with SV tags) count of SV (default: %(default)s)",
        default=4,
    )
    call_parser.add_argument(
        "--soft_len",
        action="store",
        dest="soft_len",
        type=int,
        help="Minimal soft-clipped length to be count (default: %(default)s)",
        default=5,
    )
    call_parser.add_argument(
        "-t",
        action="store",
        dest="region",
        help="Limit analysis to targets listed in the BEDPE-format FILE",
    )
    call_parser.add_argument(
        "--log_level",
        action="store",
        dest="log",
        choices=["info", "debug"],
        default="info",
        help="set log level (default: %(default)s)",
    )

    isoform_parser = sub_parsers.add_parser(
        "isoform",
        help="infer NLS isoforms using built BAM and called VCF",
        description="%(prog)s -i input_bam_file_from_ScanNLS_build -o output_vcf_file_prefix [opts]",
        epilog=textwrap.dedent(
            """Authors: Ting-You Wang and Yangyang Li, Hormel Institute, University of Minnesota, 2021"""
        ),
    )
    isoform_parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="Input BAM file",
        required=True,
    )
    isoform_parser.add_argument(
        "-o",
        "--current_output",
        action="store",
        dest="current_output",
        help="current_output file prefix",
        required=True,
    )
    isoform_parser.add_argument(
        "-c",
        "--sr",
        action="store",
        dest="sr",
        type=int,
        help="minimal observation supporting reads for SV (default: %(default)s)",
        default=4,
    )
    isoform_parser.add_argument(
        "-d",
        "--depth",
        action="store",
        dest="depth",
        type=int,
        help="minimal depth to call SV (default: %(default)s)",
        default=10,
    )
    isoform_parser.add_argument(
        "-p",
        "--pso",
        action="store",
        dest="pso",
        type=float,
        help="minimal variant allele frequency (default: %(default)s)",
        default=0.1,
    )
    isoform_parser.add_argument(
        "--log",
        action="store",
        dest="log",
        choices=["info", "debug"],
        default="info",
        help="set log level (default: %(default)s)",
    )

    return parser


def main():
    if sys.version_info < (3, 8):
        sys.exit(
            "Sorry, this code need Python 3.8 or higher. Please update. Aborting..."
        )
    parser = parse_args()

    if len(sys.argv[1:]) < 1:
        parser.print_help()
        raise SystemExit
    else:
        options = parser.parse_args()

    if options.sub_command == "draft":
        # add logger
        logger.remove()
        logger.add(sys.stdout, level=options.log.upper())

        # check external tools used
        external_tool_checking(logger=logger)

        if options.parallel > 1:
            logger.info("ScanNLS draft starts running in parallel mode")
        else:
            logger.info("ScanNLS draft starts running in normal mode")

        logger.info(f"{options.input=}")
        start = time.time()
        blat = Blat(options.two_bit, logger, options.port, options.tmp_dir)
        blat.start_server()
        blat_logfile = blat.log_file
        # CIGAR string refinement or add SV tag
        motif_required = not options.noncanonical

        # scan_bam(
        #     input_bam=options.input,
        #     mapq_cutoff=options.mapq,
        #     output=options.output,
        #     ref_genome=options.ref,
        #     gtf=options.gtf,
        #     splice_bin=options.splice_bin,
        #     blat=blat,
        #     logger=logger,
        #     motif_required=motif_required,
        #     parallel=options.parallel,
        #     max_allowed_nm=options.max_allowed_nm,
        #     min_soft_seg_len=options.min_soft_seg_len,
        #     blat_ident_pct_cutoff=options.ident_cutoff,
        # )
        scan_run(
            two_bit=options.two_bit,
            port=options.port,
            tmp_dir=options.tmp_dir,
            blat_logfile=blat_logfile,
            in_bam_path=options.input,
            mapq_cutoff=options.mapq,
            output=options.output,
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
        #
        # bam_scanner = BamScanner(
        #     input_bam=Path(options.input),
        #     mapq_cutoff=options.mapq,
        #     output=Path(options.output),
        #     ref_genome=Path(options.ref),
        #     gtf=Path(options.gtf),
        #     splice_in=options.splice_bin,
        #     blat=blat,
        #     logger=logger,
        #     motif_required=motif_required,
        #     parallel=options.parallel,
        #     max_allowed_nm=options.max_allowed_nm,
        #     min_soft_seg_len=options.min_soft_seg_len,
        #     blat_ident_pct_cutoff=options.ident_cutoff,
        # )
        # intact_series_list = bam_scanner.run()

        logger.info("ScanNLS build running done")
        end = time.time()
        logger.info(f"ScanNLS build takes {end - start} seconds.")

    elif options.sub_command == "call":
        pass
        # print(
        #     "ScanNLS calling NLS events starts running: "
        #     + time.strftime("%Y-%m-%d %H:%M:%S")
        # )
        # start = time.time()
        # event_dict = joint_call(
        #     options.input,
        #     options.current_output,
        #     options.sr,
        #     options.depth,
        #     options.pso,
        #     options.length,
        #     options.soft_len,
        #     options.region,
        #     options.mapq,
        #     options.mismatch,
        #     options.alignment_fraction,
        #     options.seed,
        # )
        # print(
        #     "ScanNLS calling NLS events running done: "
        #     + time.strftime("%Y-%m-%d %H:%M:%S")
        # )
        # end = time.time()
        # print("ScanNLS calling NLS events takes " + str(end - start) + " seconds.")
    # infer transcript forms (GTF) and the corresponding sequences (FASTA)
    elif options.sub_command == "isoform":
        pass


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        # stop gfserver
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
