# !/usr/bin/env python
"""CLi for scannls.

@Filename:    cli.py
@license:     MIT Licence
@Time:        1/11/22 4:28 PM
"""
import argparse
import sys
import time
from pathlib import Path
from typing import Union

from loguru import logger

from . import Blat
from . import CliqueFinder
from . import FastaWriter
from . import GTFWriter
from . import Options
from . import SpliceGraph
from . import SRRescuer
from . import VCFWriter
from ._class.writer import Writers
from .core.main import scanbam_run
from .utils import external_tool_checking
from .utils import sleep


def cli(options: Union[argparse.Namespace, Options]):
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
    external_tool_checking(logger)

    if options.parallel > 1:
        logger.info("scannls starts running in parallel mode")
    else:
        logger.info("scannls starts running in normal mode")

    logger.info(f"{options.input=} {options.closed=}")

    tmp_dir = Path(options.tmp_dir)
    if not tmp_dir.exists():
        tmp_dir.mkdir()
        logger.info(f"Created temporary directory: {tmp_dir.resolve()}")
    start = time.time()
    blat = Blat(options.two_bit, logger, options.port, str(tmp_dir.resolve()))
    # delay random seconds to preventing from starting multiple servers simultaneously
    if options.nsleep:
        sleep(options.input)
    blat.start_server()
    blat_info = blat.log_file_path, blat.is_start_server
    # CIGAR string refinement
    motif_required = not options.noncanonical
    try:
        intact_series_list, in_bam_io_object = scanbam_run(
            two_bit=options.two_bit,
            port=options.port,
            tmp_dir=str(tmp_dir.resolve()),
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
            long_indel_length=options.long_indel_length,
            substitutions_num=options.substitutions_num,
            substitutions_fraction=options.substitutions_fraction,
            indels_fraction=options.indels_fraction,
        )

        intact_series_list_len = len(intact_series_list)

        if intact_series_list_len == 0:
            logger.warning("No valid series found")
            raise SystemExit

        logger.info(f"Total Series: {intact_series_list_len}")
        rescuer = SRRescuer(
            in_bam_io_object,
            options.input,
            options.mapq,
            options.soft_len,
            options.mismatch,
            options.alignment_fraction,
            logger,
        )
        splice_graph = SpliceGraph(logger, rescuer, options.prune_threshold)
        clique_finder = CliqueFinder(intact_series_list, intact_series_list_len, logger)
        # cliques is generator

        fasta_writer = FastaWriter(f"{options.output}.fasta", options.ref, logger)
        gtf_writer = GTFWriter(f"{options.output}.gtf", logger)
        vcf_writer = VCFWriter(
            f"{options.output}.vcf",
            options.ref,
            in_bam_io_object.header.as_dict(),
            logger,
        )

        writers = Writers((fasta_writer, gtf_writer, vcf_writer))

        cliques = clique_finder.find_clique()

        with writers.open() as _:
            for ind, clique in enumerate(cliques, 1):
                logger.debug(f"processing clique {ind}")
                for series in splice_graph(clique, ind, is_plot=False):
                    logger.debug(f"Output Clique{ind}: {series}")
                    if len(series) == 1:
                        logger.warning(
                            f"Single Series {ind}: {series}{series[0].query_name}"
                        )
                    if series.is_all_node_sr_higher_than_threshold(
                        options.support_reads
                    ):
                        writers.write_series(series, ind)

        in_bam_io_object.close()

        logger.info("ScanNLS build running done")
        end = time.time()
        logger.info(f"ScanNLS build takes {end - start} seconds.")

    except KeyboardInterrupt:
        if options.closed and not blat.is_stop_server:
            logger.info("KeyboardInterrupt")
            blat.stop_server()
        raise
    finally:
        if options.closed and not blat.is_stop_server:
            logger.info("Program ends")
            blat.stop_server()
