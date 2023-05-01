# !/usr/bin/env python
"""CLi for scannls.

@Filename:    cli.py
@Author:      YangyangLi
@license:     MIT Licence
@Time:        1/11/22 4:28 PM
"""
import argparse
import os
import sys
import tempfile
import time
from functools import partial
from typing import Any
from typing import Union

from loguru import logger

from .. import Blat
from .. import ClusterFinder
from .. import FastaWriter
from .. import GTFWriter
from .. import LoggerType
from .. import MyLogger
from .. import ParallelWorker
from .. import SpliceGraph
from .. import VCFWriter
from .. import Writers
from ..utils import find_2bit_file
from ..utils import sleep
from .arg import DefaultOptions
from .main import scanbam_run


def get_writers(
    output_prefix: str,
    ref_path: str,
    bam_header: Any,
) -> Writers:
    """Get writers."""
    fasta_writer = FastaWriter(f"{output_prefix}.fasta", ref_path)
    gtf_writer = GTFWriter(f"{output_prefix}.gtf")
    vcf_writer = VCFWriter(
        f"{output_prefix}.vcf",
        ref_path,
        bam_header,
    )

    return Writers((fasta_writer, gtf_writer, vcf_writer))


def parse_splice_graph_for_cliques_seq(
    cliques: Any,
    writers: Writers,
    options: Union[DefaultOptions, argparse.Namespace],
    node_rescued_sr_maximum: int,
    logger: LoggerType,
    average_read_depth: Union[int, None] = None,
) -> None:
    """Parse splice graph for cliques."""
    splice_graph = SpliceGraph.create_graph(
        options.input,
        options.mapq,
        options.soft_len,
        options.mismatch,
        options.alignment_fraction,
        logger,
        options.prune_threshold,
        node_rescued_sr_maximum,
        average_read_depth,
    )

    with writers.open() as _:
        for ind, clique in enumerate(cliques, 1):
            logger.debug(f"processing clique {ind}")
            for series in splice_graph(clique, ind):
                if len(series) == 1:
                    logger.warning(
                        f"Single Series {ind}: {series}{series[0].query_name}"
                    )
                if series.is_all_node_sr_higher_than_threshold(options.support_reads):
                    logger.debug(f"Output Clique{ind}: {series}")
                    writers.write_series(series, ind)


def _parse_splice_graph_for_cliques_par(
    cliques: Any,
    options: Union[DefaultOptions, argparse.Namespace],
    node_rescued_sr_maximum: int,
    average_read_depth: Union[int, None],
):
    """Parse splice graph for cliques."""
    from loguru import logger

    logger = MyLogger(f"PID-{os.getpid()}", logger)  # type: ignore

    splice_graph = SpliceGraph.create_graph(
        options.input,
        options.mapq,
        options.soft_len,
        options.mismatch,
        options.alignment_fraction,
        logger,
        options.prune_threshold,
        node_rescued_sr_maximum,
        average_read_depth,
    )

    result_series = []
    for ind, clique in enumerate(cliques, 1):
        series_list = []
        for series in splice_graph(clique, ind):
            series_list.append(series)
        result_series.append(series_list)
    return result_series


def parse_splice_graph_for_cliques_par(
    cliques: Any,
    writers: Writers,
    options: Union[DefaultOptions, argparse.Namespace],
    node_rescued_sr_maximum: int,
    logger: LoggerType,
    average_read_depth: Union[int, None] = None,
) -> None:
    """Parse splice graph for cliques."""
    parallel_workers = ParallelWorker(
        partial(
            _parse_splice_graph_for_cliques_par,
            options=options,
            node_rescued_sr_maximum=node_rescued_sr_maximum,
            average_read_depth=average_read_depth,
        ),
        logger,
        options.parallel,
    )
    cliques = [[list(clique)] for clique in cliques]
    result = parallel_workers.map(
        cliques, chunksize=max(1, len(cliques) // parallel_workers.n_jobs)
    )
    with writers.open() as _:
        for ind, clique in enumerate(result, 1):
            for series in clique[0]:  # reduce list depth
                if len(series) == 1:
                    logger.warning(
                        f"Single Series {ind}: {series}{series[0].query_name}"
                    )
                if series.is_all_node_sr_higher_than_threshold(options.support_reads):
                    logger.debug(f"Output Clique{ind}: {series}")
                    writers.write_series(series, ind)


def cli(options: Union[argparse.Namespace, DefaultOptions]):
    """Cli function."""
    start = time.perf_counter()
    # add logger
    logger.remove()
    logger.add(
        sys.stdout,
        level=options.log.upper(),
        enqueue=True,
        colorize=True,
        backtrace=False,
        diagnose=True,
    )

    running_mode = "parallel" if options.parallel > 1 else "normal"
    logger.info(f"scannls starts running in {running_mode} mode PID-{os.getpid()}")
    logger.info(f"{options.input=} {options.closed=}")
    logger.info(f"{options.bound=}")

    tmp_dir = tempfile.TemporaryDirectory()
    # find 2bit file
    if options.two_bit is None:
        options.two_bit = find_2bit_file(options.ref)
    blat = Blat(options.two_bit, options.port, tmp_dir.name)
    # delay random seconds to preventing from starting multiple servers simultaneously
    if options.sleep:
        sleep(options.input)
    blat.start_server()
    blat_info = blat.log_file_path, blat.is_start_server
    # CIGAR string refinement
    motif_required = not options.noncanonical
    try:
        intact_series_list, in_bam_header, avg_cov = scanbam_run(
            two_bit=options.two_bit,
            port=options.port,
            tmp_dir=tmp_dir.name,
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
            indels_fraction=options.indel_fraction,
            species=options.species,
        )

        avg_cov = None if not options.bound else avg_cov

        intact_series_list_len = len(intact_series_list)

        if intact_series_list_len == 0:
            logger.warning("No valid series found")
            raise SystemExit

        logger.info(f"Total Series: {intact_series_list_len}")

        clique_finder = ClusterFinder(intact_series_list, intact_series_list_len)
        # cliques is generator
        cliques = clique_finder.find_cluster()

        writers = get_writers(options.output, options.ref, in_bam_header)

        parse_splice_graph_for_cliques = (
            parse_splice_graph_for_cliques_seq
            if options.parallel == 1
            else parse_splice_graph_for_cliques_par
        )

        node_rescued_sr_max = 100
        parse_splice_graph_for_cliques(
            cliques, writers, options, node_rescued_sr_max, logger, avg_cov
        )

        logger.info(f"ScanNLS takes {time.perf_counter() - start:.2f} seconds.")

    except KeyboardInterrupt:
        if options.closed and not blat.is_stop_server:
            logger.warning("KeyboardInterrupt")
            blat.stop_server()
            tmp_dir.cleanup()
        raise
    finally:
        if options.closed and not blat.is_stop_server:
            logger.info("Program ends")
            blat.stop_server()
            tmp_dir.cleanup()
