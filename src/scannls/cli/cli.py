"""CLi for scannls."""

from __future__ import annotations

import os
import sys
import tempfile
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any

from loguru import logger

from scannls.base import Blat
from scannls.graph import ClusterFinder, NLGraph
from scannls.utils import find_2bit_file, sleep, wait_for_aligner
from scannls.writer import FastaWriter, GTFWriter, VCFWriter, Writers
from scannls.writer.tsg_writer import TSGWriter

from .main import scanbam_run

if TYPE_CHECKING:
    import argparse

    from scannls.mtype import LoggerType

    from .arg import DefaultOptions


sys.setrecursionlimit(10000)


def get_writers(
    output_prefix: str,
    ref_path: str,
    rescue_sr: bool,
    output_sequence_choice: str,
    read_name_to_seq_dict: dict,
    bam_header: Any,
) -> Writers:
    """Get writers."""
    gtf_writer = GTFWriter(f"{output_prefix}.gtf", rescue_sr)
    vcf_writer = VCFWriter(
        f"{output_prefix}.vcf",
        rescue_sr,
        ref_path,
        bam_header,
    )
    tsg_writer = TSGWriter(f"{output_prefix}.tsg")

    if output_sequence_choice in {"reference", "consensus"}:
        fasta_writer = FastaWriter(f"{output_prefix}.fasta", ref_path, output_sequence_choice, read_name_to_seq_dict)
        return Writers((fasta_writer, gtf_writer, vcf_writer, tsg_writer))

    fasta_writer1 = FastaWriter(f"{output_prefix}.reference.fasta", ref_path, "reference", read_name_to_seq_dict)
    fasta_writer2 = FastaWriter(f"{output_prefix}.consensus.fasta", ref_path, "consensus", read_name_to_seq_dict)

    return Writers((fasta_writer1, fasta_writer2, gtf_writer, vcf_writer, tsg_writer))


def parse_nlgraph_for_cluster_seq(
    clusters: Any,
    writers: Writers,
    options: DefaultOptions | argparse.Namespace,
    node_rescued_sr_maximum: int,
    logger: LoggerType,
    output_dir: Path,
    average_read_depth: int | None = None,
) -> None:
    """Parse splice graph for cliques."""
    nlgraph = NLGraph.create_graph(
        options.input,
        options.mapq,
        options.soft_len,
        options.mismatch,
        options.alignment_fraction,
        logger,
        options.prune_threshold,
        options.support_reads,
        node_rescued_sr_maximum,
        average_read_depth,
        output_dir,
        ignore_circle=options.ignore_circle,
        rescue_sr=options.rescue_sr,
        refine=options.refine,
    )

    with writers.open():
        for ind, cluster in enumerate(clusters, 1):
            logger.debug(f"Read guided: Processing Cluster {ind=}")
            graph_id = f"TSG{ind:010d}"
            for nlpath in nlgraph(cluster, graph_id, is_plot=options.graph):
                if len(nlpath) == 1:
                    logger.warning(
                        f"Single nlpath {ind=}: {nlpath}{nlpath[0].query_name}",
                    )

                logger.debug(f"cluster {graph_id=} output {nlpath=} ")
                writers.write_path(nlpath, graph_id)

            writers.write_graph(nlgraph, graph_id)


def cli(options: argparse.Namespace | DefaultOptions):
    """Cli function."""
    start = time.perf_counter()

    working_dir = Path.cwd()

    options.input = Path(options.input).resolve().as_posix()

    logger.remove()
    if options.log.upper() == "INFO":
        info_format = "<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> | <level>{message}</level>"
        logger.add(
            sys.stdout,
            level=options.log.upper(),
            format=info_format,
            enqueue=True,
            colorize=True,
            backtrace=False,
            diagnose=True,
        )
    else:
        logger.add(
            sys.stdout,
            level=options.log.upper(),
            enqueue=True,
            colorize=True,
            backtrace=False,
            diagnose=True,
        )

    running_mode = "parallel" if options.thread > 1 else "normal"

    output_prefix_path = Path(options.output)
    if output_prefix_path.is_absolute():
        output_file_path = output_prefix_path
    else:
        out_file_path = working_dir / output_prefix_path
        output_file_path = out_file_path.resolve()

    output_dir = output_file_path.parent
    if not output_dir.is_dir():
        msg = f"Error: The directory {output_dir} does not exist."
        raise SystemExit(msg)

    logger.info(f"scannls starts running in {running_mode} mode PID-{os.getpid()}")
    logger.info(f"{options.input=} {options.blat_closed=}")
    logger.info(f"{output_file_path=}")
    logger.info(f"{options.bound=}")

    tmp_dir = tempfile.TemporaryDirectory()
    if options.aligner == "blat":
        # find 2bit file
        if options.blat_two_bit is None:
            options.two_bit = find_2bit_file(options.ref)
        aligner = Blat(options.blat_two_bit, options.blat_port, tmp_dir.name)
        # delay random seconds to preventing from starting multiple servers simultaneously
        if options.blat_sleep:
            sleep(options.input)
        aligner.start_server()
        blat_info = aligner.log_file_path, aligner.is_start_server
    else:
        blat_info, aligner = None, None

    if aligner:
        try:
            wait_for_aligner(aligner)
        except TimeoutError as e:
            logger.error(str(e))
            aligner.stop_server()
            tmp_dir.cleanup()
            raise SystemExit

    # CIGAR string refinement
    motif_required = not options.noncanonical
    try:
        intact_nlpaths, intact_read_query_name_to_sequence, in_bam_header, avg_cov = scanbam_run(
            blat_two_bit=options.blat_two_bit,
            blat_port=options.blat_port,
            tmp_dir=tmp_dir.name,
            blat_info=blat_info,
            in_bam_path=options.input,
            mapq_cutoff=options.mapq,
            ref_genome=options.ref,
            gtf=options.gtf,
            splice_bin=options.splice_bin,
            aligner=aligner,
            logger=logger,
            motif_required=motif_required,
            parallel=options.thread,
            max_allowed_nm=options.max_allowed_nm,
            min_soft_seg_len=options.min_soft_seg_len,
            blat_ident_pct_cutoff=options.ident_cutoff,
            long_indel_length=options.long_indel_length,
            substitutions_num=options.substitutions_num,
            substitutions_fraction=options.substitutions_fraction,
            indels_fraction=options.indel_fraction,
            species=options.species,
            circular_rna=options.circular_rna,
            exon_filter=options.exon_filter,
            rt_switching_filter_len=options.rt_switching_filter_len,
            prune_threshold=options.prune_threshold,
            max_allowed_ins=options.max_allowed_ins,
            read_insertion_len_threshold=options.min_required_ins,
        )

        avg_cov = None if not options.bound else avg_cov

        intact_nlpaths_len = len(intact_nlpaths)

        if intact_nlpaths_len == 0:
            logger.warning("No valid path found")
            raise SystemExit

        logger.info(f"Total nlpaths: {intact_nlpaths_len}")

        cluster_finder = ClusterFinder(intact_nlpaths, options.prune_threshold)
        # cliques is generator
        clusters = cluster_finder.merge_cluster()

        writers = get_writers(
            str(output_file_path), options.ref, options.rescue_sr, options.output_sequence_choice, intact_read_query_name_to_sequence, in_bam_header
        )
        parse_splice_graph_for_cluster = parse_nlgraph_for_cluster_seq

        node_rescued_sr_max = 100
        parse_splice_graph_for_cluster(
            clusters,
            writers,
            options,
            node_rescued_sr_max,
            logger,  # type: ignore
            output_dir,
            avg_cov,
        )

        logger.info(f"ScanNLS takes {time.perf_counter() - start:.2f} seconds.")

    except KeyboardInterrupt:
        logger.warning("KeyboardInterrupt")
        if options.aligner == "blat" and aligner and options.blat_closed and not aligner.is_stop_server:
            aligner.stop_server()
            tmp_dir.cleanup()
        raise
    finally:
        logger.info("Program ends")
        if options.aligner == "blat" and aligner and options.blat_closed and not aligner.is_stop_server:
            aligner.stop_server()
            tmp_dir.cleanup()
