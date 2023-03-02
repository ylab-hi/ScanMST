"""Module contains the main function of the draft scannls."""
from asyncio import threads
import copy
import inspect
import math
import re
from itertools import chain
from pathlib import Path
from typing import Any

import HTSeq
import pyfaidx
import pysam
from pyfaidx import Fasta, FastaNotFoundError
import rscannls
from loguru import logger

from .. import (
    Blat,
    Event,
    MyLogger,
    ParallelWorker,
    Series,
    cigarstring2cigartuple,
    detect_read_read_connections_from_cigar,
    get_longest_insertion_sequence,
    get_softclip_length,
    reverse_complement,
)

from ..base.type import LoggerType
from .helper import (
    blat2chimeric_alignment,
    extract_splice_sites,
    get_transcriptome_length,
    insertion2chimeric_alignment,
    obtain_variants_stats,
    strand_mode_checker,
)

from .nls_inference import infer_nls_from_connected_reads


def get_genome_fasta(ref_genome):
    """Get the genome fasta file."""
    try:
        return Fasta(str(ref_genome), sequence_always_upper=True)
    except FastaNotFoundError:
        raise SystemExit(
            f"Reference File {ref_genome} is Not Found!"
        ) from FastaNotFoundError


def get_cvg_gene_iv(gtf, splice_bin):
    """Get the gene coverage interval.

    :param gtf: gtf file
    :param splice_bin: splice_bin file
    """
    try:
        return extract_splice_sites(str(gtf), splice_bin)
    except OSError:
        raise SystemExit(f"Reading GTF file {gtf} error!") from OSError


def detect_sv_from_cigar(
    *,
    read: Any,
    mapq_cutoff: int,
    max_allowed_nm: int,
    splice_bin: int,
    genome_fasta: pyfaidx.Fasta,
    cvg: HTSeq.GenomicArrayOfSets,
    gene_iv: HTSeq.GenomicArrayOfSets,
    motif_required: bool,
    blat: Blat,
    logger: LoggerType,
) -> Any:
    """Detect SV from cigar string.

    :param logger: logger for logging
    :param blat: `class.Blat`
    :param read: A read from pysam.AlignedSegment
    :param mapq_cutoff: MAPQ cutoff
    :param max_allowed_nm: NM cutoff
    :param splice_bin: a small bin for splice site searching
    :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
    :param cvg: annotated splice sites (HTSeq.GenomicArrayOfSets) of reference gene
           annotation (GTF file)
    :param gene_iv: annotated gene region (HTSeq.GenomicArrayOfSets) of reference
           gene annotation (GTF file)
    :param motif_required: considering canonical splice sites only OR considering both canonical
           and noncanonical splice sites
    :return: event groups in a list, every group is also a list
    :rtype: list (list of lists)
    """
    (
        read_chains,
        reads_pair_mode_dict,
        num_added_reads,
    ) = detect_read_read_connections_from_cigar(
        read=read,
        mapq_cutoff=mapq_cutoff,
        max_allowed_nm=max_allowed_nm,
        blat=blat,
        logger=logger,
    )

    logger.trace(f"{read_chains=}")
    event_list: list[Event] = []

    if read_chains:
        # every chain is a group of connected reads
        # every chain may have a list of events
        for _lt, _rt in zip(read_chains[::1], read_chains[1::1]):
            if (_lt, _rt) in reads_pair_mode_dict:
                _lt_mode, _rt_mode = reads_pair_mode_dict[(_lt, _rt)]

            elif (_rt, _lt) in reads_pair_mode_dict:
                _rt_mode, _lt_mode = reads_pair_mode_dict[(_rt, _lt)]

            if not strand_mode_checker(_lt.strand, _rt.strand, _lt_mode, _rt_mode):
                logger.warning(
                    f"{_lt.strand=}, {_rt.strand=}, {_lt_mode=}, {_rt_mode=}"
                )

            event = Event(
                infer_nls_from_connected_reads(
                    read_lt=_lt,
                    read_rt=_rt,
                    lt_mode=_lt_mode,
                    rt_mode=_rt_mode,
                    splice_bin=splice_bin,
                    genome_fasta=genome_fasta,
                    cvg=cvg,
                    gene_iv=gene_iv,
                    motif_required=motif_required,
                    logger=logger,
                )
            )

            if not event.is_type_na():
                event_list.append(event)
                logger.trace(str(event))
            else:  # temporary solution
                logger.warning(f"Event Type is NA {event=}")

    return event_list, read_chains, num_added_reads


class Mrecord:
    def __init__(
        self,
        reference_name,
        reference_start,
        cigarstring,
        mapping_quality,
        is_reverse,
        query_sequence,
    ):
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.cigarstring = cigarstring
        self.mapping_quality = mapping_quality
        self.is_reverse = is_reverse
        self.query_sequence = query_sequence
        self.tags = {}
        self.query_qualities = None
        self.is_supplementary = None
        self.query_name = None

    def __str__(self):
        return f"{self.reference_start=} {self.cigarstring=} {self.is_reverse=} {self.is_supplementary=} {self.tags['SA']}"

    __repr__ = __str__  # for debugging

    def set_tag(self, key, value):
        self.tags[key] = value

    def get_tag(self, key):
        return self.tags[key]

    @classmethod
    def from_alignment(cls, alignment_info):
        # TTTGAGGTTTCTAAATACATTAAAGTTATTTCTTAAGAA-false-name;chr1,3847474,-,841S140M994N174M3513N127M4467N309M,60,0; chr1,3479514,+,746S77M221N102M534N13M1D46M1I606M,60,2
        logger.warning(f"{alignment_info=}")
        aln_info_list = alignment_info.split(";")

        if len(aln_info_list) < 2:
            logger.warning(f"Invalid alignment info: {alignment_info}")
            return None

        # chr1, 3479514, +, 746S77M221N102M534N13M1D46M1I606M,60,2
        read_info = aln_info_list.pop().split(",")

        temp = aln_info_list.pop(0).split("-")

        if len(temp) != 3:
            logger.warning(f"Invalid alignment info: {alignment_info}")
            return None

        (sequence, is_supplementary, read_name) = temp

        record = cls(
            reference_name=read_info[0],
            reference_start=int(read_info[1]),
            cigarstring=read_info[3],
            mapping_quality=int(read_info[4]),
            is_reverse=True if read_info[2] == "-" else False,
            query_sequence=sequence,
        )

        record.is_supplementary = False
        record.query_name = read_name

        record.set_tag("SA", ";".join(aln_info_list) + ";")
        record.set_tag("NM", int(read_info[5]))

        logger.debug(f"{record.query_sequence=}")

        return record


def detect_sv_from_cigar_wrapper(
    read,
    mapq_cutoff,
    max_allowed_nm,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    blat,
    logger,
):
    return detect_sv_from_cigar(
        read=read,
        mapq_cutoff=mapq_cutoff,
        max_allowed_nm=max_allowed_nm,
        splice_bin=splice_bin,
        genome_fasta=genome_fasta,
        cvg=cvg,
        gene_iv=gene_iv,
        motif_required=motif_required,
        blat=blat,
        logger=logger,
    )


def scanbam_run(
    two_bit,
    port,
    tmp_dir,
    blat_info,
    in_bam_path,
    mapq_cutoff,
    ref_genome,
    gtf,
    splice_bin,
    blat,
    logger,
    motif_required,
    parallel,
    max_allowed_nm,
    min_soft_seg_len,
    blat_ident_pct_cutoff,
    long_indel_length,
    substitutions_num,
    substitutions_fraction,
    indels_fraction,
    species,
):
    """Main function to run scanbam."""
    #       bam_path: &str,
    #       long_indel_threshold: usize,
    #       substitutions_threshold: usize,
    #       substitutions_fraction_threshold: f32,
    #       indels_fraction_threshold: f32

    collect_cigar_option = rscannls.CollectCigarOption(
        bam_path=in_bam_path,
        long_indel_threshold=long_indel_length,
        substitutions_threshold=substitutions_num,
        substitutions_fraction_threshold=substitutions_fraction,
        indels_fraction_threshold=indels_fraction,
    )
    # bam_path: &str,
    #       fasta_path: &str,
    #       gtf_path: &str,
    #       threads: usize,
    #       mapping_quality_threshold: usize,
    #       max_allowed_nm: usize,
    #       insertion_length_threshold: usize,
    #       insertion_alignment_diff_theshold: usize,
    #       softclip_length_threshold: usize,
    #       collect_cigar_option: PyCollectCigarOption,
    #   )
    scan_bam_option = rscannls.ScanBamOption(
        bam_path=in_bam_path,
        fasta_path=ref_genome,
        gtf_path=gtf,
        threads=parallel,
        mapping_quality_threshold=mapq_cutoff,
        max_allowed_nm=max_allowed_nm,
        insertion_length_threshold=50,
        insertion_alignment_diff_theshold=5,
        softclip_length_threshold=min_soft_seg_len,
        collect_cigar_option=collect_cigar_option,
    )

    result = rscannls.scan_bam(scan_bam_option)

    genome_fasta = get_genome_fasta(ref_genome)
    cvg, gene_iv = get_cvg_gene_iv(gtf, splice_bin)
    blat_log_file, blat_is_start_server = blat_info

    blat = Blat(
        two_bit, logger, port, tmp_dir, blat_log_file, blat_is_start_server, None
    )

    nls_src_forms_list = []

    for alignment_info in result:
        read = Mrecord.from_alignment(alignment_info)
        if read is None:
            continue

        event_lists, read_chains, num_added_reads = detect_sv_from_cigar_wrapper(
            read=read,
            mapq_cutoff=mapq_cutoff,
            max_allowed_nm=max_allowed_nm,
            splice_bin=splice_bin,
            genome_fasta=genome_fasta,
            cvg=cvg,
            gene_iv=gene_iv,
            motif_required=motif_required,
            blat=blat,
            logger=logger,
        )

        nls_event_list = []
        for event in event_lists:
            if event.sv_type in {"TDUP", "INV", "TRA", "DEL", "IDUP"}:
                nls_event_list.append(event)

        if nls_event_list:
            series = Series(blat=blat, logger=logger)
            logger.debug(f"{nls_event_list=}")
            series.init(
                nls_event_list,
                read_chains,
                splice_bin,
                genome_fasta,
                cvg,
                gene_iv,
                motif_required,
            )
            series.disable_blat_logger()
            if not series.is_all_type_del():
                nls_src_forms_list.append(series)
                logger.trace(f"{series=}")

    return nls_src_forms_list
