# !/usr/bin/env python
"""Module for nls inference."""
from typing import Any

import HTSeq
import pyfaidx

from ..base.type import LoggerType
from .helper import (
    diff_chrom_diff_strand_handler,
    diff_chrom_same_strand_handler,
    same_chrom_diff_strand_handler,
    same_chrom_same_strand_handler,
)


def infer_nls_from_connected_reads(
    read_lt,
    read_rt,
    lt_mode: int,
    rt_mode: int,
    splice_bin: int,
    genome_fasta: pyfaidx.Fasta,
    cvg: HTSeq.GenomicArrayOfSets,
    gene_iv: HTSeq.GenomicArrayOfSets,
    motif_required: bool,
    logger: LoggerType,
    microinsertion_cutoff: int = 60,
) -> Any:
    """Infer NLS event from connected reads.

    :param logger:
    :param read_lt: Read 1
    :param read_rt: Read 2
    :param lt_mode: mode of Read 1
    :param rt_mode: mode of Read 2
    :param splice_bin: a small bin for splice site searching
    :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
    :param cvg: annotated splice sites (HTSeq.GenomicArrayOfSets) of reference gene annotation
        (GTF file)
    :param gene_iv: annotated gene region (HTSeq.GenomicArrayOfSets) of reference gene annotation
        (GTF file)
    :param motif_required: considering canonical splice sites only OR considering both canonical
        and noncanonical splice sites
    :param microinsertion_cutoff: threshold of dicarding NLS event with long microinsertion (>60bp)
    :return: putative event from reads-pair

    .. note:: putative event

    examples: * 'NA', 0, 0, (), (), (), (), (), [] * 'TDUP', annotation, canonical/noncanonical,
     ('bp_chrm1:bp_pos1', 'bp_chrm2:bp_pos2', bp_mode1, bp_mode2),
     (bp_read1_ref_start, bp_read1_ref_end, bp_read1_exons), (bp_read2_ref_start, bp_read2_ref_end,
      bp_read2_exons), (lt_bp_seq, rt_bp_seq), (strand1, strand2), [gene1, gene2]

       annotation explanation:
       3(11) => both breakpoints overlap with coding exons boundary
       2(10) => one breakpoint overlap with coding exons boundary
       1(01) => one breakpoint overlap with coding exons boundary
       0(00) => none breakpoint overlap with coding exons boundary
    """
    noreturn = "NA", 0, 0, (), (), (), (), (), []  # type: ignore

    logger.trace(f"{read_lt=} {read_rt=}")

    if lt_mode == 3 or rt_mode == 3:
        return noreturn

    lt_chrm, lt_strand = (
        read_lt.chrom,
        read_lt.strand,
    )
    rt_chrm, rt_strand = (
        read_rt.chrom,
        read_rt.strand,
    )

    logger.trace(f"{lt_chrm=} {rt_chrm=}")
    logger.trace(f"{lt_strand=} {rt_strand=}")
    logger.trace(f"{lt_mode=} {rt_mode=}")

    if lt_chrm == rt_chrm:
        if lt_strand == rt_strand:  # deletion, insertion, duplication
            # If using noncanonical splice site, return NA
            if not read_lt.splice_site_checker(
                genome_fasta
            ) or not read_rt.splice_site_checker(genome_fasta):
                logger.debug(
                    f"Splice site checking[Same chroms, same strands]: "
                    f"{read_lt.query_name=}, {read_lt.cigarstring=}, {read_rt.cigarstring=}"
                )
                return noreturn
            return same_chrom_same_strand_handler(
                read_lt,
                read_rt,
                lt_mode,
                rt_mode,
                splice_bin,
                genome_fasta,
                cvg,
                gene_iv,
                motif_required,
                logger,
                microinsertion_cutoff,
            )
        else:  # lt_strand != rt_strand
            # IDUP and INV detection in this category
            return same_chrom_diff_strand_handler(
                read_lt,
                read_rt,
                lt_mode,
                rt_mode,
                splice_bin,
                genome_fasta,
                cvg,
                gene_iv,
                motif_required,
                logger,
                microinsertion_cutoff,
            )
    else:  # lt_chrm != rt_chrm
        # If using noncanonical splice site, return NA
        if not read_lt.splice_site_checker(
            genome_fasta
        ) or not read_rt.splice_site_checker(genome_fasta):
            logger.debug(
                f"Splice site checking[different chroms]: "
                f"{read_lt.query_name=}, {read_lt.cigarstring=}, {read_rt.cigarstring=}"
            )
            return noreturn
        if lt_strand == rt_strand:
            return diff_chrom_same_strand_handler(
                read_lt,
                read_rt,
                lt_mode,
                rt_mode,
                splice_bin,
                genome_fasta,
                cvg,
                gene_iv,
                motif_required,
                logger,
                microinsertion_cutoff,
            )
        else:  # lt_strand != rt_strand
            return diff_chrom_diff_strand_handler(
                read_lt,
                read_rt,
                lt_mode,
                rt_mode,
                splice_bin,
                genome_fasta,
                cvg,
                gene_iv,
                motif_required,
                logger,
                microinsertion_cutoff,
            )
