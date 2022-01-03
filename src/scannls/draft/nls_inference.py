# !/usr/bin/env python
"""Module for nls inference."""
from typing import Any
from typing import Optional

import HTSeq  # type: ignore
import pyfaidx  # type: ignore
from loguru._logger import Logger

from .helper import gene_annotation
from .helper import splicing_confirmation


def short_tdup_or_not(
    ra_mode: int, read_sa, bp_region_seq_len: int, ins_seq_in_read: str, logger: Logger
) -> Optional[bool]:
    """Judge the ins_seq_in_read is a TDUP (TDUP size < reads length).

    OR novel sequence insertion using reference sequence inferred
    from chimeric alignment start position and indel_size from 'query_offset - target_offset'

    :param ra_mode: representative alignment mode
    :param read_sa: supplementary read mode
    :param bp_region_seq_len: softclipping read size
    :param ins_seq_in_read: putative insertion sequence from the representative read
    :param logger: logger
    :type ra_mode: int
    :type read_sa: Read
    :type bp_region_seq_len: int
    :type ins_seq_in_read: str
    :type logger: Logger
    :returns: True if it is a short TDUP
    :rtype: bool
    """
    event_size = len(ins_seq_in_read)

    if bp_region_seq_len >= 0:
        soft_extension_size = bp_region_seq_len
        matched_reduced_size = 0
    else:
        soft_extension_size = 0
        matched_reduced_size = -bp_region_seq_len

    # logger.trace(f"{event_size=} {matched_reduced_size=} {soft_extension_size=}")

    read_sa_matched_segment = read_sa.query_sequence[
        read_sa.lt_soft_len : read_sa.query_length - read_sa.rt_soft_len
    ]

    diff_len = (
        len(read_sa_matched_segment)
        - matched_reduced_size
        + soft_extension_size
        - event_size
    )
    # len(ins_seq_in_read) > read_sa_matched_segment
    if diff_len < 0:
        compared_seq = (
            read_sa.query_sequence[: read_sa.lt_soft_len][
                read_sa.lt_soft_len - soft_extension_size :
            ]
            + read_sa_matched_segment[
                matched_reduced_size : len(read_sa_matched_segment)
            ]
            if ra_mode == 1
            else read_sa_matched_segment[
                : len(read_sa_matched_segment) - matched_reduced_size
            ]
            + read_sa.query_sequence[-read_sa.rt_soft_len :][:soft_extension_size]
        )
        logger.trace(f"{compared_seq=}")
        ins_seq_in_read_modified = (
            ins_seq_in_read[: event_size + diff_len]
            if ra_mode == 1
            else ins_seq_in_read[-diff_len:]
        )
        logger.trace(f"{ins_seq_in_read_modified=}")
        return compared_seq == ins_seq_in_read_modified
    else:
        compared_seq = (
            read_sa.query_sequence[: read_sa.lt_soft_len][
                read_sa.lt_soft_len - soft_extension_size :
            ]
            + read_sa_matched_segment[
                matched_reduced_size : len(read_sa_matched_segment) - diff_len
            ]
            if ra_mode == 1
            else read_sa_matched_segment[
                diff_len : len(read_sa_matched_segment) - matched_reduced_size
            ]
            + read_sa.query_sequence[-read_sa.rt_soft_len :][:soft_extension_size]
        )
        logger.trace(f"{compared_seq=}")
        logger.trace(f"{ins_seq_in_read=}")
        return compared_seq == ins_seq_in_read


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
    logger: Logger,
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

    def softclipped_length_and_event_size_checker(
        read, mode, event_size, bp_region_seq_len
    ) -> bool:
        """When read length > predicted tandem duplication size.

        check whether the softclipped length is less than the inferred event size

        :param read: a chimeric read
        :param mode: mode for the chimeric read
        :param event_size: event size inferred from 'query_offset - target_offset'

        :type read : Read
        :type mode: int
        :type event_size: int
        :return: whether event_size > softclipped_length (If it is True, it will be a TDUP event)
        :rtype: bool
        """
        flag = False
        if mode == 2:
            if read.lt_soft_len < event_size + bp_region_seq_len:
                flag = True
        else:
            if read.rt_soft_len < event_size + bp_region_seq_len:
                flag = True
        return flag

    def obtain_ins_seq_from_softclipped_part_read(
        read, mode, event_size, bp_region_seq_len
    ) -> str:
        """Obtain insertion sequence from soft-clipped part of read.

        :param read: a chimeric read
        :param mode: mode for the chimeric read
        :param event_size: event size inferred from 'query_offset - target_offset'
        :param bp_region_seq_len: bp_region_seq_len

        :type read : Read
        :type mode: int
        :type event_size: int
        :return: putative insertion sequence from the read

        """
        read_seq = read.query_sequence

        ins_seq_in_read = (
            read_seq[: read.lt_soft_len][-(event_size + bp_region_seq_len) :]
            if mode == 2
            else read_seq[-read.rt_soft_len :][: (event_size + bp_region_seq_len)]
        )

        return ins_seq_in_read

    def obtain_bp_region_seq(read, mode, bp_region_seq_len) -> Any:
        """Obtain breakpoint region sequence from read.

        :param bp_region_seq_len: the length of the breakpoint region sequence
        :param read:  the chimeirc read
        :param mode: mode for the chimeirc read
        :type mode: int
        :return: (putative insertion/microhomology sequence from the read; + means insertion,
            - means microhomology, mode)

        ..note:
            * inserted sequence:
              S-----SM---M    M---MS-----S
              SSSSSXXMMMMM    MMMMMXXSSSSS
            * microhomology:
              S-----SM---M    M---MS-----S
              SSSSSSSXXMMM    MMMXXSSSSSSS
        """
        read_seq = read.query_sequence
        bp_region_seq = ""
        # inserted sequence
        if bp_region_seq_len > 0:
            if mode == 2:  # SM
                bp_region_seq = read_seq[: read.lt_soft_len][-bp_region_seq_len:]
            elif mode == 1:  # MS
                bp_region_seq = read_seq[-read.rt_soft_len :][:bp_region_seq_len]
            bp_region_seq = "+" + bp_region_seq
        # microhomology
        elif bp_region_seq_len < 0:
            if mode == 2:  # SM
                bp_region_seq = read_seq[
                    read.lt_soft_len : read.lt_soft_len - bp_region_seq_len
                ]
            elif mode == 1:  # MS
                bp_region_seq = read_seq[
                    bp_region_seq_len - read.rt_soft_len : -read.rt_soft_len
                ]
            bp_region_seq = "-" + bp_region_seq
        else:
            bp_region_seq = ""
        return bp_region_seq

    noreturn = "NA", 0, 0, (), (), (), (), (), []  # type: ignore
    logger.trace(f"{read_lt=} {read_rt=}")
    if lt_mode == 3 or rt_mode == 3:
        return noreturn
    lt_chrm, lt_strand, lt_start, lt_end, lt_cigartuples, lt_cigarstring = (
        read_lt.chrom,
        read_lt.strand,
        read_lt.ref_start,
        read_lt.ref_end,
        read_lt.cigartuples,
        read_lt.cigarstring,
    )
    rt_chrm, rt_strand, rt_start, rt_end, rt_cigartuples, rt_cigarstring = (
        read_rt.chrom,
        read_rt.strand,
        read_rt.ref_start,
        read_rt.ref_end,
        read_rt.cigartuples,
        read_rt.cigarstring,
    )

    lt_exons, lt_introns = read_lt.get_exons_and_introns()
    rt_exons, rt_introns = read_rt.get_exons_and_introns()

    logger.trace(f"{lt_chrm=} {rt_chrm=}")
    logger.trace(f"{lt_strand=} {rt_strand=}")
    logger.trace(f"{lt_mode=} {rt_mode=}")

    if lt_chrm == rt_chrm:
        if lt_strand == rt_strand:  # deletion, insertion, duplication
            if lt_mode == 2 and rt_mode == 1:
                target_start = read_rt.ref_start
                target_end = read_lt.ref_end
                target_offset = target_end - target_start

                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.rt_soft_len
                    - read_rt.lt_soft_len
                    - read_lt.read_match_size
                    - read_rt.read_match_size
                )

                logger.trace(f"{bp_region_seq_len=}")

                if bp_region_seq_len > 0:
                    query_offset = (
                        read_lt.reference_match_size + read_rt.reference_match_size
                    )
                else:
                    query_offset = (
                        read_lt.reference_match_size
                        + read_rt.reference_match_size
                        + bp_region_seq_len
                    )

                lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
                rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)

                evt_size = query_offset - target_offset

                logger.trace(f"{evt_size=}, {query_offset=}")
                if evt_size <= 0:  # deletion

                    del_start = read_rt.ref_start + read_rt.reference_match_size
                    del_end = del_start + abs(evt_size)
                    _, _anno, _can = splicing_confirmation(
                        lt_chrm,
                        del_start,
                        lt_chrm,
                        del_end,
                        splice_bin,
                        genome_fasta,
                        cvg,
                        False,
                        motif_required,
                    )
                    _genes = gene_annotation(
                        lt_chrm, del_start, lt_chrm, del_end, gene_iv
                    )
                    return (
                        "DEL",
                        _anno,
                        _can,
                        (
                            f"{lt_chrm}:{del_start}",
                            f"{lt_chrm}:{del_end}",
                            1,
                            2,
                        ),
                        (rt_start, rt_end, rt_exons),
                        (lt_start, lt_end, lt_exons),
                        (rt_bp_seq, lt_bp_seq),
                        (rt_strand, lt_strand),
                        [*_genes],
                    )

                # reads length < tandem duplication size
                elif evt_size >= query_offset:  # large tandem duplication
                    chrm_start = lt_chrm
                    junc_start = read_lt.ref_start
                    chrm_end = lt_chrm
                    junc_end = junc_start + evt_size
                    _nls, _anno, _can = splicing_confirmation(
                        chrm_start,
                        junc_start,
                        chrm_end,
                        junc_end,
                        splice_bin,
                        genome_fasta,
                        cvg,
                        False,
                        motif_required,
                    )
                    _genes = gene_annotation(
                        chrm_start, junc_start, chrm_end, junc_end, gene_iv
                    )
                    if _nls:
                        return (
                            "TDUP",
                            _anno,
                            _can,
                            (
                                f"{lt_chrm}:{junc_start}",
                                f"{lt_chrm}:{junc_end}",
                                2,
                                1,
                            ),
                            (lt_start, lt_end, lt_exons),
                            (rt_start, rt_end, rt_exons),
                            (lt_bp_seq, rt_bp_seq),
                            (lt_strand, rt_strand),
                            [*_genes],
                        )
                    else:
                        return noreturn
                # read length > tandem duplication size
                else:
                    # softclipped length < tandem duplication size
                    if softclipped_length_and_event_size_checker(
                        read_lt, lt_mode, evt_size, bp_region_seq_len
                    ):
                        logger.trace("softclipped length < event size: TDUP")
                        is_dup = True
                    # softclipped length >= tandem duplication size
                    # Novel sequence insertion OR TDUP
                    else:
                        is_dup = False
                        ins_start = read_lt.ref_start
                        ref_allele = genome_fasta[lt_chrm][
                            ins_start : ins_start + 1
                        ].seq
                        ins_seq_in_read = obtain_ins_seq_from_softclipped_part_read(
                            read_lt, lt_mode, evt_size, bp_region_seq_len
                        )

                        if short_tdup_or_not(
                            lt_mode, read_rt, bp_region_seq_len, ins_seq_in_read, logger
                        ):
                            logger.trace("softclipped length >= event size: TDUP")
                            is_dup = True
                        else:
                            logger.trace("softclipped length >= event size: INS")
                            is_dup = False
                    if is_dup:
                        chrm_start = lt_chrm
                        junc_start = read_lt.ref_start
                        chrm_end = lt_chrm
                        junc_end = junc_start + evt_size
                        _nls, _anno, _can = splicing_confirmation(
                            chrm_start,
                            junc_start,
                            chrm_end,
                            junc_end,
                            splice_bin,
                            genome_fasta,
                            cvg,
                            False,
                            motif_required,
                        )
                        _genes = gene_annotation(
                            chrm_start, junc_start, chrm_end, junc_end, gene_iv
                        )
                        if _nls:
                            return (
                                "TDUP",
                                _anno,
                                _can,
                                (
                                    f"{lt_chrm}:{junc_start}",
                                    f"{lt_chrm}:{junc_end}",
                                    2,
                                    1,
                                ),
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_bp_seq, rt_bp_seq),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                        else:
                            return noreturn
                    else:  # it's a short insertion
                        _genes = gene_annotation(
                            lt_chrm, ins_start, lt_chrm, ins_start, gene_iv
                        )
                        return (
                            "INS",
                            ref_allele,
                            ins_seq_in_read,
                            (f"{lt_chrm}:{ins_start}", len(ins_seq_in_read), 2, 1),
                            (lt_start, lt_end, lt_exons),
                            (rt_start, rt_end, rt_exons),
                            (lt_bp_seq, rt_bp_seq),
                            (lt_strand, rt_strand),
                            [*_genes],
                        )
            elif lt_mode == 1 and rt_mode == 2:
                target_start = read_lt.ref_start
                target_end = read_rt.ref_start + read_rt.reference_match_size
                target_offset = target_end - target_start
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.lt_soft_len
                    - read_rt.rt_soft_len
                    - read_lt.read_match_size
                    - read_rt.read_match_size
                )

                logger.trace(f"{bp_region_seq_len=}")

                if bp_region_seq_len > 0:
                    query_offset = (
                        read_lt.reference_match_size + read_rt.reference_match_size
                    )
                else:
                    query_offset = (
                        read_lt.reference_match_size
                        + read_rt.reference_match_size
                        + bp_region_seq_len
                    )

                lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
                rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
                evt_size = query_offset - target_offset

                logger.trace(f"{evt_size=}, {query_offset=}")

                if evt_size <= 0:  # deletion
                    del_start = read_lt.ref_start + read_lt.reference_match_size
                    del_end = del_start + abs(evt_size)
                    _, _anno, _can = splicing_confirmation(
                        lt_chrm,
                        del_start,
                        lt_chrm,
                        del_end,
                        splice_bin,
                        genome_fasta,
                        cvg,
                        False,
                        motif_required,
                    )
                    _genes = gene_annotation(
                        lt_chrm, del_start, lt_chrm, del_end, gene_iv
                    )
                    return (
                        "DEL",
                        _anno,
                        _can,
                        (
                            f"{lt_chrm}:{del_start}",
                            f"{lt_chrm}:{del_end}",
                            1,
                            2,
                        ),
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                elif evt_size >= query_offset:
                    chrm_start = rt_chrm
                    junc_start = read_rt.ref_start
                    chrm_end = rt_chrm
                    junc_end = junc_start + evt_size
                    _nls, _anno, _can = splicing_confirmation(
                        chrm_start,
                        junc_start,
                        chrm_end,
                        junc_end,
                        splice_bin,
                        genome_fasta,
                        cvg,
                        False,
                        motif_required,
                    )
                    _genes = gene_annotation(
                        chrm_start, junc_start, chrm_end, junc_end, gene_iv
                    )
                    if _nls:
                        return (
                            "TDUP",
                            _anno,
                            _can,
                            (
                                f"{rt_chrm}:{junc_start}",
                                f"{rt_chrm}:{junc_end}",
                                2,
                                1,
                            ),
                            (rt_start, rt_end, rt_exons),
                            (lt_start, lt_end, lt_exons),
                            (rt_bp_seq, lt_bp_seq),
                            (rt_strand, lt_strand),
                            [*_genes],
                        )
                    else:
                        return noreturn
                # read length > tandem duplication size
                else:
                    # softclipped length < tandem duplication size
                    if softclipped_length_and_event_size_checker(
                        read_lt, lt_mode, evt_size, bp_region_seq_len
                    ):
                        logger.trace("softclipped length < event size: TDUP")
                        is_dup = True
                    # softclipped length >= tandem duplication size
                    # Novel sequence insertion OR TDUP
                    else:
                        is_dup = False
                        ins_start = read_lt.ref_start + read_lt.reference_match_size
                        ref_allele = genome_fasta[lt_chrm][
                            ins_start : ins_start + 1
                        ].seq
                        ins_seq_in_read = obtain_ins_seq_from_softclipped_part_read(
                            read_lt, lt_mode, evt_size, bp_region_seq_len
                        )

                        if short_tdup_or_not(
                            lt_mode, read_rt, bp_region_seq_len, ins_seq_in_read, logger
                        ):
                            logger.trace("softclipped length >= event size: TDUP")
                            is_dup = True
                        else:
                            logger.trace("softclipped length >= event size: INS")
                            is_dup = False
                    if is_dup:
                        chrm_start = rt_chrm
                        junc_start = rt_start
                        chrm_end = rt_chrm
                        junc_end = junc_start + evt_size
                        _nls, _anno, _can = splicing_confirmation(
                            chrm_start,
                            junc_start,
                            chrm_end,
                            junc_end,
                            splice_bin,
                            genome_fasta,
                            cvg,
                            False,
                            motif_required,
                        )
                        _genes = gene_annotation(
                            chrm_start, junc_start, chrm_end, junc_end, gene_iv
                        )
                        if _nls:
                            return (
                                "TDUP",
                                _anno,
                                _can,
                                (
                                    f"{rt_chrm}:{junc_start}",
                                    f"{rt_chrm}:{junc_end}",
                                    2,
                                    1,
                                ),
                                (rt_start, rt_end, rt_exons),
                                (lt_start, lt_end, lt_exons),
                                (rt_bp_seq, lt_bp_seq),
                                (rt_strand, lt_strand),
                                [*_genes],
                            )
                        else:
                            return noreturn
                    # it is a short insertion
                    else:
                        _genes = gene_annotation(
                            rt_chrm, ins_start, rt_chrm, ins_start, gene_iv
                        )
                        return (
                            "INS",
                            ref_allele,
                            ins_seq_in_read,
                            (f"{rt_chrm}:{ins_start}", len(ins_seq_in_read), 1, 2),
                            (rt_start, rt_end, rt_exons),
                            (lt_start, lt_end, lt_exons),
                            (rt_bp_seq, lt_bp_seq),
                            (rt_strand, lt_strand),
                            [*_genes],
                        )
            else:
                return noreturn
        else:  # lt_strand != rt_strand
            if lt_mode == rt_mode == 1:
                ra_bp = read_lt.ref_start + read_lt.reference_match_size
                sa_bp = read_rt.ref_start + read_rt.reference_match_size
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.lt_soft_len
                    - read_rt.lt_soft_len
                    - read_lt.read_match_size
                    - read_rt.read_match_size
                )
                logger.trace(f"{bp_region_seq_len=}")

                if ra_bp == sa_bp:  # inverted duplication (IDUP)
                    chrm_start = lt_chrm
                    junc_start = ra_bp
                    chrm_end = lt_chrm
                    junc_end = ra_bp
                    _nls, _anno, _can = splicing_confirmation(
                        chrm_start,
                        junc_start,
                        chrm_end,
                        junc_end,
                        splice_bin,
                        genome_fasta,
                        cvg,
                        True,
                        motif_required,
                    )
                    strands = (lt_strand, rt_strand)
                    lt_start_end_exons = (lt_start, lt_end, lt_exons)
                    rt_start_end_exons = (rt_start, rt_end, rt_exons)
                    lt_bp_seq = obtain_bp_region_seq(
                        read_lt, lt_mode, bp_region_seq_len
                    )
                    rt_bp_seq = obtain_bp_region_seq(
                        read_rt, rt_mode, bp_region_seq_len
                    )
                    if _nls:
                        _genes = gene_annotation(
                            chrm_start, junc_start, chrm_end, junc_end, gene_iv
                        )
                        return (
                            "IDUP",
                            _anno,
                            _can,
                            (
                                f"{lt_chrm}:{junc_start}",
                                f"{lt_chrm}:{junc_end}",
                                1,
                                1,
                            ),
                            lt_start_end_exons,
                            rt_start_end_exons,
                            (lt_bp_seq, rt_bp_seq),
                            tuple([*strands]),
                            [*_genes],
                        )
                    else:
                        return noreturn
                else:
                    chrm_start = lt_chrm
                    junc_start = min(ra_bp, sa_bp)
                    chrm_end = lt_chrm
                    junc_end = junc_start + abs(ra_bp - sa_bp)
                    _nls, _anno, _can = splicing_confirmation(
                        chrm_start,
                        junc_start,
                        chrm_end,
                        junc_end,
                        splice_bin,
                        genome_fasta,
                        cvg,
                        True,
                        motif_required,
                    )

                    if junc_start == ra_bp:
                        strands = (lt_strand, rt_strand)
                        lt_start_end_exons = (lt_start, lt_end, lt_exons)
                        rt_start_end_exons = (rt_start, rt_end, rt_exons)
                        lt_bp_seq = obtain_bp_region_seq(
                            read_lt, lt_mode, bp_region_seq_len
                        )
                        rt_bp_seq = obtain_bp_region_seq(
                            read_rt, rt_mode, bp_region_seq_len
                        )
                    elif junc_start == sa_bp:
                        strands = (rt_strand, lt_strand)
                        lt_start_end_exons = (rt_start, rt_end, rt_exons)
                        rt_start_end_exons = (lt_start, lt_end, lt_exons)
                        lt_bp_seq = obtain_bp_region_seq(
                            read_rt, rt_mode, bp_region_seq_len
                        )
                        rt_bp_seq = obtain_bp_region_seq(
                            read_lt, lt_mode, bp_region_seq_len
                        )
                    _genes = gene_annotation(
                        chrm_start, junc_start, chrm_end, junc_end, gene_iv
                    )
                    if _nls:
                        # check whether the chimeric alignments uses
                        # canonical splice sites or not (60% fraction by default)
                        if read_lt.splice_site_checker(
                            genome_fasta
                        ) and read_rt.splice_site_checker(genome_fasta):
                            return (
                                "INV",
                                _anno,
                                _can,
                                (
                                    f"{lt_chrm}:{junc_start}",
                                    f"{lt_chrm}:{junc_end}",
                                    1,
                                    1,
                                ),
                                lt_start_end_exons,
                                rt_start_end_exons,
                                (lt_bp_seq, rt_bp_seq),
                                tuple([*strands]),
                                [*_genes],
                            )
                        else:
                            return noreturn
                    else:
                        return noreturn
            elif lt_mode == rt_mode == 2:  # inversion
                ra_bp = read_lt.ref_start
                sa_bp = read_rt.ref_start
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.rt_soft_len
                    - read_rt.rt_soft_len
                    - read_lt.read_match_size
                    - read_rt.read_match_size
                )

                logger.trace(f"{bp_region_seq_len=}")

                if ra_bp == sa_bp:  # inverted duplication (IDUP)
                    chrm_start = lt_chrm
                    junc_start = ra_bp
                    chrm_end = lt_chrm
                    junc_end = ra_bp
                    _nls, _anno, _can = splicing_confirmation(
                        chrm_start,
                        junc_start,
                        chrm_end,
                        junc_end,
                        splice_bin,
                        genome_fasta,
                        cvg,
                        True,
                        motif_required,
                    )
                    strands = (lt_strand, rt_strand)
                    lt_start_end_exons = (lt_start, lt_end, lt_exons)
                    rt_start_end_exons = (rt_start, rt_end, rt_exons)
                    lt_bp_seq = obtain_bp_region_seq(
                        read_lt, lt_mode, bp_region_seq_len
                    )
                    rt_bp_seq = obtain_bp_region_seq(
                        read_rt, rt_mode, bp_region_seq_len
                    )
                    if _nls:
                        _genes = gene_annotation(
                            chrm_start, junc_start, chrm_end, junc_end, gene_iv
                        )
                        return (
                            "IDUP",
                            _anno,
                            _can,
                            (
                                f"{lt_chrm}:{junc_start}",
                                f"{lt_chrm}:{junc_end}",
                                2,
                                2,
                            ),
                            lt_start_end_exons,
                            rt_start_end_exons,
                            (lt_bp_seq, rt_bp_seq),
                            tuple([*strands]),
                            [*_genes],
                        )
                    else:
                        return noreturn
                else:
                    chrm_start = lt_chrm
                    junc_start = min(ra_bp, sa_bp)
                    chrm_end = lt_chrm
                    junc_end = junc_start + abs(ra_bp - sa_bp)
                    _nls, _anno, _can = splicing_confirmation(
                        chrm_start,
                        junc_start,
                        chrm_end,
                        junc_end,
                        splice_bin,
                        genome_fasta,
                        cvg,
                        True,
                        motif_required,
                    )
                    if junc_start == ra_bp:
                        strands = (lt_strand, rt_strand)
                        lt_start_end_exons = (lt_start, lt_end, lt_exons)
                        rt_start_end_exons = (rt_start, rt_end, rt_exons)
                        lt_bp_seq = obtain_bp_region_seq(
                            read_lt, lt_mode, bp_region_seq_len
                        )
                        rt_bp_seq = obtain_bp_region_seq(
                            read_rt, rt_mode, bp_region_seq_len
                        )
                    elif junc_start == sa_bp:
                        strands = (rt_strand, lt_strand)
                        lt_start_end_exons = (rt_start, rt_end, rt_exons)
                        rt_start_end_exons = (lt_start, lt_end, lt_exons)
                        lt_bp_seq = obtain_bp_region_seq(
                            read_rt, rt_mode, bp_region_seq_len
                        )
                        rt_bp_seq = obtain_bp_region_seq(
                            read_lt, lt_mode, bp_region_seq_len
                        )
                    _genes = gene_annotation(
                        chrm_start, junc_start, chrm_end, junc_end, gene_iv
                    )
                    if _nls:
                        # check whether the chimeric alignments
                        # uses canonical splice sites or not (60% fraction by default)
                        if read_lt.splice_site_checker(
                            genome_fasta
                        ) and read_rt.splice_site_checker(genome_fasta):
                            return (
                                "INV",
                                _anno,
                                _can,
                                (
                                    f"{lt_chrm}:{junc_start}",
                                    f"{lt_chrm}:{junc_end}",
                                    2,
                                    2,
                                ),
                                lt_start_end_exons,
                                rt_start_end_exons,
                                (lt_bp_seq, rt_bp_seq),
                                tuple([*strands]),
                                [*_genes],
                            )
                        else:
                            return noreturn
                    else:
                        return noreturn
            else:
                return noreturn
    else:  # lt_chrm != rt_chrm
        if lt_strand == rt_strand:
            if lt_mode == 1 and rt_mode == 2:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start + read_lt.reference_match_size
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.lt_soft_len
                    - read_rt.rt_soft_len
                    - read_lt.read_match_size
                    - read_rt.read_match_size
                )

                logger.trace(f"{bp_region_seq_len=}")

                lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
                rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
                _nls, _anno, _can = splicing_confirmation(
                    chrm_start,
                    junc_start,
                    chrm_end,
                    junc_end,
                    splice_bin,
                    genome_fasta,
                    cvg,
                    False,
                    motif_required,
                )
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                if _nls:
                    return (
                        "TRA",
                        _anno,
                        _can,
                        (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 1, 2),
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                else:
                    return noreturn
            elif lt_mode == 2 and rt_mode == 1:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start + read_rt.reference_match_size
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.rt_soft_len
                    - read_rt.lt_soft_len
                    - read_lt.read_match_size
                    - read_rt.read_match_size
                )
                logger.trace(f"{bp_region_seq_len=}")
                lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
                rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
                _nls, _anno, _can = splicing_confirmation(
                    chrm_start,
                    junc_start,
                    chrm_end,
                    junc_end,
                    splice_bin,
                    genome_fasta,
                    cvg,
                    False,
                    motif_required,
                )
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                if _nls:
                    return (
                        "TRA",
                        _anno,
                        _can,
                        (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 2, 1),
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                else:
                    return noreturn
            else:
                return noreturn
        else:  # lt_strand != rt_strand
            if lt_mode == rt_mode == 1:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start + read_lt.reference_match_size
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start + read_rt.reference_match_size
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.lt_soft_len
                    - read_rt.lt_soft_len
                    - read_lt.read_match_size
                    - read_rt.read_match_size
                )
                logger.trace(f"{bp_region_seq_len=}")
                lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
                rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
                _nls, _anno, _can = splicing_confirmation(
                    chrm_start,
                    junc_start,
                    chrm_end,
                    junc_end,
                    splice_bin,
                    genome_fasta,
                    cvg,
                    True,
                    motif_required,
                )
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                if _nls:
                    return (
                        "TRA",
                        _anno,
                        _can,
                        (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 1, 1),
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                else:
                    return noreturn
            elif lt_mode == rt_mode == 2:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.rt_soft_len
                    - read_rt.rt_soft_len
                    - read_lt.read_match_size
                    - read_rt.read_match_size
                )

                logger.trace(f"{bp_region_seq_len=}")
                lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
                rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
                _nls, _anno, _can = splicing_confirmation(
                    chrm_start,
                    junc_start,
                    chrm_end,
                    junc_end,
                    splice_bin,
                    genome_fasta,
                    cvg,
                    True,
                    motif_required,
                )
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                if _nls:
                    return (
                        "TRA",
                        _anno,
                        _can,
                        (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 2, 2),
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                else:
                    return noreturn
            else:
                return noreturn
