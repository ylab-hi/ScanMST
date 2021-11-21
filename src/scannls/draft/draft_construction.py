#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ===========================================================
import argparse
import os
import re
import subprocess
import sys
import textwrap
import time
from collections import defaultdict

from loguru import logger
from pyfaidx import Fasta
from pyfaidx import FastaNotFoundError

from .. import __version__
from ..classes import Blat
from ..classes import LengthAction
from ..classes import Read
from ..classes import ReadsConnecter
from ..classes import Series
from ..common import get_softclip_length
from ..externals import blat2chimeric_alignment
from ..externals import external_tool_checking
from ..utils import reverse_complement
from .helper import extract_splice_sites

try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)

__funcs__ = {"detect_sv_from_cigar", "scan_bam"}


def detect_sv_from_cigar(
    *,
    read,
    mapq_cutoff,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    blat,
    logger,
    update_bps=False,
) -> list:
    """
    :param logger: logger for logging
    :param blat: `class.Blat`
    :param update_bps:
    :param read: A read from pysam.AlignedSegment
    :param mapq_cutoff: MAPQ cutoff
    :param splice_bin: a small bin for splice site searching
    :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
    :param cvg: annotated splice sites (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param gene_iv: annotated gene region (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param motif_required: considering canonical splice sites only OR considering both canonical and noncanonical splice sites
    :type read: pysam.AlignedSegment
    :type mapq_cutoff: int
    :type splice_bin: int
    :type genome_fasta: pyfaidx.Fasta
    :type cvg: HTSeq.GenomicArrayOfSets
    :type gene_iv: HTSeq.GenomicArrayOfSets
    :type motif_required: bool
    :return: event groups in a list, every group is also a list
    :rtype: list (list of lists)
    """
    (
        read_to_read_chains,
        reads_pair_mode_dict,
        insertion_dict,
    ) = detect_read_read_connections_from_cigar(
        read=read,
        mapq_cutoff=mapq_cutoff,
        blat=blat,
        logger=logger,
    )

    read_to_read_chains = [read_to_read_chains]

    logger.debug("Read-to-Read chain: ", read_to_read_chains)
    logger.debug("Read-to-Read pair modes: ", reads_pair_mode_dict)

    event_groups = []
    if read_to_read_chains:
        # every chain is a group of connected reads
        # every chain may have a list of events
        for chain in read_to_read_chains:
            event_list = []
            for _lt, _rt in zip(chain[::1], chain[1::1]):
                # print(_lt, _rt)
                if (_lt, _rt) in reads_pair_mode_dict:
                    _lt_mode, _rt_mode = reads_pair_mode_dict[(_lt, _rt)]
                elif (_rt, _lt) in reads_pair_mode_dict:
                    _rt_mode, _lt_mode = reads_pair_mode_dict[(_rt, _lt)]
                # print(_lt_mode, _rt_mode)
                (
                    nls_type,
                    _anno,
                    _canonical,
                    positions,
                    lt_info,
                    rt_info,
                    strands,
                    genes,
                ) = infer_nls_from_connected_reads(
                    _lt,
                    _rt,
                    _lt_mode,
                    _rt_mode,
                    splice_bin,
                    genome_fasta,
                    cvg,
                    gene_iv,
                    motif_required,
                    update_bps,
                )

                if nls_type != "NA":
                    event_list.append(
                        (
                            nls_type,
                            _anno,
                            _canonical,
                            positions,
                            lt_info,
                            rt_info,
                            strands,
                            genes,
                        )
                    )
            if event_list:
                event_groups.append(event_list)
    for group in event_groups:
        for i in group:
            print(i)
    return event_groups


def scan_bam(
    input_bam,
    mapq_cutoff,
    output,
    ref_genome,
    gtf,
    splice_bin,
    blat,
    logger,
    motif_required,
    blat_ident_pct_cutoff=0.9,
    max_allowed_nm=50,
    min_soft_seg_len=200,
):
    """(1) update CIGAR strings of supplementary alignments in the primary alignment SA tag.
       (2) add SA tag for reads with long length of softclipped segment using BLAT
       (3) identify putative regions of NLS events using connected chimeric reads
       (4) add putative regions of NLS events to SV tag of primary alignment
       (5) output regions of NLS events in BEDPE file

    :param blat:
    :param logger:
    :param input_bam: Transcriptomic long-read sorted BAM file
    :param mapq_cutoff: MAPQ cutoff
    :param output: file full name for output rebuild BAM file
    :param ref_genome: reference genome (FASTA file)
    :param gtf: reference gene annotations (GTF file)
    :param splice_bin: bin size for splice site searching
    :param ref_2bit: reference 2bit file for BLAT
    :param motif_required: canonical splice sites required; if True: considering canonical splice sites only; else: considering canonical and noncanonical splice sites both
    :param port: BLAT server port
    :param output_dir: BLAT output directory for psl files
    :param blat_ident_pct_cutoff: BLAT HSP identity cutoff
    :param max_allowed_nm: mismatches cutoff used for discarding supplementary alignments
    :param min_soft_seg_len: minium softclipped segement length to trigger BLAT for reads with softcliping but no SA tag
    :type input_bam: str (BAM filename)
    :type mapq_cutoff: int
    :type output: str
    :type ref_genome: str
    :type gtf: str
    :type splice_bin: int
    :type ref_2bit: str
    :type motif_required: bool
    :type port: int
    :type output_dir: str
    :type blat_ident_pct_cutoff: float
    :type max_allowed_nm: int
    :type min_soft_seg_len: int
    :return: No returns
    :rtype: None
    ..note ::
        SV tag uses the same genomic corrdinate as SA tag,
        So position should be always add 1
    """
    logger.debug("scan_bam")

    in_bam = pysam.AlignmentFile(input_bam, "rb")
    output_bam = pysam.AlignmentFile(f"{output}", "wb", template=in_bam)

    candidate_ao_dict = defaultdict(int)
    nls_src_forms_list = []

    if "~" in ref_genome:
        ref_genome = os.path.expanduser(ref_genome)
    if "~" in gtf:
        gtf = os.path.expanduser(gtf)
    try:
        genome_fasta = Fasta(ref_genome, sequence_always_upper=True)
    except FastaNotFoundError as e:
        print("read reference genome " + ref_genome + " error!", e)
        sys.exit(1)
    try:
        cvg, gene_iv = extract_splice_sites(gtf, splice_bin)
    except IOError as e:
        print("read GTF file " + gtf + " error!", e)
        sys.exit(1)

    # supplementary alignment cigarstring extraction
    # key: read.query_name + left S + right S
    # For minimap2, "-Y" need to be used, use soft clipping for supplementary alignments
    representative_alignments_new_cigar = {}
    pat_left_S = re.compile(r"^(\d+)S")
    pat_right_S = re.compile(r"(\d+)S$")
    try:
        for read in in_bam.fetch(until_eof=False):
            if read.is_supplementary:
                sup_aln_cigar = read.cigarstring
                left_mat = pat_left_S.search(sup_aln_cigar)
                right_mat = pat_right_S.search(sup_aln_cigar)
                if left_mat:
                    l_S_len = left_mat.group(1)
                else:
                    l_S_len = ""
                if right_mat:
                    r_S_len = right_mat.group(1)
                else:
                    r_S_len = ""
                representative_alignments_new_cigar[
                    "{}\t{}\t{}".format(read.qname, l_S_len, r_S_len)
                ] = sup_aln_cigar
    except ValueError as e:
        print(
            "BAM index file is not found in supplementary alignments!\n",
            e,
            file=sys.stderr,
        )
        sys.exit(1)

    # update SA tags and iterate the BAM file
    # in_bam = pysam.AlignmentFile(input_bam, "rb")
    try:
        for read in in_bam.fetch(until_eof=False):
            if (
                read.mapq >= mapq_cutoff
                and not read.is_secondary
                and not read.has_tag("XA")
                and not read.is_unmapped
            ):
                chrm = read.reference_name
                # update SA tag of representative alignments (START)
                if read.has_tag("SA") and not read.is_supplementary:
                    updated_chimeric_alns = []
                    chimeric_alns = read.get_tag("SA")[:-1].split(";")
                    # print(read.get_tag('SA'))
                    # one representative alignment could have multiple corresponding supplementary alignments
                    for _aln in chimeric_alns:
                        (
                            __chr_sa,
                            __pos_sa,
                            __strand_sa,
                            __cigar_sa,
                            __mapq_sa,
                            __nm_sa,
                        ) = _aln.split(",")
                        left_mat = pat_left_S.search(__cigar_sa)
                        right_mat = pat_right_S.search(__cigar_sa)

                        l_S_len = left_mat.group(1) if left_mat else ""
                        r_S_len = right_mat.group(1) if right_mat else ""

                        tgt_key = "{}\t{}\t{}".format(read.qname, l_S_len, r_S_len)
                        if tgt_key in representative_alignments_new_cigar:
                            __updated_cigar = representative_alignments_new_cigar[
                                tgt_key
                            ]
                            # discard supplementary alignments with too many mismatches or lower MAPQ
                            if not (
                                int(__nm_sa) > max_allowed_nm
                                or int(__mapq_sa) < mapq_cutoff
                            ):
                                updated_chimeric_alns.append(
                                    "{},{},{},{},{},{}".format(
                                        __chr_sa,
                                        __pos_sa,
                                        __strand_sa,
                                        __updated_cigar,
                                        __mapq_sa,
                                        __nm_sa,
                                    )
                                )
                    if len(updated_chimeric_alns) == 0:
                        read.set_tag("SA", None)
                    else:
                        read.set_tag(
                            "SA", "{};".format(";".join(updated_chimeric_alns))
                        )
                    # remove SA tags of representative alignments with too much mismatches
                # update SA tag of representative alignments (END)

                # Detect novel chimeric alignments for reads with long softclipped segment but without SA tags using BLAT
                if not read.has_tag("SA") and not read.is_supplementary:
                    read_strand = "-" if read.is_reverse else "+"
                    read_length = int(read.query_length)
                    # assert read.cigarstring, f"{read.query_name}" # TEST
                    _, _soft_seq, _, read_mode = get_softclip_length(read)

                    if read.is_reverse:
                        soft_seq_ori = reverse_complement(_soft_seq)
                    else:
                        soft_seq_ori = _soft_seq

                    if (
                        read_mode in {1, 2}
                        and soft_seq_ori
                        and len(soft_seq_ori) >= min_soft_seg_len
                    ):
                        chimeric_aln_str = blat2chimeric_alignment(
                            soft_seq_ori,
                            read_length,
                            read_strand,
                            read_mode,
                            blat,
                            mapq_cutoff,
                            max_allowed_nm,
                            blat_ident_pct_cutoff,
                        )
                        if chimeric_aln_str:
                            read.set_tag("SA", chimeric_aln_str)
                # _anno:annotated exon boundary (0/1/2); _can: canonical_or_not(1/0);
                # newpos=[pos,size/pos2_of_translocation, rep_aln_mode, sup_aln_mode]

                # select reads with SA tags (original or newly-added), ignore supplementary alignment
                if read.has_tag("SA") and not read.is_supplementary:
                    event_groups = detect_sv_from_cigar(
                        read=read,
                        mapq_cutoff=mapq_cutoff,
                        splice_bin=splice_bin,
                        genome_fasta=genome_fasta,
                        cvg=cvg,
                        gene_iv=gene_iv,
                        motif_required=motif_required,
                        blat=blat,
                        logger=logger,
                    )
                    sv_tag_list = []
                    ot_tag_list = []
                    for group in event_groups:
                        nls_event_list = []
                        for event in group:
                            (
                                _type,
                                _anno,
                                _canonical,
                                _positions,
                                read1_info,
                                read2_info,
                                strands,
                                genes,
                            ) = event
                            _bp1, _bp2, _mode1, _mode2 = _positions
                            _strand1, _strand2 = strands
                            _gene1, _gene2 = genes
                            if _type in {"TDUP", "INV", "TRA"}:
                                _chrm1, _pos1 = _bp1.split(":")
                                _chrm2, _pos2 = _bp2.split(":")
                                # SV tag uses SA tag corrdinate system (start with 1)
                                # So, position should always add 1
                                sv_tag_list.append(
                                    f"{_type},{_anno}|{_canonical},{_chrm1}:{int(_pos1) + 1},{_chrm2}:{int(_pos2) + 1},{_mode1}{_mode2},{_strand1}{_strand2},{_gene1}|{_gene2};"
                                )
                                nls_event_list.append(event)

                                event_key = f"{_type}\t{_canonical}\t{_chrm1}:{int(_pos1) + 1}\t{_chrm2}:{int(_pos2) + 1}\t{_strand1}{_strand2}"
                                reversed_event_key = f"{_type}\t{_canonical}\t{_chrm2}:{int(_pos2) + 1}\t{_chrm1}:{int(_pos1) + 1}\t{_strand2}{_strand1}"
                                if event_key in candidate_ao_dict:
                                    candidate_ao_dict[event_key] += 1
                                elif reversed_event_key in candidate_ao_dict:
                                    candidate_ao_dict[reversed_event_key] += 1

                                # candidate_group_dict[
                                #    f"{_type}\t{_canonical}\t{_chrm1}:{int(_pos1)+1}\t{_chrm2}:{int(_pos2)+1}\t{_strand1}{_strand2}"
                                # ] = group_counter
                            elif _type in {"INS"}:
                                _end_pos = int(_bp1) + int(_bp2)
                                # 'INS', ref_allele, ins_seq_in_read, [ins_start, len(ins_seq_in_read), 1, 2], [lt_strand, rt_strand], [*_genes]
                                ot_tag_list.append(
                                    f"{_type},{_anno}|{_canonical},{_bp1},{_end_pos},{_mode1}{_mode2},{_strand1}{_strand2},{_gene1}|{_gene2};"
                                )

                                candidate_ao_dict[
                                    f"{_type}\t{_canonical}\t{chrm}:{_bp1}\t{chrm}:{_end_pos}\t{_strand1}{_strand2}"
                                ] += 1

                        if nls_event_list:
                            # print("nls_event_list: ", len(nls_event_list))
                            bp_series = Series()
                            bp_series.init(nls_event_list)
                            nls_src_forms_list.append(bp_series)

                    if sv_tag_list:
                        read.set_tag("SV", "".join(sv_tag_list))
                    if ot_tag_list:
                        read.set_tag("OT", "".join(ot_tag_list))

            output_bam.write(read)
    except ValueError as e:
        print("BAM index file is not found!", e, file=sys.stderr)
        # sys.exit(1)
    in_bam.close()
    output_bam.close()

    subprocess.check_call("samtools index {}".format(output), shell=True)
    print("NLS Src forms: ", nls_src_forms_list)
    # output_candidates = aggregate_candidates(candidate_ao_dict, len_cutoff=0)
    prefix = output.split(".")[0]
    # output_bedpe_file(candidate_ao_dict, candidate_group_dict, prefix, splice_bin)
    return None
