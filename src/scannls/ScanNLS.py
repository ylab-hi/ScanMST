#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ===========================================================
"""
2021-10-01:
detect_sv_from_cigar output a list of putative NLS events
modify SV tag endswith ";", SV:Z:XXX;YYY;ZZZ;

"""
import argparse
import copy
import glob
import logging
import math
import os
import random
import re
import subprocess
import sys
import textwrap
import time
from collections import defaultdict
from collections import OrderedDict
from typing import Iterable

from align import aligner
from Bio import SearchIO
from Bio.Seq import Seq
from pyfaidx import Fasta

from . import __version__
from .classes import LengthAction
from .classes import Path
from .classes import Read
from .classes import Series
from .common import remove
from .common import remove_files
from .common import status_message
from .externals import blat_mapq_calculator
from .externals import checkIfProcessRunning
from .externals import external_tool_checking
from .externals import gfClient_query
from .externals import gfserver_tester
from .externals import psl2sam
from .externals import softclipped_seq2SA_tag
from .externals import start_gfServer
from .externals import stop_gfServer
from .utils import aggregate_candidates
from .utils import chimeric_aln_order_finder
from .utils import extract_splice_sites
from .utils import gene_annotation
from .utils import infer_sv_from_connected_reads
from .utils import output_bedpe_file
from .utils import short_TDUP_or_not
from .utils import similar_hit
from .utils import splicing_confirmation
from .utils import test_is_connected
from .utils import update_breakpoints

# from .call import sv_scan

try:
    import pysam
except:
    sys.exit("pysam module not found.\nPlease install it before.")
try:
    import numpy as np
except:
    sys.exit("numpy module not found.\nPlease install it before.")
try:
    import HTSeq
except:
    sys.exit("HTSeq module not found.\nPlease install it before.")
try:
    import skbio
except:
    sys.exit("scikit-bio module not found.\nPlease install it before.")


def detect_read_read_connections_from_cigar(
    chrm, read, mapq_cutoff, allowed_difference
) -> tuple:
    """Detecting read-read connections with chimeric alignments CIGAR string

    :param chrm: chromosome of input read
    :param read: read of pysam.AlignedSegment object, expecting representative read with SA tag
    :param mapq_cutoff: MAPQ cutoff
    :param allowed_difference: the difference of read_match_size (Read1) and softclipped length (Read2) to determine the S-M match
    :type chrm: str
    :type read: pysam.AlignedSegment object
    :type mapq_cutoff: int
    :type allowed_difference: int
    :return: Read-to-Read chain (a list of lists), a dictionary of Read-pair(Read1, Read2) => mode-of-Read1, mode-of-Read2
    :rtype: tuple
    .. note::
        Read-to-Read chain scenarios
        * [[Read1, Read2, Read3]]
        * [[Read1, Read2, Read3],[Read4,Read5]]

        Dictionary of Read-pair scenarios
        * (Read1, Read2) => mode-of-Read1, mode-of-Read2
        * (Read2, Read1) => mode-of-Read2, mode-of-Read1

    .. important::
        If no 'SA' tag is found in this read, read-to-read chain and the read-pair => mode dictionary will become empty.

    #return: NLS_type(TDUP/INV), exon_boundary(0/1/2/3), canonical_or_not (1/0), [position, size, rep_aln_mode, sup_aln_mode], [++]
    #        TRA, canonical_or_not (1/0), [position, sup_position, rep_aln_mode, sup_aln_mode], [+-]
    #        e.g., INV,1,43947377,181934993,1,1,++
    #              TRA,1,160289623,chr17:17189212,1,1,+-
    """

    def format_sa_tag(in_str):
        """
        To keep read.reference_start and start position of SA alignment consistent, start position of SA alignment need to substract 1
        :param in_str: string of supplementary read item in the SA tag
        :type in_str: str
        :return: chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa
        :rtype: tuple
        .. note::
             pos_sa, mapq_sa and nm_sa are integral variables now.
        """
        chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = in_str.split(",")
        pos_sa = int(pos_sa) - 1
        mapq_sa = int(mapq_sa)
        nm_sa = int(nm_sa)
        return chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa

    def obtain_sa_query_seq_from_ra(query_seq_ra, strand_ra, strand_sa):
        """a helper function to define query_seq for the supplementary alignment
        :param query_seq_ra: query sequence of representative alignment
        :type query_seq_ra: str
        :param strand_ra: direction of representative read (-|+)
        :type strand_ra: str
        :param strand_sa: direction of supplementary read (-|+)
        :type strand_sa: str
        :return: query sequence of supplementary alignment
        :rtype: str
        """
        if strand_ra == strand_sa:
            return query_seq_ra
        else:
            __seq = Seq(query_seq_ra)
            return str(__seq.reverse_complement())

    if read.has_tag("SV"):
        return [], {}

    if read.is_supplementary:
        return [], {}

    # if no 'SA' tag was found, read-to-read chain will be empty
    try:
        chimeric_aln = read.get_tag("SA")[:-1].split(";")
    except KeyError:
        return [], {}

    # chimeric alignments for a chimeric read
    # a chimeric read can have multiple chimeric alignments
    chimeric_aln_list = []

    chrm_ra = read.reference_name
    pos_ra = read.reference_start
    if read.is_reverse:
        strand_ra = "-"
    else:
        strand_ra = "+"
    cigar_ra = read.cigarstring
    mapq_ra = read.mapping_quality
    nm_ra = read.get_tag("NM")
    seq_ra = read.query_sequence

    if mapq_ra > mapq_cutoff:
        chimeric_aln_list.append(
            Read.init(chrm_ra, pos_ra, strand_ra, cigar_ra, mapq_ra, nm_ra, seq_ra)
        )

    for sa_string in chimeric_aln:
        chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = format_sa_tag(sa_string)
        seq_sa = obtain_sa_query_seq_from_ra(seq_ra, strand_ra, strand_sa)
        if mapq_sa > mapq_cutoff:
            chimeric_aln_list.append(
                Read.init(chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa, seq_sa)
            )

    # for i in chimeric_aln_list:
    #    print(i, '*', i.reference_match_size,'|', i.sms)

    if not chimeric_aln_list:
        return [], {}
    else:
        read_to_read_chains, reads_pair_mode_dict = chimeric_aln_order_finder(
            chimeric_aln_list, allowed_difference
        )

        # print(read_to_read_chains, reads_pair_mode_dict)
        # print('Read-to-Read chain: ',read_to_read_chains)
        return read_to_read_chains, reads_pair_mode_dict


def detect_sv_from_cigar(
    chrm,
    read,
    mapq_cutoff,
    splice_bin,
    allowed_difference,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    update_bps=False,
) -> list:
    """
    :param chrm: chromosome
    :param read: A read from pysam.AlignedSegment
    :param mapq_cutoff: MAPQ cutoff
    :param splice_bin: a small bin for splice site searching
    :param allowed_difference: the difference of read_match_size (Read1) and softclipped length (Read2) to determine the S-M match
    :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
    :param cvg: annotated splice sites (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param gene_iv: annotated gene region (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param motif_required: considering canonical splice sites only OR considering both canonical and noncanonical splice sites
    :type chrm: str
    :type read: pysam.AlignedSegment
    :type mapq_cutoff: int
    :type splice_bin: int
    :type allowed_difference: int
    :type genome_fasta: pyfaidx.Fasta
    :type cvg: HTSeq.GenomicArrayOfSets
    :type gene_iv: HTSeq.GenomicArrayOfSets
    :type motif_required: bool
    :return: event groups in a list, every group is also a list
    :rtype: list (list of lists)
    """
    read_to_read_chains, reads_pair_mode_dict = detect_read_read_connections_from_cigar(
        chrm, read, mapq_cutoff, allowed_difference
    )

    print("Read-to-Read chain: ", read_to_read_chains)
    print("Read-to-Read pair modes: ", reads_pair_mode_dict)

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
                ) = infer_sv_from_connected_reads(
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


def softclipping_realignment(
    input_bam,
    mapq_cutoff,
    output,
    ref_genome,
    gtf,
    splice_bin,
    ref_2bit,
    motif_required=True,
    blat=False,
    port=88888,
    output_dir="/tmp",
    blat_ident_pct_cutoff=0.9,
    max_allowed_nm=50,
    min_soft_seg_len=200,
    allowed_difference=30,
):
    """(1) update CIGAR strings of supplementary alignments in the primary alignment SA tag.
       (2) add SA tag for reads with long length of softclipped segment using BLAT (Optional)
       (3) identify putative regions of NLS events using connected chimeric reads
       (4) add putative regions of NLS events to SV tag of primary alignment
       (5) output regions of NLS events in BEDPE file

    :param input_bam: Transcriptomic long-read sorted BAM file
    :param mapq_cutoff: MAPQ cutoff
    :param output: file full name for output rebuild BAM file
    :param ref_genome: reference genome (FASTA file)
    :param gtf: reference gene annotations (GTF file)
    :param splice_bin: bin size for splice site searching
    :param ref_2bit: reference 2bit file for BLAT
    :param motif_required: canonical splice sites required; if True: considering canonical splice sites only; else: considering canonical and noncanonical splice sites both
    :param blat: use BLAT OR not
    :param port: BLAT server port
    :param output_dir: BLAT output directory for psl files
    :param blat_ident_pct_cutoff: BLAT HSP identity cutoff
    :param max_allowed_nm: mismatches cutoff used for discarding supplementary alignments
    :param min_soft_seg_len: minium softclipped segement length to trigger BLAT for reads with softcliping but no SA tag
    :param allowed_difference: the difference of read_match_size (Read1) and softclipped length (Read2) to determine the S-M match
    :type input_bam: str (BAM filename)
    :type mapq_cutoff: int
    :type output: str
    :type ref_genome: str
    :type gtf: str
    :type splice_bin: int
    :type ref_2bit: str
    :type motif_required: bool
    :type blat: bool
    :type port: int
    :type output_dir: str
    :type blat_ident_pct_cutoff: float
    :type max_allowed_nm: int
    :type min_soft_seg_len: int
    :type allowed_difference: int
    :return: No returns
    :rtype: None
    ..note ::
        SV tag uses the same genomic corrdinate as SA tag,
        So position should be always add 1
    """
    in_bam = pysam.AlignmentFile(input_bam, "rb")
    output_bam = pysam.AlignmentFile(f"{output}", "wb", template=in_bam)

    candidate_regions = set()
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

    ## supplementary alignment cigarstring extraction
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
    # print(representative_alignments_new_cigar)

    ## update SA tags and iterate the BAM file
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
                        if left_mat:
                            l_S_len = left_mat.group(1)
                        else:
                            l_S_len = ""
                        if right_mat:
                            r_S_len = right_mat.group(1)
                        else:
                            r_S_len = ""
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
                if blat:
                    if not read.has_tag("SA") and not read.is_supplementary:
                        if read.is_reverse:
                            read_strand = "-"
                        else:
                            read_strand = "+"
                        read_length = int(read.query_length)
                        # assert read.cigarstring, f"{read.query_name}" # TEST
                        _, _soft_seq, _, read_mode = get_softclip_length(read)
                        soft_seq = Seq(_soft_seq)
                        if read.is_reverse:
                            soft_seq_ori = str(__soft_seq.reverse_complement())
                        else:
                            soft_seq_ori = str(__soft_seq)
                        if (
                            read_mode in {1, 2}
                            and soft_seq_ori
                            and len(soft_seq_ori) >= min_soft_seg_len
                        ):
                            chimeric_aln_str = softclipped_seq2SA_tag(
                                soft_seq_ori,
                                read_length,
                                read_strand,
                                read_mode,
                                ref_2bit,
                                port,
                                mapq_cutoff,
                                max_allowed_nm,
                                output_dir,
                                blat_ident_pct_cutoff,
                            )
                            if chimeric_aln_str:
                                read.set_tag("SA", chimeric_aln_str)
                # _anno:annotated exon boundary (0/1/2); _can: canonical_or_not(1/0);newpos=[pos, size/pos2_of_translocation, rep_aln_mode, sup_aln_mode]

                # select reads with SA tags (original or newly-added), ignore supplementary alignment
                if read.has_tag("SA") and not read.is_supplementary:
                    event_groups = detect_sv_from_cigar(
                        chrm,
                        read,
                        mapq_cutoff,
                        splice_bin,
                        allowed_difference,
                        genome_fasta,
                        cvg,
                        gene_iv,
                        motif_required,
                    )
                    # TODO
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
                                    f"{_type},{_anno}|{_canonical},{_chrm1}:{int(_pos1)+1},{_chrm2}:{int(_pos2)+1},{_mode1}{_mode2},{_strand1}{_strand2},{_gene1}|{_gene2};"
                                )
                                nls_event_list.append(event)

                                event_key = f"{_type}\t{_canonical}\t{_chrm1}:{int(_pos1)+1}\t{_chrm2}:{int(_pos2)+1}\t{_strand1}{_strand2}"
                                reversed_event_key = f"{_type}\t{_canonical}\t{_chrm2}:{int(_pos2)+1}\t{_chrm1}:{int(_pos1)+1}\t{_strand2}{_strand1}"
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="ScanNLS: Nonlinear splicing (NLS) events identification using transcriptomic long-reads data",
        epilog=textwrap.dedent(
            """Author: Ting-You Wang <tywang@umn.edu>, Hormel Institute, University of Minnesota, 2021"""
        ),
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s {}".format(__version__)
    )
    sub_parsers = parser.add_subparsers(help="sub-command help", dest="sub_command")

    build_parser = sub_parsers.add_parser(
        "build",
        help="add additional tags to build.BAM",
        description="%(prog)s -i input_bam_file -o output_bam_file -r ref_genome_fasta -g gtf_file [opts]",
        epilog=textwrap.dedent(
            """Author: Ting-You Wang <tywang@umn.edu>, Hormel Institute, University of Minnesota, 2021"""
        ),
    )

    build_parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="Input BAM file",
        required=True,
    )
    build_parser.add_argument(
        "-r",
        "--ref",
        action="store",
        dest="ref",
        help="reference genome in FASTA format (with fai index)",
        required=True,
    )
    build_parser.add_argument(
        "-g",
        "--gtf",
        action="store",
        dest="gtf",
        help="gene annotations in GTF format",
        required=True,
    )
    build_parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        help="output BAM file",
        required=True,
    )
    build_parser.add_argument(
        "-s",
        "--splice_bin",
        action="store",
        dest="splice_bin",
        type=int,
        help="minimal observation count for ITD (default: %(default)s)",
        default=5,
    )
    build_parser.add_argument(
        "-m",
        "--mapq",
        action="store",
        dest="mapq",
        type=int,
        help="minimal MAPQ in BAM for calling NLS (default: %(default)s)",
        default=15,
    )
    build_parser.add_argument(
        "-n",
        "--noncanonical",
        action="store_true",
        dest="noncanonical",
        default=False,
        help="Considering Non-canonical spliced sites",
    )
    build_parser.add_argument(
        "-b",
        "--blat",
        action="store_true",
        dest="blat",
        default=False,
        help="Using BLAT to remap softclipped reads (default: %(default)s)",
    )
    build_parser.add_argument(
        "--2bit", action="store", dest="two_bit", help="reference genome in 2bit format"
    )
    build_parser.add_argument(
        "-p",
        "--port",
        action="store",
        dest="port",
        type=int,
        help="port for BLAT server (default: %(default)s)",
        default=88888,
    )
    build_parser.add_argument(
        "--min_soft_seg_len",
        action="store",
        dest="min_soft_seg_len",
        type=int,
        help="minimum softclipped segement length to trigger BLAT alignment (default: %(default)s)",
        default=200,
    )
    build_parser.add_argument(
        "--max_allowed_nm",
        action="store",
        dest="max_allowed_nm",
        type=int,
        help="Maximum allowed NM to keep AS tag (default: %(default)s)",
        default=100,
    )
    build_parser.add_argument(
        "--allowed_difference",
        action="store",
        dest="allowed_difference",
        type=int,
        help="Maximum allowed difference in length between one read matched part and other read softclipped part (default: %(default)s)",
        default=80,
    )
    build_parser.add_argument(
        "--identity",
        action="store",
        dest="ident_cutoff",
        type=float,
        help="blat_ident_pct_cutoff (default: %(default)s)",
        default=0.95,
    )
    build_parser.add_argument(
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
            """Author: Ting-You Wang <tywang@umn.edu>, Hormel Institute, University of Minnesota, 2021"""
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

    infer_parser = sub_parsers.add_parser(
        "infer",
        help="infer NLS isoforms using built BAM and called VCF",
        description="%(prog)s -i input_bam_file_from_ScanNLS_build -o output_vcf_file_prefix [opts]",
        epilog=textwrap.dedent(
            """Author: Ting-You Wang <tywang@umn.edu>, Hormel Institute, University of Minnesota, 2021"""
        ),
    )

    return parser


def main():
    if sys.version_info < (3, 7):
        sys.exit(
            "Sorry, this code need Python 3.7 or higher. Please update. Aborting..."
        )
    parser = parse_args()

    if len(sys.argv[1:]) < 1:
        parser.print_help()
        sys.exit(1)
    else:
        options = parser.parse_args()

    if options.sub_command == "build":
        use_blat = options.blat
        # check external tools used
        external_tool_checking(blat=use_blat)

        print("ScanNLS build starts running: " + time.strftime("%Y-%m-%d %H:%M:%S"))
        start = time.time()

        if use_blat:
            # start gfserver
            print("starting gfserver\n")
            if not checkIfProcessRunning("gfServer"):
                start_gfServer(
                    ref_2bit=options.two_bit,
                    port=options.port,
                    output_dir=options.tmp_dir,
                )

        # CIGAR string refinement or add SV tag
        motif_required = not options.noncanonical
        softclipping_realignment(
            input_bam=options.input,
            mapq_cutoff=options.mapq,
            output=options.output,
            ref_genome=options.ref,
            gtf=options.gtf,
            splice_bin=options.splice_bin,
            ref_2bit=options.two_bit,
            motif_required=motif_required,
            blat=use_blat,
            port=options.port,
            output_dir=options.tmp_dir,
            blat_ident_pct_cutoff=options.ident_cutoff,
            max_allowed_nm=options.max_allowed_nm,
            min_soft_seg_len=options.min_soft_seg_len,
            allowed_difference=options.allowed_difference,
        )

        if use_blat:
            # stop gfserver
            stop_gfServer(port=options.port, output_dir=options.tmp_dir)

        print("ScanNLS build running done: " + time.strftime("%Y-%m-%d %H:%M:%S"))
        end = time.time()
        print("ScanNLS build takes " + str(end - start) + " seconds.")

    elif options.sub_command == "call":
        print(
            "ScanNLS calling NLS events starts running: "
            + time.strftime("%Y-%m-%d %H:%M:%S")
        )
        start = time.time()
        event_dict = sv_scan(
            options.input,
            options.output,
            options.sr,
            options.depth,
            options.pso,
            options.length,
            options.soft_len,
            options.region,
            options.mapq,
            options.mismatch,
            options.alignment_fraction,
            options.seed,
        )
        print(
            "ScanNLS calling NLS events running done: "
            + time.strftime("%Y-%m-%d %H:%M:%S")
        )
        end = time.time()
        print("ScanNLS calling NLS events takes " + str(end - start) + " seconds.")
    # infer transcript forms (GTF) and the corresponding sequences (FASTA)
    elif options.sub_command == "infer":
        pass


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        # stop gfserver
        stop_gfServer(port=options.port, output_dir=options.tmp_dir)
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
