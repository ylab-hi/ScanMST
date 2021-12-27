#!/usr/bin/env python
# ===============================================================================
import argparse
import copy
import glob
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from collections import defaultdict
from collections import OrderedDict
from io import BytesIO
from io import StringIO

import HTSeq
import networkx as nx
import numpy as np
import pandas as pd
import pysam
import skbio
import vcf
from align import aligner
from Bio import AlignIO
from Bio.Application import ApplicationError
from Bio.Seq import Seq
from pyfaidx import Fasta
from varname import nameof

# from Bio.Align.Applications import MuscleCommandline, MafftCommandline, ClustalwCommandline
# from Bio.Align import AlignInfo


def reverse_complement(in_str):
    """
    obtain reverse complement sequence
    """
    my_dna = Seq(in_str)
    return str(my_dna.reverse_complement())


def mismatch_filter(
    anchor_seq, anchor_soft_pos, in_list, alignment_frac, mismatch_cutoff
):
    out_list = []
    anchor_seq = skbio.DNA(anchor_seq)
    for _mode, _seq, _pos, _read in in_list:
        mismatch = 1e6
        if int(_pos) != anchor_soft_pos:
            continue
        each_seq = skbio.DNA(_seq)
        # if seq1 or each_seq is an empty string, ignore it
        if not anchor_seq or not each_seq:
            continue
        try:
            alignment, score, start_end_pos = skbio.alignment.local_pairwise_align_ssw(
                anchor_seq, each_seq
            )
        # raise IndexError if SSW cannot work
        except IndexError:
            continue
        except ValueError:
            continue
        if (
            len(alignment[0]) / float(len(anchor_seq)) < alignment_frac
            and len(alignment[1]) / float(len(each_seq)) < alignment_frac
        ):
            continue
        # Align left
        if start_end_pos[0][0] == 0 and start_end_pos[1][0] == 0:
            if sum(alignment[0].mismatches(alignment[1])) < mismatch:
                mismatch = sum(alignment[0].mismatches(alignment[1]))
        # Align right
        elif (
            start_end_pos[0][1] == len(anchor_seq) - 1
            and start_end_pos[1][1] == len(each_seq) - 1
        ):
            if sum(alignment[0].mismatches(alignment[1])) < mismatch:
                mismatch = sum(alignment[0].mismatches(alignment[1]))
        else:
            continue
        if mismatch < mismatch_cutoff:
            out_list.append((_mode, _seq, _pos, _read))
    return out_list


def get_softclip_length(read):
    """0 => M
    1 => I
    2 => D
    3 => N
    4 => S
    5 => H
    Name changes:
    reference_start  == pos
    reference_end    == aend
    query_length     == rlen
    reference_length == alen
    query_sequence   == seq
    return: length of soft-clipped part, sequence of soft-clipped part, the connection point of soft-clipped part (left/right), left/right soft-clipped part: 0:other; 2:left[SM]; 1:right[MS]
    """
    if read.cigar[0][0] == 4:
        # there are soft-clipped segments in left and right both
        if read.cigar[-1][0] == 4:
            # length of left soft-clipped segment is bigger
            if read.cigar[0][1] > read.cigar[-1][1]:
                return (
                    read.cigar[0][1],
                    read.query_sequence[: read.cigar[0][1]],
                    read.reference_start,
                    2,
                )
            # length of right soft-clipped segment is bigger
            else:
                return (
                    read.cigar[-1][1],
                    read.query_sequence[read.query_length - read.cigar[-1][1] :],
                    read.reference_end - 1,
                    1,
                )
        # there are soft-clipped segments in left only
        else:
            return (
                read.cigar[0][1],
                read.query_sequence[: read.cigar[0][1]],
                read.reference_start,
                2,
            )
    # there are soft-clipped segments in right only
    elif read.cigar[-1][0] == 4:
        return (
            read.cigar[-1][1],
            read.query_sequence[read.query_length - read.cigar[-1][1] :],
            read.reference_end - 1,
            1,
        )
    else:
        return 0, "", -1, 0


def get_softclip_length_from_cigar(
    cigar_str_sa, ref_start_sa, strand_sa, strand_ra, query_sequence_ra
):
    """
    return: length of soft-clipped part, sequence of soft-clipped part, the connection point of soft-clipped part (left/right), left/right soft-clipped part: 0:other; 2:left[SM]; 1:right[MS]
    """
    cigartuple = re.findall(r"(\d+)(\w)", cigar_str_sa)

    soft_mode = cigarstring2mode(cigar_str_sa)
    matched_size, _, _ = measure_size(cigartuple)
    # SM
    if soft_mode == 2:
        soft_bp_pos = ref_start_sa
        soft_len = int(cigartuple[0][0])
        if strand_sa == strand_ra:
            soft_seq = query_sequence_ra[:soft_len]
        else:
            soft_seq = reverse_complement(query_sequence_ra)[:soft_len]
    # MS
    elif soft_mode == 1:
        soft_bp_pos = ref_start_sa + matched_size
        soft_len = int(cigartuple[-1][0])
        if strand_sa == strand_ra:
            soft_seq = query_sequence_ra[-soft_len:]
        else:
            soft_seq = reverse_complement(query_sequence_ra)[-soft_len:]
    return soft_len, soft_seq, soft_bp_pos, soft_mode


def measure_size(cigartuple):
    """
    [('2859', 'S'),
    ('60', 'M'),
    ('1', 'I'),
    ('15', 'M'),
    ('1', 'D')]
    """
    indel_size = 0
    match_size = 0
    cigar = []
    for each_type in cigartuple:
        if each_type[1] == "I":
            indel_size += -1 * int(each_type[0])
            # match_size += -1*int(each_type[0])
            cigar.append([1, int(each_type[0])])
        elif each_type[1] == "D" or each_type[1] == "N":
            indel_size += int(each_type[0])
            match_size += int(each_type[0])
            cigar.append([2, int(each_type[0])]) if each_type[
                1
            ] == "D" else cigar.append([3, int(each_type[0])])
        elif each_type[1] == "M":
            match_size += int(each_type[0])
            cigar.append([0, int(each_type[0])])
    return match_size, indel_size, cigar


def reference_span_from_cigarstring(cigarstring):
    """
    :param cigarstring: a cigar string
    :type cigarstring: str
    :return: reference_span
    :rtype: int
    """
    cigartuple = re.findall(r"(\d+)(\w)", cigarstring)
    reference_span, _, _ = measure_size(cigartuple)
    return reference_span


def multiple_sa_tag_selector(sa_tag_list):
    """1) select one SA tag with highest MAPQ
    2) if the MAPQ is the same, select the one with bigger value of (reference span) / (number of mismatches)
    """
    # chr22,23179918,-,2859S2543M107954D1417S,60,381

    mapq_dict = {}
    span_nm_dict = {}
    for i in sa_tag_list:
        cigar_str, mapq, nm = i.split(",")[3:]
        ref_span = reference_span_from_cigarstring(cigar_str)
        mapq_dict[i] = int(mapq)
        span_nm_dict[i] = float(ref_span / (int(nm) + 1))
    max_mapq = max(mapq_dict.values())
    span_nm_max_mapq_dict = {
        j: span_nm_dict[j] for j in mapq_dict if mapq_dict[j] == max_mapq
    }
    optimal_sa_tag = max(span_nm_max_mapq_dict, key=lambda k: span_nm_max_mapq_dict[k])

    return optimal_sa_tag


def min_and_max_softclip_segment(soft_list):
    """
    soft_mode, soft_seq, soft_pos, read.alignment
    """
    soft_dict = {}
    for i, j, k, n in soft_list:
        soft_dict[(i, j, k, n)] = len(j)

    max_len_soft = max(soft_dict, key=lambda k: soft_dict[k])
    min_len_soft = min(soft_dict, key=lambda k: soft_dict[k])
    return min_len_soft, max_len_soft


def min_and_max_softclip_segment2(soft_list):
    """
    soft_mode, soft_seq, soft_pos, read.alignment
    """
    soft_dict = {}
    for i, j, k, n in soft_list:
        soft_dict[(i, j, k, n)] = len(j)

    if len(soft_list) >= 2:
        second_max = sorted(list(soft_dict.values()), reverse=True)[1]
        second_min = sorted(list(soft_dict.values()), reverse=False)[1]
        for m in soft_dict:
            if soft_dict[m] == second_max:
                max_len_soft = m
            if soft_dict[m] == second_min:
                min_len_soft = m
    else:
        only_key = list(soft_dict.keys())[0]
        max_len_soft = only_key
        min_len_soft = only_key

    return min_len_soft, max_len_soft


def sequence_extracter(cigar, chrm, strand, pos, genome_fasta):
    # strand - or +
    # cigar_str is the SA tag removing heading and tailing 'S'
    cigartuple = re.findall(r"(\d+)(\w)", cigar)
    exon_positions = []
    current_pos = int(pos)
    start_pos = int(pos)
    for _len, _type in cigartuple:
        # adding corrodinates: D, N, M
        # pass corrodinates: I,S
        length = int(_len)
        if _type in {"M", "D"}:
            current_pos = current_pos + length
        elif _type == "N":
            exon_positions.append(
                HTSeq.GenomicInterval(chrm, start_pos, current_pos, strand)
            )
            current_pos = current_pos + length
            start_pos = current_pos
    # exon_positions.append((start_pos,  current_pos))
    exon_positions.append(HTSeq.GenomicInterval(chrm, start_pos, current_pos, strand))

    out_seq = ""
    if strand == "+":
        for __item in exon_positions:
            seq = genome_fasta[__item.chrom][__item.start : __item.end].seq
            out_seq += seq
    elif strand == "-":
        for __item in exon_positions[::-1]:
            seq = genome_fasta[__item.chrom][
                __item.start : __item.end
            ].reverse.complement.seq
            out_seq += seq
    return out_seq, tuple(exon_positions)


def obtain_seq_from_interval_list(in_list, genome_fasta):
    """
    in_list must be a genomic coordinate ordered tuple
    """
    strand = in_list[0].strand
    out_seq = ""
    if strand == "+":
        for _item in in_list:
            seq = genome_fasta[_item.chrom][_item.start : _item.end].seq
            out_seq += seq
    elif strand == "-":
        for _item in in_list[::-1]:
            seq = genome_fasta[_item.chrom][
                _item.start : _item.end
            ].reverse.complement.seq
            out_seq += seq
    return out_seq


def obtain_transcript_seq_from_up_down_tuple(in_tuple, genome_fasta):
    upstream = in_tuple[0]
    downstream = in_tuple[1]
    upstream_seq = obtain_seq_from_interval_list(upstream, genome_fasta)
    downstream_seq = obtain_seq_from_interval_list(downstream, genome_fasta)
    return upstream_seq + downstream_seq


def choose_longer_one(seq1, seq2):
    if len(seq1) > len(seq2):
        print("seq from left BP")
        return seq1
    else:
        print("seq from right BP")
        return seq2


def candidates_dict_to_GTF_contents(
    candidates_gtf_count_dict, gene_id, event_type, bp_pos1, bp_pos2, genome_fasta
):
    """"""
    gtf_list = []
    count = 1
    fasta_dict = {}
    for candidate in candidates_gtf_count_dict:
        fasta_id = ""
        fasta_seq = ""
        supporting_reads = candidates_gtf_count_dict[candidate]
        up_stream = candidate[0]
        down_stream = candidate[1]
        up_total_len = 0
        down_total_len = 0
        if up_stream and down_stream:
            gtf_list.append(
                f'{up_stream[0].chrom}\tScanNLS\ttranscript\t{up_stream[0].start+1}\t{up_stream[-1].end}\t.\t{up_stream[0].strand}\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}.{count}"; event_type "{event_type}"; segment "upstream"; support_reads "{supporting_reads}";'
            )
            for exon in up_stream:
                gtf_list.append(
                    f'{exon.chrom}\tScanNLS\texon\t{exon.start+1}\t{exon.end}\t.\t{exon.strand}\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}.{count}"; event_type "{event_type}"; segment "upstream";'
                )
                up_total_len += exon.end - exon.start
            gtf_list.append(
                f'{down_stream[0].chrom}\tScanNLS\ttranscript\t{down_stream[0].start+1}\t{down_stream[-1].end}\t.\t{down_stream[0].strand}\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}.{count}"; event_type "{event_type}"; segment "downstream"; support_reads "{supporting_reads}";'
            )
            for exon in down_stream:
                gtf_list.append(
                    f'{exon.chrom}\tScanNLS\texon\t{exon.start+1}\t{exon.end}\t.\t{exon.strand}\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}.{count}"; event_type "{event_type}"; segment "downstream";'
                )
                down_total_len += exon.end - exon.start
            total_len = up_total_len + down_total_len
            fasta_id = f">{gene_id}.{count} {event_type} {bp_pos1}-{bp_pos2} {up_total_len}|{total_len}"
            fasta_seq = obtain_transcript_seq_from_up_down_tuple(
                candidate, genome_fasta
            )
            fasta_dict[fasta_id] = fasta_seq
            count += 1
    return gtf_list, fasta_dict


def exon_positions_to_GTF_content(event_exons, gene_id, event_type):
    # chr1»■■■PacBio»■transcript»■36322970»■■■36324146»■■■.»■■-»■■.»■■gene_id "PBfusion.1"; transcript_id "PBfusion.1.1";
    # chr1»■■■PacBio»■exon»■■■36322970»■■■36323067»■■■.»■■-»■■.»■■gene_id "PBfusion.1"; transcript_id "PBfusion.1.1";
    # print(event_exons)
    if not event_exons:
        return []
    # elif event_exons == ([], []):
    #    return []
    elif len(event_exons) == 2 and (not event_exons[0] or not event_exons[1]):
        return []
    up_stream = event_exons[0]
    down_stream = event_exons[1]

    up_total_len = 0
    down_total_len = 0

    out_list = []
    out_list.append(
        f'{up_stream[0].chrom}\tScanNLS\ttranscript\t{up_stream[0].start+1}\t{up_stream[-1].end}\t.\t{up_stream[0].strand}\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}.1"; event_type "{event_type}";'
    )
    for exon in up_stream:
        out_list.append(
            f'{exon.chrom}\tScanNLS\texon\t{exon.start+1}\t{exon.end}\t.\t{exon.strand}\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}.1"; event_type "{event_type}";'
        )
        up_total_len += exon.end - exon.start
    out_list.append(
        f'{down_stream[0].chrom}\tScanNLS\ttranscript\t{down_stream[0].start+1}\t{down_stream[-1].end}\t.\t{down_stream[0].strand}\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}.2"; event_type "{event_type}";'
    )
    for exon in down_stream:
        out_list.append(
            f'{exon.chrom}\tScanNLS\texon\t{exon.start+1}\t{exon.end}\t.\t{exon.strand}\t.\tgene_id "{gene_id}"; transcript_id "{gene_id}.2"; event_type "{event_type}";'
        )
        down_total_len += exon.end - exon.start
    total_len = up_total_len + down_total_len
    return out_list, up_total_len, total_len


def cigartuple2cigartuple_in_reads(cigartuple):
    return list(map(tuple, cigartuple))


def cigartuple2cigarstring(cigartuple):
    """[[4,40],[0,100],[2,2],[0,10],[4,10]] => 40S100M2D10M10S"""
    cigar_char_dict = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H"}
    cigarstring = ""
    for i, j in cigartuple:
        _code = int(i)
        _len = int(j)
        if _code in cigar_char_dict:
            __char = cigar_char_dict[_code]
        else:
            __char = "X"
        cigarstring += f"{_len}{__char}"
    return cigarstring


def cigarstring2mode(cigar_str):
    """MS => 1
    SM => 2
    """
    cigartuple = re.findall(r"(\d+)(\w)", cigar_str)

    if cigartuple[0][1] == "S":
        # there are soft-clipped segments in left and right both
        if cigartuple[-1][1] == "S":
            # length of left soft-clipped segment is bigger
            if int(cigartuple[0][0]) > int(cigartuple[-1][0]):
                soft_mode = 2
            # length of right soft-clipped segment is bigger
            else:
                soft_mode = 1
        # there are soft-clipped segments in left only
        else:
            soft_mode = 2
    # there are soft-clipped segments in right only
    elif cigartuple[-1][1] == "S":
        soft_mode = 1
    return soft_mode


def update_cigar(
    chrm, ref_start, strand, cigar_str, tgt_breakpoint, genome_fasta, is_can=True
):
    """update the cigar string using the breakpoint
    increase the length of softclipped part, decrease the length of matched part
    both appliable for primary alignment and supplementary alignment

    return: updated_ref_start, updated_cigartuple
    """
    cigar_dict = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5}
    # if cigarstring is 40M25N5M then cigartuple is [('40', 'M'), ('25', 'N'), ('5', 'M')]
    cigartuple_list = []
    cigartuple = re.findall(r"(\d+)(\w)", cigar_str)
    for i, j in cigartuple:
        _len = int(i)
        _type = cigar_dict[j]
        cigartuple_list.append([_type, _len])

    soft_mode = cigarstring2mode(cigar_str)
    matched_size, _, _ = measure_size(cigartuple)

    # CANONICAL splice sites
    if is_can:
        shift_len = 0
        tgt_motif = ""
        motifs = set()
        # SM
        if soft_mode == 2:
            soft_bp_pos = ref_start
            updated_ref_start = soft_bp_pos
            updated_ref_end = soft_bp_pos + matched_size
            if strand == "+":
                motifs = {"AG", "AC"}
            elif strand == "-":
                motifs = {"AC", "GC", "AT"}
            boundary_seq = genome_fasta[chrm][soft_bp_pos - 2 : soft_bp_pos + 8].seq
            # left to right searching
            for i in range(len(boundary_seq) - 1):
                _motif = boundary_seq[i : i + 2]
                if _motif in motifs:
                    shift_len = i
                    tgt_motif = _motif
                    break
            cigartuple_list[0][1] = cigartuple_list[0][1] + shift_len
            cigartuple_list[1][1] = cigartuple_list[1][1] - shift_len
            updated_ref_start = updated_ref_start + shift_len
        # MS
        elif soft_mode == 1:
            soft_bp_pos = ref_start + matched_size
            updated_ref_end = soft_bp_pos
            updated_ref_start = ref_start
            if strand == "+":
                motifs = {"GT", "GC", "AT"}
            elif strand == "-":
                motifs = {"CT", "GT"}
            boundary_seq = genome_fasta[chrm][soft_bp_pos - 8 : soft_bp_pos + 2].seq
            # right to left searching
            boundary_seq = boundary_seq[::-1]
            for i in range(len(boundary_seq) - 1):
                _motif = boundary_seq[i : i + 2][::-1]
                if _motif in motifs:
                    shift_len = i
                    tgt_motif = _motif
                    break
            cigartuple_list[-1][1] = cigartuple_list[-1][1] + shift_len
            cigartuple_list[-2][1] = cigartuple_list[-2][1] - shift_len
            updated_ref_end = updated_ref_end - shift_len
    # NONCANONICAL splice sites
    else:
        if soft_mode == 2:
            soft_bp_pos = ref_start
            updated_ref_start = soft_bp_pos
            updated_ref_end = soft_bp_pos + matched_size

            # shift_len could be negative
            shift_len = tgt_breakpoint - soft_bp_pos
            if shift_len > 10:
                return False, False, False

            cigartuple_list[0][1] = cigartuple_list[0][1] + shift_len
            cigartuple_list[1][1] = cigartuple_list[1][1] - shift_len
            if cigartuple_list[0][1] <= 0 or cigartuple_list[1][1] <= 0:
                return False, False, False
            updated_ref_start = updated_ref_start + shift_len
        elif soft_mode == 1:
            soft_bp_pos = ref_start + matched_size
            updated_ref_end = soft_bp_pos
            updated_ref_start = ref_start
            # shift_len could be negative
            shift_len = soft_bp_pos - tgt_breakpoint

            if shift_len > 10:
                return False, False, False
            cigartuple_list[-1][1] = cigartuple_list[-1][1] + shift_len
            cigartuple_list[-2][1] = cigartuple_list[-2][1] - shift_len
            if cigartuple_list[-1][1] <= 0 or cigartuple_list[-2][1] <= 0:
                return False, False, False
            updated_ref_end = updated_ref_end - shift_len
        # cigarstring, cigar, reference_start, reference_end
        # print(soft_mode, shift_len)
    return (
        cigartuple2cigarstring(cigartuple_list),
        cigartuple2cigartuple_in_reads(cigartuple_list),
        updated_ref_start,
    )
    # print(boundary_seq, shift_len, tgt_motif)


def aln_to_position(short_seq, long_seq):
    """Searching for short_seq in long_seq to find the matched position
    A-AG -TAG
    ATAG ATAG
    """
    # alignment_result = aligner(short_seq, long_seq, method='global_cfe')[0]
    alignment_result = aligner(short_seq, long_seq, method="glocal")[0]
    search_seq = alignment_result.seq1.decode("utf-8")
    target_seq = alignment_result.seq2.decode("utf-8")
    search_seq_len = len(short_seq)
    target_seq_len = len(long_seq)
    search_start, search_end = alignment_result.start1, alignment_result.end1 - 1
    target_start, target_end = alignment_result.start2, alignment_result.end2 - 1
    aln_len = target_end - target_start + 1
    total_mismatches = len(search_seq) - aln_len + alignment_result.n_mismatches

    shift_len_L = 100
    shift_len_R = 100
    if total_mismatches <= 2:
        # aln_left
        # MS, left shift length
        # if target_start <= target_seq_len - target_end - 1:
        # shift_len_L = target_start - search_start
        # aln_right
        # SM, right shift length
        # elif target_start > target_seq_len - target_end - 1:
        #    shift_len = target_seq_len - target_end - 1
        shift_len_L = target_start - search_start
        shift_len_R = target_seq_len - target_end - 1
    return shift_len_L, shift_len_R, aln_len


def correct_breakpoint(
    chrm, ref_start, cigar_str, strand_sa, strand_ra, soft_seq, genome_fasta
):
    """calculate the correct supplementary alignment breakpoint
    (alternative soft-clipped boundary point)
    """

    SOFT_LEN = len(soft_seq)

    cigar_dict = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5}
    # if cigarstring is 40M25N5M then cigartuple is [('40', 'M'), ('25', 'N'), ('5', 'M')]
    cigartuple_list = []
    cigartuple = re.findall(r"(\d+)(\w)", cigar_str)
    for i, j in cigartuple:
        _len = int(i)
        _type = cigar_dict[j]
        cigartuple_list.append([_type, _len])

    soft_mode = cigarstring2mode(cigar_str)
    matched_size, _, _ = measure_size(cigartuple)

    if strand_sa == strand_ra:
        searching_seq = soft_seq
    else:
        searching_seq = reverse_complement(soft_seq)
    print("############# BEGIN ################################")
    ## print('Searching sequence: ',searching_seq)
    new_breakpoint = 0
    shift_len = 0
    # SM, right to left finding
    if soft_mode == 2:
        cigartuple_list.pop(0)
        soft_bp_pos = ref_start
        updated_ref_start = soft_bp_pos
        updated_ref_end = soft_bp_pos + matched_size

        _total_len = 0
        _matched_size = 0
        shorter_flag = False
        for _t, _l in reversed(cigartuple_list):
            # M
            if _t == 0:
                _total_len += _l
                _matched_size += _l
            # S+I
            elif _t in {1, 4}:
                _total_len += _l
            # D+N
            elif _t in {2, 3}:
                _matched_size += _l
            if _total_len >= SOFT_LEN:
                shorter_flag = True
                break
        ## print(nameof(_total_len), _total_len)
        if shorter_flag:
            ## print(nameof(_matched_size), _matched_size)
            shift_len = -(
                updated_ref_end
                - _matched_size
                + (_total_len - SOFT_LEN)
                - updated_ref_start
            )
            ## print('IF soft-part is shorter than the matched-part, shift_len:', shift_len)
        # SOFT_LEN > _total_len
        else:
            diff_len = SOFT_LEN - _total_len
            ## print('IF soft-part is longer than the matched-part, diff_len:', diff_len)
            _soft_seq = searching_seq[: diff_len + 5]
            _ref_seq = genome_fasta[chrm][
                soft_bp_pos - (diff_len + 10) : soft_bp_pos + 10
            ].seq
            ## print('_soft_seq:', _soft_seq)
            ## print('_ref_seq:', _ref_seq)
            __shift_len_L, __shift_len_R, __aln_len = aln_to_position(
                _soft_seq, _ref_seq
            )
            ## print('__shift_len:', __shift_len_R)
            if __shift_len_R != 100:
                shift_len = __shift_len_R - 10 + __aln_len
            else:
                shift_len = 0
            ## print('shift_len:', shift_len)

        updated_ref_start = updated_ref_start - shift_len
        new_breakpoint = updated_ref_start

    # MS
    elif soft_mode == 1:
        cigartuple_list.pop(-1)
        soft_bp_pos = ref_start + matched_size
        print("soft_bp_pos:", soft_bp_pos)
        updated_ref_end = soft_bp_pos
        updated_ref_start = ref_start

        _total_len = 0
        _matched_size = 0
        shorter_flag = False
        for _t, _l in cigartuple_list:
            # M
            if _t == 0:
                _total_len += _l
                _matched_size += _l
            # S+I
            elif _t in {1, 4}:
                _total_len += _l
            # D+N
            elif _t in {2, 3}:
                _matched_size += _l
            if _total_len >= SOFT_LEN:
                shorter_flag = True
                break
        if shorter_flag:
            print(nameof(_matched_size), _matched_size)
            shift_len = -(
                updated_ref_end
                - _matched_size
                + (_total_len - SOFT_LEN)
                - updated_ref_start
            )
            print(
                "IF soft-part is shorter than the matched-part, shift_len:", shift_len
            )
        else:
            diff_len = SOFT_LEN - _total_len
            print("IF soft-part is longer than the matched-part, diff_len:", diff_len)
            _soft_seq = searching_seq[-(diff_len + 5) : len(searching_seq)]
            print(soft_bp_pos)
            _ref_seq = genome_fasta[chrm][
                soft_bp_pos - 10 : soft_bp_pos + diff_len + 10
            ].seq
            __shift_len_L, __shift_len_R, __aln_len = aln_to_position(
                _soft_seq, _ref_seq
            )
            print("_soft_seq:", _soft_seq)
            print("_ref_seq:", _ref_seq)
            print("__shift_len:", __shift_len_L, __shift_len_R)
            if __shift_len_L != 100:
                shift_len = __shift_len_L - 10 + __aln_len
            else:
                shift_len = 0
            print("shift_len:", shift_len)

        updated_ref_end = updated_ref_end + shift_len
        new_breakpoint = updated_ref_end
    print("############# END ######################")
    return new_breakpoint


def string2mode_num(in_str):
    if in_str == "SM":
        mode = 2
    elif in_str == "MS":
        mode = 1
    else:
        mode = 0
    return mode


def obtain_upstream_breakpoint_seq(
    chrm1, chrm2, pos1, pos2, strand1, strand2, mode1, mode2, genome_fasta
):
    """
    breakpoint upstream 10 bp sequence
    """
    if strand1 == "+" and strand2 == "-":
        if mode1 == 1 and mode2 == 1:
            seq = genome_fasta[chrm1][pos1 - 1 - 10 : pos1 - 1].seq
        elif mode1 == 2 and mode2 == 2:
            seq = genome_fasta[chrm2][pos2 - 1 : pos2 - 1 + 10].reverse.complement.seq
    elif strand1 == "-" and strand2 == "+":
        if mode1 == 1 and mode2 == 1:
            seq = genome_fasta[chrm2][pos2 - 1 - 10 : pos2 - 1].seq
        elif mode1 == 2 and mode2 == 2:
            seq = genome_fasta[chrm1][pos1 - 1 : pos1 - 1 + 10].reverse.complement.seq
    elif strand1 == "+" and strand2 == "+":
        if mode1 == 1 and mode2 == 2:
            seq = genome_fasta[chrm1][pos1 - 1 - 10 : pos1 - 1].seq
        elif mode1 == 2 and mode2 == 1:
            seq = genome_fasta[chrm2][pos2 - 1 - 10 : pos2 - 1].seq
    elif strand1 == "-" and strand2 == "-":
        if mode1 == 1 and mode2 == 2:
            seq = genome_fasta[chrm2][pos2 - 1 : pos2 - 1 + 10].reverse.complement.seq
        elif mode1 == 2 and mode2 == 1:
            seq = genome_fasta[chrm1][pos1 - 1 : pos1 - 1 + 10].reverse.complement.seq
    return seq


def refine_upstream_downstream_breakpoint(
    in_up_dict, in_down_dict, chr1, chr2, pos1, pos2, strand1, strand2, mode1, mode2
):
    """"""
    up_chrm = ""
    up_pos = 0
    down_chrm = ""
    down_pos = 0

    if strand1 == "+" and strand2 == "-":
        if mode1 == 1 and mode2 == 1:
            up_chrm = chr1
            up_pos = pos1
            down_chrm = chr2
            down_pos = pos2
        elif mode1 == 2 and mode2 == 2:
            up_chrm = chr2
            up_pos = pos2
            down_chrm = chr1
            down_pos = pos1
    elif strand1 == "-" and strand2 == "+":
        if mode1 == 1 and mode2 == 1:
            up_chrm = chr2
            up_pos = pos2
            down_chrm = chr1
            down_pos = pos1
        elif mode1 == 2 and mode2 == 2:
            up_chrm = chr1
            up_pos = pos1
            down_chrm = chr2
            down_pos = pos2
    elif strand1 == "+" and strand2 == "+":
        if mode1 == 1 and mode2 == 2:
            up_chrm = chr1
            up_pos = pos1
            down_chrm = chr2
            down_pos = pos2
        elif mode1 == 2 and mode2 == 1:
            up_chrm = chr2
            up_pos = pos2
            down_chrm = chr1
            down_pos = pos1
    elif strand1 == "-" and strand2 == "-":
        if mode1 == 1 and mode2 == 2:
            up_chrm = chr2
            up_pos = pos2
            down_chrm = chr1
            down_pos = pos1
        elif mode1 == 2 and mode2 == 1:
            up_chrm = chr1
            up_pos = pos1
            down_chrm = chr2
            down_pos = pos2
    up_strand = list(in_up_dict.keys())[0].strand
    down_strand = list(in_down_dict.keys())[0].strand

    if up_strand == "+":
        pass
    else:
        pass


# TODO
def vcf_construction(in_vcf, in_bam, out_gtf, ref_genome, output_prefix):
    is_can = False
    output_fa = open(f"{output_prefix}.fasta", "w")
    try:
        genome_fasta = Fasta(ref_genome, sequence_always_upper=True)
    except FastaNotFoundError as e:
        print("read reference genome " + ref_genome + " error!", e)
        sys.exit(1)
    bam_file = pysam.AlignmentFile(in_bam, "rb")
    vcf_reader = vcf.Reader(open(in_vcf), "r")
    output_gtf_file = open(out_gtf, "w")
    isoform_dict = OrderedDict()
    position_dict = {}

    for record in vcf_reader:
        chr1 = record.CHROM
        pos1 = int(record.POS)
        chr2 = record.INFO["CHR2"]
        pos2 = int(record.INFO["END"])
        _bp1 = f"{chr1}:{pos1}"
        _bp2 = f"{chr2}:{pos2}"
        sv_len = int(record.INFO["SVLEN"])
        strand1 = str(record.INFO["STRAND1"])
        strand2 = str(record.INFO["STRAND2"])
        _mode1 = str(record.INFO["MODE1"])
        _mode2 = str(record.INFO["MODE2"])
        mode1 = string2mode_num(_mode1)
        mode2 = string2mode_num(_mode2)
        _sv_type = record.INFO["SVTYPE"]
        gene_id = str(record.ID)

        region_l = f"{chr1}:{pos1 - 1}-{pos1}"
        region_r = f"{chr2}:{pos2 - 1}-{pos2}"

        print("DEALING WITH: ", gene_id, region_l, region_r)
        isoform_candidates = []
        GTF_candidates = []

        # left_breakpoint
        new_pos1 = 0
        new_pos2 = 0
        left_exons = []
        right_exons = []
        left_seq = ""
        right_seq = ""
        transcript_seq_l = ""
        transcript_exons_l = ()
        for pileupcolumn in bam_file.pileup(
            region=region_l, truncate=True, stepper="all", min_base_quality=0
        ):
            dp = pileupcolumn.nsegments
            chrm = pileupcolumn.reference_name
            sv_soft_segments = []
            soft_segments = []
            for read in pileupcolumn.pileups:
                read_ra = read.alignment
                if read_ra.has_tag("SV"):
                    (
                        sv_type,
                        _can,
                        position,
                        size_or_sup_position,
                        rep_aln_mode,
                        sup_aln_mode,
                        strands,
                        genes,
                    ) = read_ra.get_tag("SV").split(",")
                    bp1 = f"{chrm}:{position}"
                    if ":" in size_or_sup_position:
                        bp2 = size_or_sup_position
                    else:
                        bp2 = "{}:{}".format(
                            chrm, int(position) + int(size_or_sup_position)
                        )
                    if (
                        sv_type == _sv_type
                        and (bp1 == _bp1 and bp2 == _bp2)
                        or (bp1 == _bp2 and bp2 == _bp1)
                    ):
                        soft_len, soft_seq, soft_pos, soft_mode = get_softclip_length(
                            read_ra
                        )
                        new_pos1 = soft_pos
                        if soft_mode == 1:
                            new_pos1 = new_pos1 + 1
                        search_seq = soft_seq

                        if read_ra.is_reverse:
                            strand_ra = "-"
                        else:
                            strand_ra = "+"

                        chimeric_aln = read_ra.get_tag("SA")[:-1].split(";")
                        if len(chimeric_aln) == 1:
                            (
                                chrm_sa,
                                pos_sa,
                                strand_sa,
                                cigar_sa,
                                mapq_sa,
                                nm_sa,
                            ) = chimeric_aln[0].split(",")
                        else:
                            (
                                chrm_sa,
                                pos_sa,
                                strand_sa,
                                cigar_sa,
                                mapq_sa,
                                nm_sa,
                            ) = multiple_sa_tag_selector(chimeric_aln).split(",")

                        soft_mode = cigarstring2mode(cigar_sa)
                        pos_sa = int(pos_sa) - 1

                        new_pos2 = correct_breakpoint(
                            chrm_sa,
                            pos_sa,
                            cigar_sa,
                            strand_sa,
                            strand_ra,
                            search_seq,
                            genome_fasta,
                        )

                        if not new_pos2:
                            (
                                soft_len_sa,
                                soft_seq_sa,
                                soft_pos_sa,
                                soft_mode_sa,
                            ) = get_softclip_length_from_cigar(
                                cigar_sa,
                                pos_sa,
                                strand_sa,
                                strand_ra,
                                read_ra.query_sequence,
                            )
                            new_pos2 = soft_pos_sa
                            new_pos2 = new_pos2 - 1
                            search_seq = soft_seq_sa
                            if soft_mode_sa == 1:
                                new_pos2 = new_pos2 + 1
                            new_pos1 = correct_breakpoint(
                                read_ra.reference_name,
                                read_ra.reference_start,
                                read_ra.cigarstring,
                                strand_sa,
                                strand_ra,
                                search_seq,
                                genome_fasta,
                            )
                        print("new left positions: ", new_pos1, new_pos2)
                        #############################################################################################
                        (
                            updated_cigarstring,
                            updated_cigar,
                            updated_reference_start,
                        ) = update_cigar(
                            chrm,
                            read.alignment.reference_start,
                            strand1,
                            read.alignment.cigarstring,
                            new_pos1,
                            genome_fasta,
                            is_can,
                        )
                        if updated_reference_start:
                            (
                                read.alignment.cigarstring,
                                read.alignment.cigar,
                                read.alignment.reference_start,
                            ) = (
                                updated_cigarstring,
                                updated_cigar,
                                updated_reference_start,
                            )
                            left_seq, left_exons = sequence_extracter(
                                read_ra.cigarstring,
                                read_ra.reference_name,
                                strand_ra,
                                read_ra.reference_start,
                                genome_fasta,
                            )
                        (
                            updated_cigarstring_sa,
                            _,
                            updated_reference_start_sa,
                        ) = update_cigar(
                            chrm_sa,
                            pos_sa,
                            strand_sa,
                            cigar_sa,
                            new_pos2,
                            genome_fasta,
                            is_can,
                        )

                        if updated_reference_start_sa:
                            new_cigar_sa, new_pos_sa = (
                                updated_cigarstring_sa,
                                updated_reference_start_sa,
                            )
                            new_pos_sa = int(new_pos_sa)
                            right_seq, right_exons = sequence_extracter(
                                new_cigar_sa,
                                chrm_sa,
                                strand_sa,
                                new_pos_sa,
                                genome_fasta,
                            )

                        if strand1 == "+" and strand2 == "-":
                            if mode1 == 1 and mode2 == 1:
                                transcript_seq_l = left_seq + right_seq
                                transcript_exons_l = (left_exons, right_exons)
                            elif mode1 == 2 and mode2 == 2:
                                transcript_seq_l = right_seq + left_seq
                                transcript_exons_l = (right_exons, left_exons)
                        elif strand1 == "-" and strand2 == "+":
                            if mode1 == 1 and mode2 == 1:
                                transcript_seq_l = right_seq + left_seq
                                transcript_exons_l = (right_exons, left_exons)
                            elif mode1 == 2 and mode2 == 2:
                                transcript_seq_l = left_seq + right_seq
                                transcript_exons_l = (left_exons, right_exons)
                        elif strand1 == "+" and strand2 == "+":
                            if mode1 == 1 and mode2 == 2:
                                transcript_seq_l = left_seq + right_seq
                                transcript_exons_l = (left_exons, right_exons)
                            elif mode1 == 2 and mode2 == 1:
                                transcript_seq_l = right_seq + left_seq
                                transcript_exons_l = (right_exons, left_exons)
                        elif strand1 == "-" and strand2 == "-":
                            if mode1 == 1 and mode2 == 2:
                                transcript_seq_l = right_seq + left_seq
                                transcript_exons_l = (right_exons, left_exons)
                            elif mode1 == 2 and mode2 == 1:
                                transcript_seq_l = left_seq + right_seq
                                transcript_exons_l = (left_exons, right_exons)
                        if len(transcript_seq_l) > 0:
                            isoform_candidates.append(transcript_seq_l)
                            GTF_candidates.append(transcript_exons_l)
        ############################################################################
        new_pos1 = 0
        new_pos2 = 0
        left_exons = []
        right_exons = []
        left_seq = ""
        right_seq = ""
        transcript_seq_r = ""
        transcript_exons_r = ()
        # right_breakpoint
        for pileupcolumn in bam_file.pileup(
            region=region_r, truncate=True, stepper="all", min_base_quality=0
        ):
            dp = pileupcolumn.nsegments
            chrm = pileupcolumn.reference_name
            sv_soft_segments = []
            soft_segments = []
            for read in pileupcolumn.pileups:
                read_ra = read.alignment
                if read_ra.has_tag("SV"):
                    (
                        sv_type,
                        _can,
                        position,
                        size_or_sup_position,
                        rep_aln_mode,
                        sup_aln_mode,
                        strands,
                        genes,
                    ) = read_ra.get_tag("SV").split(",")
                    bp1 = f"{chrm}:{position}"
                    if ":" in size_or_sup_position:
                        bp2 = size_or_sup_position
                    else:
                        bp2 = "{}:{}".format(
                            chrm, int(position) + int(size_or_sup_position)
                        )
                    if (
                        sv_type == _sv_type
                        and (bp1 == _bp1 and bp2 == _bp2)
                        or (bp1 == _bp2 and bp2 == _bp1)
                    ):
                        soft_len, soft_seq, soft_pos, soft_mode = get_softclip_length(
                            read_ra
                        )
                        new_pos2 = soft_pos
                        if soft_mode == 1:
                            new_pos2 = new_pos2 + 1
                        search_seq = soft_seq

                        if read_ra.is_reverse:
                            strand_ra = "-"
                        else:
                            strand_ra = "+"

                        chimeric_aln = read_ra.get_tag("SA")[:-1].split(";")
                        if len(chimeric_aln) == 1:
                            (
                                chrm_sa,
                                pos_sa,
                                strand_sa,
                                cigar_sa,
                                mapq_sa,
                                nm_sa,
                            ) = chimeric_aln[0].split(",")
                        else:
                            (
                                chrm_sa,
                                pos_sa,
                                strand_sa,
                                cigar_sa,
                                mapq_sa,
                                nm_sa,
                            ) = multiple_sa_tag_selector(chimeric_aln).split(",")

                        soft_mode = cigarstring2mode(cigar_sa)
                        pos_sa = int(pos_sa) - 1

                        new_pos1 = correct_breakpoint(
                            chrm_sa,
                            pos_sa,
                            cigar_sa,
                            strand_sa,
                            strand_ra,
                            search_seq,
                            genome_fasta,
                        )

                        if not new_pos1:
                            (
                                soft_len_sa,
                                soft_seq_sa,
                                soft_pos_sa,
                                soft_mode_sa,
                            ) = get_softclip_length_from_cigar(
                                cigar_sa,
                                pos_sa,
                                strand_sa,
                                strand_ra,
                                read_ra.query_sequence,
                            )
                            new_pos1 = soft_pos_sa
                            new_pos1 = new_pos1 - 1
                            search_seq = soft_seq_sa
                            if soft_mode_sa == 1:
                                new_pos1 = new_pos1 + 1
                            new_pos2 = correct_breakpoint(
                                read_ra.reference_name,
                                read_ra.reference_start,
                                read_ra.cigarstring,
                                strand_sa,
                                strand_ra,
                                search_seq,
                                genome_fasta,
                            )
                        print("new right positions: ", new_pos1, new_pos2)
                        ########################################################################
                        (
                            updated_cigarstring,
                            updated_cigar,
                            updated_reference_start,
                        ) = update_cigar(
                            chrm,
                            read.alignment.reference_start,
                            strand2,
                            read.alignment.cigarstring,
                            new_pos2,
                            genome_fasta,
                            is_can,
                        )
                        if updated_reference_start:
                            (
                                read.alignment.cigarstring,
                                read.alignment.cigar,
                                read.alignment.reference_start,
                            ) = (
                                updated_cigarstring,
                                updated_cigar,
                                updated_reference_start,
                            )
                            right_seq, right_exons = sequence_extracter(
                                read_ra.cigarstring,
                                read_ra.reference_name,
                                strand_ra,
                                read_ra.reference_start,
                                genome_fasta,
                            )

                        (
                            updated_cigarstring_sa,
                            _,
                            updated_reference_start_sa,
                        ) = update_cigar(
                            chrm_sa,
                            pos_sa,
                            strand_sa,
                            cigar_sa,
                            new_pos1,
                            genome_fasta,
                            is_can,
                        )
                        if updated_reference_start_sa:
                            new_cigar_sa, new_pos_sa = (
                                updated_cigarstring_sa,
                                updated_reference_start_sa,
                            )
                            new_pos_sa = int(new_pos_sa)
                            left_seq, left_exons = sequence_extracter(
                                new_cigar_sa,
                                chrm_sa,
                                strand_sa,
                                new_pos_sa,
                                genome_fasta,
                            )

                        if strand1 == "+" and strand2 == "-":
                            if mode1 == 1 and mode2 == 1:
                                transcript_seq_r = left_seq + right_seq
                                transcript_exons_r = (left_exons, right_exons)
                            elif mode1 == 2 and mode2 == 2:
                                transcript_seq_r = right_seq + left_seq
                                transcript_exons_r = (right_exons, left_exons)
                        elif strand1 == "-" and strand2 == "+":
                            if mode1 == 1 and mode2 == 1:
                                transcript_seq_r = right_seq + left_seq
                                transcript_exons_r = (right_exons, left_exons)
                            elif mode1 == 2 and mode2 == 2:
                                transcript_seq_r = left_seq + right_seq
                                transcript_exons_r = (left_exons, right_exons)
                        elif strand1 == "+" and strand2 == "+":
                            if mode1 == 1 and mode2 == 2:
                                transcript_seq_r = left_seq + right_seq
                                transcript_exons_r = (left_exons, right_exons)
                            elif mode1 == 2 and mode2 == 1:
                                transcript_seq_r = right_seq + left_seq
                                transcript_exons_r = (right_exons, left_exons)
                        elif strand1 == "-" and strand2 == "-":
                            if mode1 == 1 and mode2 == 2:
                                transcript_seq_r = right_seq + left_seq
                                transcript_exons_r = (right_exons, left_exons)
                            elif mode1 == 2 and mode2 == 1:
                                transcript_seq_r = left_seq + right_seq
                                transcript_exons_r = (left_exons, right_exons)

                        if len(transcript_seq_r) > 0:
                            isoform_candidates.append(transcript_seq_r)
                            GTF_candidates.append(transcript_exons_r)

        ### candidate_seqs_fasta = seqs2tempfasta(isoform_candidates)

        # __seqs_fasta = seqs2fasta(isoform_candidates)

        ### consensus_seq = MSA_to_consensus(candidate_seqs_fasta)

        # print(consensus_seq)

        ### remove(os.path.dirname(candidate_seqs_fasta))

        ################################################################
        # obtain_upstream_breakpoint_seqfor _l_e, _r_e in GTF_candidates:
        #    print('{}\t{}'.format(_l_e, _r_e))

        # TEST
        # candidate_up_dict, candidate_down_dict = candidates_GTF_clustering(GTF_candidates)
        # up_num, down_num = candidates_GTF_clustering(GTF_candidates)
        # if up_num == 1 and down_num == 1:
        #    stats['1-1'] += 1
        # elif up_num >= 2 and down_num >= 2:
        #    stats['2-2'] += 1
        # else:
        #    stats['1-2'] += 1
        # up_dict, down_dict = refine_upstream_downstream_breakpoint(candidate_up_dict, candidate_down_dict, chr1, chr2, pos1, pos2, strand1, strand2, mode1, mode2)

        # GTF_consensus = candidates_to_consensus_GTF(GTF_candidates, True)
        # GTF_candidates = None
        # print('GTF_consensus: ', GTF_consensus)
        candidates_gtf_count_dict = candidates_GTF_clustering(GTF_candidates)
        if candidates_gtf_count_dict:
            gtf_record_list, fasta_dict = candidates_dict_to_GTF_contents(
                candidates_gtf_count_dict, gene_id, _sv_type, _bp1, _bp2, genome_fasta
            )
            # Output GTF
            if len(gtf_record_list) > 0:
                for j in gtf_record_list:
                    output_gtf_file.write(f"{j}\n")
            else:
                print("No current_output GTF!\n")
            for _id in fasta_dict:
                output_fa.write(f"{_id}\n")
                output_fa.write(f"{fasta_dict[_id]}\n")

    output_gtf_file.close()
    output_fa.close()
    return True


def find_up_len(upstream_seq, consensus_seq, possible_up_len):
    matches = []
    for match in re.finditer(rf"{upstream_seq}", consensus_seq):
        matches.append(int(match.end()))
    up_len = 0
    if len(matches) > 0:
        for pos in matches:
            if abs(pos - possible_up_len) < 10:
                up_len = pos
                break
        matches = None
    else:
        new_matches = []
        _seq = upstream_seq[:-1]
        for match in re.finditer(rf"{_seq}", consensus_seq):
            new_matches.append(int(match.end()))
        if len(new_matches) > 0:
            for pos in new_matches:
                if abs(pos - possible_up_len) < 10:
                    up_len = pos
                    break
        new_matches = None
    if up_len:
        return up_len
    else:
        return possible_up_len


def get_intron_positions(inlist):
    tmp_list = []
    for i in inlist:
        tmp_list.append(i.start)
        tmp_list.append(i.end)
    tmp_list.pop(0)
    tmp_list.pop()
    return tuple(tmp_list)


def gtf_distance(a, b):
    """"""
    if len(a) == len(b):
        if len(a) == 1:
            if a[0].start == b[0].start or a[0].end == b[0].end:
                return 0
            else:
                return 1
        elif len(a) > 1:
            if get_intron_positions(a) == get_intron_positions(b):
                return 0
            else:
                return 1
        else:
            return 1
    else:
        return 1


def gtf_clustering(inlist):
    g = nx.Graph()
    length = len(inlist)
    for x in range(length):
        node_x = inlist[x]
        g.add_similar_node(node_x)
        for y in range(length):
            if x < y:
                node_y = inlist[y]
                distance = gtf_distance(node_x, node_y)
                if distance == 0:
                    g.add_edge(node_x, node_y)
    cliques = list(nx.find_cliques(g))
    return cliques


def merge_cluster(inlist):
    """Input: a list of intervals in the same clique
    Merge the list of intervals according to smallest start point and largest end point
    """
    in_list = copy.deepcopy(inlist)
    cluster_support_reads = len(inlist)
    if cluster_support_reads == 1:
        return inlist[0]
    else:
        merged_list = []
        first_one = in_list[0]
        start = first_one[0].start
        end = first_one[-1].end
        for i in in_list:
            if i[0].start < start:
                start = i[0].start
            if i[-1].end > end:
                end = i[-1].end
        first_one[0].start = start
        first_one[-1].end = end
        return first_one


def candidates_GTF_clustering(candidates_list):
    """Input: list of upstream exons(tuple) and downstream exons (tuple)
    Output: Counter dictionary of upstream exons and downstream exons
    """
    upstreams = []
    downstreams = []
    for _up, _down in candidates_list:
        if _up and _down:
            upstreams.append(_up)
            downstreams.append(_down)

    if len(upstreams) == 0 or len(downstreams) == 0:
        return False

    up_cluster = gtf_clustering(upstreams)
    down_cluster = gtf_clustering(downstreams)

    up_dict = {}
    for i in up_cluster:
        meta_cluster = merge_cluster(i)
        for _i in i:
            up_dict[_i] = meta_cluster

    down_dict = {}
    for j in down_cluster:
        meta_cluster = merge_cluster(j)
        for _j in j:
            down_dict[_j] = meta_cluster

    out_candidates = Counter()
    for _up, _down in candidates_list:
        merge_up = up_dict[_up]
        merge_down = down_dict[_down]
        out_candidates.update([(merge_up, merge_down)])

    return out_candidates


def interval_merger(list_a, list_b):
    tmp_list = [list_a, list_b]
    sorted_list = sorted(tmp_list, key=lambda x: len(x), reverse=True)
    list_long = sorted_list[0]
    list_short = sorted_list[1]
    output_list = []
    for m in list_long:
        flag = True
        for n in list_short:
            if m.contains(n):
                output_list.append(m)
                flag = False
            elif n.contains(m):
                output_list.append(n)
                flag = False
        if flag:
            output_list.append(m)
    tmp_list = None
    sorted_list = None
    list_long = None
    list_short = None
    return output_list


def interval_list_merger(interval_list, shrink=False):
    sorted_list = sorted(interval_list, key=lambda x: len(x), reverse=True)
    tot_len = len(sorted_list)
    if tot_len == 1:
        return interval_list[0]
    elif tot_len == 2:
        return interval_merger(sorted_list[0], sorted_list[1])
    else:
        if tot_len > 30 and shrink:
            mid = int(tot_len / 2)
            indices = [
                0,
                1,
                2,
                mid,
                mid + 1,
                mid + 2,
                tot_len - 3,
                tot_len - 2,
                tot_len - 1,
            ]
            shrinked_list = [sorted_list[index] for index in indices]
            deal_list = shrinked_list
            shrinked_list = None
        else:
            deal_list = sorted_list
        first_two_list = deal_list[:2]
        remaining_list = deal_list[2:]
        output_list = interval_merger(first_two_list[0], first_two_list[1])
        for j in remaining_list:
            output_list = interval_merger(output_list, j)
        return output_list


def dumb_consensus(align_io_object):
    """Output a fast consensus sequence of the alignment.
    Arguments:
    """
    # Iddo Friedberg, 1-JUL-2004: changed ambiguous default to "X"
    consensus = ""

    # find the length of the consensus we are creating
    con_len = align_io_object.get_alignment_length()
    num_records = len(align_io_object)
    # print(con_len)
    pos_dict = OrderedDict()
    depth_dict = OrderedDict()
    # go through each seq item
    for n in range(con_len):
        # keep track of the counts of the different atoms we get
        atom_dict = {}
        num_atoms = 0

        for record in align_io_object:
            # make sure we haven't run past the end of any sequences
            # if they are of different lengths
            if n < len(record.seq):
                # print(record.seq)
                if record.seq[n] != "-" and record.seq[n] != ".":
                    if record.seq[n] not in atom_dict:
                        atom_dict[record.seq[n]] = 1
                    else:
                        atom_dict[record.seq[n]] += 1
                    num_atoms = num_atoms + 1

        pos_dict[n] = atom_dict
        depth_dict[n] = sum(atom_dict.values())

    first_N = 0
    for n in depth_dict:
        if depth_dict[n] / num_records <= 0.5:
            first_N += 1
        if depth_dict[n] / num_records > 0.5:
            break
    last_N = 0
    for n in reversed(depth_dict):
        if depth_dict[n] / num_records <= 0.5:
            last_N += 1
        if depth_dict[n] / num_records > 0.5:
            break
    # print(first_N, last_N)

    for n in pos_dict:
        atom_dict = pos_dict[n]
        max_atoms = []
        max_size = 0

        for atom in atom_dict:
            if atom_dict[atom] > max_size:
                max_atoms = [atom]
                max_size = atom_dict[atom]
            elif atom_dict[atom] == max_size:
                max_atoms.append(atom)
        if len(max_atoms) == 1:
            if n < first_N or n >= con_len - last_N:
                consensus += max_atoms[0]
            else:
                if atom_dict[max_atoms[0]] / len(align_io_object) >= 0.5:
                    consensus += max_atoms[0]
    return consensus


def seqs2fasta(in_list):

    # tmp_dir = tempfile.mkdtemp()
    tmp_file = "tmp.fasta"
    out = open(tmp_file, "w")
    count = 1
    for i in in_list:
        out.write(f">{count}\n")
        out.write(f"{i}\n")
        count += 1
    out.close()
    return tmp_file


def MSA_to_consensus(in_fasta, program="Mafft"):
    if program == "Muscle":
        cline = MuscleCommandline(input=in_fasta)
    elif program == "Mafft":
        out_fasta = os.path.join(os.path.dirname(in_fasta), "mafft.fasta")
        cline = f"mafft {in_fasta} > {out_fasta}"
        # cline = MafftCommandline(input=in_fasta)
    try:
        # print(cline)
        return_code = subprocess.check_call(
            cline, shell=True, stderr=subprocess.STDOUT, stdin=subprocess.PIPE
        )
        # print(stdout)
        # stdout = child.stdout.read().decode('utf8')
        # print(type(stdout))
    except CalledProcessError as e:
        print("Mafft error at: ", e)
        consensus_seq = ""
    else:
        if not return_code:
            alignment = AlignIO.read(out_fasta, "fasta")
            consensus_seq = dumb_consensus(alignment)
    # summary_align = AlignInfo.SummaryInfo(alignment)
    # consensus_seq = str(summary_align.dumb_consensus())
    consensus_seq = consensus_seq.upper()
    return consensus_seq


def output_fasta(position_dict, seq_dict, prefix):
    # TDUP chr1:89114016-chr1:89181143 950|1243
    output = open(f"{prefix}.fasta", "w")
    for _id in position_dict:
        print("########")
        print(position_dict[_id])
        print("########")
        sv_type, bp1, bp2, up_len, total_len = position_dict[_id].split("\t")
        output.write(f">{_id} {sv_type} {bp1}-{bp2} {up_len}|{total_len}\n")
        output.write(f"{seq_dict[_id]}\n")
    output.close()


def vcf_parser(
    in_vcf, pos_dict, in_bam, out_gtf, ref_genome, alignment_frac, mismatch_cutoff
):
    try:
        genome_fasta = Fasta(ref_genome, sequence_always_upper=True)
    except FastaNotFoundError as e:
        print("read reference genome " + ref_genome + " error!", e)
        sys.exit(1)
    bam_file = pysam.AlignmentFile(in_bam, "rb")
    vcf_reader = vcf.Reader(open(in_vcf), "r")
    # count = 0
    for record in vcf_reader:
        chrm1 = record.CHROM
        pos1 = int(record.POS)
        chrm2 = record.INFO["CHR2"]
        pos2 = int(record.INFO["END"])
        _bp1 = f"{chrm1}:{pos1}"
        _bp2 = f"{chrm2}:{pos2}"
        sv_len = int(record.INFO["SVLEN"])
        strand1 = str(record.INFO["STRAND1"])
        strand2 = str(record.INFO["STRAND2"])
        _mode1 = str(record.INFO["MODE1"])
        _mode2 = str(record.INFO["MODE2"])
        mode1 = string2mode_num(_mode1)
        mode2 = string2mode_num(_mode2)
        _sv_type = record.INFO["SVTYPE"]
        gene_id = str(record.ID)

        # if 'CANONICAL' in record.INFO:
        #    is_can = True
        # else:

        new_chrm_pos1, new_chrm_pos2 = pos_dict[gene_id].split("\t")
        new_pos1 = int(new_chrm_pos1.split(":")[1])
        new_pos2 = int(new_chrm_pos2.split(":")[1])
        is_can = False
        if new_pos1 > 0 and new_pos2 > 0:
            region_l = f"{chrm1}:{pos1 - 1}-{pos1}"
            region_r = f"{chrm2}:{pos2 - 1}-{pos2}"
            # count += 1
            print(region_l, region_r)
            # Left breakpoint
            transcript_seq_l, transcript_exons_l = transcript_contruction(
                bam_file,
                chrm1,
                chrm2,
                pos1,
                pos2,
                new_pos1,
                new_pos2,
                strand1,
                strand2,
                mode1,
                mode2,
                _sv_type,
                genome_fasta,
                alignment_frac,
                mismatch_cutoff,
                left=True,
            )
            # Right breakpoint
            transcript_seq_r, transcript_exons_r = transcript_contruction(
                bam_file,
                chrm1,
                chrm2,
                pos1,
                pos2,
                new_pos1,
                new_pos2,
                strand1,
                strand2,
                mode1,
                mode2,
                _sv_type,
                genome_fasta,
                alignment_frac,
                mismatch_cutoff,
                left=False,
            )

            transcript_seq = choose_longer_one(transcript_seq_l, transcript_seq_r)
            if transcript_seq == transcript_seq_l:
                transcript_exons = transcript_exons_l
            elif transcript_seq == transcript_seq_r:
                transcript_exons = transcript_exons_r
            # print(transcript_seq_l)
            # print(transcript_seq_r)
            # print(transcript_exons)
            gtf_record_list = exon_positions_to_GTF_content(
                transcript_exons, gene_id, _sv_type
            )
            # print(gtf_record_list)
            ####################BUG######################################
            if len(gtf_record_list) > 0:
                for j in gtf_record_list:
                    gtf_file.write(f"{j}\n")
            else:
                print("No current_output GTF!\n")
    gtf_file.close()


protein_coding = {"protein_coding"}


def GTF_reader(in_gtf, field="exon"):
    """transcript => exons in genomicinterval"""
    trx_to_exon = defaultdict(list)
    gtf_file = HTSeq.GFF_Reader(in_gtf)
    gene_positions = []
    gene_positions_dict = {}
    available_chroms = set()
    status_message("Loading GTF file!")
    for feature in gtf_file:
        biotype = feature.attr["gene_type"]
        gene_name = feature.attr["gene_name"]
        gene_id = feature.attr["gene_id"]
        if feature.type == field:
            if "transcript_type" in feature.attr:
                trx_biotype = feature.attr["transcript_type"]
                if biotype in protein_coding and trx_biotype in protein_coding:
                    trx_id = feature.attr["transcript_id"]
                    trx_to_exon[trx_id].append(feature.iv)

    sorted_trx_to_exon = {}
    for trx_id in trx_to_exon:
        tmpList = trx_to_exon[trx_id]
        tmpList.sort(key=lambda x: x.start)
        sorted_trx_to_exon[trx_id] = tmpList
    trx_to_exon = None
    status_message("GTF load finished!")
    return sorted_trx_to_exon


def update_gtf(in_gtf, out_gtf, out_fasta, ref_genome):
    """
    :param in_gtf: reads-inferred GTF
    :type in_gtf: str
    :param out_gtf: updated GTF
    :type out_gtf: str
    :param out_fasta: transcript sequence (5' to 3')
    :type out_fasta: str
    :param ref_genome: reference genome FASTA
    :type ref_genome: str
    :return: None
    :rtype: None
    """
    try:
        genome_fasta = Fasta(ref_genome, sequence_always_upper=True)
    except FastaNotFoundError as e:
        print("read reference genome " + ref_genome + " error!", e)
        sys.exit(1)

    trx_to_exons = defaultdict(list)
    gene_ids = set()
    event_type_dict = {}
    gtf_file = HTSeq.GFF_Reader(in_gtf)
    for feature in gtf_file:
        gene_id = feature.attr["gene_id"]
        event_type = feature.attr["event_type"]
        event_type_dict[gene_id] = event_type
        gene_ids.add(gene_id)
        trx_id = feature.attr["transcript_id"]
        if feature.type == "exon":
            trx_to_exons[trx_id].append(feature.iv)
    gene_id_list = list(gene_ids)
    gene_id_list.sort(key=lambda x: int(x.split("_")[1]))

    remove(in_gtf)
    out_fasta_file = open(out_fasta, "w")
    out_gtf_file = open(out_gtf, "w")
    for gene_id in gene_id_list:
        trx_upstream = f"{gene_id}.1"
        trx_downstream = f"{gene_id}.2"
        upstream_seq = ""
        downstream_seq = ""
        upstream_titles = []
        downstream_titles = []
        event_type = event_type_dict[gene_id]
        _up = trx_to_exons[trx_upstream]
        _down = trx_to_exons[trx_downstream]
        # current_output to GTF file
        out_gtf_file.write(
            f'{_up[0].chrom}\tScanNLS\ttranscript\t{_up[0].start+1}\t{_up[-1].end}\t.\t{_up[0].strand}\t.\tgene_id "{gene_id}"; transcript_id "{trx_upstream}"; event_type "{event_type}";\n'
        )
        for exon in _up:
            out_gtf_file.write(
                f'{exon.chrom}\tScanNLS\texon\t{exon.start+1}\t{exon.end}\t.\t{exon.strand}\t.\tgene_id "{gene_id}"; transcript_id "{trx_upstream}"; event_type "{event_type}";\n'
            )
            upstream_titles.append(
                f"{exon.chrom}:{exon.start+1}-{exon.end}({exon.strand})"
            )
        out_gtf_file.write(
            f'{_down[0].chrom}\tScanNLS\ttranscript\t{_down[0].start+1}\t{_down[-1].end}\t.\t{_down[0].strand}\t.\tgene_id "{gene_id}"; transcript_id "{trx_downstream}"; event_type "{event_type}";\n'
        )
        for exon in _down:
            out_gtf_file.write(
                f'{exon.chrom}\tScanNLS\texon\t{exon.start+1}\t{exon.end}\t.\t{exon.strand}\t.\tgene_id "{gene_id}"; transcript_id "{trx_downstream}"; event_type "{event_type}";\n'
            )
            downstream_titles.append(
                f"{exon.chrom}:{exon.start+1}-{exon.end}({exon.strand})"
            )
        bp_up = ""
        bp_down = ""
        # current_output upstream sequences
        if _up[0].strand == "+":
            chrm, start_end_pos = upstream_titles[-1].split("(")[0].split(":")
            bp_up = "{}:{}".format(chrm, int(start_end_pos.split("-")[1]) + 1)
            for exon in _up:
                upstream_seq += genome_fasta[exon.chrom][exon.start : exon.end].seq
        elif _up[0].strand == "-":
            chrm, start_end_pos = upstream_titles[0].split("(")[0].split(":")
            bp_up = "{}:{}".format(chrm, start_end_pos.split("-")[0])
            for exon in _up[::-1]:
                upstream_seq += genome_fasta[exon.chrom][
                    exon.start : exon.end
                ].reverse.complement.seq
        # current_output downstream sequences
        if _down[0].strand == "+":
            chrm, start_end_pos = downstream_titles[0].split("(")[0].split(":")
            bp_down = "{}:{}".format(chrm, start_end_pos.split("-")[0])
            for exon in _down:
                downstream_seq += genome_fasta[exon.chrom][exon.start : exon.end].seq
        elif _down[0].strand == "-":
            chrm, start_end_pos = downstream_titles[-1].split("(")[0].split(":")
            bp_down = "{}:{}".format(chrm, int(start_end_pos.split("-")[1]) + 1)
            for exon in _down[::-1]:
                downstream_seq += genome_fasta[exon.chrom][
                    exon.start : exon.end
                ].reverse.complement.seq

        seq = upstream_seq + downstream_seq
        title = ">{} {}|{}|{} {}-{} {}|{}".format(
            gene_id,
            "&".join(upstream_titles),
            "&".join(downstream_titles),
            event_type,
            bp_up,
            bp_down,
            len(upstream_seq),
            len(seq),
        )
        out_fasta_file.write(f"{title}\n")
        out_fasta_file.write(f"{seq}\n")
    out_gtf_file.close()


def ref_guided_update(in_exons, ref_exons, shift=50):
    """
    :param shift:
    :param in_exons: reads-inferred upstream transcript exons form
    :type in_exons: list
    :param ref_exons: overlapped reference transcript exons form
    :type ref_exons: list
    :return: out_exons
    :rtype: list
    """
    in_exons.sort(key=lambda x: x.start)
    strand = in_exons[0].strand
    out_exons = []
    flag = False
    if strand == "+":
        tgt_point = in_exons[-1].end
        for i in ref_exons:
            out_exons.append(i)
            if abs(i.end - tgt_point) <= shift:
                flag = True
                break
        # no break triggered
        # if out_exons[-1].end != tgt_point:
        #    out_exons = in_exons
        if not flag:
            out_exons = in_exons
    elif strand == "-":
        tgt_point = in_exons[0].start
        for i in ref_exons[::-1]:
            out_exons.append(i)
            if abs(i.start - tgt_point) <= shift:
                flag = True
                break
        out_exons = out_exons[::-1]
        if not flag:
            out_exons = in_exons
        # if out_exons[0].start != tgt_point:
        #    out_exons = in_exons
    return out_exons


def both_nonzero_positions(in_str):
    """"""
    _a, _b = in_str.split("\t")
    pos_a = int(_a.split(":")[1])
    pos_b = int(_b.split(":")[1])
    if not pos_a and not pos_b:
        return True
    else:
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="Input VCF file",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--bam",
        action="store",
        dest="bam",
        help="SV built BAM file",
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
        "-a",
        "--alignment_fraction",
        action="store",
        dest="alignment_fraction",
        type=float,
        help="minimal fraction of aligned part for smith-waterman local alignment (default: %(default)s)",
        default=0.8,
    )
    parser.add_argument(
        "-n",
        action="store",
        dest="mismatch",
        type=int,
        help="maximum mismatch bases of pairwise local alignment (default: %(default)s)",
        default=3,
    )
    # parser.add_argument('-g', '--gtf', action='store', dest='gtf', help="Gene annotation GTF file", required=False)
    # parser.add_argument('--ao', action='store', dest='ao', type=int, help="minimal AO (default: %(default)s)", default=2)
    # parser.add_argument('--dp', action='store', dest='dp', type=int, help="minimal DP (default: %(default)s)", default=5)
    # parser.add_argument('--af', action='store', dest='af', type=float, help="minimal AF (default: %(default)s)", default=0.01)
    # parser.add_argument('-t', '--type', action='store', dest='type', help="current_output type (default: %(default)s)", choices=['gene', 'transcript'], default='gene')
    parser.add_argument(
        "-o",
        "--current_output",
        action="store",
        dest="current_output",
        help="current_output file prefix (default: %(default)s)",
        default="current_output",
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")
    args = parser.parse_args()

    in_vcf = args.input
    in_bam = args.bam
    # ref_gtf = args.gtf
    # interim_gtf = f'{os.path.splitext(os.path.basename(in_vcf))[0]}.tmp.gtf'
    ref = args.ref
    align_frac = args.alignment_fraction
    nm = args.mismatch
    output_prefix = args.output
    out_gtf = f"{output_prefix}.gtf"
    vcf_construction(in_vcf, in_bam, out_gtf, ref, output_prefix)
