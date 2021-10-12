#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
import time
from collections import defaultdict
from collections import OrderedDict
from typing import Iterable

import psutil
from align import aligner
from Bio import SearchIO
from Bio.Seq import Seq
from pyfaidx import Fasta

from . import __version__
from .classes import Path

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


def extract_splice_sites(file, bin):
    gtf_file = HTSeq.GFF_Reader(file)
    cvg = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gene_iv = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    trx_to_exon = defaultdict(list)
    for feature in gtf_file:
        biotype = feature.attr["gene_type"]
        gene_name = feature.attr["gene_name"]
        if biotype == "protein_coding":
            if feature.type == "exon":
                trx_id = feature.attr["transcript_id"]
                trx_to_exon[trx_id].append(feature.iv)
            if feature.type == "gene":
                gene_iv[
                    HTSeq.GenomicInterval(
                        feature.iv.chrom,
                        feature.iv.start - bin,
                        feature.iv.end + bin,
                        ".",
                    )
                ] += str(gene_name)
        # if feature.type == 'gene':
        #    gene_iv[HTSeq.GenomicInterval(feature.iv.chrom, feature.iv.start-bin, feature.iv.end+bin, '.')] += f'{gene_name}:{biotype}'

    for trx_id in trx_to_exon:
        exonList = trx_to_exon[trx_id]
        exonList.sort(key=lambda x: x.start)
        exon_num = len(exonList)
        first_exon = exonList[0]
        last_exon = exonList[-1]
        strand = first_exon.strand
        if exon_num == 1:
            cvg[
                HTSeq.GenomicInterval(
                    first_exon.chrom,
                    first_exon.start - bin,
                    first_exon.start + bin,
                    ".",
                )
            ] += "XX"
            cvg[
                HTSeq.GenomicInterval(
                    first_exon.chrom, first_exon.end - bin, first_exon.end + bin, "."
                )
            ] += "XX"
        elif exon_num == 2:
            if strand == "+":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin,
                        first_exon.start + bin,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin,
                        first_exon.end + bin,
                        ".",
                    )
                ] += "GT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin,
                        last_exon.start + bin,
                        ".",
                    )
                ] += "AG"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom, last_exon.end - bin, last_exon.end + bin, "."
                    )
                ] += "XX"
            elif strand == "-":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin,
                        first_exon.start + bin,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin,
                        first_exon.end + bin,
                        ".",
                    )
                ] += "CT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin,
                        last_exon.start + bin,
                        ".",
                    )
                ] += "AC"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom, last_exon.end - bin, last_exon.end + bin, "."
                    )
                ] += "XX"
        else:
            if strand == "+":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin,
                        first_exon.start + bin,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin,
                        first_exon.end + bin,
                        ".",
                    )
                ] += "GT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin,
                        last_exon.start + bin,
                        ".",
                    )
                ] += "AG"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom, last_exon.end - bin, last_exon.end + bin, "."
                    )
                ] += "XX"
            elif strand == "-":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin,
                        first_exon.start + bin,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin,
                        first_exon.end + bin,
                        ".",
                    )
                ] += "CT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin,
                        last_exon.start + bin,
                        ".",
                    )
                ] += "AC"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom, last_exon.end - bin, last_exon.end + bin, "."
                    )
                ] += "XX"
            for _exon in exonList[1:-1]:
                iv1 = HTSeq.GenomicInterval(
                    _exon.chrom, _exon.start - bin, _exon.start + bin, "."
                )
                iv2 = HTSeq.GenomicInterval(
                    _exon.chrom, _exon.end - bin, _exon.end + bin, "."
                )
                if strand == "+":
                    cvg[iv1] += "AG"
                    cvg[iv2] += "GT"
                elif strand == "-":
                    cvg[iv1] += "AC"
                    cvg[iv2] += "CT"
    sys.stdout.write("{} is fully loaded!\n".format(file))
    return cvg, gene_iv


def junc_site_checker(junc_seq):
    sites = ["GT", "AG", "CT", "AC"]
    matches = list((i in junc_seq for i in sites))
    if True in matches:
        posInA = matches.index(True)
        tgt_site = sites[posInA]
        return tgt_site
    else:
        return False


def gene_annotation(chrm1, pos1, chrm2, pos2, gene_iv):
    """breakpoints gene annotations"""
    # print(chrm1, pos1, type(chrm1), type(pos1))
    # print(chrm2, pos2, type(chrm2), type(pos2))
    # print(gene_iv[HTSeq.GenomicPosition(chrm1, pos1)])
    # print(gene_iv[HTSeq.GenomicPosition(chrm2, pos2)])
    try:
        gene1 = "&".join(list(gene_iv[HTSeq.GenomicPosition(chrm1, pos1)]))
    except IndexError:
        gene1 = ""
    except TypeError:
        print(chrm1, pos1)
    try:
        gene2 = "&".join(list(gene_iv[HTSeq.GenomicPosition(chrm2, pos2)]))
    except IndexError:
        gene2 = ""
    except TypeError:
        print(chrm2, pos2)

    if not gene1:
        gene1 = "INTERGENIC"
    if not gene2:
        gene2 = "INTERGENIC"
    return (gene1, gene2)


def splicing_confirmation(
    chrm1,
    pos1,
    chrm2,
    pos2,
    splice_bin,
    fastafile,
    cvg,
    strand_changed,
    motif_required=True,
):

    """motif_required = True  => considering canonical splice sites only
    motif_required = False => considering canonical and noncanonical splice sites both

    True,  3(11), 1 => reported, both breakpoints overlap with coding exons boundary, using canonical splice motif
    True,  2(10), 0 => reported, one breakpoint overlap with coding exons boundary, using noncanonical splice motif
    True,  1(01), 0 => reported, one breakpoint overlap with coding exons boundary, using noncanonical splice motif
    False, 0(00), 0 => not reported, none breakpoint overlap with coding exons boundary, using noncanonical splice motif
    """
    if strand_changed:
        splice_motif_dict = {"GT": "CT", "AG": "AC", "CT": "GT", "AC": "AG"}
    else:
        splice_motif_dict = {"GT": "AG", "AG": "GT", "CT": "AC", "AC": "CT"}
    try:
        junc1 = list(cvg[HTSeq.GenomicPosition(chrm1, pos1)])[0]
    except IndexError:
        junc1 = ""
    try:
        junc2 = list(cvg[HTSeq.GenomicPosition(chrm2, pos2)])[0]
    except IndexError:
        junc2 = ""
    # Non-annotated coding exon boundary
    if junc1 not in splice_motif_dict and junc2 not in splice_motif_dict:
        junc_seq1 = fastafile[chrm1][pos1 - splice_bin : pos1 + splice_bin].seq
        junc_seq2 = fastafile[chrm2][pos2 - splice_bin : pos2 + splice_bin].seq
        _junc1 = junc_site_checker(junc_seq1)
        _junc2 = junc_site_checker(junc_seq2)
        if _junc1 and splice_motif_dict[_junc1] in junc_seq2:
            if motif_required:
                return True, 0, 1
            else:
                return True, 0, 1
        elif _junc2 and splice_motif_dict[_junc2] in junc_seq1:
            if motif_required:
                return True, 0, 1
            else:
                return True, 0, 1
        else:
            if motif_required:
                return False, 0, 0
            else:
                return True, 0, 0
    # pos1 in annotated coding exon boundary, pos2 not.
    elif junc1 in splice_motif_dict and junc2 not in splice_motif_dict:
        junc_seq = fastafile[chrm2][pos2 - splice_bin : pos2 + splice_bin].seq
        if splice_motif_dict[junc1] in junc_seq:
            if motif_required:
                return True, 2, 1
            else:
                return True, 2, 1
        else:
            if motif_required:
                return False, 2, 0
            else:
                return True, 2, 0

    # pos2 in annotated coding exon boundary, pos1 not.
    elif junc1 not in splice_motif_dict and junc2 in splice_motif_dict:
        junc_seq = fastafile[chrm1][pos1 - splice_bin : pos1 + splice_bin].seq
        if splice_motif_dict[junc2] in junc_seq:
            if motif_required:
                return True, 1, 1
            else:
                return True, 1, 1
        else:
            if motif_required:
                return False, 1, 0
            else:
                return True, 1, 0

    # pos1 and pos2 both in annotated coding exon boundary
    # junc1 in splice_motif_dict and junc2 in splice_motif_dict
    else:
        if splice_motif_dict[junc1] == junc2:
            if motif_required:
                return True, 3, 1
            else:
                return True, 3, 1
        else:
            if motif_required:
                return False, 3, 0
            else:
                return True, 3, 0


def update_breakpoints(
    bp1_chrm,
    bp1_pos,
    bp2_chrm,
    bp2_pos,
    bp1_strand,
    bp2_strand,
    bp1_mode,
    bp2_mode,
    splice_bin,
    genome_fasta,
):
    """
    update the breakpoints of NLS events with canonical splice sites
    keep NLS events with noncanonical splice sites unchanged
    GT-AG;GC-AG;AT-AC
    """
    tgt_motifs = {
        1: {"+": {"GT", "GC", "AT"}, "-": {"AG", "AC"}},
        2: {"+": {"AG", "AC"}, "-": {"GT", "GC", "AT"}},
    }

    # SM: breakpoint move to right
    # MS: breakpoint move to left
    # boundary_seq1
    # print(bp1_chrm, bp1_pos, splice_bin, bp1_mode, bp1_strand)
    # print(bp2_chrm, bp2_pos, splice_bin, bp2_mode, bp2_strand)
    bp1_pos_dict = {}
    if bp1_mode == 2:
        if bp1_strand == "+":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].seq
            # left to right: AG, AC
            bp1_pos_dict["AG"] = splice_site_search(boundary_seq1, "AG", splice_bin)
            bp1_pos_dict["AC"] = splice_site_search(boundary_seq1, "AC", splice_bin)
        elif bp1_strand == "-":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].reverse.complement.seq
            # right to left: GT, GC, AT
            # print(boundary_seq1)
            bp1_pos_dict["GT"] = splice_site_search(
                boundary_seq1, "GT", splice_bin, False
            )
            # print(bp1_pos_dict)

            bp1_pos_dict["GC"] = splice_site_search(
                boundary_seq1, "GC", splice_bin, False
            )
            bp1_pos_dict["AT"] = splice_site_search(
                boundary_seq1, "AT", splice_bin, False
            )
    elif bp1_mode == 1:
        if bp1_strand == "+":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].seq
            # right to left: GT, GC, AT
            bp1_pos_dict["GT"] = splice_site_search(
                boundary_seq1, "GT", splice_bin, False
            )
            bp1_pos_dict["GC"] = splice_site_search(
                boundary_seq1, "GC", splice_bin, False
            )
            bp1_pos_dict["AT"] = splice_site_search(
                boundary_seq1, "AT", splice_bin, False
            )
        elif bp1_strand == "-":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].reverse.complement.seq
            # left to right: AG, AC
            bp1_pos_dict["AG"] = splice_site_search(boundary_seq1, "AG", splice_bin)
            bp1_pos_dict["AC"] = splice_site_search(boundary_seq1, "AC", splice_bin)
    # boundary_seq2
    bp2_pos_dict = {}
    if bp2_mode == 2:
        if bp2_strand == "+":
            # if TRA:
            #    boundary_seq2 = genome_fasta[bp2_chrm][bp2_pos-1-splice_bin:bp2_pos-1+splice_bin].seq
            # else:
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].seq
            # left to right: AG, AC
            bp2_pos_dict["AG"] = splice_site_search(boundary_seq2, "AG", splice_bin)
            bp2_pos_dict["AC"] = splice_site_search(boundary_seq2, "AC", splice_bin)
        elif bp2_strand == "-":
            # if TRA:
            #    boundary_seq2 = genome_fasta[bp2_chrm][bp2_pos-1-splice_bin:bp2_pos-1+splice_bin].reverse.complement.seq
            # else:
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].reverse.complement.seq
            # right to left: GT, GC, AT
            bp2_pos_dict["GT"] = splice_site_search(
                boundary_seq2, "GT", splice_bin, False
            )
            bp2_pos_dict["GC"] = splice_site_search(
                boundary_seq2, "GC", splice_bin, False
            )
            bp2_pos_dict["AT"] = splice_site_search(
                boundary_seq2, "AT", splice_bin, False
            )
    elif bp2_mode == 1:
        if bp2_strand == "+":
            # if TRA:
            #    boundary_seq2 = genome_fasta[bp2_chrm][bp2_pos-1-splice_bin:bp2_pos-1+splice_bin].seq
            # else:
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].seq
            # right to left: GT, GC, AT
            bp2_pos_dict["GT"] = splice_site_search(
                boundary_seq2, "GT", splice_bin, False
            )
            bp2_pos_dict["GC"] = splice_site_search(
                boundary_seq2, "GC", splice_bin, False
            )
            bp2_pos_dict["AT"] = splice_site_search(
                boundary_seq2, "AT", splice_bin, False
            )
        elif bp2_strand == "-":
            # if TRA:
            #    boundary_seq2 = genome_fasta[bp2_chrm][bp2_pos-1-splice_bin:bp2_pos-1+splice_bin].reverse.complement.seq
            # else:
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].reverse.complement.seq
            # left to right: AG, AC
            bp2_pos_dict["AG"] = splice_site_search(boundary_seq2, "AG", splice_bin)
            bp2_pos_dict["AC"] = splice_site_search(boundary_seq2, "AC", splice_bin)

    # print(f'bp1_pos_dict: {bp1_pos_dict}', f'bp2_pos_dict: {bp2_pos_dict}')
    # print('bp1_mode: {}, bp2_mode: {}'.format(bp1_mode, bp2_mode))
    # print('bp1_seq: {}, bp2_seq: {}, bp1: {}:{}, bp2: {}:{}'.format(boundary_seq1, boundary_seq2, bp1_chrm, bp1_pos, bp2_chrm, bp2_pos))
    ####
    shift1 = 0
    shift2 = 0
    if bp1_strand == "+" and bp2_strand == "-":
        if bp1_mode == 1 and bp2_mode == 1:
            shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, True)
        elif bp1_mode == 2 and bp2_mode == 2:
            shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, False)
    elif bp1_strand == "-" and bp2_strand == "+":
        if bp1_mode == 1 and bp2_mode == 1:
            shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, False)
        elif bp1_mode == 2 and bp2_mode == 2:
            shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, True)

    elif bp1_strand == "+" and bp2_strand == "+":
        if bp1_mode == 1 and bp2_mode == 2:
            shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, True)
        elif bp1_mode == 2 and bp2_mode == 1:
            shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, False)
    elif bp1_strand == "-" and bp2_strand == "-":
        if bp1_mode == 1 and bp2_mode == 2:
            shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, False)
        elif bp1_mode == 2 and bp2_mode == 1:
            shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, True)

    # print(bp1_pos, bp2_pos, shift1, shift2)

    # No canonical splice sites found, so change it to noncanonical label
    if shift1 == True and shift2 == True:
        return 0, 0

    new_bp1_pos = bp1_pos
    new_bp2_pos = bp2_pos

    if bp1_mode == 2:
        new_bp1_pos = bp1_pos + (shift1 - splice_bin + 2) - 1
    elif bp1_mode == 1:
        new_bp1_pos = bp1_pos - (shift1 - splice_bin + 2) + 1

    if bp2_mode == 2:
        new_bp2_pos = bp2_pos + (shift2 - splice_bin + 2) - 1
    elif bp2_mode == 1:
        new_bp2_pos = bp2_pos - (shift2 - splice_bin + 2) + 1

    # print('new_bp1: ',new_bp1_pos, 'new_bp2: ',new_bp2_pos)
    return new_bp1_pos, new_bp2_pos


def obtain_bps_shift_len(bp1_dict, bp2_dict, bp1_upstream=True):
    """bp1_upstream: breakpoint 1 locates at the upstream half of the transcript"""
    if bp1_upstream:
        if (
            "GT" in bp1_dict
            and "AG" in bp2_dict
            and bp1_dict["GT"] != -1
            and bp2_dict["AG"] != -1
        ):
            return bp1_dict["GT"], bp2_dict["AG"]
        elif (
            "GC" in bp1_dict
            and "AG" in bp2_dict
            and bp1_dict["GC"] != -1
            and bp2_dict["AG"] != -1
        ):
            return bp1_dict["GC"], bp2_dict["AG"]
        elif (
            "AT" in bp1_dict
            and "AC" in bp2_dict
            and bp1_dict["AT"] != -1
            and bp2_dict["AC"] != -1
        ):
            return bp1_dict["AT"], bp2_dict["AC"]
        else:
            # No splice site found
            return True, True
    else:
        if (
            "GT" in bp2_dict
            and "AG" in bp1_dict
            and bp2_dict["GT"] != -1
            and bp1_dict["AG"] != -1
        ):
            return bp1_dict["AG"], bp2_dict["GT"]
        elif (
            "GC" in bp2_dict
            and "AG" in bp1_dict
            and bp2_dict["GC"] != -1
            and bp1_dict["AG"] != -1
        ):
            return bp1_dict["AG"], bp2_dict["GC"]
        elif (
            "AT" in bp2_dict
            and "AC" in bp1_dict
            and bp2_dict["AT"] != -1
            and bp1_dict["AC"] != -1
        ):
            return bp1_dict["AC"], bp2_dict["AT"]
        else:
            # No splice site found
            return True, True


def splice_site_search(in_str, s_site, splice_bin, to_right=True):
    """search for 's_site' in 'in_str'
    search for donor site: from right to left
    search for acceptor site: from left to right
    """
    shift_list = []
    # from left to right searching
    if to_right:
        for i in range(len(in_str) - 1):
            _motif = in_str[i : i + 2]
            if _motif == s_site:
                shift_list.append(i)

    # from right to left searching
    else:
        for i in range(len(in_str) - 1):
            if i == 0:
                _motif = in_str[-i - 2 :]
            else:
                _motif = in_str[-i - 2 : -i]
            if _motif == s_site:
                shift_list.append(i)
    if len(shift_list) == 0:
        return -1
    else:
        ordered_shift_list = sorted(shift_list, key=lambda k: abs(k - (splice_bin - 1)))
        return ordered_shift_list[0] + 1


def cigar_validity(cigar_str):
    """
    40M25N5M then cigartuple is [('40', 'M'), ('25', 'N'), ('5', 'M')]
    """
    cigartuple = list(map(list, re.findall(r"(\d+)(\w)", cigar_str)))
    if cigartuple[0][1] == cigartuple[1][1]:
        cigartuple[1][0] = str(int(cigartuple[0][0]) + int(cigartuple[1][0]))
        del cigartuple[0]

    elif cigartuple[-1][1] == cigartuple[-2][1]:
        cigartuple[-2][0] = str(int(cigartuple[-1][0]) + int(cigartuple[-2][0]))
        del cigartuple[-1]

    valid_cigar = ""
    for i in cigartuple:
        valid_cigar = valid_cigar + i[0] + i[1]
    return valid_cigar


def closest(a, b, tgt):
    """return the closest number to tgt from two numbers (a, b)"""
    if abs(a - tgt) <= abs(b - tgt):
        return a
    else:
        return b


def detect_fusion_from_cigar(
    chrm,
    read,
    mapq_cutoff,
    splice_bin,
    fastafile,
    cvg,
    gene_iv,
    min_dist_between_loci=10000,
    motif_required=True,
):
    """
    Detect linear splicing events
    """

    def strict_inequality(genes_tuple):
        flag = True
        gene1, gene2 = genes_tuple
        for i in gene1.split("&"):
            for j in gene2.split("&"):
                if i == j:
                    flag = False
                    break
        return flag

    strand = "-" if read.is_reverse else "+"

    if not "N" in read.cigarstring:
        return []
    if read.is_supplementary:
        return []
    if read.mapping_quality < mapq_cutoff:
        return []

    cigartuple = re.findall(r"(\d+)(\w)", read.cigarstring)
    pos = read.reference_start
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
            exon_positions.append((start_pos, current_pos))
            current_pos = current_pos + length
            start_pos = current_pos
    exon_positions.append((start_pos, current_pos))

    # No junctions in the representative alignment
    if len(exon_positions) == 1:
        return []

    hit_positions = []
    for i, j in exon_positions:
        hit_positions.extend([i, j])
    hit_positions.sort()
    hit_positions.pop(0)
    hit_positions.pop(-1)
    intron_positions = list(zip(hit_positions[::2], hit_positions[1::2]))

    fusion_list = []
    can_sites = {"GT-AG", "GC-AG", "AT-AC"}
    for start, end in intron_positions:
        if end - start >= min_dist_between_loci:
            if strand == "-":
                left_site = fastafile[chrm][
                    end - 2 - 1 : end - 1
                ].reverse.complement.seq
                right_site = fastafile[chrm][
                    start - 1 : start + 2 - 1
                ].reverse.complement.seq
            else:
                left_site = fastafile[chrm][start - 1 : start + 2 - 1].seq
                right_site = fastafile[chrm][end - 2 - 1 : end - 1].seq
            if f"{left_site}-{right_site}" in can_sites:
                _can = 1
            else:
                _can = 0
            _nls, _anno, __can = splicing_confirmation(
                chrm,
                start,
                chrm,
                end,
                splice_bin,
                fastafile,
                cvg,
                False,
                motif_required,
            )
            _genes = gene_annotation(chrm, start, chrm, end, gene_iv)
            if _nls and strict_inequality(_genes):
                fusion_list.append(
                    (
                        "FUS",
                        _anno,
                        __can,
                        [start, end - start, 0, 0],
                        [strand, strand],
                        [*_genes],
                    )
                )
    return fusion_list


def short_TDUP_or_not(
    chrm, ra_mode, sa_start, sa_end, ins_seq_in_read, fastafile
) -> bool:
    """judge the ins_seq_in_read is a TDUP (TDUP size < reads length) OR novel sequence insertion using
    reference sequence infered from chimeirc alignment start position and indel_size from 'query_offset - target_offset'
    :param chrm: the chromosome
    :param ra_mode: representative alignment mode
    :param sa_start: supplementary alignment reference start position
    :param sa_end: supplementary alignment reference end position
    :param ins_seq_in_read: putative inertion sequence from the read
    :param fastafile: pyfaidx.Fasta object of reference genome (FASTA file)
    :type chrm: str
    :type ra_mode: int
    :type sa_start: int
    :type sa_end: int
    :type ins_seq_in_read: str
    :type fastafile: pyfaidx.Fasta object
    :returns: True if it is a short TDUP
    :rtype: bool
    """
    indel_size = len(ins_seq_in_read)
    if ra_mode == 1:
        ref_seq = fastafile[chrm][sa_start - 10 : sa_start + indel_size].seq
    elif ra_mode == 2:
        ref_seq = fastafile[chrm][sa_end - indel_size : sa_end + 10].seq

    alignment_result = aligner(ins_seq_in_read, ref_seq, method="glocal")[0]
    search_seq = alignment_result.seq1.decode("utf-8")
    target_seq = alignment_result.seq2.decode("utf-8")
    search_seq_len = len(search_seq)
    target_seq_len = len(target_seq)
    search_start, search_end = alignment_result.start1, alignment_result.end1 - 1
    target_start, target_end = alignment_result.start2, alignment_result.end2 - 1
    aln_len = search_end - search_start + 1
    total_mismatches = len(search_seq) - aln_len + alignment_result.n_mismatches
    if total_mismatches <= 3:
        return True
    else:
        return False


def test_is_connected(sms_read1, sms_read2, len_cutoff=30) -> tuple:
    """Test whether two SMS tuples of chimeric reads can be connected or not.
    :param sms_read1: triple tuple for (left soft-clipped length, middle read matched size, right softclipped length) of read1 OR path
    :type sms_read1: tuple
    :param sms_read2: triple tuple for (left soft-clipped length, middle read matched size, right softclipped length) of read2 OR path
    :type sms_read2: tuple
    :param len_cutoff: length cutoff to determine the S-M match
    :type len_cutoff: int
    :return: is_connected flag, summed 'SMS' value, mode for read1 and read2
    :rtype: tuple
    """

    def propinquity(a, b, len_cutoff) -> bool:
        if abs(a - b) <= len_cutoff:
            return True
        else:
            return False

    def init_mode_judge(sms) -> int:
        _lt, _read_match, _rt = sms
        # SM
        if _lt > _rt:
            return 2
        # MS
        else:
            return 1

    _lt_len_r1, _read_match_r1, _rt_len_r1 = sms_read1
    _lt_len_r2, _read_match_r2, _rt_len_r2 = sms_read2

    mode_r1 = init_mode_judge(sms_read1)
    mode_r2 = init_mode_judge(sms_read2)

    is_connected = False
    out_lt_len = 0
    out_rt_len = 0
    out_read_match = 0
    if propinquity(_read_match_r1, _lt_len_r2, len_cutoff):
        out_lt_len = 0
        out_read_match = _lt_len_r2 + _read_match_r2
        out_rt_len = _rt_len_r2
        mode_r2 = 2
        is_connected = True

    elif propinquity(_read_match_r1, _rt_len_r2, len_cutoff):
        out_lt_len = _lt_len_r2
        out_read_match = _read_match_r2 + _rt_len_r2
        out_rt_len = 0
        mode_r2 = 1
        is_connected = True

    elif propinquity(_read_match_r2, _lt_len_r1, len_cutoff):
        out_lt_len = 0
        out_read_match = _lt_len_r1 + _read_match_r1
        out_rt_len = _rt_len_r1
        mode_r1 = 2
        is_connected = True

    elif propinquity(_read_match_r2, _rt_len_r1, len_cutoff):
        out_lt_len = _rt_len_r1
        out_read_match = _read_match_r1 + _rt_len_r1
        out_rt_len = 0
        mode_r1 = 1
        is_connected = True

    else:
        out_lt_len = 0
        out_read_match = 0
        out_rt_len = 0
        is_connected = False

    return is_connected, (out_lt_len, out_read_match, out_rt_len), (mode_r1, mode_r2)


def chimeric_aln_order_finder(aln_list, soft_len_cutoff=30) -> tuple:
    """
    :param aln_list: list of Read object
    :type aln_list: list
    :return: Read-to-Read chain, a dictionary of Read-pair(Read1, Read2) => mode-of-Read1, mode-of-Read2
    :rtype: tuple
    .. note::
        Read-to-Read chain scenarios
        * [[Read1, Read2, Read3]]
        * [[Read1, Read2, Read3],[Read4,Read5]]

        Dictionary of Read-pair scenarios
        * (Read1, Read2) => mode-of-Read1, mode-of-Read2
        * (Read2, Read1) => mode-of-Read2, mode-of-Read1
    """

    def two_overlapped_lists(a, b) -> list:
        """
        :param a: list of Reads
        :param b: list of Reads
        :type a: list
        :type b: list
        :return: merged list (union)
        :rtype: list
        """
        count = 0
        overlap_len = 0
        first_hit_index_a = 0
        merged_list = []
        if len(a) == 0 or len(b) == 0:
            return []
        for i in range(len(a)):
            if a[i] == b[overlap_len]:
                overlap_len += 1
                if count == 0:
                    first_hit_index_a = i
                count += 1
        if len(a) - first_hit_index_a == overlap_len:
            merged_list.extend(a)
            merged_list.extend(b[overlap_len:])
        return merged_list

    def select_path_for_one_node(node, paths) -> tuple:
        """
        :param node: start node
        :param paths: linked paths (list of Path) for node
        :type node: Read
        :type paths: list
        :return: end-to-end chain for one node and the corresponding reads-pair mode dict
        :rtype: tuple
        .. note::
            * 1) select one path with highest length
            * 2) if the length is the same, select the one with smaller value of (number of mismatches)
        """
        if len(paths) == 0:
            return [], {}
        len_dict = {}
        for i in paths:
            len_dict[i] = len(i)
        nm_max_len_dict = dict(
            (j, j.nm) for j in len_dict if len_dict[j] == max(len_dict.values())
        )
        min_nm_max_len_path = min(nm_max_len_dict, key=lambda k: k.nm)
        if min_nm_max_len_path:
            return [node] + min_nm_max_len_path.nodes, min_nm_max_len_path.mode

    def optimal_path_selector(lt_node, rt_node, lt_paths, rt_paths) -> list:
        """
        :param lt_node: start node (left side)
        :param rt_node: start node (right side)
        :param lt_paths: linked paths (list of Path) for lt_node
        :param rt_paths: linked paths (list of Path) for rt_node
        :type lt_node: Read
        :type rt_node: Read
        :type lt_paths: list
        :type rt_paths: list
        :return: end-to-end chains
        :rtype: list (list of lists)
        """
        end_to_end_chain = []
        reads_pair_mode_dict = {}
        if len(lt_paths) > 0 and len(rt_paths) > 0:
            # check if lt_paths can overlap with rt_paths
            # common boundary nodes or common boundary edges
            overlap_path_flag = False
            overlap_path = []
            for i in lt_paths:
                for j in rt_paths:
                    j_rev = j.nodes[::-1]
                    _overlap_path = two_overlapped_lists(i.nodes, j_rev)
                    if _overlap_path:
                        overlap_path_flag = True
                        overlap_path = _overlap_path
                        reads_pair_mode_dict.update(i.mode)
                        reads_pair_mode_dict.update(j.mode)
                        break
            # overlapped path can be found
            if overlap_path_flag:
                end_to_end_chain.append([lt_node] + overlap_path + [rt_node])
            # if no overlapped path can be found
            # select the longest path in lt_paths and rt_paths, respectively
            else:
                lt_chain, lt_mode_dict = select_path_for_one_node(lt_node, lt_paths)
                rt_chain, rt_mode_dict = select_path_for_one_node(rt_node, rt_paths)
                if lt_chain:
                    end_to_end_chain.append(lt_chain)
                    reads_pair_mode_dict.update(lt_mode_dict)
                if rt_chain:
                    end_to_end_chain.append(rt_chain)
                    reads_pair_mode_dict.update(rt_mode_dict)
        # either lt_paths or rt_paths is empty
        else:
            # select the longest path in lt_paths and rt_paths, respectively
            lt_chain, lt_mode_dict = select_path_for_one_node(lt_node, lt_paths)
            rt_chain, rt_mode_dict = select_path_for_one_node(rt_node, rt_paths)
            if lt_chain:
                end_to_end_chain.append(lt_chain)
                reads_pair_mode_dict.update(lt_mode_dict)
            if rt_chain:
                end_to_end_chain.append(rt_chain)
                reads_pair_mode_dict.update(rt_mode_dict)
        return end_to_end_chain, reads_pair_mode_dict

    ######################################################################
    start_nodes = []
    candidate_nodes = []
    reads_pair_mode_dict = {}
    end_to_end_chain = []

    for read in aln_list:
        if read.lt_soft_len < soft_len_cutoff or read.rt_soft_len < soft_len_cutoff:
            start_nodes.append(read)
        else:
            candidate_nodes.append(read)

    # print('start_nodes: ', start_nodes)
    # print('candidate_nodes: ', candidate_nodes)

    if len(start_nodes) == 2:
        # print('len(start_nodes) == 2')
        if len(candidate_nodes) == 0:
            _is_connected, _, _mode = test_is_connected(
                start_nodes[0].sms, start_nodes[1].sms
            )
            if _is_connected:
                reads_pair_mode_dict[(start_nodes[0], start_nodes[1])] = _mode
                end_to_end_chain.append(start_nodes)
        elif len(candidate_nodes) == 1:
            for _node in start_nodes:
                _is_connected, _, _mode = test_is_connected(
                    _node.sms, candidate_nodes[0].sms
                )
                if _is_connected:
                    reads_pair_mode_dict[(_node, candidate_nodes[0])] = _mode
            tmpList = []
            for i, j in reads_pair_mode_dict:
                if not tmpList:
                    tmpList.append(i)
                    tmpList.append(j)
                else:
                    last_item = tmpList[-1]
                    if i == last_item:
                        tmpList.append(j)
                    elif j == last_item:
                        tmpList.append(i)
            end_to_end_chain.append(tmpList)
            tmpList = None
        elif len(candidate_nodes) > 1:
            # searching linked paths from start node 1
            # start_nodes[0]
            tgt_node = start_nodes[0]
            stop_signal = False
            count = 0
            while not stop_signal:
                stop_signal = True
                for _node in candidate_nodes:
                    if count < len(candidate_nodes):
                        _is_connected, _sum_sms, _mode = test_is_connected(
                            tgt_node.sms, _node.sms
                        )
                        if _is_connected:
                            stop_signal = False
                            # reads_pair_mode_dict[(tgt_node, _node)] = _mode
                            if not tgt_node.linked_paths:
                                tmp_path = Path()
                                tmp_path.add(_node)
                                tmp_path.add_mode({(tgt_node, _node): _mode})
                                tmp_path.sms = _sum_sms
                                tgt_node.add_path(tmp_path)
                            else:
                                # check if _node was already the last node in any of the existing path
                                flag = False
                                for j in tgt_node.linked_paths:
                                    if _node == j.nodes[-1]:
                                        flag = True
                                # add the _node to a newly created path, then add this path to the tgt_node (linked_paths)
                                if not flag:
                                    tmp_path = Path()
                                    tmp_path.add(_node)
                                    tmp_path.add_mode({(tgt_node, _node): _mode})
                                    tmp_path.sms = _sum_sms
                                    tgt_node.add_path(tmp_path)
                    else:
                        # print(tgt_node.linked_paths)
                        for path in tgt_node.linked_paths:
                            # print(path.sms, _node.sms)
                            _is_connected, _sum_sms, _mode = test_is_connected(
                                path.sms, _node.sms
                            )
                            if _is_connected:
                                if path.nodes[-1] != _node:
                                    stop_signal = False
                                    path.add(_node)
                                    path.add_mode({(path.nodes[-1], _node): _mode})
                                    path.sms = _sum_sms
                                    # reads_pair_mode_dict[(path.nodes[-1], _node)] = _mode
                    count += 1
                if stop_signal:
                    break

            lt_candidate_paths = copy.deepcopy(tgt_node.linked_paths)
            lt_node = tgt_node
            # searching linked paths from start node 2
            # start_nodes[1]
            tgt_node = start_nodes[1]
            stop_signal = False
            count = 0
            while not stop_signal:
                stop_signal = True
                for _node in candidate_nodes:
                    if count < len(candidate_nodes):
                        _is_connected, _sum_sms, _mode = test_is_connected(
                            tgt_node.sms, _node.sms
                        )
                        if _is_connected:
                            stop_signal = False
                            # reads_pair_mode_dict[(tgt_node, _node)] = _mode
                            if not tgt_node.linked_paths:
                                tmp_path = Path()
                                tmp_path.add(_node)
                                tmp_path.add_mode({(tgt_node, _node): _mode})
                                tmp_path.sms = _sum_sms
                                tgt_node.add_path(tmp_path)
                            else:
                                # check if _node was already the last node in any of the existing path
                                flag = False
                                for j in tgt_node.linked_paths:
                                    if _node == j.nodes[-1]:
                                        flag = True
                                # add the _node to a newly created path, then add this path to the tgt_node (linked_paths)
                                if not flag:
                                    tmp_path = Path()
                                    tmp_path.add(_node)
                                    tmp_path.add_mode({(tgt_node, _node): _mode})
                                    tmp_path.sms = _sum_sms
                                    tgt_node.add_path(tmp_path)

                    else:
                        for path in tgt_node.linked_paths:
                            _is_connected, _sum_sms, _mode = test_is_connected(
                                path.sms, _node.sms
                            )
                            if _is_connected:
                                if path.nodes[-1] != _node:
                                    stop_signal = False
                                    path.add(_node)
                                    path.add_mode({(path.nodes[-1], _node): _mode})
                                    path.sms = _sum_sms
                                    # reads_pair_mode_dict[(path.nodes[-1], _node)] = _mode

                    count += 1
                if stop_signal:
                    break
            rt_candidate_paths = copy.deepcopy(tgt_node.linked_paths)
            rt_node = tgt_node
            # print(lt_node)
            # print(rt_node)
            # print(lt_candidate_paths)
            # print(rt_candidate_paths)
            end_to_end_chain, reads_pair_mode_dict = optimal_path_selector(
                lt_node, rt_node, lt_candidate_paths, rt_candidate_paths
            )

    return end_to_end_chain, reads_pair_mode_dict


def output_bedpe_file(sr_dict, group_dict, prefix, splice_bin):
    """
    :param sr_dict: sv candidate to number of supporting reads(SR) dictionary
    :param group_dict: sv candidate to group of events dictionary, connected chimeirc reads are included in one group
    :param prefix: output file prefix
    :param splice_bin: bin size for splice site searching
    :type sr_dict: dict
    :type group_dict: dict
    :type prefix: str
    :type splice_bin: int
    :return: output BEDPE file
    :rtype: str
    .. note::
    """
    output = open("{}.sv.bedpe".format(prefix), "w")
    count = 0
    for key in sr_dict:
        sr = sr_dict[key]
        num_of_group = group_dict[key]
        _type, _can, bp1, bp2, strands = key.split("\t")
        chrm1, pos1 = bp1.split(":")
        chrm2, pos2 = bp2.split(":")
        strand1 = strands[0]
        strand2 = strands[1]
        pos1 = int(pos1)
        pos2 = int(pos2)
        if pos1 - splice_bin > 0 and pos2 - splice_bin > 0:
            output.write(
                f"{chrm1}\t{pos1-splice_bin}\t{pos1+splice_bin}\t{chrm2}\t{pos2-splice_bin}\t{pos2+splice_bin}\tgroup_{num_of_group}\t{sr}\t{strand1}\t{strand2}\n"
            )
    output.close()
    return output


def aggregate_candidates(in_dict, len_cutoff=10):
    if len_cutoff == 0:
        return in_dict
    else:
        discarded_items = set([])
        L = len(in_dict)
        items = list(in_dict.keys())
        for i in range(L):
            for r2 in items[i + 1 :]:
                r1 = items[i]
                if similar_hit(r1, r2, len_cutoff):
                    can_1 = r1.split("\t")[1]
                    can_2 = r2.split("\t")[1]
                    ao_1 = in_dict[r1]
                    ao_2 = in_dict[r2]
                    new_ao = ao_1 + ao_2
                    if ao_1 > ao_2:
                        in_dict[r1] = new_ao
                        discarded_items.add(r2)
                    elif ao_1 < ao_2:
                        in_dict[r2] = new_ao
                        discarded_items.add(r1)
                    elif ao_1 == ao_2:
                        if can_1 == 1:
                            in_dict[r1] = new_ao
                            discarded_items.add(r2)
                        elif can_2 == 1:
                            in_dict[r2] = new_ao
                            discarded_items.add(r1)
        out_dict = {}
        for m in in_dict:
            if not m in discarded_items:
                out_dict[m] = in_dict[m]
        return out_dict


def similar_hit(r1, r2, len_cutoff=10):
    r1_type, r1_can, A1, A2, strand_1 = r1.split("\t")
    r2_type, r2_can, B1, B2, strand_2 = r2.split("\t")

    chrm_a, pos_a = A1.split(":")
    chrm_b, pos_b = A2.split(":")

    chrm_A, pos_A = B1.split(":")
    chrm_B, pos_B = B2.split(":")

    if r1_type != r2_type:
        return False
    else:
        if chrm_a == chrm_A and chrm_b == chrm_B:
            if (
                abs(int(pos_a) - int(pos_A)) <= len_cutoff
                and abs(int(pos_b) - int(pos_B)) <= len_cutoff
            ):
                return True
            else:
                return False

        elif chrm_a == chrm_B and chrm_b == chrm_A:
            if (
                abs(int(pos_a) - int(pos_B)) <= len_cutoff
                and abs(int(pos_b) - int(pos_A)) <= len_cutoff
            ):
                return True
            else:
                return False
        else:
            return False


def vcf_header(output_prefix, bam_header):
    """output VCF header using information from BAM header
    :param output_prefix: file prefix for VCF file, usually uses sample name
    :param bam_header: header of BAM file
    :type output_prefix: str
    :type bam_header: str
    :return: header of VCF file
    :rtype: str
    .. note::
        GROUP field is important,
        multiple NLS events with the same 'group' information will form one transcript
    """

    _aligners = {
        "CLC",
        "ContextMap2",
        "CRAC",
        "GSNAP",
        "HISAT",
        "HISAT2",
        "MapSplice2",
        "Novoalign",
        "OLego",
        "RUM",
        "SOAPsplice",
        "STAR",
        "Subread",
        "TopHat",
        "TopHap2",
        "bwa",
        "bowtie",
        "bowtie2",
        "minimap2",
    }
    avail_aligners = set([x.upper() for x in _aligners])
    if "PG" in bam_header:
        for j in bam_header["PG"]:
            if j["ID"].upper() in avail_aligners:
                aln_cmd = j["CL"]
                break
            else:
                aln_cmd = "Unknown"
    else:
        aln_cmd = "Unknown"

    header = ["##fileformat=VCFv4.1"]
    header.append("##source=ScanNLS " + __version__)
    header.append(
        '##reference=<CMD={},Description="Alignment parameters">'.format(aln_cmd)
    )
    header.append('##ALT=<ID=TDUP,Description="Tandem duplication">')
    header.append('##ALT=<ID=INV,Description="Inversion">')
    header.append('##ALT=<ID=TRA,Description="Translocation">')
    header.append(
        '##INFO=<ID=CANONICAL,Number=0,Type=Flag,Description="Canonical splice site">'
    )
    header.append(
        '##INFO=<ID=NONCANONICAL,Number=0,Type=Flag,Description="Noncanonical splice site">'
    )
    header.append(
        '##INFO=<ID=BOUNDARY,Number=1,Type=String,Description="The coding exon boundary type of event, BOTH, LEFT, RIGHT, NEITHER.">'
    )
    header.append(
        '##INFO=<ID=DP1,Number=1,Type=Integer,Description="Total read depth at the breakpoint1">'
    )
    header.append(
        '##INFO=<ID=DP2,Number=1,Type=Integer,Description="Total read depth at the breakpoint2">'
    )
    header.append(
        '##INFO=<ID=SR,Number=1,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">'
    )
    header.append(
        '##INFO=<ID=PSO,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1], representing the ratio of reads showing the alternative allele to all reads">'
    )
    header.append(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The type of event, INS, DEL, TDUP, INV, TRA.">'
    )
    header.append(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'
    )
    header.append(
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">'
    )
    header.append(
        '##INFO=<ID=GROUP,Number=1,Type=String,Description="events in the group from the same transcript">'
    )
    header.append(
        '##INFO=<ID=GENE1,Number=1,Type=String,Description="Overlapped coding gene for breakpoint1">'
    )
    header.append(
        '##INFO=<ID=GENE2,Number=1,Type=String,Description="Overlapped coding gene for breakpoint2">'
    )
    header.append(
        '##INFO=<ID=STRAND1,Number=1,Type=String,Description="Strand for breakpoint1">'
    )
    header.append(
        '##INFO=<ID=STRAND2,Number=1,Type=String,Description="Strand for breakpoint2">'
    )
    header.append(
        '##INFO=<ID=MODE1,Number=1,Type=String,Description="mode for softclipped reads at breakpoint1">'
    )
    header.append(
        '##INFO=<ID=MODE2,Number=1,Type=String,Description="mode for softclipped reads at breakpoint2">'
    )
    header.append(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="nd position of the structural variant">'
    )
    header.append(
        '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">'
    )
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + output_prefix
    )
    return "\n".join(header)


def sv_checker(read, ref_site, mapq_cutoff, sv_len_cutoff):
    """"""
    chimeric_aln = read.alignment.get_tag("SA")[:-1].split(";")
    if len(chimeric_aln) == 1:
        chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = chimeric_aln[0].split(",")
    else:
        chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = multiple_sa_tag_selector(
            chimeric_aln
        ).split(",")

    if int(mapq_sa) < mapq_cutoff:
        return False
    (
        sv_type,
        sv_anno_and_can,
        sv_pos,
        sv_length,
        read_mode,
        sa_mode,
        strands,
        genes,
    ) = read.alignment.get_tag("SV").split(",")

    if sv_type == "TRA":
        if int(sv_pos) - 2 + (int(read_mode) - 1) == ref_site:
            return True
        else:
            return False
    # TDUP and INV
    else:
        if int(sv_length) < sv_len_cutoff:
            return False
        if int(sv_pos) - 2 + (int(read_mode) - 1) == ref_site:
            return True
        elif int(sv_pos) - 2 + int(sv_length) + (int(read_mode) - 1) == ref_site:
            return True
        else:
            return False
