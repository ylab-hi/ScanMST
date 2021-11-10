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

from align import aligner
from Bio import SearchIO
from Bio.Seq import Seq
from pyfaidx import Fasta

from . import __version__
from .classes import Path
from .common import transcript_upstream_part_determiner

try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)

__funcs__ = {"extract_splice_sites", "junc_site_checker"}



def extract_splice_sites(in_file, bin_size) -> tuple:
    """Extract splice sites and gene regions from input GTF file
    :param in_file: gene annotation file (GTF file)
    :param bin_size: bin size to search splice site
    :type in_file: str
    :type bin_size: int
    :return: annotated splice sites (HTSeq.GenomicArrayOfSets) and annotated gene regions (HTSeq.GenomicArrayOfSets)
    :rtype: tuple
    """
    # todo using real splice sites from reference genome
    gtf_file = HTSeq.GFF_Reader(in_file)
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
                        feature.iv.start - bin_size,
                        feature.iv.end + bin_size,
                        ".",
                    )
                ] += str(gene_name)

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
                    first_exon.start - bin_size,
                    first_exon.start + bin_size,
                    ".",
                )
            ] += "XX"
            cvg[
                HTSeq.GenomicInterval(
                    first_exon.chrom,
                    first_exon.end - bin_size,
                    first_exon.end + bin_size,
                    ".",
                )
            ] += "XX"
        elif exon_num == 2:
            if strand == "+":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin_size,
                        first_exon.start + bin_size,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin_size,
                        first_exon.end + bin_size,
                        ".",
                    )
                ] += "GT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin_size,
                        last_exon.start + bin_size,
                        ".",
                    )
                ] += "AG"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.end - bin_size,
                        last_exon.end + bin_size,
                        ".",
                    )
                ] += "XX"
            elif strand == "-":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin_size,
                        first_exon.start + bin_size,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin_size,
                        first_exon.end + bin_size,
                        ".",
                    )
                ] += "CT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin_size,
                        last_exon.start + bin_size,
                        ".",
                    )
                ] += "AC"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.end - bin_size,
                        last_exon.end + bin_size,
                        ".",
                    )
                ] += "XX"
        # exon_num > 2
        else:
            if strand == "+":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin_size,
                        first_exon.start + bin_size,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin_size,
                        first_exon.end + bin_size,
                        ".",
                    )
                ] += "GT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin_size,
                        last_exon.start + bin_size,
                        ".",
                    )
                ] += "AG"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.end - bin_size,
                        last_exon.end + bin_size,
                        ".",
                    )
                ] += "XX"
            elif strand == "-":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin_size,
                        first_exon.start + bin_size,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin_size,
                        first_exon.end + bin_size,
                        ".",
                    )
                ] += "CT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin_size,
                        last_exon.start + bin_size,
                        ".",
                    )
                ] += "AC"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.end - bin_size,
                        last_exon.end + bin_size,
                        ".",
                    )
                ] += "XX"
            for _exon in exonList[1:-1]:
                iv1 = HTSeq.GenomicInterval(
                    _exon.chrom, _exon.start - bin_size, _exon.start + bin_size, "."
                )
                iv2 = HTSeq.GenomicInterval(
                    _exon.chrom, _exon.end - bin_size, _exon.end + bin_size, "."
                )
                if strand == "+":
                    cvg[iv1] += "AG"
                    cvg[iv2] += "GT"
                elif strand == "-":
                    cvg[iv1] += "AC"
                    cvg[iv2] += "CT"
    sys.stdout.write("{} is fully loaded!\n".format(in_file))
    return cvg, gene_iv


def gene_annotation(chrm1, pos1, chrm2, pos2, gene_iv) -> tuple:
    """obtain gene annotations for breakpoints
    :param chrm1: chromosome for breakpoint1
    :param chrm2: chromosome for breakpoint2
    :param pos1: position for breakpoint1
    :param pos2: position for breakpoint2
    :param gene_iv: gene annotations in HTSeq.GenomicArrayOfSets
    :type chrm1: str
    :type chrm2: str
    :type pos1: int
    :type pos2: int
    :type gene_iv: HTSeq.GenomicArrayOfSets
    :return: overlapped genes for breakpoints
    :rtype: tuple
    """
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
    genome_fasta,
    cvg,
    strand_changed,
    motif_required=True,
) -> tuple:
    """Judge whether the breakpoints are NLS events or not
    if motif_required is ON: it will only report NLS events with 'canonical splice sites';
    otherwise: it will report NLS events whatever the splice sites they used
    :param chrm1: chromosome for breakpoint1
    :param chrm2: chromosome for breakpoint2
    :param pos1: position for breakpoint1
    :param pos2: position for breakpoint2
    :param splice_bin: bin size for splice sites searching
    :param genome_fasta: reference genome (pyfaidx.Fasta object)
    :param cvg: splice site annotations (HTSeq.GenomicArrayOfSets)
    :param strand_changed: whether breakpoint1 and breakpoint2 use the same strand or not
    :param motif_required: canonical splice sites required; if True: considering canonical splice sites only; else: considering canonical and noncanonical splice sites both
    :type chrm1: str
    :type chrm2: str
    :type pos1: int
    :type pos2: int
    :type splice_bin: int
    :type genome_fasta: pyfaidx.Fasta
    :type cvg: HTSeq.GenomicArrayOfSets
    :type strand_changed: bool
    :type motif_required: bool
    :return: report/not report, overlapping boundary in bits, canonical splice site/noncanonical splice site
    :rtype: tuple

    .. note::
        Possible output scenarios
        * True,  3(11), 1 => reported, both breakpoints overlap with annotated coding exons boundary, using canonical splice motif
        * True,  2(10), 0 => reported, one breakpoint overlap with annotated coding exons boundary, using noncanonical splice motif
        * True,  1(01), 0 => reported, one breakpoint overlap with annotated coding exons boundary, using noncanonical splice motif
        * False, 0(00), 0 => not reported, none breakpoint overlap with annotated coding exons boundary, using noncanonical splice motif
    """

    def canonical_site_finder(in_seq) -> list:
        """find canonical splice sites in the input sequence
        :param in_seq: input sequence (usually sequence nearby the breakpoints)
        :type in_seq: str
        :return: a list of canonical splice sites in the input sequence, it can be a empty list
        :rtype: list
        """
        candidate_sites = ["GT", "AG", "CT", "AC"]
        matches = (i in in_seq for i in candidate_sites)
        hit_sites = [j for i, j in zip(matches, candidate_sites) if i]
        return hit_sites

    def splice_paired_checker(hit_sites, pair_seq, splice_motif_dict) -> bool:
        """find canonical splice sites in the input sequence
        :param hit_sites: canonical splice site at one end
        :param pair_seq: sequence at the other pair end
        :param splice_motif_dict: paired splice sites (same strand or different strand)
        :type hit_sites: list
        :type pair_seq: str
        :type splice_motif_dict: dict
        :return: canonical splice sites are paired or not
        :rtype: bool
        """
        paired = False
        for site in hit_sites:
            if site in splice_motif_dict and splice_motif_dict[site] in pair_seq:
                paired = True
                break
        return paired

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
        junc_seq1 = genome_fasta[chrm1][pos1 - splice_bin : pos1 + splice_bin].seq
        junc_seq2 = genome_fasta[chrm2][pos2 - splice_bin : pos2 + splice_bin].seq
        _junc1 = canonical_site_finder(junc_seq1)
        _junc2 = canonical_site_finder(junc_seq2)
        if splice_paired_checker(_junc1, junc_seq2, splice_motif_dict):
            if motif_required:
                return True, 0, 1
            else:
                return True, 0, 1
        elif splice_paired_checker(_junc2, junc_seq1, splice_motif_dict):
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
        junc_seq = genome_fasta[chrm2][pos2 - splice_bin : pos2 + splice_bin].seq
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
        junc_seq = genome_fasta[chrm1][pos1 - splice_bin : pos1 + splice_bin].seq
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
) -> tuple:
    """
    Update the breakpoints of NLS events with canonical splice sites
    keep NLS events with noncanonical splice sites unchanged (GT-AG;GC-AG;AT-AC)
    :param bp1_chrm: chromosome for breakpoint1
    :param bp1_pos: position for breakpoint1
    :param bp2_chrm: chromosome for breakpoint2
    :param bp2_pos: position for breakpoint2
    :param bp1_strand: strand for breakpoint1
    :param bp2_strand: strand for breakpoint2
    :param bp1_mode: mode for breakpoint1
    :param bp2_mode: mode for breakpoint2
    :param splice_bin: bin size for splice site searching
    :param genome_fasta: reference genome (pyfaidx.Fasta)
    :type bp1_chrm: str
    :type bp1_pos: int
    :type bp2_chrm: str
    :type bp2_pos: int
    :type bp1_strand: str (-/+)
    :type bp2_strand: str (-/+)
    :type bp1_mode: int (1/2)
    :type bp2_mode: int (1/2)
    :type splice_bin: int
    :type genome_fasta: pyfaidx.Fasta object
    :return: updated position for breakpoint1 and updated position for breakpoint2
    :rtype: tuple
    ..note :
        mode moving rules:
        * mode (SM) [2]: breakpoint move to right
        * mode (MS) [1]: breakpoint move to left
    """

    def splice_site_search(in_str, s_site, splice_bin, for_acceptor=True):
        """search for 's_site' in 'in_str';
        searching for acceptor site (searching from left to right [==>> ...XXX])
        searching for donor site (searching from right to left [XXX... <<==])
        :param in_str: nearby sequence of breakpoints
        :param s_site: splice site
        :param splice_bin: bin size for splice site searching
        :param for_acceptor: searching for acceptor site [True] OR donor site[False]
        :type in_str: str
        :type s_site: str
        :type splice_bin: int
        :type for_acceptor: bool
        :return: optimal splice site position in 'in_str'; not found(-1)
        :rtype: int
        ..note ::
             For donor, splice position will be ZXXXXGTYYY
                                                |||||^
             For acceptor, splice position will be ZXXXAGTYYY
                                                        ^||||
             assert splice_site_search('AGXXAGTXPX', 'AG', 5, True) == 5
             assert splice_site_search('AGXXAGTXPX', 'GT', 5, False) == 4
        """
        shift_positions = []
        if for_acceptor:  # ==>>
            for i in range(len(in_str) - 1):
                _motif = in_str[i : i + 2]
                if _motif == s_site:
                    shift_positions.append(i)
        else:  # <<==
            for i in range(len(in_str) - 1):
                if i == 0:
                    _motif = in_str[-i - 2 :]
                else:
                    _motif = in_str[-i - 2 : -i]
                if _motif == s_site:
                    shift_positions.append(i)
        # No found canonical splice site
        if len(shift_positions) == 0:
            return -1
        else:
            # select position closest to the center point of the 'in_str'
            ordered_shift_positions = sorted(
                shift_positions, key=lambda k: abs(k - (splice_bin - 1))
            )
            return ordered_shift_positions[0] + 1

    def obtain_bps_shift_len(bp1_dict, bp2_dict, bp1_is_upstream=True) -> tuple:
        """obtain shift length for breakpoint1 and breakpoint2 according to putative canonical splice
        site positions in 'bp1_dict' and 'bp2_dict'
        In order to find the correct breakpoint pair, canonical splice site matches is needed (GT-AG, GC-AG, AT-AC)
        :param bp1_dict: breakpoint1 splice site positions
        :param bp2_dict: breakpoint2 splice site positions
        :param bp1_is_upstream: wheather sequence at breakpoint1 is the upstream segment of the transcript
        :type bp1_dict: dict
        :type bp2_dict: dict
        :type bp1_is_upstream: bool
        :return: shift length for breakpoint1, shift length for breakpoint2
        :rtype: tuple
        """
        if bp1_is_upstream:
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
                return None, None

    bp1_pos_dict = {}
    if bp1_mode == 2:
        if bp1_strand == "+":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].seq
            # acceptor site: AG/AC
            bp1_pos_dict["AG"] = splice_site_search(boundary_seq1, "AG", splice_bin)
            bp1_pos_dict["AC"] = splice_site_search(boundary_seq1, "AC", splice_bin)
        elif bp1_strand == "-":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].reverse.complement.seq
            # donor site: GT/GC/AT
            bp1_pos_dict["GT"] = splice_site_search(
                boundary_seq1, "GT", splice_bin, False
            )
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
            # donor site: GT/GC/AT
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
            # acceptor site: AG/AC
            bp1_pos_dict["AG"] = splice_site_search(boundary_seq1, "AG", splice_bin)
            bp1_pos_dict["AC"] = splice_site_search(boundary_seq1, "AC", splice_bin)

    bp2_pos_dict = {}
    if bp2_mode == 2:
        if bp2_strand == "+":
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].seq
            # acceptor site: AG/AC
            bp2_pos_dict["AG"] = splice_site_search(boundary_seq2, "AG", splice_bin)
            bp2_pos_dict["AC"] = splice_site_search(boundary_seq2, "AC", splice_bin)
        elif bp2_strand == "-":
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].reverse.complement.seq
            # donor site: GT/GC/AT
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
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].seq
            # donor site: GT/GC/AT
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
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].reverse.complement.seq
            # acceptor site: AG/AC
            bp2_pos_dict["AG"] = splice_site_search(boundary_seq2, "AG", splice_bin)
            bp2_pos_dict["AC"] = splice_site_search(boundary_seq2, "AC", splice_bin)

    shift1 = 0
    shift2 = 0
    bp1_is_upstream = transcript_upstream_part_determiner(
        bp1_strand, bp2_strand, bp1_mode, bp2_mode
    )
    shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, bp1_is_upstream)

    # No canonical splice sites found, so change the annotation to noncanonical splice site
    if shift1 == None and shift2 == None:
        return 0, 0

    new_bp1_pos = bp1_pos
    new_bp2_pos = bp2_pos
    # [SM:2] breakpoint move to right
    # [MS:1] breakpoint move to left
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
    reference sequence infered from chimeric alignment start position and indel_size from 'query_offset - target_offset'
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


def test_is_connected(
    sms_read1, sms_read2, seq_read1, seq_read2, allowed_difference=80, SM_align=False
) -> tuple:
    """Test whether two SMS tuples of chimeric reads can be connected or not.
    :param sms_read1: triple tuple for (left soft-clipped length, middle read matched size, right softclipped length) of read1 OR path
    :type sms_read1: tuple
    :param sms_read2: triple tuple for (left soft-clipped length, middle read matched size, right softclipped length) of read2 OR path
    :type sms_read2: tuple
    :param seq_read1: reads sequence (as it is stored in the BAM file) of read2 OR path
    :type seq_read1: str
    :param seq_read2: reads sequence (as it is stored in the BAM file) of read2 OR path
    :type seq_read2: str
    :param allowed_difference: the difference of read_match_size (Read1) and softclipped length (Read2) to determine the S-M match
    :type allowed_difference: int
    :param SM_align: whether use the alignment of softclipped segment of one read and matched segment of another read
    :type SM_align: bool
    :return: is_connected flag, summed 'SMS' value, summed 'Reads' sequence, mode for read1 and read2
    :rtype: tuple
    """

    def propinquity(a, b, allowed_difference) -> bool:
        if abs(a - b) <= allowed_difference:
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

    def SM_alignment(query_seq, target_seq, same_strand=True) -> bool:
        """matched segment of read1 align to softclipped segment of read2
        solve the problem that the length of the left softclipped segment and the right softclipped segment may be quite similar, eliminate the ambiguity of connected reads
        :param query_seq: matched segment of read1
        :param target_seq: softclipped segment of read2
        :param same_strand: read1 and read2 one the same strand or not
        :type query_seq: str
        :type target_seq: str
        :type same_strand: bool
        """
        if not same_strand:
            target_seq = str(Seq(target_seq).reverse_complement())
        allowed_mismatches = abs(len(target_seq) - len(query_seq))
        # print('allowed_mismatches:', allowed_mismatches)
        alignment_result = aligner(query_seq, target_seq, method="glocal")[0]
        _query_seq = alignment_result.seq1.decode("utf-8")
        _target_seq = alignment_result.seq2.decode("utf-8")
        _query_seq_len = len(_query_seq)
        _target_seq_len = len(_target_seq)
        _query_start, _query_end = alignment_result.start1, alignment_result.end1 - 1
        _target_start, _target_end = alignment_result.start2, alignment_result.end2 - 1
        aln_len = _query_end - _query_start + 1
        total_mismatches = len(_query_seq) - aln_len + alignment_result.n_mismatches

        # print('query_seq :', query_seq,  'target_seq :', target_seq)
        # print('total_mismatches :', total_mismatches, _query_seq, _target_seq)

        if total_mismatches <= allowed_mismatches:
            return True
        else:
            return False

    _lt_len_r1, _read_match_r1, _rt_len_r1 = sms_read1
    _lt_len_r2, _read_match_r2, _rt_len_r2 = sms_read2

    same_strand = None
    if seq_read1 == seq_read2:
        same_strand = True
    else:
        same_strand = False

    mode_r1 = init_mode_judge(sms_read1)
    mode_r2 = init_mode_judge(sms_read2)

    is_connected = False
    out_lt_len = 0
    out_rt_len = 0
    out_read_match = 0
    out_seq = ""

    if propinquity(_read_match_r1, _lt_len_r2, allowed_difference):
        if SM_align:
            if SM_alignment(
                seq_read1[_lt_len_r1 : _lt_len_r1 + _read_match_r1],
                seq_read2[:_lt_len_r2],
                same_strand,
            ):
                is_connected = True
                out_lt_len = 0
                out_read_match = _lt_len_r2 + _read_match_r2
                out_rt_len = _rt_len_r2
                mode_r2 = 2
                out_seq = seq_read2
        else:
            is_connected = True
            out_lt_len = 0
            out_read_match = _lt_len_r2 + _read_match_r2
            out_rt_len = _rt_len_r2
            mode_r2 = 2
            out_seq = seq_read2

    if propinquity(_read_match_r1, _rt_len_r2, allowed_difference):
        if SM_align:
            if SM_alignment(
                seq_read1[_lt_len_r1 : _lt_len_r1 + _read_match_r1],
                seq_read2[-_rt_len_r2:],
                same_strand,
            ):
                is_connected = True
                out_lt_len = _lt_len_r2
                out_read_match = _read_match_r2 + _rt_len_r2
                out_rt_len = 0
                mode_r2 = 1
                out_seq = seq_read2
        else:
            is_connected = True
            out_lt_len = _lt_len_r2
            out_read_match = _read_match_r2 + _rt_len_r2
            out_rt_len = 0
            mode_r2 = 1
            out_seq = seq_read2

    if propinquity(_read_match_r2, _lt_len_r1, allowed_difference):
        if SM_align:
            if SM_alignment(
                seq_read2[_lt_len_r2 : _lt_len_r2 + _read_match_r2],
                seq_read1[:_lt_len_r1],
                same_strand,
            ):
                is_connected = True
                out_lt_len = 0
                out_read_match = _lt_len_r1 + _read_match_r1
                out_rt_len = _rt_len_r1
                mode_r1 = 2
                out_seq = seq_read1
        else:
            is_connected = True
            out_lt_len = 0
            out_read_match = _lt_len_r1 + _read_match_r1
            out_rt_len = _rt_len_r1
            mode_r1 = 2
            out_seq = seq_read1

    if propinquity(_read_match_r2, _rt_len_r1, allowed_difference):
        if SM_align:
            if SM_alignment(
                seq_read2[_lt_len_r2 : _lt_len_r2 + _read_match_r2],
                seq_read1[-_rt_len_r1:],
                same_strand,
            ):
                is_connected = True
                out_lt_len = _rt_len_r1
                out_read_match = _read_match_r1 + _rt_len_r1
                out_rt_len = 0
                mode_r1 = 1
                out_seq = seq_read1
        else:
            is_connected = True
            out_lt_len = _rt_len_r1
            out_read_match = _read_match_r1 + _rt_len_r1
            out_rt_len = 0
            mode_r1 = 1
            out_seq = seq_read1

    return (
        is_connected,
        (out_lt_len, out_read_match, out_rt_len),
        out_seq,
        (mode_r1, mode_r2),
    )


def test_is_connected2(
    sms_read1, sms_read2, seq_read1, seq_read2, strand1, strand2, SM_align=False
) -> tuple:
    """Test whether two SMS tuples of chimeric reads can be connected or not.
    :param sms_read1: triple tuple for (left soft-clipped length, middle read matched size, right softclipped length) of read1 OR path
    :type sms_read1: tuple
    :param sms_read2: triple tuple for (left soft-clipped length, middle read matched size, right softclipped length) of read2 OR path
    :type sms_read2: tuple
    :param seq_read1: reads sequence (as it is stored in the BAM file) of read2 OR path
    :type seq_read1: str
    :param seq_read2: reads sequence (as it is stored in the BAM file) of read2 OR path
    :type seq_read2: str
    :param allowed_difference: the difference of read_match_size (Read1) and softclipped length (Read2) to determine the S-M match
    :type allowed_difference: int
    :param SM_align: whether use the alignment of softclipped segment of one read and matched segment of another read
    :type SM_align: bool
    :return: is_connected flag, summed 'SMS' value, summed 'Reads' sequence, mode for read1 and read2
    :rtype: tuple
    ..note:
       ------>    ---------->     ---------->    ---------->
       MMMMMM VS. SSSSSSSSSS   => MMMMMM      OR     MMMMMM
                                  SSSSSSSSSS     SSSSSSSSSS
       <------    <----------     <----------    <----------
        MMMMMM VS. SSSSSSSSSS   =>     MMMMMM OR  MMMMMM
                                   SSSSSSSSSS     SSSSSSSSSS
    """

    def propinquity(a, b, allowed_difference) -> bool:
        if abs(a - b) <= allowed_difference:
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

    def SM_alignment(query_seq, target_seq, same_strand=True) -> bool:
        """matched segment of read1 align to softclipped segment of read2
        solve the problem that the length of the left softclipped segment and the right softclipped segment may be quite similar, eliminate the ambiguity of connected reads
        :param query_seq: matched segment of read1 (M)
        :param target_seq: softclipped segment of read2 (S)
        :param same_strand: read1 and read2 one the same strand or not
        :type query_seq: str
        :type target_seq: str
        :type same_strand: bool
        """
        if not same_strand:
            target_seq = str(Seq(target_seq).reverse_complement())

        allowed_mismatches = abs(len(target_seq) - len(query_seq))
        # print('allowed_mismatches:', allowed_mismatches)
        alignment_result = aligner(query_seq, target_seq, method="glocal")[0]
        _query_seq = alignment_result.seq1.decode("utf-8")
        _target_seq = alignment_result.seq2.decode("utf-8")
        _query_seq_len = len(_query_seq)
        _target_seq_len = len(_target_seq)
        _query_start, _query_end = alignment_result.start1, alignment_result.end1 - 1
        _target_start, _target_end = alignment_result.start2, alignment_result.end2 - 1
        aln_len = _query_end - _query_start + 1
        total_mismatches = len(_query_seq) - aln_len + alignment_result.n_mismatches

        # print('query_seq :', query_seq,  'target_seq :', target_seq)
        # print('total_mismatches :', total_mismatches, _query_seq, _target_seq)

        if total_mismatches <= allowed_mismatches:
            return True
        else:
            return False

    _lt_len_r1, _read_match_r1, _rt_len_r1 = sms_read1
    _lt_len_r2, _read_match_r2, _rt_len_r2 = sms_read2

    same_strand = None
    if seq_read1 == seq_read2:
        same_strand = True
    else:
        same_strand = False

    mode_r1 = init_mode_judge(sms_read1)
    mode_r2 = init_mode_judge(sms_read2)

    is_connected = False
    out_lt_len = 0
    out_rt_len = 0
    out_read_match = 0
    out_seq = ""

    if SM_align:
        if SM_alignment(
            seq_read1[_lt_len_r1 : _lt_len_r1 + _read_match_r1],
            seq_read2[:_lt_len_r2],
            same_strand,
        ):
            is_connected = True
            out_lt_len = 0
            out_read_match = _lt_len_r2 + _read_match_r2
            out_rt_len = _rt_len_r2
            mode_r2 = 2
            out_seq = seq_read2
        else:
            is_connected = True
            out_lt_len = 0
            out_read_match = _lt_len_r2 + _read_match_r2
            out_rt_len = _rt_len_r2
            mode_r2 = 2
            out_seq = seq_read2

    if propinquity(_read_match_r1, _rt_len_r2, allowed_difference):
        if SM_align:
            if SM_alignment(
                seq_read1[_lt_len_r1 : _lt_len_r1 + _read_match_r1],
                seq_read2[-_rt_len_r2:],
                same_strand,
            ):
                is_connected = True
                out_lt_len = _lt_len_r2
                out_read_match = _read_match_r2 + _rt_len_r2
                out_rt_len = 0
                mode_r2 = 1
                out_seq = seq_read2
        else:
            is_connected = True
            out_lt_len = _lt_len_r2
            out_read_match = _read_match_r2 + _rt_len_r2
            out_rt_len = 0
            mode_r2 = 1
            out_seq = seq_read2

    if propinquity(_read_match_r2, _lt_len_r1, allowed_difference):
        if SM_align:
            if SM_alignment(
                seq_read2[_lt_len_r2 : _lt_len_r2 + _read_match_r2],
                seq_read1[:_lt_len_r1],
                same_strand,
            ):
                is_connected = True
                out_lt_len = 0
                out_read_match = _lt_len_r1 + _read_match_r1
                out_rt_len = _rt_len_r1
                mode_r1 = 2
                out_seq = seq_read1
        else:
            is_connected = True
            out_lt_len = 0
            out_read_match = _lt_len_r1 + _read_match_r1
            out_rt_len = _rt_len_r1
            mode_r1 = 2
            out_seq = seq_read1

    if propinquity(_read_match_r2, _rt_len_r1, allowed_difference):
        if SM_align:
            if SM_alignment(
                seq_read2[_lt_len_r2 : _lt_len_r2 + _read_match_r2],
                seq_read1[-_rt_len_r1:],
                same_strand,
            ):
                is_connected = True
                out_lt_len = _rt_len_r1
                out_read_match = _read_match_r1 + _rt_len_r1
                out_rt_len = 0
                mode_r1 = 1
                out_seq = seq_read1
        else:
            is_connected = True
            out_lt_len = _rt_len_r1
            out_read_match = _read_match_r1 + _rt_len_r1
            out_rt_len = 0
            mode_r1 = 1
            out_seq = seq_read1

    return (
        is_connected,
        (out_lt_len, out_read_match, out_rt_len),
        out_seq,
        (mode_r1, mode_r2),
    )


def chimeric_aln_order_finder(
    aln_list, allowed_difference=80, soft_len_cutoff=30
) -> tuple:
    """Find the best connected paths for a list of chimeric alignments

    :param aln_list: list of Read object
    :type aln_list: list
    :param soft_len_cutoff: softclipped segment length cutoff to determine 'two starting reads'
    :type soft_len_cutoff: int
    :param allowed_difference: the difference of read_match_size (Read1) and softclipped length (Read2) to determine the S-M match
    :type allowed_difference: int
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
            if overlap_len < len(b) and a[i] == b[overlap_len]:
                overlap_len += 1
                if count == 0:
                    first_hit_index_a = i
                count += 1
        if len(a) - first_hit_index_a == overlap_len:
            merged_list.extend(a)
            merged_list.extend(b[overlap_len:])
        return merged_list

    def select_path_for_one_node(node, paths) -> tuple:
        """If one read has multiple paths, select one with higest length and lowest number of mismatches
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

    def optimal_path_selector(lt_node, rt_node, lt_paths, rt_paths) -> tuple:
        """
        :param lt_node: start node (left side)
        :param rt_node: start node (right side)
        :param lt_paths: linked paths (list of Path) for lt_node
        :param rt_paths: linked paths (list of Path) for rt_node
        :type lt_node: Read
        :type rt_node: Read
        :type lt_paths: list
        :type rt_paths: list
        :return: end-to-end chains (list of lists) and a dict of reads-pair mode
        :rtype: tuple
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
            _is_connected, _sum_sms, _sum_seq, _mode = test_is_connected(
                start_nodes[0].sms,
                start_nodes[1].sms,
                start_nodes[0].query_sequence,
                start_nodes[1].query_sequence,
                allowed_difference,
                SM_align=False,
            )
            if _is_connected:
                reads_pair_mode_dict[(start_nodes[0], start_nodes[1])] = _mode
                end_to_end_chain.append(start_nodes)
        elif len(candidate_nodes) == 1:
            for _node in start_nodes:
                _is_connected, _sum_sms, _sum_seq, _mode = test_is_connected(
                    _node.sms,
                    candidate_nodes[0].sms,
                    _node.query_sequence,
                    candidate_nodes[0].query_sequence,
                    allowed_difference,
                    SM_align=True,
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
                        _is_connected, _sum_sms, _sum_seq, _mode = test_is_connected(
                            tgt_node.sms,
                            _node.sms,
                            tgt_node.query_sequence,
                            _node.query_sequence,
                            allowed_difference,
                            SM_align=True,
                        )
                        if _is_connected:
                            stop_signal = False
                            # reads_pair_mode_dict[(tgt_node, _node)] = _mode
                            if not tgt_node.linked_paths:
                                tmp_path = Path()
                                tmp_path.add(_node)
                                tmp_path.add_mode({(tgt_node, _node): _mode})
                                tmp_path.sms = _sum_sms
                                tmp_path.sequence = _sum_seq
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
                                    tmp_path.sequence = _sum_seq
                                    tgt_node.add_path(tmp_path)
                    else:
                        # print(tgt_node.linked_paths)
                        for path in tgt_node.linked_paths:
                            # print(path.sms, _node.sms)
                            (
                                _is_connected,
                                _sum_sms,
                                _sum_seq,
                                _mode,
                            ) = test_is_connected(
                                path.sms,
                                _node.sms,
                                path.sequence,
                                _node.query_sequence,
                                allowed_difference,
                                SM_align=True,
                            )
                            if _is_connected:
                                if path.nodes[-1] != _node:
                                    stop_signal = False
                                    path.add_mode({(path.nodes[-1], _node): _mode})
                                    path.add(_node)
                                    path.sms = _sum_sms
                                    path.sequence = _sum_seq
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
                        _is_connected, _sum_sms, _sum_seq, _mode = test_is_connected(
                            tgt_node.sms,
                            _node.sms,
                            tgt_node.query_sequence,
                            _node.query_sequence,
                            allowed_difference,
                            SM_align=True,
                        )
                        if _is_connected:
                            stop_signal = False
                            # reads_pair_mode_dict[(tgt_node, _node)] = _mode
                            if not tgt_node.linked_paths:
                                tmp_path = Path()
                                tmp_path.add(_node)
                                tmp_path.add_mode({(tgt_node, _node): _mode})
                                tmp_path.sms = _sum_sms
                                tmp_path.sequence = _sum_seq
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
                                    tmp_path.sequence = _sum_seq
                                    tgt_node.add_path(tmp_path)

                    else:
                        for path in tgt_node.linked_paths:
                            (
                                _is_connected,
                                _sum_sms,
                                _sum_seq,
                                _mode,
                            ) = test_is_connected(
                                path.sms,
                                _node.sms,
                                path.sequence,
                                _node.query_sequence,
                                allowed_difference,
                                SM_align=True,
                            )
                            if _is_connected:
                                if path.nodes[-1] != _node:
                                    stop_signal = False
                                    path.add_mode({(path.nodes[-1], _node): _mode})
                                    path.add(_node)
                                    path.sms = _sum_sms
                                    path.sequence = _sum_seq

                    count += 1
                if stop_signal:
                    break
            rt_candidate_paths = copy.deepcopy(tgt_node.linked_paths)
            rt_node = tgt_node
            # print(lt_node)
            # print(rt_node)
            # print('lt_candidate_paths: ', lt_candidate_paths)
            # print('rt_candidate_paths: ',rt_candidate_paths)
            # print('lt_candidate_paths: modes: ',lt_candidate_paths[0].mode)
            # print('rt_candidate_paths: modes: ',rt_candidate_paths[0].mode)
            # print("######################")
            end_to_end_chain, reads_pair_mode_dict = optimal_path_selector(
                lt_node, rt_node, lt_candidate_paths, rt_candidate_paths
            )

    return end_to_end_chain, reads_pair_mode_dict


def infer_sv_from_connected_reads(
    read_lt,
    read_rt,
    lt_mode,
    rt_mode,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
) -> tuple:
    """
    :param read_lt: Read 1
    :param read_rt: Read 2
    :param lt_mode: mode of Read 1
    :param rt_mode: mode of Read 2
    :param splice_bin: a small bin for splice site searching
    :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
    :param cvg: annotated splice sites (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param gene_iv: annotated gene region (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param motif_required: considering canonical splice sites only OR considering both canonical and noncanonical splice sites
    :type read_lt: Read
    :type read_rt: Read
    :type lt_mode: int
    :type rt_mode: int
    :type splice_bin : int
    :type genome_fasta: pyfaidx.Fasta
    :type cvg: HTSeq.GenomicArrayOfSets
    :type gene_iv: HTSeq.GenomicArrayOfSets
    :type motif_required: bool
    :return: putative event from reads-pair
    :rtype: tuple
    .. note::
        putative event examples:
            * 'NA', 0, 0, (), (), (), (), []
            * 'TDUP', annotation, canonical/noncanonical, ('chrm1:pos1', 'chrm2:pos2', mode1, mode2),
            (read1_ref_start, read1_ref_end, read1_exons), (read2_ref_start, read2_ref_end, read2_exons), (strand1, strand2), [gene1, gene2]
       annotation explanation:
       3(11) => both breakpoints overlap with coding exons boundary
       2(10) => one breakpoint overlap with coding exons boundary
       1(01) => one breakpoint overlap with coding exons boundary
       0(00) => none breakpoint overlap with coding exons boundary
    """

    def obtain_ins_seq_from_softclipped_part_read(read, mode, indel_size) -> str:
        """
        :param read: a chimeirc read
        :param mode: mode for the chimeirc read
        :param indel_size: indel size infered from 'query_offset - target_offset'
        :type read : Read
        :type mode: int
        :type indel_size: int
        :return: putative insertion sequence from the read
        :rtype: str
        """
        read_seq = read.query_sequence
        ins_seq_in_read = ""
        if mode == 2:  # SM
            ins_seq_in_read = read_seq[: read.lt_soft_len][-indel_size:]
        elif mode == 1:  # MS
            ins_seq_in_read = read_seq[-read.rt_soft_len :][:indel_size]
        return ins_seq_in_read

    if lt_mode == 3 or rt_mode == 3:
        return "NA", 0, 0, (), (), (), (), []

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

    target_start = 0
    target_end = 0
    if lt_chrm == rt_chrm:
        if lt_strand == rt_strand:  # deletion, insertion, duplication
            if lt_mode == 2 and rt_mode == 1:
                target_start = read_rt.ref_start
                target_end = read_lt.ref_end
                target_offset = target_end - target_start
                query_offset = (
                    read_lt.query_length
                    - read_rt.lt_soft_len
                    - read_lt.rt_soft_len
                    + read_lt.indel_size
                    + read_rt.indel_size
                )
                indel_size = query_offset - target_offset
                if indel_size == 0:  # micro-inversion
                    return "NA", 0, 0, (), (), (), (), []
                elif indel_size < 0:  # deletion
                    return "NA", 0, 0, (), (), (), (), []
                elif indel_size >= query_offset:  # large tandem duplication
                    chrm_start = lt_chrm
                    junc_start = read_lt.ref_start
                    chrm_end = lt_chrm
                    junc_end = junc_start + indel_size
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
                    # print(read.reference_start, junc_start)
                    if _nls:
                        if _can == 1:
                            new_junc_start, new_junc_end = update_breakpoints(
                                lt_chrm,
                                junc_start,
                                lt_chrm,
                                junc_end,
                                lt_strand,
                                rt_strand,
                                2,
                                1,
                                splice_bin,
                                genome_fasta,
                            )
                            # print('OLD: ', junc_start, junc_end)
                            # print('NEW: ', new_junc_start, new_junc_end)
                            if not new_junc_start and not new_junc_end:
                                return (
                                    "TDUP",
                                    _anno,
                                    0,
                                    (
                                        f"{lt_chrm}:{junc_start}",
                                        f"{lt_chrm}:{junc_end}",
                                        2,
                                        1,
                                    ),
                                    (lt_start, lt_end, lt_exons),
                                    (rt_start, rt_end, rt_exons),
                                    (lt_strand, rt_strand),
                                    [*_genes],
                                )
                            else:
                                return (
                                    "TDUP",
                                    _anno,
                                    1,
                                    (
                                        f"{lt_chrm}:{new_junc_start}",
                                        f"{lt_chrm}:{new_junc_end}",
                                        2,
                                        1,
                                    ),
                                    (lt_start, lt_end, lt_exons),
                                    (rt_start, rt_end, rt_exons),
                                    (lt_strand, rt_strand),
                                    [*_genes],
                                )
                        else:
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
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                    else:
                        return "NA", 0, 0, (), (), (), (), []
                else:  # read length > tandem duplication size
                    ins_start = read_lt.ref_start
                    ref_allele = genome_fasta[lt_chrm][ins_start : ins_start + 1].seq
                    ins_seq_in_read = obtain_ins_seq_from_softclipped_part_read(
                        read_lt, lt_mode, indel_size
                    )

                    is_DUP = None
                    if short_TDUP_or_not(
                        lt_chrm,
                        lt_mode,
                        rt_start,
                        rt_end,
                        ins_seq_in_read,
                        genome_fasta,
                    ):
                        is_DUP = True
                    else:
                        is_DUP = False
                    if is_DUP:
                        chrm_start = lt_chrm
                        junc_start = read_lt.ref_start
                        chrm_end = lt_chrm
                        junc_end = junc_start + indel_size
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
                            if _can == 1:
                                new_junc_start, new_junc_end = update_breakpoints(
                                    lt_chrm,
                                    junc_start,
                                    lt_chrm,
                                    junc_end,
                                    lt_strand,
                                    rt_strand,
                                    2,
                                    1,
                                    splice_bin,
                                    genome_fasta,
                                )
                                if not new_junc_start and not new_junc_end:
                                    return (
                                        "TDUP",
                                        _anno,
                                        0,
                                        (
                                            f"{lt_chrm}:{junc_start}",
                                            f"{lt_chrm}:{junc_end}",
                                            2,
                                            1,
                                        ),
                                        (lt_start, lt_end, lt_exons),
                                        (rt_start, rt_end, rt_exons),
                                        (lt_strand, rt_strand),
                                        [*_genes],
                                    )
                                else:
                                    return (
                                        "TDUP",
                                        _anno,
                                        1,
                                        (
                                            f"{lt_chrm}:{new_junc_start}",
                                            f"{lt_chrm}:{new_junc_end}",
                                            2,
                                            1,
                                        ),
                                        (lt_start, lt_end, lt_exons),
                                        (rt_start, rt_end, rt_exons),
                                        (lt_strand, rt_strand),
                                        [*_genes],
                                    )
                            else:
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
                                    (lt_strand, rt_strand),
                                    [*_genes],
                                )
                        else:
                            return "NA", 0, 0, (), (), (), (), []
                    else:  # it's a short insertion
                        _genes = gene_annotation(
                            lt_chrm, ins_start, lt_chrm, ins_start, gene_iv
                        )
                        return (
                            "INS",
                            ref_allele,
                            ins_seq_in_read,
                            (ins_start, len(ins_seq_in_read), 2, 1),
                            (lt_start, lt_end, lt_exons),
                            (rt_start, rt_end, rt_exons),
                            (lt_strand, rt_strand),
                            [*_genes],
                        )
            elif lt_mode == 1 and rt_mode == 2:
                target_start = read_lt.ref_start
                target_end = read_rt.ref_start + read_rt.reference_match_size
                target_offset = target_end - target_start
                query_offset = (
                    read_lt.query_length
                    - read_rt.rt_soft_len
                    - read_lt.lt_soft_len
                    + read_lt.indel_size
                    + read_rt.indel_size
                )
                indel_size = query_offset - target_offset
                # print('indel_size: ', indel_size)
                # print('query_offset: ', query_offset)
                if indel_size == 0:  # micro-inversion
                    return "NA", 0, 0, (), (), (), (), []
                elif indel_size < 0:  # deletion
                    return "NA", 0, 0, (), (), (), (), []
                elif indel_size >= query_offset:
                    chrm_start = rt_chrm
                    junc_start = read_rt.ref_start
                    chrm_end = rt_chrm
                    junc_end = junc_start + indel_size
                    # print('junc_start: ', junc_start)
                    # print('junc_end: ', junc_end)
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
                        if _can == 1:
                            new_junc_start, new_junc_end = update_breakpoints(
                                lt_chrm,
                                junc_start,
                                lt_chrm,
                                junc_end,
                                rt_strand,
                                lt_strand,
                                2,
                                1,
                                splice_bin,
                                genome_fasta,
                            )
                            # print('OLD: ', junc_start, junc_end)
                            # print('NEW: ', new_junc_start, new_junc_end)

                            if not new_junc_start and not new_junc_end:
                                return (
                                    "TDUP",
                                    _anno,
                                    0,
                                    (
                                        f"{rt_chrm}:{junc_start}",
                                        f"{rt_chrm}:{junc_end}",
                                        2,
                                        1,
                                    ),
                                    (rt_start, rt_end, rt_exons),
                                    (lt_start, lt_end, lt_exons),
                                    (rt_strand, lt_strand),
                                    [*_genes],
                                )
                            else:
                                return (
                                    "TDUP",
                                    _anno,
                                    1,
                                    (
                                        f"{rt_chrm}:{new_junc_start}",
                                        f"{rt_chrm}:{new_junc_end}",
                                        2,
                                        1,
                                    ),
                                    (rt_start, rt_end, rt_exons),
                                    (lt_start, lt_end, lt_exons),
                                    (rt_strand, lt_strand),
                                    [*_genes],
                                )
                        else:
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
                                (rt_strand, lt_strand),
                                [*_genes],
                            )
                    else:
                        return "NA", 0, 0, (), (), (), (), []
                # indel_size < query_offset
                else:
                    ins_start = read_lt.ref_start + read_lt.reference_match_size
                    ref_allele = genome_fasta[lt_chrm][ins_start : ins_start + 1].seq
                    ins_seq_in_read = obtain_ins_seq_from_softclipped_part_read(
                        read_lt, lt_mode, indel_size
                    )

                    is_DUP = None

                    if short_TDUP_or_not(
                        lt_chrm,
                        lt_mode,
                        rt_start,
                        rt_end,
                        ins_seq_in_read,
                        genome_fasta,
                    ):
                        is_DUP = True
                    else:
                        is_DUP = False
                    if is_DUP:
                        chrm_start = rt_chrm
                        junc_start = rt_start
                        chrm_end = rt_chrm
                        junc_end = junc_start + indel_size
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
                            if _can == 1:
                                new_junc_start, new_junc_end = update_breakpoints(
                                    rt_chrm,
                                    junc_start,
                                    rt_chrm,
                                    junc_end,
                                    rt_strand,
                                    lt_strand,
                                    2,
                                    1,
                                    splice_bin,
                                    genome_fasta,
                                )
                                if not new_junc_start and not new_junc_end:
                                    return (
                                        "TDUP",
                                        _anno,
                                        0,
                                        (
                                            f"{rt_chrm}:{junc_start}",
                                            f"{rt_chrm}:{junc_end}",
                                            2,
                                            1,
                                        ),
                                        (rt_start, rt_end, rt_exons),
                                        (lt_start, lt_end, lt_exons),
                                        (rt_strand, lt_strand),
                                        [*_genes],
                                    )
                                else:
                                    return (
                                        "TDUP",
                                        _anno,
                                        1,
                                        (
                                            f"{rt_chrm}:{new_junc_start}",
                                            f"{rt_chrm}:{new_junc_end}",
                                            2,
                                            1,
                                        ),
                                        (rt_start, rt_end, rt_exons),
                                        (lt_start, lt_end, lt_exons),
                                        (rt_strand, lt_strand),
                                        [*_genes],
                                    )
                            else:
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
                                    (rt_strand, lt_strand),
                                    [*_genes],
                                )
                        else:
                            return "NA", 0, 0, (), (), (), (), []
                    # it is a short insertion
                    else:
                        _genes = gene_annotation(
                            rt_chrm, ins_start, rt_chrm, ins_start, gene_iv
                        )
                        return (
                            "INS",
                            ref_allele,
                            ins_seq_in_read,
                            (ins_start, len(ins_seq_in_read), 1, 2),
                            (rt_start, rt_end, rt_exons),
                            (lt_start, lt_end, lt_exons),
                            (rt_strand, lt_strand),
                            [*_genes],
                        )
            else:
                return "NA", 0, 0, (), (), (), (), []
        else:  # lt_strand != rt_strand
            if lt_mode == rt_mode == 1:
                ra_bp = read_lt.ref_start + read_lt.reference_match_size
                sa_bp = read_rt.ref_start + read_rt.reference_match_size
                if ra_bp == sa_bp:
                    return "NA", 0, 0, (), (), (), (), []
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
                    elif junc_start == sa_bp:
                        strands = (rt_strand, lt_strand)
                        lt_start_end_exons = (rt_start, rt_end, rt_exons)
                        rt_start_end_exons = (lt_start, lt_end, lt_exons)
                    _genes = gene_annotation(
                        chrm_start, junc_start, chrm_end, junc_end, gene_iv
                    )
                    if _nls:
                        # check whether the chimeric alignments uses canonical splice sites or not (60% fraction by default)
                        if read_lt.splice_site_checker(
                            genome_fasta
                        ) and read_rt.splice_site_checker(genome_fasta):
                            if _can == 1:
                                strand_l, strand_r = strands
                                new_junc_start, new_junc_end = update_breakpoints(
                                    lt_chrm,
                                    junc_start,
                                    lt_chrm,
                                    junc_end,
                                    strand_l,
                                    strand_r,
                                    1,
                                    1,
                                    splice_bin,
                                    genome_fasta,
                                )
                                if not new_junc_start and not new_junc_end:
                                    return (
                                        "INV",
                                        _anno,
                                        0,
                                        (
                                            f"{lt_chrm}:{junc_start}",
                                            f"{lt_chrm}:{junc_end}",
                                            1,
                                            1,
                                        ),
                                        lt_start_end_exons,
                                        rt_start_end_exons,
                                        tuple([*strands]),
                                        [*_genes],
                                    )
                                else:
                                    return (
                                        "INV",
                                        _anno,
                                        1,
                                        (
                                            f"{lt_chrm}:{new_junc_start}",
                                            f"{lt_chrm}:{new_junc_end}",
                                            1,
                                            1,
                                        ),
                                        lt_start_end_exons,
                                        rt_start_end_exons,
                                        tuple([*strands]),
                                        [*_genes],
                                    )
                            else:
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
                                    tuple([*strands]),
                                    [*_genes],
                                )
                        else:
                            return "NA", 0, 0, (), (), (), (), []
                    else:
                        return "NA", 0, 0, (), (), (), (), []
            elif lt_mode == rt_mode == 2:  # inversion
                ra_bp = read_lt.ref_start
                sa_bp = read_rt.ref_start
                if ra_bp == sa_bp:
                    return "NA", 0, 0, (), (), (), (), []
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
                    elif junc_start == sa_bp:
                        strands = (rt_strand, lt_strand)
                        lt_start_end_exons = (rt_start, rt_end, rt_exons)
                        rt_start_end_exons = (lt_start, lt_end, lt_exons)
                    _genes = gene_annotation(
                        chrm_start, junc_start, chrm_end, junc_end, gene_iv
                    )
                    if _nls:
                        # check whether the chimeric alignments uses canonical splice sites or not (60% fraction by default)
                        if read_lt.splice_site_checker(
                            genome_fasta
                        ) and read_rt.splice_site_checker(genome_fasta):
                            if _can == 1:
                                strand_l, strand_r = strands
                                new_junc_start, new_junc_end = update_breakpoints(
                                    lt_chrm,
                                    junc_start,
                                    lt_chrm,
                                    junc_end,
                                    strand_l,
                                    strand_r,
                                    2,
                                    2,
                                    splice_bin,
                                    genome_fasta,
                                )
                                if not new_junc_start and not new_junc_end:
                                    return (
                                        "INV",
                                        _anno,
                                        0,
                                        (
                                            f"{lt_chrm}:{junc_start}",
                                            f"{lt_chrm}:{junc_end}",
                                            2,
                                            2,
                                        ),
                                        lt_start_end_exons,
                                        rt_start_end_exons,
                                        tuple([*strands]),
                                        [*_genes],
                                    )
                                else:
                                    return (
                                        "INV",
                                        _anno,
                                        1,
                                        (
                                            f"{lt_chrm}:{new_junc_start}",
                                            f"{lt_chrm}:{new_junc_end}",
                                            2,
                                            2,
                                        ),
                                        lt_start_end_exons,
                                        rt_start_end_exons,
                                        tuple([*strands]),
                                        [*_genes],
                                    )
                            else:
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
                                    tuple([*strands]),
                                    [*_genes],
                                )
                        else:
                            return "NA", 0, 0, (), (), (), (), []
                    else:
                        return "NA", 0, 0, (), (), (), (), []
            else:
                return "NA", 0, 0, (), (), (), (), []
    else:  # lt_chrm != rt_chrm
        if lt_strand == rt_strand:
            if lt_mode == 1 and rt_mode == 2:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start + read_lt.reference_match_size
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start
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
                    if _can == 1:
                        new_junc_start, new_junc_end = update_breakpoints(
                            chrm_start,
                            junc_start,
                            chrm_end,
                            junc_end,
                            lt_strand,
                            rt_strand,
                            1,
                            2,
                            splice_bin,
                            genome_fasta,
                        )
                        if not new_junc_start and not new_junc_end:
                            return (
                                "TRA",
                                _anno,
                                0,
                                (
                                    f"{lt_chrm}:{junc_start}",
                                    f"{rt_chrm}:{junc_end}",
                                    1,
                                    2,
                                ),
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                        else:
                            return (
                                "TRA",
                                _anno,
                                1,
                                (
                                    f"{lt_chrm}:{new_junc_start}",
                                    f"{rt_chrm}:{new_junc_end}",
                                    1,
                                    2,
                                ),
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                    else:
                        return (
                            "TRA",
                            _anno,
                            _can,
                            (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 1, 2),
                            (lt_start, lt_end, lt_exons),
                            (rt_start, rt_end, rt_exons),
                            (lt_strand, rt_strand),
                            [*_genes],
                        )
                else:
                    return "NA", 0, 0, (), (), (), (), []
            elif lt_mode == 2 and rt_mode == 1:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start + read_rt.reference_match_size
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
                    if _can == 1:
                        new_junc_start, new_junc_end = update_breakpoints(
                            chrm_start,
                            junc_start,
                            chrm_end,
                            junc_end,
                            lt_strand,
                            rt_strand,
                            2,
                            1,
                            splice_bin,
                            genome_fasta,
                        )
                        if not new_junc_start and not new_junc_end:
                            return (
                                "TRA",
                                _anno,
                                0,
                                (
                                    f"{lt_chrm}:{junc_start}",
                                    f"{rt_chrm}:{junc_end}",
                                    2,
                                    1,
                                ),
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                        else:
                            return (
                                "TRA",
                                _anno,
                                1,
                                (
                                    f"{lt_chrm}:{new_junc_start}",
                                    f"{rt_chrm}:{new_junc_end}",
                                    2,
                                    1,
                                ),
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                    else:
                        return (
                            "TRA",
                            _anno,
                            _can,
                            (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 2, 1),
                            (lt_start, lt_end, lt_exons),
                            (rt_start, rt_end, rt_exons),
                            (lt_strand, rt_strand),
                            [*_genes],
                        )
                else:
                    return "NA", 0, 0, (), (), (), (), []
            else:
                return "NA", 0, 0, (), (), (), (), []
        else:  # lt_strand != rt_strand
            if lt_mode == rt_mode == 1:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start + read_lt.reference_match_size
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start + read_rt.reference_match_size
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
                    if _can == 1:
                        new_junc_start, new_junc_end = update_breakpoints(
                            chrm_start,
                            junc_start,
                            chrm_end,
                            junc_end,
                            lt_strand,
                            rt_strand,
                            1,
                            1,
                            splice_bin,
                            genome_fasta,
                        )
                        if not new_junc_start and not new_junc_end:
                            return (
                                "TRA",
                                _anno,
                                0,
                                (
                                    f"{lt_chrm}:{junc_start}",
                                    f"{rt_chrm}:{junc_end}",
                                    1,
                                    1,
                                ),
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                        else:
                            return (
                                "TRA",
                                _anno,
                                1,
                                (
                                    f"{lt_chrm}:{new_junc_start}",
                                    f"{rt_chrm}:{new_junc_end}",
                                    1,
                                    1,
                                ),
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                    else:
                        return (
                            "TRA",
                            _anno,
                            _can,
                            (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 1, 1),
                            (lt_start, lt_end, lt_exons),
                            (rt_start, rt_end, rt_exons),
                            (lt_strand, rt_strand),
                            [*_genes],
                        )
                else:
                    return "NA", 0, 0, (), (), (), (), []
            elif lt_mode == rt_mode == 2:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start
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
                    if _can == 1:
                        new_junc_start, new_junc_end = update_breakpoints(
                            chrm_start,
                            junc_start,
                            chrm_end,
                            junc_end,
                            lt_strand,
                            rt_strand,
                            2,
                            2,
                            splice_bin,
                            genome_fasta,
                        )
                        if not new_junc_start and not new_junc_end:
                            return (
                                "TRA",
                                _anno,
                                0,
                                (
                                    f"{lt_chrm}:{junc_start}",
                                    f"{rt_chrm}:{junc_end}",
                                    2,
                                    2,
                                ),
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                        else:
                            return (
                                "TRA",
                                _anno,
                                1,
                                (
                                    f"{lt_chrm}:{new_junc_start}",
                                    f"{rt_chrm}:{new_junc_end}",
                                    2,
                                    2,
                                ),
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                    else:
                        return (
                            "TRA",
                            _anno,
                            _can,
                            (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 2, 2),
                            (lt_start, lt_end, lt_exons),
                            (rt_start, rt_end, rt_exons),
                            (lt_strand, rt_strand),
                            [*_genes],
                        )
                else:
                    return "NA", 0, 0, (), (), (), (), []
            else:
                return "NA", 0, 0, (), (), (), (), []


def nls_series_assembly(nls_series_list, overlap_len_cutoff):
    """assemble overlapped NLS series
    A->-B->-C and C->-D => A->-B->-C->-D
    """

    def is_connected_series(series_A, series_B, overlap_len_cutoff):
        """find two series with one overlapping node and return the merged series
        series_A: A->-B; series_B: B->-C
        merged_series: A->-B->-C
        """
        merged_series = copy.deepcopy(series_A)
        last_of_A = series_A[-1]
        first_of_B = series_B[0]
        if (
            last_of_A.strand == first_of_B.strand
            and last_of_A.ref_start <= first_of_B.ref_start
            and last_of_A.ref_end <= first_of_B.ref_end
            and last_of_A.ref_end - first_of_B.ref_start > overlap_len_cutoff
            and last_of_A.introns == first_of_B.introns
        ):
            merged_series[-1].next_breakpoint = first_of_B.next_breakpoint
            merged_series[-1].exons[0][0] = last_of_A.ref_start
            merged_series[-1].exons[-1][1] = first_of_B.ref_end
            for i in series_B[1:]:
                merged_series.append(i)
            return merged_series
        else:
            return None


def merge_nls_forms(nls_src_forms):
    """
    :param nls_src_forms: NLS candidate forms (transcript forms)
    :type nls_src_forms: list of Series
    :return: a dictionary of supporting reads, a list of unique Series, a dictionary of breakpoint pair
    :rtype: tuple
    .. note::
        *
    """
    breakpoint_pair_dict = defaultdict(list)
    sr_dict = defaultdict(int)
    sorted_nls_src_forms = sorted(nls_src_forms, reverse=False)
    nls_forms = []
    for series in sorted_nls_src_forms:
        pairs = series.decompose()
        for pair in pairs:
            if not breakpoint_pair_dict:
                pass
            # if not pair in breakpoint_pair_dict:
            #    pass

    return None


def output_bedpe_file(sr_dict, group_dict, prefix, splice_bin):
    """
    :param sr_dict: sv candidate to number of supporting reads(SR) dictionary
    :param group_dict: sv candidate to group of events dictionary, connected chimeric reads are included in one group
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
