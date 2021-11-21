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

from .. import __version__
from ..classes import Path
from ..utils import get_softclip_length
from .helper import vcf_header

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


def scan_region(bam_object, in_region, mapq_cutoff, soft_len_cutoff):
    """Scan candidate region to find soft-clipped segment that can aligned to the anchor soft-clipped segment in reads with SV tag.

    :param bam_object: pysam parsed SV-tag added BAM
    :param in_region: one candidate region from BEDPE file
    :param mapq_cutoff: MAPQ cutoff
    :param soft_len_cutoff: minimum length of soft-clipped segments considered to be align to anchor soft-clipped segment in reads with SV tag
    :type bam_object: pysam.libcalignmentfile.AlignmentFile object
    :type in_region: str (chrm:start-end)
    :type mapq_cutoff: int
    :type soft_len_cutoff: int
    :return: dictionary of SV => SR
    :rtype: dict

    ..note ::

    """
    for col in bam_object.pileup(
        region=in_region, truncate=True, stepper="nofilter", min_base_quality=0
    ):
        chrm = col.reference_name
        ao_dict = {}
        sv_seq_dict = defaultdict(set)
        sr_list = {1: [], 2: []}
        sv_soft_mode = 0
        for read in col.pileups:
            # read is an instance of pysam.PileupRead
            aln_read = read.alignment
            query_position = read.query_position
            if aln_read.mapq >= mapq_cutoff and query_position != None:
                # the read has soft-clipped part but no SV tag
                if not aln_read.has_tag("SV"):
                    if "S" in aln_read.cigarstring:
                        # print(read.alignment.query_name)
                        soft_len, soft_seq, soft_pos, soft_mode = get_softclip_length(
                            aln_read
                        )
                        # print(soft_mode, read.alignment.query_name)
                        # only consider soft_seq without undetermined (N) bases
                        if "N" not in soft_seq:
                            if soft_len >= soft_len_cutoff:
                                if soft_mode == 1:
                                    if (
                                        abs(
                                            len(
                                                aln_read.query_sequence[
                                                    query_position + 1 :
                                                ]
                                            )
                                            - soft_len
                                        )
                                        < 20
                                    ):
                                        sr_list[1].append(
                                            aln_read.query_sequence[
                                                query_position + 1 :
                                            ]
                                        )
                                elif soft_mode == 2:
                                    if (
                                        abs(
                                            len(
                                                aln_read.query_sequence[:query_position]
                                            )
                                            - soft_len
                                        )
                                        < 20
                                    ):
                                        sr_list[2].append(
                                            aln_read.query_sequence[:query_position]
                                        )
                # the read has SV tag, it must be a representive alignment
                else:
                    if sv_checker(read, col.reference_pos, mapq_cutoff, sv_len_cutoff):
                        _, _, _, sv_soft_mode = get_softclip_length(read.alignment)
                        (
                            sv_type,
                            _anno_can,
                            _position,
                            size_or_sup_position,
                            rep_aln_mode,
                            sup_aln_mode,
                            strands,
                            genes,
                        ) = read.alignment.get_tag("SV")[:-1].split(",")
                        # store gene information
                        _gene1, _gene2 = genes.split("|")
                        # bp1 and bp2 modes determination
                        if ":" in size_or_sup_position:
                            _chrm2, _bp2_pos = size_or_sup_position.split(":")
                        else:
                            _bp2_pos = int(size_or_sup_position) + int(_position)
                            _chrm2 = chrm
                        if "_" in _chrm2 or "_" in chrm:
                            continue
                        # breakpoint with small chrm as the bp1, breakpoint with small position as the bp1
                        if is_bps_order_changed(
                            f"{chrm}:{_position}", f"{_chrm2}:{_bp2_pos}"
                        ):
                            if _anno_can[0] == "1":
                                _anno_can = f"2{_anno_can[1]}"
                            elif _anno_can[0] == "2":
                                _anno_can = f"1{_anno_can[1]}"
                                sv_id = "{},{},{},{},{}:{}".format(
                                    sv_type,
                                    _anno_can,
                                    _chrm2,
                                    _bp2_pos,
                                    chrm,
                                    _position,
                                )
                                strand_dict[sv_id] = (strands[1], strands[0])
                                gene_dict[sv_id] = (_gene2, _gene1)
                            if abs(int(_position) - col.reference_pos) < abs(
                                int(_bp2_pos) - col.reference_pos
                            ):
                                mode_dict[sv_id] = (
                                    int(sup_aln_mode),
                                    int(rep_aln_mode),
                                )
                            else:
                                mode_dict[sv_id] = (
                                    int(rep_aln_mode),
                                    int(sup_aln_mode),
                                )
                        else:
                            # sv_id = '{},{},{},{},{}'.format(sv_type, _anno_can, chrm, _position, size_or_sup_position)
                            sv_id = "{},{},{},{},{}:{}".format(
                                sv_type, _anno_can, chrm, _position, _chrm2, _bp2_pos
                            )
                            strand_dict[sv_id] = (strands[0], strands[1])
                            gene_dict[sv_id] = (_gene1, _gene2)
                            if abs(int(_position) - col.reference_pos) < abs(
                                int(_bp2_pos) - col.reference_pos
                            ):
                                mode_dict[sv_id] = (
                                    int(rep_aln_mode),
                                    int(sup_aln_mode),
                                )
                            else:
                                mode_dict[sv_id] = (
                                    int(sup_aln_mode),
                                    int(rep_aln_mode),
                                )

                                # MS mode (representive alignment)
                            if int(rep_aln_mode) == 1:
                                sv_seq_dict[sv_id].add(
                                    read.alignment.query_sequence[
                                        read.query_position + 1 :
                                    ]
                                )
                            # SM mode (representive alignment)
                            else:
                                sv_seq_dict[sv_id].add(
                                    read.alignment.query_sequence[: read.query_position]
                                )

                        if not sv_id in ao_dict:
                            ao_dict[sv_id] = 1
                        else:
                            ao_dict[sv_id] += 1
                    for (
                        sv_id
                    ) in ao_dict:  # only consider one event per genomic position
                        if sv_id in sv_set:
                            continue
                        ao = ao_dict[sv_id]
                        if ao < seed_ao_cutoff:
                            continue
                        # print(sr_list)
                        if sv_soft_mode in sr_list:
                            for sr_read in sr_list[sv_soft_mode]:
                                # print(sr_read, sv_seq_dict[sv_id])
                                if (
                                    mismatch_count(
                                        sr_read,
                                        sv_seq_dict[sv_id],
                                        alignment_frac,
                                        sv_soft_mode,
                                    )
                                    <= mismatch_cutoff
                                ):
                                    ao += 1
                        # print('Left AO:', ao)
                        ao_total_dict[sv_id] = ao


def joint_call(
    input_bam,
    target,
    output_prefix,
    sr_cutoff,
    dp_cutoff,
    pso_cutoff,
    sv_len_cutoff,
    soft_len_cutoff,
    mapq_cutoff,
    mismatch_cutoff,
    alignment_frac,
    seed_cutoff,
):
    """joint calling NLS events using BAM file with SV tags and BEDPE file
    :param input_bam: BAM file with SV tags !important
    :param target: BEDPE file !important
    :param output_prefix: description
    :param sr_cutoff: supporting reads (SR) cutoff
    :param dp_cutoff: depth (DP) cutoff
    :param pso_cutoff: PSO (cutoff)
    :param sv_len_cutoff: description
    :param soft_len_cutoff: description
    :param mapq_cutoff: description
    :param mismatch_cutoff: description
    :param alignment_frac: description
    :param seed_cutoff: description
    :type output_prefix: str
    :type sr_cutoff: int
    :type dp_cutoff: int
    :type pso_cutoff: float
    :type sv_len_cutoff: int
    :type soft_len_cutoff: int
    :type mapq_cutoff: int
    :type mismatch_cutoff: int
    :type alignment_frac: float
    :type seed_cutoff: int
    :return: no returned value
    :rtype: None
    """
    # store the called NLS events
    sv_set = set()
    bam_file = pysam.AlignmentFile(input_bam, "rb")
    vcf_file = open(f"{output_prefix}.sv.vcf", "w")
    bam_header = bam_file.header
    vcf_file.write(f"{vcf_header(output_prefix, bam_header)}\n")

    regions = OrderedDict()
    try:
        with open(target) as ifp:
            for line in ifp:
                (
                    chrm1,
                    s1,
                    e1,
                    chrm2,
                    s2,
                    e2,
                    group_id,
                    sv_ao,
                    strand1,
                    strand2,
                ) = line.rstrip("\n").split("\t")
                if not "_" in chrm1 and not "_" in chrm2:
                    regions[(f"{chrm1}:{s1}-{e1}", f"{chrm2}:{s2}-{e2}")] = {
                        "ao": int(sv_ao),
                        "group": group_id,
                    }
    except TypeError as e:
        sys.stderr.write("Error! BEDPE file {}\n".format(e))

    # output_NLS_events = OrderedDict()
    ID_count = 1
    INS_count = 1
    for region_l, region_r in regions:
        # Insertion (novel sequence)
        if region_l == region_r:
            ins_ao_dict = {}
            ins_strand_dict = {}
            ins_gene_dict = {}
            try:
                for col in bam_file.pileup(
                    region=region_l,
                    truncate=True,
                    stepper="nofilter",
                    min_base_quality=0,
                ):
                    # print(col.reference_name, col.reference_pos)
                    chrm = col.reference_name
                    dp = col.nsegments
                    for read in col.pileups:
                        # read is an instance of pysam.PileupRead
                        if (
                            read.alignment.mapq >= mapq_cutoff
                            and read.query_position != None
                        ):
                            if read.alignment.has_tag("OT"):
                                (
                                    sv_type,
                                    _anno_can,
                                    _position,
                                    size_or_sup_position,
                                    rep_aln_mode,
                                    sup_aln_mode,
                                    strands,
                                    genes,
                                ) = read.alignment.get_tag("OT")[:-1].split(",")
                                if sv_type == "INS":
                                    ref_allele, alt_allele = _anno_can.split("|")
                                    sv_id = "{},{},{},{}".format(
                                        sv_type, _position, ref_allele, alt_allele
                                    )

                                    # store strands information: rep_strand, sup_strand
                                    ins_strand_dict[sv_id] = (strands[0], strands[1])
                                    # store gene information
                                    ins_gene_dict[sv_id] = genes.split("|")

                                    if not sv_id in ins_ao_dict:
                                        ins_ao_dict[sv_id] = 1
                                    else:
                                        ins_ao_dict[sv_id] += 1

                    for (
                        sv_id
                    ) in ins_ao_dict:  # only consider one event per genomic position
                        if sv_id in sv_set:
                            continue
                        ao = ins_ao_dict[sv_id]
                        if ao < seed_ao_cutoff:
                            continue
                        vaf = ao / dp
                        if dp >= dp_cutoff and ao >= sr_cutoff and vaf >= pso_cutoff:
                            sv_set.add(sv_id)
                            sv_type, _position, ref_allele, alt_allele = sv_id.split(
                                ","
                            )
                            bp_distance = len(alt_allele)
                            bp1_strand, bp2_strand = ins_strand_dict[sv_id]
                            bp1_gene, bp2_gene = ins_gene_dict[sv_id]

                            can_field = "NA"
                            anno_field = "NA"
                            mode1 = "NA"
                            mode2 = "NA"

                            vcf_field_gt = "GT\t0/1"
                            vcf_field_info = "{};AO={};DP={};VAF={:.2g};SVLEN={};SVTYPE={};SVMETHOD=ScanNLS;GENE1={};GENE2={};STRAND1={};STRAND2={};END={}".format(
                                can_field,
                                ao,
                                dp,
                                vaf,
                                bp_distance,
                                sv_type,
                                bp1_gene,
                                bp2_gene,
                                bp1_strand,
                                bp2_strand,
                                _position,
                            )
                            vcf_file.write(
                                "{}\t{}\tINS_{}\t{}\t{}\t.\t.\t{}\t{}\n".format(
                                    chrm,
                                    _position,
                                    INS_count,
                                    ref_allele,
                                    alt_allele,
                                    vcf_field_info,
                                    vcf_field_gt,
                                )
                            )
                            INS_count += 1

            except ValueError as e:
                sys.stderr.write("Error in INS! {} at {}\n".format(e, region_l))
                continue
        else:  # DUP/INV/TRA
            sv_ao = regions[(region_l, region_r)]["ao"]
            group_id = regions[(region_l, region_r)]["group"]

            print("DUP/INV/TRA", region_l, region_r)

            ao_total_dict = {}
            strand_dict = {}
            gene_dict = {}
            mode_dict = {}
            # print('Left region:', region_l)
            ################left breakpoint##############################
            try:
                for col in bam_file.pileup(
                    region=region_l,
                    truncate=True,
                    stepper="nofilter",
                    min_base_quality=0,
                ):
                    # print(col.reference_name, col.reference_pos)
                    chrm = col.reference_name
                    ao_dict = {}
                    sv_seq_dict = defaultdict(set)
                    sr_list = {1: [], 2: []}
                    sv_soft_mode = 0
                    for read in col.pileups:
                        # read is an instance of pysam.PileupRead
                        if (
                            read.alignment.mapq >= mapq_cutoff
                            and read.query_position != None
                        ):
                            # the read has soft-clipped part but no SV tag
                            # print(read.alignment.query_name)
                            if not read.alignment.has_tag("SV"):
                                if "S" in read.alignment.cigarstring:
                                    # print(read.alignment.query_name)
                                    (
                                        soft_len,
                                        soft_seq,
                                        soft_pos,
                                        soft_mode,
                                    ) = get_softclip_length(read.alignment)
                                    # print(soft_mode, read.alignment.query_name)
                                    if "N" not in soft_seq:
                                        if soft_len >= soft_len_cutoff:
                                            if soft_mode == 1:
                                                if (
                                                    abs(
                                                        len(
                                                            read.alignment.query_sequence[
                                                                read.query_position
                                                                + 1 :
                                                            ]
                                                        )
                                                        - soft_len
                                                    )
                                                    < 20
                                                ):
                                                    sr_list[1].append(
                                                        read.alignment.query_sequence[
                                                            read.query_position + 1 :
                                                        ]
                                                    )
                                            elif soft_mode == 2:
                                                if (
                                                    abs(
                                                        len(
                                                            read.alignment.query_sequence[
                                                                : read.query_position
                                                            ]
                                                        )
                                                        - soft_len
                                                    )
                                                    < 20
                                                ):
                                                    sr_list[2].append(
                                                        read.alignment.query_sequence[
                                                            : read.query_position
                                                        ]
                                                    )
                            else:
                                if sv_checker(
                                    read, col.reference_pos, mapq_cutoff, sv_len_cutoff
                                ):
                                    _, _, _, sv_soft_mode = get_softclip_length(
                                        read.alignment
                                    )
                                    (
                                        sv_type,
                                        _anno_can,
                                        _position,
                                        size_or_sup_position,
                                        rep_aln_mode,
                                        sup_aln_mode,
                                        strands,
                                        genes,
                                    ) = read.alignment.get_tag("SV")[:-1].split(",")
                                    # store gene information
                                    _gene1, _gene2 = genes.split("|")
                                    # bp1 and bp2 modes determination
                                    if ":" in size_or_sup_position:
                                        _chrm2, _bp2_pos = size_or_sup_position.split(
                                            ":"
                                        )
                                    else:
                                        _bp2_pos = int(size_or_sup_position) + int(
                                            _position
                                        )
                                        _chrm2 = chrm
                                    if "_" in _chrm2 or "_" in chrm:
                                        continue
                                    # breakpoint with small chrm as the bp1, breakpoint with small position as the bp1
                                    if is_bps_order_changed(
                                        f"{chrm}:{_position}", f"{_chrm2}:{_bp2_pos}"
                                    ):
                                        if _anno_can[0] == "1":
                                            _anno_can = f"2{_anno_can[1]}"
                                        elif _anno_can[0] == "2":
                                            _anno_can = f"1{_anno_can[1]}"
                                        sv_id = "{},{},{},{},{}:{}".format(
                                            sv_type,
                                            _anno_can,
                                            _chrm2,
                                            _bp2_pos,
                                            chrm,
                                            _position,
                                        )
                                        strand_dict[sv_id] = (strands[1], strands[0])
                                        gene_dict[sv_id] = (_gene2, _gene1)
                                        if abs(
                                            int(_position) - col.reference_pos
                                        ) < abs(int(_bp2_pos) - col.reference_pos):
                                            mode_dict[sv_id] = (
                                                int(sup_aln_mode),
                                                int(rep_aln_mode),
                                            )
                                        else:
                                            mode_dict[sv_id] = (
                                                int(rep_aln_mode),
                                                int(sup_aln_mode),
                                            )
                                    else:
                                        # sv_id = '{},{},{},{},{}'.format(sv_type, _anno_can, chrm, _position, size_or_sup_position)
                                        sv_id = "{},{},{},{},{}:{}".format(
                                            sv_type,
                                            _anno_can,
                                            chrm,
                                            _position,
                                            _chrm2,
                                            _bp2_pos,
                                        )
                                        strand_dict[sv_id] = (strands[0], strands[1])
                                        gene_dict[sv_id] = (_gene1, _gene2)
                                        if abs(
                                            int(_position) - col.reference_pos
                                        ) < abs(int(_bp2_pos) - col.reference_pos):
                                            mode_dict[sv_id] = (
                                                int(rep_aln_mode),
                                                int(sup_aln_mode),
                                            )
                                        else:
                                            mode_dict[sv_id] = (
                                                int(sup_aln_mode),
                                                int(rep_aln_mode),
                                            )

                                    # MS mode (representive alignment)
                                    if int(rep_aln_mode) == 1:
                                        sv_seq_dict[sv_id].add(
                                            read.alignment.query_sequence[
                                                read.query_position + 1 :
                                            ]
                                        )
                                    # SM mode (representive alignment)
                                    else:
                                        sv_seq_dict[sv_id].add(
                                            read.alignment.query_sequence[
                                                : read.query_position
                                            ]
                                        )

                                    if not sv_id in ao_dict:
                                        ao_dict[sv_id] = 1
                                    else:
                                        ao_dict[sv_id] += 1
                    for (
                        sv_id
                    ) in ao_dict:  # only consider one event per genomic position
                        if sv_id in sv_set:
                            continue
                        ao = ao_dict[sv_id]
                        if ao < seed_ao_cutoff:
                            continue
                        # print(sr_list)
                        if sv_soft_mode in sr_list:
                            for sr_read in sr_list[sv_soft_mode]:
                                # print(sr_read, sv_seq_dict[sv_id])
                                if (
                                    mismatch_count(
                                        sr_read,
                                        sv_seq_dict[sv_id],
                                        alignment_frac,
                                        sv_soft_mode,
                                    )
                                    <= mismatch_cutoff
                                ):
                                    ao += 1
                        # print('Left AO:', ao)
                        ao_total_dict[sv_id] = ao
            except ValueError as e:
                sys.stderr.write("Error! Left:{} at {}\n".format(e, region_l))
                continue
            ##########right breakpoint######################################
            try:
                for col in bam_file.pileup(
                    region=region_r,
                    truncate=True,
                    stepper="nofilter",
                    min_base_quality=0,
                ):
                    dp_r = col.nsegments
                    chrm = col.reference_name
                    ao_dict = {}
                    sv_seq_dict = defaultdict(set)
                    sr_list = {1: [], 2: []}
                    sv_soft_mode = 3
                    for read in col.pileups:
                        # read is an instance of pysam.PileupRead
                        if (
                            read.alignment.mapq >= mapq_cutoff
                            and read.query_position != None
                        ):
                            # the read has soft-clipped part but no SV tag
                            if not read.alignment.has_tag("SV"):
                                if "S" in read.alignment.cigarstring:
                                    (
                                        soft_len,
                                        soft_seq,
                                        soft_pos,
                                        soft_mode,
                                    ) = get_softclip_length(read.alignment)
                                    if "N" not in soft_seq:
                                        if soft_len >= soft_len_cutoff:
                                            if soft_mode == 1:
                                                if (
                                                    abs(
                                                        len(
                                                            read.alignment.query_sequence[
                                                                read.query_position
                                                                + 1 :
                                                            ]
                                                        )
                                                        - soft_len
                                                    )
                                                    < 20
                                                ):
                                                    sr_list[1].append(
                                                        read.alignment.query_sequence[
                                                            read.query_position + 1 :
                                                        ]
                                                    )
                                            elif soft_mode == 2:
                                                if (
                                                    abs(
                                                        len(
                                                            read.alignment.query_sequence[
                                                                : read.query_position
                                                            ]
                                                        )
                                                        - soft_len
                                                    )
                                                    < 20
                                                ):
                                                    sr_list[2].append(
                                                        read.alignment.query_sequence[
                                                            : read.query_position
                                                        ]
                                                    )
                            else:
                                if sv_checker(
                                    read, col.reference_pos, mapq_cutoff, sv_len_cutoff
                                ):
                                    _, _, _, sv_soft_mode = get_softclip_length(
                                        read.alignment
                                    )
                                    (
                                        sv_type,
                                        _anno_can,
                                        _position,
                                        size_or_sup_position,
                                        rep_aln_mode,
                                        sup_aln_mode,
                                        strands,
                                        genes,
                                    ) = read.alignment.get_tag("SV")[:-1].split(",")
                                    # store gene information
                                    _gene1, _gene2 = genes.split("|")
                                    # bp1 and bp2 modes determination
                                    if ":" in size_or_sup_position:
                                        _chrm2, _bp2_pos = size_or_sup_position.split(
                                            ":"
                                        )
                                    else:
                                        _bp2_pos = int(size_or_sup_position) + int(
                                            _position
                                        )
                                        _chrm2 = chrm
                                    if "_" in _chrm2 or "_" in chrm:
                                        continue
                                    if is_bps_order_changed(
                                        f"{chrm}:{_position}", f"{_chrm2}:{_bp2_pos}"
                                    ):
                                        if _anno_can[0] == "1":
                                            _anno_can = f"2{_anno_can[1]}"
                                        elif _anno_can[0] == "2":
                                            _anno_can = f"1{_anno_can[1]}"
                                        sv_id = "{},{},{},{},{}:{}".format(
                                            sv_type,
                                            _anno_can,
                                            _chrm2,
                                            _bp2_pos,
                                            chrm,
                                            _position,
                                        )
                                        strand_dict[sv_id] = (strands[1], strands[0])
                                        gene_dict[sv_id] = (_gene2, _gene1)
                                        if abs(
                                            int(_position) - col.reference_pos
                                        ) < abs(int(_bp2_pos) - col.reference_pos):
                                            mode_dict[sv_id] = (
                                                int(sup_aln_mode),
                                                int(rep_aln_mode),
                                            )
                                        else:
                                            mode_dict[sv_id] = (
                                                int(rep_aln_mode),
                                                int(sup_aln_mode),
                                            )
                                    else:
                                        sv_id = "{},{},{},{},{}:{}".format(
                                            sv_type,
                                            _anno_can,
                                            chrm,
                                            _position,
                                            _chrm2,
                                            _bp2_pos,
                                        )
                                        strand_dict[sv_id] = (strands[0], strands[1])
                                        gene_dict[sv_id] = (_gene1, _gene2)
                                        if abs(
                                            int(_position) - col.reference_pos
                                        ) < abs(int(_bp2_pos) - col.reference_pos):
                                            mode_dict[sv_id] = (
                                                int(rep_aln_mode),
                                                int(sup_aln_mode),
                                            )
                                        else:
                                            mode_dict[sv_id] = (
                                                int(sup_aln_mode),
                                                int(rep_aln_mode),
                                            )

                                    # MS mode (representive alignment)
                                    if int(rep_aln_mode) == 1:
                                        sv_seq_dict[sv_id].add(
                                            read.alignment.query_sequence[
                                                read.query_position + 1 :
                                            ]
                                        )
                                    # SM mode (representive alignment)
                                    else:
                                        sv_seq_dict[sv_id].add(
                                            read.alignment.query_sequence[
                                                : read.query_position
                                            ]
                                        )

                                    if not sv_id in ao_dict:
                                        ao_dict[sv_id] = 1
                                    else:
                                        ao_dict[sv_id] += 1
                    for (
                        sv_id
                    ) in ao_dict:  # only consider one event per genomic position
                        if sv_id in sv_set:
                            continue
                        ao = ao_dict[sv_id]
                        if ao < seed_ao_cutoff:
                            continue
                        if sv_soft_mode in sr_list:
                            for sr_read in sr_list[sv_soft_mode]:
                                if (
                                    mismatch_count(
                                        sr_read,
                                        sv_seq_dict[sv_id],
                                        alignment_frac,
                                        sv_soft_mode,
                                    )
                                    <= mismatch_cutoff
                                ):
                                    ao += 1
                        # print('Right AO:', ao)
                        if sv_id in ao_total_dict:
                            ao_total_dict[sv_id] = ao_total_dict[sv_id] + ao
                            # print(ao_total_dict[sv_id])
                        else:
                            ao_total_dict[sv_id] = ao
                            # print(ao_total_dict[sv_id])
            except ValueError as e:
                sys.stderr.write("Error! Right:{} at {}\n".format(e, region_r))
                continue
            ################################################################################
            for sv_id in ao_total_dict:
                sv_type, _anno_can, bp1_chrm, _position, _bp2 = sv_id.split(",")
                bp1_pos = int(_position)
                bp2_chrm, bp2_pos = _bp2.split(":")
                bp2_pos = int(bp2_pos)

                sr = ao_total_dict[sv_id]
                try:
                    dp1 = bam_file.count(
                        region="{}:{}-{}".format(bp1_chrm, bp1_pos - 1, bp1_pos)
                    )
                except ValueError as e:
                    print(f"{e} at {bp1_chrm}:{bp1_pos}")
                    dp1 = 0
                try:
                    dp2 = bam_file.count(
                        region="{}:{}-{}".format(bp2_chrm, bp2_pos - 1, bp2_pos)
                    )
                except ValueError as e:
                    print(f"{e} at {bp2_chrm}:{bp2_pos}")
                    dp2 = 0
                pso = sr / (sr + (dp1 + dp2) / 2)
                if (
                    (dp1 + dp2) / 2 >= dp_cutoff
                    and sr >= sr_cutoff
                    and pso >= pso_cutoff
                ):
                    sv_set.add(sv_id)
                    # size_or_sup_position: sup_chrom:sup_position
                    if sv_type == "TRA":
                        bp_distance = 0
                    # size_or_sup_position: length of TDUP/INV
                    else:
                        bp2_chrm = bp1_chrm
                        bp_distance = bp2_pos - bp1_pos

                    bp1_strand, bp2_strand = strand_dict[sv_id]
                    bp1_mode, bp2_mode = mode_dict[sv_id]
                    bp1_gene, bp2_gene = gene_dict[sv_id]

                    anno_flag = _anno_can[0]
                    can_flag = _anno_can[1]

                    if can_flag == "1":
                        can_field = "CANONICAL"
                    else:
                        can_field = "NONCANONICAL"

                    if anno_flag == "0":
                        anno_field = "NEITHER"
                    elif anno_flag == "1":
                        anno_field = "RIGHT"
                    elif anno_flag == "2":
                        anno_field = "LEFT"
                    else:
                        anno_field = "BOTH"

                    mode1 = mode2string(bp1_mode)
                    mode2 = mode2string(bp2_mode)

                    vcf_field_gt = "GT\t0/1"
                    vcf_field_info = "{};SR={};DP1={};DP2={};PSO={:.2g};SVLEN={};SVTYPE={};BOUNDARY={};SVMETHOD=ScanNLS;GENE1={};GENE2={};STRAND1={};STRAND2={};MODE1={};MODE2={};CHR2={};END={}".format(
                        can_field,
                        sr,
                        dp1,
                        dp2,
                        pso,
                        bp_distance,
                        sv_type,
                        anno_field,
                        bp1_gene,
                        bp2_gene,
                        bp1_strand,
                        bp2_strand,
                        mode1,
                        mode2,
                        bp2_chrm,
                        bp2_pos,
                    )
                    vcf_file.write(
                        "{}\t{}\tNLS_{}\t.\t<{}>\t.\t.\t{}\t{}\n".format(
                            bp1_chrm,
                            bp1_pos,
                            ID_count,
                            sv_type,
                            vcf_field_info,
                            vcf_field_gt,
                        )
                    )
                    ID_count += 1
    bam_file.close()
    return True
