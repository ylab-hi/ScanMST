#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ===========================================================
"""
2021-10-01:
detect_sv_from_cigar output a list of putative NLS events
modify SV tag endswith ";", SV:Z:XXX;YYY;ZZZ;

"""
import argparse
import os
import re
import subprocess
import sys
import textwrap
import time
from collections import defaultdict

from Bio.Seq import Seq
from loguru import logger
from pyfaidx import Fasta
from pyfaidx import FastaNotFoundError

from .. import __version__
from ..classes import Blat
from ..classes import LengthAction
from ..classes import Read
from ..classes import ReadsConnecter
from ..classes import Series

try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)

__funcs__ = {"output_bedpe_file", "aggregate_candidates", "similar_hit"}


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
                f"{chrm1}\t{pos1 - splice_bin}\t{pos1 + splice_bin}\t{chrm2}\t{pos2 - splice_bin}\t{pos2 + splice_bin}\tgroup_{num_of_group}\t{sr}\t{strand1}\t{strand2}\n"
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
