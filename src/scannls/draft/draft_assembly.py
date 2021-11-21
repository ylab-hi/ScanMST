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

from . import __version__
from .classes import Blat
from .classes import LengthAction
from .classes import Read
from .classes import ReadsConnecter
from .classes import Series
from .common import get_softclip_length
from .externals import external_tool_checking
from .externals import softclipped_seq2SA_tag
from .utils import extract_splice_sites
from .utils import infer_sv_from_connected_reads

try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)


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
