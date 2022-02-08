#!/usr/bin/env python
"""Helper functions."""
# ===============================================================================
from itertools import tee


def pairwise(iterable):
    """The same as itertools.pairwise from python 3.10.

    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def to_wt_and_mt_fasta(input_dict: dict, output_prefix: str):
    """Output WT and MT fasta."""
    with open(f"{output_prefix}.WT.fa", "w") as wt_fa, open(
        f"{output_prefix}.MT.fa", "w"
    ) as mt_fa:
        for trx_id in input_dict:
            _metaexon_list = input_dict[trx_id]
            mt_seq = ""
            for _idx, _metaexon in enumerate(_metaexon_list, 1):
                wt_seq = _metaexon.wt_seq
                mt_seq += _metaexon.mt_seq
                if wt_seq:
                    wt_fa.write(f">{trx_id}_{_idx}\n{wt_seq}\n")
            if mt_seq:
                mt_fa.write(f">{trx_id}\n{mt_seq}\n")
