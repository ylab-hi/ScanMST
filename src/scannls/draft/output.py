#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Output module for draft scannls.

2021-10-01:
detect_sv_from_cigar current_output a list of putative NLS events
modify SV tag endswith ";", SV:Z:XXX;YYY;ZZZ;
"""
from typing import Dict
from typing import TextIO

__funcs__ = {"output_bedpe_file", "aggregate_candidates", "similar_hit"}


def output_bedpe_file(
    sr_dict: Dict, group_dict: Dict, prefix: str, splice_bin: int
) -> TextIO:
    """OUTPUT BEDPE FILE.

    :param sr_dict: sv candidate to number of supporting reads(SR) dictionary
    :param group_dict: sv candidate to group of events dictionary, connected
        chimeric reads are included in one group
    :param prefix: current_output file prefix
    :param splice_bin: bin size for splice site searching
    :type sr_dict: dict
    :type group_dict: dict
    :type prefix: str
    :type splice_bin: int
    :return: current_output BEDPE file
    :rtype: str
    """
    output = open("{}.sv.bedpe".format(prefix), "w")
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
                f"{chrm1}\t{pos1 - splice_bin}\t{pos1 + splice_bin}\t{chrm2}"
                f"\t{pos2 - splice_bin}\t{pos2 + splice_bin}\tgroup_{num_of_group}"
                f"\t{sr}\t{strand1}\t{strand2}\n"
            )
    output.close()
    return output


def aggregate_candidates(in_dict: Dict, len_cutoff: int = 10) -> Dict:
    """Aggregate candidates."""
    if len_cutoff == 0:
        return in_dict
    else:
        discarded_items = set()
        in_dict_len = len(in_dict)
        items = list(in_dict.keys())
        for i in range(in_dict_len):
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
            if m not in discarded_items:
                out_dict[m] = in_dict[m]
        return out_dict


def similar_hit(r1: str, r2: str, len_cutoff: int = 10) -> bool:
    """Similar hit.

    :param r1: read 1
    :param r2: read 2
    :param len_cutoff: length cutoff
    :return: bool value
    """
    r1_type, r1_can, a1, a2, strand_1 = r1.split("\t")
    r2_type, r2_can, b1, b2, strand_2 = r2.split("\t")

    chrm_a1, pos_a1 = a1.split(":")
    chrm_b1, pos_b1 = a2.split(":")

    chrm_a2, pos_a2 = b1.split(":")
    chrm_b2, pos_b2 = b2.split(":")

    if r1_type != r2_type:
        return False
    else:
        if chrm_a1 == chrm_a2 and chrm_b1 == chrm_b2:
            if (
                abs(int(pos_a1) - int(pos_a2)) <= len_cutoff
                and abs(int(pos_b1) - int(pos_b2)) <= len_cutoff
            ):
                return True
            else:
                return False

        elif chrm_a1 == chrm_b2 and chrm_b1 == chrm_a2:
            if (
                abs(int(pos_a1) - int(pos_b2)) <= len_cutoff
                and abs(int(pos_b1) - int(pos_a2)) <= len_cutoff
            ):
                return True
            else:
                return False
        else:
            return False
