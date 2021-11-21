#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import os
import subprocess
import sys

try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)


def remove(infile):
    if os.path.isfile(infile):
        os.remove(infile)


def remove_files(pattern):
    for j in glob.iglob(pattern):
        remove(j)


def status_message(msg):
    print(msg)
    sys.stdout.flush()


def run_cmd(cmd, msg=None):
    status_message(cmd)
    if msg:
        if "," in msg:
            begin, finish = msg.split(",")
            status_message(begin)
        else:
            finish = msg
    try:
        subprocess.check_call(
            cmd, shell=True, stderr=subprocess.STDOUT, stdin=subprocess.PIPE
        )
    except subprocess.CalledProcessError as err:
        error_msg = "Error happend!: {}\n{}".format(err, err.output)
    else:
        error_msg = ""
    if not error_msg:
        if msg:
            status_message(finish)
        return True
    else:
        status_message(error_msg)
        return False


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
    if read.cigartuples[0][0] == 4:
        # there are soft-clipped segments in left and right both
        if read.cigartuples[-1][0] == 4:
            # length of left soft-clipped segment is bigger
            if read.cigartuples[0][1] > read.cigartuples[-1][1]:
                return (
                    read.cigartuples[0][1],
                    read.query_sequence[: read.cigartuples[0][1]],
                    read.ref_start,
                    2,
                )
            # length of right soft-clipped segment is bigger
            else:
                return (
                    read.cigartuples[-1][1],
                    read.query_sequence[read.query_length - read.cigartuples[-1][1] :],
                    read.ref_end - 1,
                    1,
                )
        # there are soft-clipped segments in left only
        else:
            return (
                read.cigartuples[0][1],
                read.query_sequence[: read.cigartuples[0][1]],
                read.ref_start,
                2,
            )
    # there are soft-clipped segments in right only
    elif read.cigartuples[-1][0] == 4:
        return (
            read.cigartuples[-1][1],
            read.query_sequence[read.query_length - read.cigartuples[-1][1] :],
            read.reference_end - 1,
            1,
        )
    else:
        return 0, "", -1, 0


'''
def transcript_upstream_part_determiner(strand1, strand2, mode1, mode2) -> bool:
    """Determine the transcript upstream part using strand and mode information
    :param strand1: strand for breakpoint1
    :param strand2: strand for breakpoint2
    :param mode1: mode for breakpoint1
    :param mode2: mode for breakpoint2
    :type strand1: str(-/+)
    :type strand2: str(-/+)
    :type mode1: int(1/2)
    :type mode2: int(1/2)
    :return: wheather breakpoint1 is the upstream segment of the transcript or not
    :rtype: bool
    """
    bp1_is_upstream = None
    if strand1 == "+" and strand2 == "-":
        if mode1 == 1 and mode2 == 1:
            is_bp1_upstream = True
        elif mode1 == 2 and mode2 == 2:
            is_bp1_upstream = False
    elif strand1 == "-" and strand2 == "+":
        if mode1 == 1 and mode2 == 1:
            is_bp1_upstream = False
        elif mode1 == 2 and mode2 == 2:
            is_bp1_upstream = True
    elif strand1 == "+" and strand2 == "+":
        if mode1 == 1 and mode2 == 2:
            is_bp1_upstream = True
        elif mode1 == 2 and mode2 == 1:
            is_bp1_upstream = False
    elif strand1 == "-" and strand2 == "-":
        if mode1 == 1 and mode2 == 2:
            is_bp1_upstream = False
        elif mode1 == 2 and mode2 == 1:
            is_bp1_upstream = True
    return is_bp1_upstream
'''
