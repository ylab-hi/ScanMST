#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import time
import subprocess
import psutil
import argparse
import re
import glob
import random
import math
import logging
import copy
from typing import Iterable
from Bio import SearchIO
from Bio.Seq import Seq
from pyfaidx import Fasta
from align import aligner
from collections import OrderedDict,defaultdict
import psutil


try:
    import pysam
except:
    sys.exit('pysam module not found.\nPlease install it before.')
try:
    import numpy as np
except:
    sys.exit('numpy module not found.\nPlease install it before.')
try:
    import HTSeq
except:
    sys.exit('HTSeq module not found.\nPlease install it before.')


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
        if ',' in msg:
            begin, finish = msg.split(',')
            status_message(begin)
        else:
            finish = msg
    try:
        subprocess.check_call(cmd, shell=True, stderr=subprocess.STDOUT, stdin=subprocess.PIPE)
    except subprocess.CalledProcessError as err:
        error_msg = 'Error happend!: {}\n{}'.format(err, err.output)
    else:
        error_msg = ''
    if not error_msg:
        if msg:
            status_message(finish)
        return True
    else:
        status_message(error_msg)
        return False

def get_softclip_length(read):
    ''' 0 => M
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
    '''
    if read.cigartuples[0][0] == 4:
        # there are soft-clipped segments in left and right both
        if read.cigartuples[-1][0] == 4:
            # length of left soft-clipped segment is bigger
            if read.cigartuples[0][1] > read.cigartuples[-1][1]:
                return read.cigartuples[0][1], read.query_sequence[:read.cigartuples[0][1]], read.ref_start, 2
            # length of right soft-clipped segment is bigger
            else:
                return read.cigartuples[-1][1], read.query_sequence[read.query_length - read.cigartuples[-1][1]:], read.ref_end - 1, 1
        # there are soft-clipped segments in left only
        else:
            return read.cigartuples[0][1], read.query_sequence[:read.cigartuples[0][1]], read.ref_start, 2
    # there are soft-clipped segments in right only
    elif read.cigartuples[-1][0] == 4:
        return read.cigartuples[-1][1], read.query_sequence[read.query_length - read.cigartuples[-1][1]:], read.ref_end - 1, 1
    else:
        return 0, '', -1, 0

