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
from Bio import SearchIO

from . import __version__

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


def checkIfProcessRunning(processName):
    """
    Check if there is any running process that contains the given name processName.
    """
    # Iterate over the all the running process
    for proc in psutil.process_iter():
        try:
            # Check if process name contains the given name string.
            if processName.lower() in proc.name().lower():
                return True
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return False


def start_gfServer(ref_2bit, timeout=300, port=88888, output_dir="/tmp"):
    """gfServer should run at the directory where gfServer, gfClient and hg38.2bit located"""
    sys.stderr.write("Starting BLAT gfServer\n")
    if ref_2bit.startswith("~"):
        abs_2bit = os.path.join(os.path.expanduser("~"), ref_2bit.replace("~/", ""))
        ref_dir = os.path.dirname(abs_2bit)
        base_2bit = os.path.basename(abs_2bit)
    else:
        abs_2bit = os.path.abspath(ref_2bit)
        ref_dir = os.path.dirname(abs_2bit)
        base_2bit = os.path.basename(abs_2bit)
    cwd = os.path.abspath(os.getcwd())
    sys.stderr.write("BLAT 2bit file location: {}\n".format(ref_dir))
    sys.stderr.write("Current directory: {}\n".format(cwd))
    if os.path.isabs(output_dir):
        log_dir = output_dir
    else:
        log_dir = os.path.join(cwd, output_dir)
    # change to blat directory
    os.chdir(ref_dir)
    sys.stderr.write("Current directory: {}\n".format(os.getcwd()))
    try:
        cmd = "gfServer -canStop -log={0}/gfserver.temp.log -stepSize=5 start localhost {1} {2} &".format(
            log_dir, port, base_2bit
        )
        print(cmd)
        start = time.time()
        ret = subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as err:
        print(
            "Execution failed for starting BLAT gfServer:",
            err,
            err.output,
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        if not ret:
            interval = 10
            while not gfserver_tester(log_dir) and time.time() - start < timeout:
                time.sleep(interval)
            if gfserver_tester(log_dir):
                sys.stdout.write("gfServer is ready for use.\n")
                remove("{}/gfserver.temp.log".format(output_dir))
            else:
                sys.stdout.write("Timeout!\n")
                stop_gfServer(port, output_dir)
                sys.exit(1)
        else:
            sys.stderr.write(
                "Something wrong in {0}/gfserver.temp.log\n".format(output_dir)
            )
    os.chdir(cwd)
    sys.stderr.write("Current directory: {}\n".format(os.getcwd()))


def stop_gfServer(port=88888, output_dir="/tmp"):
    try:
        subprocess.check_call(
            "gfServer stop localhost {0}".format(port),
            stderr=subprocess.STDOUT,
            shell=True,
        )
        print("gfServer stop localhost {0}".format(port))
    except subprocess.CalledProcessError as e:
        print("Execution failed for stoping BLAT gfServer:", e.output, file=sys.stderr)
        sys.exit(1)
    remove_files("{}/*.temp.log".format(output_dir))


def gfClient_query(in_seq, ref_2bit, port=88888, output_dir="/tmp"):
    ran_id = random.getrandbits(30)
    in_fasta = os.path.join(output_dir, "{}.fasta".format(ran_id))

    fasta_file = open(in_fasta, "w", buffering=1)
    fasta_file.write(">{}\n".format(ran_id))
    fasta_file.write("{}\n".format(in_seq))
    fasta_file.close()

    out_psl = os.path.join(output_dir, "{}.psl".format(ran_id))

    if ref_2bit.startswith("~"):
        abs_2bit = os.path.join(os.path.expanduser("~"), ref_2bit.replace("~/", ""))
        ref_dir = os.path.dirname(abs_2bit)
        base_2bit = os.path.basename(abs_2bit)
    else:
        abs_2bit = os.path.abspath(ref_2bit)
        ref_dir = os.path.dirname(abs_2bit)
        base_2bit = os.path.basename(abs_2bit)
    cwd = os.path.abspath(os.getcwd())

    if os.path.isabs(output_dir):
        log_dir = output_dir
    else:
        log_dir = os.path.join(cwd, output_dir)

    os.chdir(ref_dir)
    # print(os.getcwd())
    try:
        cmd = "gfClient -minScore=20 -minIdentity=0 localhost {} {} {} {} > /dev/null".format(
            port, ref_dir, in_fasta, out_psl
        )
        # print(cmd)
        print(in_seq, in_fasta)
        ret = subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as err:
        print(
            "Execution failed for starting BLAT gfServer:",
            err,
            err.output,
            file=sys.stderr,
        )
        sys.exit(1)
    os.chdir(cwd)
    # print(os.getcwd())
    remove(in_fasta)
    # remove(out_psl)
    return out_psl


def psl2sam(hsp, query_seq_len):
    """psl2sam try to implement the psl2sam.pl script and return the cigar and mapping position estimated from psl file"""

    cigar = ""
    query_start = hsp.query_start
    query_end = hsp.query_end

    _strand = hsp.query_strand_all[0]  # may need replace by qery_strand
    ref_start, ref_end = hsp.hit_range
    ref_chrom = hsp.hit_id
    num_of_mismatch = hsp.mismatch_num

    soft_len = 0
    if _strand == -1:
        query_start = query_seq_len - hsp.query_end
        query_end = query_seq_len - hsp.query_start
    if query_start:
        # 5'-end clipping
        soft_len = query_start
        cigar += str(query_start) + "S"
    x = hsp.query_span_all
    if _strand == -1:
        y = [
            query_seq_len - item[1] for item in hsp.query_range_all
        ]  # may need replace by query_start_all when the bug is fixed in Biopython
    else:
        y = [
            item[0] for item in hsp.query_range_all
        ]  # may need replace by query_start_all when the bug is fixed in Biopython
    z = hsp.hit_start_all
    y0, z0 = y[0], z[0]
    for i in range(1, len(hsp)):
        ly = y[i] - y[i - 1] - x[i - 1]
        lz = z[i] - z[i - 1] - x[i - 1]
        if ly < lz:
            # del: the reference gap is longer
            cigar += str(y[i] - y0) + "M"
            if lz - ly >= 10:
                cigar += str(lz - ly) + "N"
            else:
                cigar += str(lz - ly) + "D"
            y0, z0 = y[i], z[i]
        elif lz < ly:
            # ins: the query gap is longer
            cigar += str(z[i] - z0) + "M"
            cigar += str(ly - lz) + "I"
            y0, z0 = y[i], z[i]

    cigar += str(query_end - y0) + "M"
    # print(cigar)
    # return cigar, soft_len
    if query_seq_len != query_end:
        # 3'-end clipping
        end3 = query_seq_len - query_end
        if end3 > soft_len:
            soft_len = end3
        cigar += str(end3) + "S"
    # return cigar, soft_len
    if _strand == 1:
        strand = "+"
    elif _strand == -1:
        strand = "-"
    return ref_chrom, ref_start + 1, strand, cigar, num_of_mismatch


def blat_mapq_calculator(hsps, query_len, blat_ident_pct_cutoff=0.95):
    """
    60 = Uniquely mapping
    3  = Maps to 2 locations in the target
    2  = Maps to 3 locations in the target
    1  = Maps to 4-9 locations in the target
    0  = Maps to 10 or more locations in the target
    """
    num_of_locations = 0
    for hsp in hsps:
        if (
            hsp.ident_pct / 100 >= blat_ident_pct_cutoff
            and hsp.query_span / query_len >= blat_ident_pct_cutoff
        ):
            num_of_locations += 1
    if num_of_locations == 1:
        mapq = 60
    elif num_of_locations == 2:
        mapq = 3
    elif num_of_locations == 3:
        mapq = 2
    elif num_of_locations >= 4 and num_of_locations <= 9:
        mapq = 1
    else:
        mapq = 0
    return mapq


def cigar_validity(cigar_str):
    """
    40M25N5M then cigartuple is [('40', 'M'), ('25', 'N'), ('5', 'M')]
    :param cigar_str: BLAT generated cigarstring from 'softclipped_seq2SA_tag'
    :type cigar_str: str
    :return: valid cigarstring
    :rtype: str
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


def softclipped_seq2SA_tag(
    in_seq,
    read_length,
    read_strand,
    read_mode,
    ref_2bit,
    port,
    mapq_cutoff,
    max_allowed_nm,
    output_dir="/tmp",
    blat_ident_pct_cutoff=0.95,
):
    """
    find new supplementary alignment for alignment which has longer softclipped segment

    in_seq: softclipped segment of alignment
    read_mode: read mode of alignment
    read_strand: read strand of alignment
    return: putative supplementary alignment of the alignment (input softclipped segment)
    """
    in_seq_len = len(in_seq)
    psl_file = gfClient_query(in_seq, ref_2bit, port, output_dir)
    chimeric_aln_str = ""
    try:
        blat = SearchIO.read(psl_file, "blat-psl")
    except ValueError as err:
        print("No BLAT hit! {}".format(in_seq), err, file=sys.stderr)
    else:
        hsps = blat.hsps
        hsps.sort(key=lambda k: k.score, reverse=True)
        top_hsp = hsps[0]
        remove(psl_file)
        __mapq = blat_mapq_calculator(hsps, in_seq_len, blat_ident_pct_cutoff)
        if (
            top_hsp.ident_pct / 100 >= blat_ident_pct_cutoff
            and top_hsp.query_span / in_seq_len >= blat_ident_pct_cutoff
        ):
            __chrm_sa, __pos_sa, __strand_sa, __cigar_sa_partial, __nm_sa = psl2sam(
                top_hsp, in_seq_len
            )
            if read_strand == __strand_sa:
                # same strand: different reads mode
                # MS(1) ~ SM(2) or SM(2) ~ MS(1)
                if read_mode == 1:
                    __cigar_sa = "{0}S{1}".format(
                        read_length - in_seq_len, __cigar_sa_partial
                    )  # SM
                else:
                    __cigar_sa = "{1}{0}S".format(
                        read_length - in_seq_len, __cigar_sa_partial
                    )  # MS
            else:
                # opposite strand: same reads mode
                # MS(1) ~ MS(1) or SM(2) ~ SM(2)
                if read_mode == 1:
                    __cigar_sa = "{1}{0}S".format(
                        read_length - in_seq_len, __cigar_sa_partial
                    )  # MS
                else:
                    __cigar_sa = "{0}S{1}".format(
                        read_length - in_seq_len, __cigar_sa_partial
                    )  # SM
            valid_cigar_sa = cigar_validity(__cigar_sa)
            if __mapq < mapq_cutoff and int(__nm_sa) < max_allowed_nm:
                chimeric_aln_str = "{},{},{},{},{},{};".format(
                    __chrm_sa, __pos_sa, __strand_sa, valid_cigar_sa, __mapq, __nm_sa
                )
            else:
                chimeric_aln_str = ""

    return chimeric_aln_str


def external_tool_checking(blat=False):
    """checking dependencies are installed"""
    if blat:
        software = ["sambamba", "gfClient", "gfServer"]
    else:
        software = ["sambamba"]
    cmd = "which"
    for each in software:
        try:
            path = subprocess.check_output([cmd, each], stderr=subprocess.STDOUT)
            path = str(path, "utf-8")
        except subprocess.CalledProcessError:
            print(
                "Checking for '" + each + "': ERROR - could not find '" + each + "'",
                file=sys.stderr,
            )
            print("Exiting.", file=sys.stderr)
            sys.exit(0)
        print("Checking for '" + each + "': found " + path)


def gfserver_tester(log_dir):
    log_file = "{0}/gfserver.temp.log".format(log_dir)
    if os.path.exists(log_file):
        with open(log_file) as f:
            for line in f:
                if "Server ready" in line:
                    return True
                elif "gfServer aborted" in line or "error" in line:
                    return False
    return False
