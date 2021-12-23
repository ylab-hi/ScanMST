# !/usr/bin/env python
# -*- coding:utf-8 -*-
"""Module for BLAT.

@Filename:    blat.py
@license:     MIT Licence
@Time:        12/15/21 2:00 PM
"""
import os
import random
import subprocess
import time
from multiprocessing import Process
from typing import Any
from typing import List
from typing import Tuple

import psutil  # type: ignore
from Bio import SearchIO  # type: ignore
from loguru import logger
from loguru._logger import Logger

from .basicClass import Insertion
from .basicClass import NovelInsertion


class Blat(object):
    """Blat class is used to integrate the blat service.

     (gfServer and gfClient) so that we can query certain sequences from the genome shamelessly

    :param ref_2bit: the path of reference for blat alignment
    :param logger: the logger for logging
    :param port: the port of server service for blat alignment
    :param output_dir: the path for storing alignment result

    :Example:

    >>> from  loguru import   logger
    >>> blat = Blat(ref_2bit='reference.2bit', logger=logger, output_dir='/tmp')
    >>> blat.is_running()
    True
    >>> blat.stop_server()
    >>> blat.is_running()
    False
    >>> blat.start_server()
    >>> blat.is_running()
    True
    >>> blat.query(in_seq='ATCGTCC')
    /tmp/tmp_gfClient_in_seq_out.psl
    >>> blat.query_insertion(insert_seq='ATCGTCC')
    True, Insertion(chr1:1-9:+,ATCGTCC, None, TPA, chr1:1, chr1:9 )
    >>> blat.query_insertion(insert_seq='ATCCATCC')
    False, NovelInsertion(ATCCATCC:0)
    >>> blat.query_insertion(insert_seq="ATCG")
    False, NovelInsertion(ATCG:10)
    """

    def __init__(
        self,
        ref_2bit: str,
        logger: Logger,
        port: int,
        output_dir: str,
        fix_log_file=None,
        is_start_server=True,
    ) -> None:
        """Initialize the blat class."""
        self.port, self.ref_2bit = port, ref_2bit
        self.output_dir = output_dir
        self.ran_id = random.getrandbits(30)
        self.is_start_server = is_start_server
        self.logger = logger
        self.fix_log_file = fix_log_file

    @property
    def ref_dir(self) -> str:
        """Property for ref_dir, which is the path of reference for blat.

        :return: the absolute path of reference dir
        """
        if self.ref_2bit.startswith("~"):
            abs_2bit = os.path.join(
                os.path.expanduser("~"), self.ref_2bit.replace("~/", "")
            )
            ref_dir = os.path.dirname(abs_2bit)
        else:
            abs_2bit = os.path.abspath(self.ref_2bit)
            ref_dir = os.path.dirname(abs_2bit)
        return ref_dir

    @property
    def log_file(self) -> str:
        """Property for log_file, which is the path of log file for blat."""
        return (
            f"{self.ref_dir}/gfserver.temp.{self.ran_id}.log"
            if self.fix_log_file is None
            else self.fix_log_file
        )

    def is_ready(self) -> bool:
        """Function for checking whether the blat server is ready or not.

        after starting the server service

        :return: the boolean value of whether the server is ready or not
        """
        flag = False
        self.logger.debug("check if the server starts")
        if os.path.exists(self.log_file):
            with open(self.log_file) as f:
                for line in f:
                    if "Server ready" in line:
                        flag = True
        return flag

    def is_running(self) -> bool:
        """Function for checking whether the blat server is running or not.

        :return: the boolean value of whether the server is running or not
        """
        return True if self._search_processing() else False

    def _search_processing(self) -> List:
        """Function for searching the process of blat server.

        in current system

        :return: the list of process of blat server
        """
        result = []
        self.logger.debug("searching server service")
        for proc in psutil.process_iter(["pid", "name"]):
            if "gfServer".lower() == proc.name().lower() and proc.cmdline():
                result.append(proc)
        return result

    def _run_cmd(self, cmd: str) -> None:
        """Function is used to run the command in the system.

        :param cmd: the command to be run
        """
        subprocess.check_call(cmd.split())

    def _start_server(self) -> Process:
        """The GfServer should run at the directory.

        where gfServer gfClient and hg38.2bit located.
        """
        cwd = os.path.abspath(os.getcwd())
        logger.debug(os.getcwd())

        # change to use_blat directory
        os.chdir(self.ref_dir)
        logger.trace(f"{self.ref_dir=}")
        logger.trace(f"{os.getcwd()}")

        if os.path.exists(self.log_file):
            os.remove(self.log_file)

        cmd = (
            f"gfServer -canStop -log={self.log_file} -stepSize=5 start "
            f"localhost {self.port} {os.path.basename(self.ref_2bit)}"
        )
        logger.trace(f"{cmd=}")
        process = Process(target=self._run_cmd, args=[cmd])  # type: ignore
        process.start()
        self.logger.debug("starting server service")
        os.chdir(cwd)
        logger.trace(f"{os.getcwd()}")
        return process

    def start_server(self) -> None:
        """Function for starting the server service, if the server is not running.

        we will start the server service
        """
        running_flag = self.is_running()
        if not running_flag:
            self._start_server()
        else:
            self.is_start_server = False

    def stop_server(self) -> None:
        """Function for stopping the server service, if the server is running."""
        procs = self._search_processing()
        self.logger.debug("stopping server service")
        for proc in procs:
            proc.kill()

    def _query(self, in_seq: str, mini_identity: int = 90) -> str:
        """Function is help function in order to using gfClient.

        to query 'in_seq' to generate alignment file (in PSL format).

        :param mini_identity: the threshold of the identity for aligning
        :param in_seq: sequence of softclipped segment
        :return: the path for PSL file
        """
        self.logger.debug("querying the sequence")
        ran_id = random.getrandbits(30)
        in_fasta = os.path.join(self.output_dir, "{}.fasta".format(ran_id))
        with open(in_fasta, "w", buffering=1) as fasta_file:
            fasta_file.write(">{}\n".format(ran_id))
            fasta_file.write("{}\n".format(in_seq))

        out_psl = os.path.join(self.output_dir, "{}.psl".format(ran_id))

        cwd = os.path.abspath(os.getcwd())
        logger.trace(os.getcwd())

        os.chdir(self.ref_dir)
        logger.trace(f"{self.ref_dir=}")
        logger.trace(os.getcwd())
        cmd = "gfClient -minScore=20 -minIdentity={} localhost {} . {} {} > /dev/null".format(
            mini_identity, self.port, in_fasta, out_psl
        )
        logger.trace(f"{cmd=}")
        try:
            subprocess.check_call(cmd.split(), stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as err:
            raise SystemExit(f"{err} {err.output}") from err

        os.chdir(cwd)
        logger.trace(os.getcwd())
        self._remove(in_fasta)

        return out_psl

    def _wait_ready(self, interval: int = 30) -> None:
        """Function for waiting the server service to be ready.

        :param interval: the interval time for checking the server service
        """
        while not self.is_ready():
            time.sleep(interval)

    def query(self, in_seq: str, mini_identity: int = 90) -> str:
        """Function for querying the sequence to the server service.

        :param in_seq: the sequence of input sequence
        :param mini_identity: the threshold of the identity for aligning
        :return: the path for PSL file
        """
        if self.is_start_server:
            if self.is_ready():
                out_psl = self._query(in_seq, mini_identity)
            else:
                self._wait_ready()
                out_psl = self._query(in_seq, mini_identity)
        else:
            out_psl = self._query(in_seq, mini_identity)

        return out_psl

    @staticmethod
    def _query_insertion(
        blat_result: Any, insert_seq: str, threshold_identity: float, top: int
    ) -> Any:
        """Helper function for querying insertion sequence."""
        hsps = blat_result.hsps
        hsps.sort(key=lambda x: x.score, reverse=True)
        hsps = hsps[:top]
        keep_hsp = []
        for hsp in hsps:
            if (sum(hsp.hit_span_all) - hsp.mismatch_num) / len(
                insert_seq
            ) > threshold_identity:
                keep_hsp.append(hsp)
        hit = len(keep_hsp)

        return hit, keep_hsp

    def query_insertion(
        self,
        insert_seq: str,
        threshold_identity: float = 0.99,
        top: int = 3,
        align_len_threshold: int = 20,
    ) -> Any:
        """Function for querying the insertion sequence to the server service.

        the function is a specific version of the :func: 'Blat.query'.

        :param insert_seq: insertion sequence
        :param threshold_identity: the threshold of the identity for aligning
        :param top: the top number of the alignments
        :param align_len_threshold: the threshold of the insertion sequence length
        :return: insertion sequence alignment in NamedTuple format
        """
        flag = False  # flag for checking the insertion  if its hit is only one

        if len(insert_seq) < align_len_threshold:
            return flag, NovelInsertion(hit_num=0, query_sequence=insert_seq)

        out_blat = self.query(in_seq=insert_seq)
        try:
            blat_result = SearchIO.read(out_blat, "blat-psl")
        except ValueError:
            return flag, NovelInsertion(hit_num=0, query_sequence=insert_seq)

        hit, keep_hsp = Blat._query_insertion(
            blat_result, insert_seq, threshold_identity, top
        )

        if hit == 1:
            top_hsp = keep_hsp[0]
            flag = True

            ref_chrom, position, strand, cigar, num_of_mismatch = self.psl2sam(
                top_hsp, in_seq_len=len(insert_seq)
            )
            return flag, Insertion(
                hit_num=1,
                chrom=ref_chrom,
                ref_start=position,
                strand=strand,
                cigarstring=cigar,
                mapq=60,
                nm=num_of_mismatch,
                query_sequence=insert_seq,
            )
        else:
            return flag, NovelInsertion(hit_num=hit, query_sequence=insert_seq)

    @staticmethod
    def _remove(file):
        """Function for removing the file.

        :param file: the path of the file
        """
        if os.path.exists(file):
            os.remove(file)

    @staticmethod
    def _calculate_mapq(hsps: Any, in_seq_len: int, threshold_identity: float) -> int:
        """Function is used to calculate map quality of the insertion.

        :param hsps: the list of hsp after aligning the insertion sequence
        :param in_seq_len: the length of the input sequence
        :param threshold_identity: the threshold of the identity for aligning
        :return: the map quality of the insertion
        """
        num_of_locations = 0

        for hsp in hsps:
            if (
                hsp.ident_pct / 100 >= threshold_identity
                and hsp.query_span / in_seq_len >= threshold_identity
            ):
                num_of_locations += 1
        if num_of_locations == 1:
            mapq = 60
        elif num_of_locations == 2:
            mapq = 3
        elif num_of_locations == 3:
            mapq = 2
        elif 4 <= num_of_locations <= 9:
            mapq = 1
        else:
            mapq = 0
        return mapq

    def fetch_mapq(self, in_seq: str, threshold_identity: float) -> Any:
        """Function is used to fetch the map quality of the insertion.

        :param in_seq: the input sequence
        :param threshold_identity: the threshold of the identity for aligning
        :return: the top hit of the insertion sequence, and the map quality of the insertion
        """
        psl_file = self.query(in_seq=in_seq)

        try:
            blat = SearchIO.read(psl_file, "blat-psl")
        except ValueError:
            self.logger.warning(f"No Blat hit found {in_seq}")
            return None, None
        else:
            hsps = blat.hsps
            hsps.sort(key=lambda k: k.score, reverse=True)
            top_hsp = hsps[0]
            Blat._remove(psl_file)
            mapq = Blat._calculate_mapq(hsps, len(in_seq), threshold_identity)
        return top_hsp, mapq

    @staticmethod
    def psl2sam(hsp: Any, in_seq_len: int) -> Tuple[str, int, str, str, int]:
        """Convert the top HSP in PSL file to SAM fields chrom, reference_start.

        strand, cigarstring, num_of_mismatch. The function try to implement
        the psl2sam.pl script and return the cigar and mapping position
        estimated from psl file

        :param hsp: the selected HSP form BLAT
        :param in_seq_len: the length of the input sequence
        :return: chrom, reference_start, strand, cigarstring, num_of_mismatch

        .. note::
            hsp is 0-based, as same as python [ )
        """
        cigar = ""
        query_start = hsp.query_start
        query_end = hsp.query_end

        _strand = hsp.query_strand_all[0]  # may need replace by query_strand
        ref_start, ref_end = hsp.hit_range
        ref_chrom = hsp.hit_id
        num_of_mismatch = hsp.mismatch_num

        if _strand == -1:
            query_start = in_seq_len - hsp.query_end
            query_end = in_seq_len - hsp.query_start
        if query_start:
            # 5'-end clipping
            cigar += str(query_start) + "S"
        x = hsp.query_span_all
        if _strand == -1:
            y = [
                in_seq_len - item[1] for item in hsp.query_range_all
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

        if in_seq_len != query_end:
            # 3'-end clipping
            end3 = in_seq_len - query_end
            cigar += str(end3) + "S"
        strand = "+" if _strand == 1 else "-"

        return ref_chrom, ref_start, strand, cigar, num_of_mismatch
