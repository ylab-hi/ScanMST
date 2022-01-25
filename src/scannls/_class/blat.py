# !/usr/bin/env python
"""Module for BLAT.

@Filename:    blat.py
@license:     MIT Licence
@Time:        12/15/21 2:00 PM
"""
import os
import platform
import random
import subprocess
import time
from multiprocessing import Process
from pathlib import Path
from typing import Any
from typing import List
from typing import Tuple

import psutil  # type: ignore
from Bio import SearchIO  # type: ignore
from loguru import logger

from ..type import LoggerType
from .basicClass import Insertion
from .basicClass import NovelInsertion


class Blat:
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

    ENV_SCANNLS_SERVER_IS_READY = "SCANNLS_SERVER_IS_READY"
    ENV_DIR = Path.home() / ".scannls"

    def __init__(
        self,
        ref_2bit: str,
        logger: LoggerType,
        port: int,
        output_dir: str,
        fix_log_file=None,
        is_start_server=False,
    ) -> None:
        """Initialize the blat class."""
        self.port, self.ref_2bit = port, ref_2bit
        self.output_dir = output_dir
        self.ran_id = random.getrandbits(30)
        self.is_start_server = is_start_server
        self.is_stop_server = False
        self.logger = logger
        self.fix_log_file = fix_log_file
        self.env_file = Blat.ENV_DIR / f"env_{platform.node()}.conf"
        self.handle_process = None

    def set_env(self, is_ready: bool = False) -> None:
        """Set the environment variable for blat."""
        with open(self.env_file, "w") as f:
            f.write(f"{Blat.ENV_SCANNLS_SERVER_IS_READY}={is_ready}\n")

    @property
    def env_is_ready(self) -> bool:
        """Check if the blat server is ready."""
        flag = False
        if not self.env_file.exists():
            return flag

        with open(self.env_file) as f:
            content_list = [line.strip() for line in f.readlines()]
            flag = content_list[0].split("=")[1] == "True"
        return flag

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
    def log_file_path(self) -> str:
        """Property for log_file_path, which is the path of log file for blat."""
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
        self.logger.debug("check if the server starts")
        # self open check self log file
        flag = False
        if os.path.exists(self.log_file_path) and self.is_start_server:
            with open(self.log_file_path) as f:
                for line in f:
                    if "Server ready" in line:
                        # set env variable
                        self.set_env(is_ready=True)
                        flag = True
                        break
        else:
            # when do not start server check env variable
            flag = self.env_is_ready and not self.is_start_server
        return flag

    def is_running(self) -> bool:
        """Function for checking whether the blat server is running or not.

        :return: the boolean value of whether the server is running or not
        """
        return bool(self._search_processing())

    def _search_processing(self) -> List[psutil.Process]:
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
        try:
            subprocess.check_call(cmd, shell=True)
        except KeyboardInterrupt:
            self.stop_server()

    def _start_server(self) -> None:
        """The GfServer should run at the directory.

        where gfServer gfClient and hg38.2bit located.
        """
        cwd = os.path.abspath(os.getcwd())
        logger.debug(os.getcwd())

        # change to use_blat directory
        os.chdir(self.ref_dir)
        logger.trace(f"{self.ref_dir=}")
        logger.trace(f"{os.getcwd()}")

        if os.path.exists(self.log_file_path):
            os.remove(self.log_file_path)

        cmd = (
            f"gfServer -canStop -log={self.log_file_path} -stepSize=5 start "
            f"localhost {self.port} {os.path.basename(self.ref_2bit)}"
        )
        logger.trace(f"{cmd=}")
        self.handle_process = Process(target=self._run_cmd, args=[cmd])  # type: ignore
        self.handle_process.start()
        self.logger.debug("starting server service")
        os.chdir(cwd)
        logger.trace(f"{os.getcwd()}")

    def start_server(self) -> None:
        """Function for starting the server service, if the server is not running.

        we will start the server service
        """
        self.logger.debug(f"is running {self.is_running()}")
        running_flag = self.is_running()
        if not running_flag:
            self._start_server()
            self.is_start_server = True
        else:
            self.is_start_server = False

    def stop_server(self) -> None:
        """Function for stopping the server service, if the server is running."""
        # self open then self close
        if self.is_start_server:
            self.logger.debug("stopping server service")
            if self.handle_process is None:
                raise SystemExit("Not Start Server Want to Stop")
            self.handle_process.close()
            self._remove(str(self.env_file))
            self._remove(self.log_file_path)  # remove temp log file
            self.is_stop_server = True

    def _query(self, in_seq: str, mini_identity: int = 90) -> str:
        """Function is help function in order to using gfClient.

        to query 'in_seq' to generate alignment file (in PSL format).

        :param mini_identity: the threshold of the identity for aligning
        :param in_seq: sequence of softclipped segment
        :return: the path for PSL file
        """
        self.logger.debug("querying the sequence")
        ran_id = random.getrandbits(30)
        in_fasta = os.path.join(self.output_dir, f"{ran_id}.fasta")
        with open(in_fasta, "w", buffering=1) as fasta_file:
            fasta_file.write(f">{ran_id}\n")
            fasta_file.write(f"{in_seq}\n")

        out_psl = os.path.join(self.output_dir, f"{ran_id}.psl")

        cwd = os.path.abspath(os.getcwd())
        logger.trace(os.getcwd())

        os.chdir(self.ref_dir)
        logger.trace(f"{self.ref_dir=}")
        logger.trace(os.getcwd())
        cmd = "gfClient -minScore=20 -minIdentity={} localhost {} . {} {} 2> /dev/null".format(
            mini_identity, self.port, in_fasta, out_psl
        )
        logger.trace(f"{cmd=}")
        subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True)
        os.chdir(cwd)
        logger.trace(os.getcwd())
        self._remove(in_fasta)

        return out_psl

    def _wait_ready(self, interval: int = 30) -> None:
        """Function for waiting the server service to be ready.

        :param interval: the interval time for checking the server service
        """
        if self.is_ready() and self.is_running():
            return

        if not self.is_ready() and self.is_running():
            time.sleep(interval)
        elif not self.is_ready() and not self.is_running():
            self.start_server()

        self._wait_ready()

    def query(self, in_seq: str, mini_identity: int = 90) -> str:
        """Function for querying the sequence to the server service.

        :param in_seq: the sequence of input sequence
        :param mini_identity: the threshold of the identity for aligning
        :return: the path for PSL file
        """
        # check if need to start server service
        self.start_server()

        out_psl = ""

        if self.is_ready():
            out_psl = self._query(in_seq, mini_identity)
        #  not ready as no others is running
        elif not self.is_running():
            # self open check if ready
            if self.is_ready():
                out_psl = self._query(in_seq, mini_identity)
            else:
                self._wait_ready()
                out_psl = self._query(in_seq, mini_identity)
        else:  # not ready but others or self is running
            while self.is_running():
                try:
                    out_psl = self._query(in_seq, mini_identity)
                except subprocess.CalledProcessError:
                    time.sleep(30)
                else:
                    break

            if not self.is_running():
                self._wait_ready()
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
    def _remove(file: str) -> None:
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
        mapq_dict = {
            1: 60,
            2: 3,
            3: 2,
        }

        for hsp in hsps:
            if (
                hsp.ident_pct / 100 >= threshold_identity
                and hsp.query_span / in_seq_len >= threshold_identity
            ):
                num_of_locations += 1

        if 4 <= num_of_locations <= 9:
            return 1
        return mapq_dict.get(num_of_locations, 0)

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
