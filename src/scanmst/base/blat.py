"""Module for BLAT."""

import array
import contextlib
import os
import secrets
import subprocess
import time
from multiprocessing import Process
from pathlib import Path
from typing import Any

import psutil
from Bio import SearchIO
from loguru import logger

from scanmst.blat import load_gfclient, load_gfserver

from .basic_class import Insertion, NovelInsertion


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
    /tmp/tmp_gfClient_in_seq_out.pslx
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
        port: int,
        output_dir: str,
        *,
        fix_log_file=None,
        is_start_server=False,
        lock=None,
    ) -> None:
        """Initialize the blat class."""
        self.port, self.ref_2bit = port, ref_2bit
        self.output_dir = output_dir
        self.ran_id = secrets.randbits(42)
        self.is_start_server = is_start_server
        self.is_stop_server = False
        self.fix_log_file = fix_log_file
        self.handle_process = None
        self.gfserver = load_gfserver()
        self.gfclient = load_gfclient()
        self.lock = lock

        if self.ref_2bit.startswith("~"):
            abs_2bit = os.path.join(
                os.path.expanduser("~"),
                self.ref_2bit.replace("~/", ""),
            )
            self.ref_dir = os.path.dirname(abs_2bit)
        else:
            abs_2bit = os.path.abspath(self.ref_2bit)
            self.ref_dir = os.path.dirname(abs_2bit)

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

        :except: raise when the process start server but log file is not created
        :return: the boolean value of whether the server is ready or not
        """
        if not os.path.exists(self.log_file_path):
            raise RuntimeError(
                f"the BLAT has started server but the log file does not exist: {self.log_file_path}",
            )

        this_lock = self.lock if self.lock is not None else contextlib.nullcontext()
        with this_lock, open(self.log_file_path) as f:
            return any("Server ready" in line for line in f)

    def is_running(self) -> bool:
        """Function for checking whether the blat server is running or not.

        :return: the boolean value of whether the server is running or not
        """
        flag = False
        for proc in self._search_processing():
            if proc.status() in ("running", "sleeping"):
                flag = True
            elif proc.status() == "stopped":
                proc.kill()
        return flag

    def _search_processing(self) -> list[psutil.Process]:
        """Function for searching the process of blat server.

        in current system

        :return: the list of process of blat server
        """
        result = []
        logger.debug("searching server service")
        for proc in psutil.process_iter(["pid", "name"]):
            with contextlib.suppress(psutil.NoSuchProcess):
                if "gfServer".lower() == proc.name().lower() and proc.cmdline():
                    result.append(proc)
        return result

    def _run_cmd(self, cmd: str) -> None:
        """Function is used to run the command in the system.

        :param cmd: the command to be run
        """
        try:
            subprocess.check_call(cmd, shell=True)
        except (KeyboardInterrupt, subprocess.CalledProcessError) as e:
            if isinstance(e, KeyboardInterrupt) and not self.is_stop_server:
                self.stop_server()

    def _start_server(self) -> None:
        """The GfServer should run at the directory.

        where gfServer gfClient and hg38.2bit located.
        """
        self.is_start_server = True
        logger.debug(f"start server service{self.is_start_server=}")
        cwd = Path.cwd().absolute()
        logger.debug(Path.cwd().as_posix())

        # change to use_blat directory
        os.chdir(self.ref_dir)
        logger.trace(f"{self.ref_dir=}")
        logger.trace(f"{Path().cwd()}")

        cmd = (
            f"{self.gfserver} -canStop -log={self.log_file_path} -stepSize=5 start "
            f"localhost {self.port} {os.path.basename(self.ref_2bit)}"
        )
        logger.trace(f"{cmd=}")
        self.handle_process = Process(target=self._run_cmd, args=[cmd])  # type: ignore
        if self.handle_process is None:
            raise ValueError("handle process is None")
        self.handle_process.start()
        logger.debug("starting server service")
        os.chdir(cwd)
        logger.trace(f"{Path().cwd()}")

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
        # self open then self close
        logger.trace(f"{self.is_start_server=}")
        if self.is_start_server:
            logger.info("Stopping  server service")

            for proc in self._search_processing():
                proc.kill()

            self.is_stop_server = True

    def _query(self, in_seq: str, mini_identity: int = 90) -> str:
        """Function is help function in order to using gfClient.

        to query 'in_seq' to generate alignment file (in PSL format).

        :param mini_identity: the threshold of the identity for aligning
        :param in_seq: sequence of softclipped segment
        :return: the path for PSL file
        """
        logger.debug("querying the sequence")
        ran_id = secrets.randbits(42)
        in_fasta = os.path.join(self.output_dir, f"{ran_id}.fasta")
        with open(in_fasta, "w", buffering=1) as fasta_file:
            fasta_file.write(f">{ran_id}\n")
            fasta_file.write(f"{in_seq}\n")

        out_pslx = os.path.join(self.output_dir, f"{ran_id}.pslx")

        cwd = Path.cwd().absolute()
        logger.trace(f"{Path().cwd()}")

        os.chdir(self.ref_dir)
        logger.trace(f"{self.ref_dir=}")
        logger.trace(f"{Path().cwd()}")
        cmd = (
            f"{self.gfclient} -minScore=20 -out=pslx -minIdentity={mini_identity} localhost {self.port} . "
            f"{in_fasta} {out_pslx}"
        )
        logger.trace(f"{cmd=}")
        subprocess.check_call(
            cmd,
            stderr=subprocess.STDOUT,
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        os.chdir(cwd)
        logger.trace(f"{Path().cwd()}")
        self._remove(in_fasta)

        return out_pslx

    def _check_if_self_ready(self, interval: int = 60 * 2) -> None:
        """Function for waiting the server service to be ready.

        :param interval: the interval time for checking the server service
        """
        # check self start server and sever is running
        if self.is_start_server:
            # check log file
            # it will take some time from start the server to create the log file
            if os.path.exists(self.log_file_path) and self.is_ready():
                return
            #  not ready yet, wait for a while
            time.sleep(interval)
            self._check_if_self_ready()

    def query(self, in_seq: str, mini_identity: int = 90) -> str:
        """Function for querying the sequence to the server service.

        :param in_seq: the sequence of input sequence
        :param mini_identity: the threshold of the identity for aligning
        :return: the path for PSLX file
        """
        while (
            self.is_start_server or self.is_running()
        ):  # self or other is running service
            try:
                self._check_if_self_ready()  # if self start blocking, then wait for the server service to be ready
                out_pslx = self._query(in_seq, mini_identity)
            except subprocess.CalledProcessError:
                time.sleep(60 * 2)  # wait for other's service to be ready
            else:
                return out_pslx

        # other kill service and self start
        self.start_server()
        return self.query(in_seq, mini_identity)

    @staticmethod
    def hsp_matched_len(hsp: Any) -> int:
        """Helper function for HSP matched length calculation."""
        gap_num = hsp.hit_gap_num
        # `N` splicing junction
        if gap_num >= 10:
            hsp_matched_length = sum(hsp.hit_span_all) - hsp.mismatch_num
        # `D` deletion
        else:
            hsp_matched_length = sum(hsp.hit_span_all) - hsp.mismatch_num - gap_num
        return hsp_matched_length

    @staticmethod
    def _query_insertion(
        blat_result: Any,
        insert_seq: str,
        threshold_identity: float,
        top: int,
    ) -> Any:
        """Helper function for querying insertion sequence."""
        hsps = blat_result.hsps
        hsps.sort(key=lambda x: x.score, reverse=True)
        hsps = hsps[:top]
        keep_hsps = []
        for hsp in hsps:
            if (
                Blat.hsp_matched_len(hsp) / len(insert_seq) > threshold_identity
                and "_" not in hsp.hit_id
            ):
                keep_hsps.append(hsp)
        hit = len(keep_hsps)

        mapq = Blat._calculate_mapq(keep_hsps)
        top_hsp = keep_hsps[0] if hit >= 1 else ""

        return hit, top_hsp, mapq

    def query_insertion(
        self,
        insert_seq: str,
        threshold_identity: float = 0.90,
        top: int = 5,
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

        # https://biopython.org/docs/1.84/api/Bio.SearchIO.BlatIO.html
        # HSPFragment carris query sequence (query) and hit sequence(hit)
        try:
            blat_result = SearchIO.read(out_blat, "blat-psl", pslx=True)
        except ValueError:
            return flag, NovelInsertion(hit_num=0, query_sequence=insert_seq)

        hit, top_hsp, mapq = Blat._query_insertion(
            blat_result,
            insert_seq,
            threshold_identity,
            top,
        )

        # keep the top one hit
        if hit >= 1:
            flag = True

            ref_chrom, position, strand, cigar, num_of_mismatch = self.psl2sam(
                top_hsp,
                in_seq_len=len(insert_seq),
            )

            dummy_qualities = array.array("B", [40] * len(insert_seq))
            return flag, Insertion(
                hit_num=1,
                chrom=ref_chrom,
                ref_start=position,
                strand=strand,
                cigarstring=cigar,
                mapq=mapq,
                nm=num_of_mismatch,
                query_sequence=insert_seq,
                query_qualities=list(dummy_qualities),
            )
        return flag, NovelInsertion(hit_num=hit, query_sequence=insert_seq)

    @staticmethod
    def _remove(file: str) -> None:
        """Function for removing the file.

        :param file: the path of the file
        """
        if os.path.exists(file):
            os.remove(file)

    @staticmethod
    def _calculate_mapq(hsps: Any) -> int:
        """Function is used to calculate map quality of the insertion.
        We adapted the way of calculation in Minimap2.
            MAPQ=min([60 × (S1 - S2)/S1], 60)
            S1: best alignment score
            S2: the second-best alignment score
            They both passed the query identity threshold.

        :param hsps: the list of hsp ordered by score after aligning the insertion sequence
        :return: the map quality of the insertion
        """
        if len(hsps) == 0:
            return 0
        elif len(hsps) == 1:
            return 60
        else:
            return int(min(60 * (hsps[0].score - hsps[1].score) / hsps[0].score, 60))

    @staticmethod
    def psl2sam(hsp: Any, in_seq_len: int) -> tuple[str, int, str, str, int]:
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
        # NM tag in SAM is edit distance which includes SNV, INS and DEL
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
                    num_of_mismatch += lz - ly
                y0, z0 = y[i], z[i]
            elif lz < ly:
                # ins: the query gap is longer
                cigar += str(z[i] - z0) + "M"
                cigar += str(ly - lz) + "I"
                num_of_mismatch += ly - lz
                y0, z0 = y[i], z[i]

        cigar += str(query_end - y0) + "M"

        if in_seq_len != query_end:
            # 3'-end clipping
            end3 = in_seq_len - query_end
            cigar += str(end3) + "S"
        strand = "+" if _strand == 1 else "-"

        return ref_chrom, ref_start, strand, cigar, num_of_mismatch

    @staticmethod
    def obtain_variants_stats(
        hsp: Any,
        in_seq_len: int,
        large_indel_len_threshold: int = 4,
    ) -> tuple[float, float, float]:
        """Obtain variants stats from HSP."""

        _strand = hsp.query_strand_all[0]
        strand = "+" if _strand == 1 else "-"
        if strand == "-":
            query_ranges = [
                (in_seq_len - end, in_seq_len - start)
                for start, end in hsp.query_range_all
            ]
        else:
            query_ranges = hsp.query_range_all
        # reference ranges
        hit_ranges = hsp.hit_range_all

        # calculate the reference alignment length
        reference_length = 0
        for start, end in hit_ranges:
            reference_length += end - start

        deletion_length_total = 0
        substitution_num = hsp.mismatch_num
        del_outlier_num, ins_outlier_num = 0, 0
        insertion_num = 0
        deletion_num = 0
        for idx, _ in enumerate(query_ranges[:-1]):
            # Gaps between consecutive query blocks
            query_gap = query_ranges[idx + 1][0] - query_ranges[idx][1]
            hit_gap = hit_ranges[idx + 1][0] - hit_ranges[idx][1]

            # If query gap > hit gap, it's an insertion
            if query_gap > hit_gap:
                insertion_length = query_gap - hit_gap
                insertion_num += 1
                if insertion_length >= large_indel_len_threshold:
                    ins_outlier_num += 1
            # If hit gap > query gap, it's a deletion
            elif hit_gap > query_gap:
                deletion_length = hit_gap - query_gap
                # treat it as intron when deletion_length >= 10
                if deletion_length < 10:
                    deletion_length_total += deletion_length
                    deletion_num += 1
                    if deletion_length >= large_indel_len_threshold:
                        del_outlier_num += 1

        # Total reference length includes aligned blocks + deletions
        total_reference_length = reference_length + deletion_length_total

        ins_fraction = ins_outlier_num / total_reference_length
        del_fraction = del_outlier_num / total_reference_length
        subs_fraction = substitution_num / total_reference_length

        return subs_fraction, ins_fraction, del_fraction
