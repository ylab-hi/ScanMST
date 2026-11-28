"""Module for BLAT."""

import array
import contextlib
import os
import secrets
import signal
import subprocess
import time
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
    :param min_aligned_len: minimum aligned query length (bp) of a hit; shorter hits count as unaligned

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

    # class-level default so instances created without __init__ (e.g. test doubles) still have it
    min_aligned_len: int = 30

    def __init__(
        self,
        ref_2bit: str,
        port: int,
        output_dir: str,
        *,
        fix_log_file=None,
        is_start_server=False,
        lock=None,
        min_aligned_len: int = 30,
        max_wait: int = 30 * 60,
        retry_interval: int = 30,
    ) -> None:
        """Initialize the blat class.

        :param min_aligned_len: minimum aligned query length (bp) of a BLAT hit;
            shorter hits are treated as unaligned
        :param max_wait: seconds to wait for a usable server on ``port`` before giving up
        :param retry_interval: seconds between retries while the server is not answering
        """
        self.port, self.ref_2bit = port, ref_2bit
        self.min_aligned_len = min_aligned_len
        self.max_wait = max_wait
        self.retry_interval = retry_interval
        self._confirmed_ready = False
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
        """Function for checking whether a blat server for ``self.port`` is running or not.

        A server counts when a gfServer process started on our port is alive (it may still be
        building its index) or when ``gfServer status`` answers on our port (covers servers
        started by other users that psutil cannot inspect).

        :return: the boolean value of whether the server is running or not
        """
        return self._port_process_alive() or self._server_alive()

    def _port_process_alive(self) -> bool:
        """Whether a gfServer process serving ``self.port`` exists; kill stopped ones on our port."""
        flag = False
        for proc in self._search_processing():
            with contextlib.suppress(psutil.NoSuchProcess, psutil.AccessDenied):
                status = proc.status()
                if status in ("running", "sleeping", "disk-sleep"):
                    flag = True
                elif status == "stopped":
                    proc.kill()
        return flag

    def _serves_our_port(self, cmdline: list[str]) -> bool:
        """Whether a gfServer argv is ``... start <host> <port> ...`` with ``self.port``."""
        if "start" not in cmdline:
            return False
        i = cmdline.index("start")
        return len(cmdline) > i + 2 and cmdline[i + 2] == str(self.port)

    def _search_processing(self) -> list[psutil.Process]:
        """Function for searching the gfServer processes that serve ``self.port``.

        Servers on other ports belong to other jobs and are ignored, so they are never reused or killed.

        :return: the list of process of blat server
        """
        result = []
        logger.debug("searching server service")
        for proc in psutil.process_iter(["pid", "name"]):
            with contextlib.suppress(psutil.NoSuchProcess, psutil.AccessDenied):
                if proc.name().lower() == "gfserver" and self._serves_our_port(proc.cmdline()):
                    result.append(proc)
        return result

    def _gfserver_cmd(self, action: str) -> subprocess.CompletedProcess | None:
        """Run ``gfServer <action> localhost <port>``; return None if it cannot run or times out."""
        try:
            return subprocess.run(  # noqa: S603
                [str(self.gfserver), action, "localhost", str(self.port)],
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None

    def _server_alive(self) -> bool:
        """Whether a gfServer answers on ``localhost:self.port`` (``gfServer status`` exits 0)."""
        res = self._gfserver_cmd("status")
        return res is not None and res.returncode == 0

    def _served_files(self) -> list[str]:
        """Basenames of the files served on ``self.port`` (``gfServer files``); empty if unknown."""
        res = self._gfserver_cmd("files")
        if res is None or res.returncode != 0:
            return []
        return [Path(line.strip()).name for line in res.stdout.splitlines() if line.strip()]

    def _starting_server_files(self) -> list[str]:
        """Basenames of the files given to gfServer processes on our port (``start <host> <port> <files>``)."""
        files = []
        for proc in self._search_processing():
            with contextlib.suppress(psutil.NoSuchProcess, psutil.AccessDenied):
                cmdline = proc.cmdline()
                files += [Path(f).name for f in cmdline[cmdline.index("start") + 3 :]]
        return files

    def _check_served_reference(self, served: list[str] | None = None) -> None:
        """Refuse to reuse a server on our port that serves a different reference."""
        expected = Path(self.ref_2bit).name
        served = self._served_files() if served is None else served
        if served and expected not in served:
            msg = f"BLAT server on localhost:{self.port} serves {served}, expected {expected}; use a different --blat-port"
            raise RuntimeError(msg)

    def _start_server(self) -> None:
        """Start gfServer on ``self.port`` from the directory holding the 2bit file.

        The server runs in its own session through ``subprocess.Popen``, not as a non-daemon
        ``multiprocessing.Process`` child: a pool worker that restarts the server would otherwise
        wait for that child at exit forever and hang the whole ProcessPoolExecutor.
        """
        self.is_start_server = True
        logger.debug(f"start server service{self.is_start_server=}")
        cmd = [
            str(self.gfserver),
            "-canStop",
            f"-log={self.log_file_path}",
            "-stepSize=5",
            "start",
            "localhost",
            str(self.port),
            Path(self.ref_2bit).name,
        ]
        logger.trace(f"{cmd=} {self.ref_dir=}")
        # drop a stale "Server ready" line left by a previous server using the same log file
        with contextlib.suppress(OSError):
            Path(self.log_file_path).write_text("")
        self.handle_process = subprocess.Popen(  # noqa: S603
            cmd,
            cwd=self.ref_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
        logger.debug("starting server service")

    def start_server(self) -> None:
        """Reuse a server already on ``self.port``, otherwise start our own.

        :raise RuntimeError: when the server on our port serves a different reference
        """
        if self._server_alive():
            self._check_served_reference()
            self.is_start_server = False
            logger.info(f"reusing BLAT server on localhost:{self.port}")
        elif self._port_process_alive():
            # a server on our port is still building its index
            self._check_served_reference(self._starting_server_files())
            self.is_start_server = False
            logger.info(f"waiting for BLAT server on localhost:{self.port} that is starting up")
        else:
            self._start_server()

    def _own_server_gone(self) -> bool:
        """Whether the server this instance started has exited."""
        if self.handle_process is not None:
            return self.handle_process.poll() is not None
        # worker copies have no process handle: look for any process on our port
        return not self._search_processing()

    def stop_server(self) -> None:
        """Stop the server this job started on ``self.port``; servers of other jobs are left alone."""
        # self open then self close
        logger.trace(f"{self.is_start_server=}")
        if not self.is_start_server:
            return
        if self.handle_process is not None and self.handle_process.poll() is not None:
            # ours exited early (e.g. lost a bind race); whatever serves the port now is not ours
            logger.info("own BLAT server already exited; leaving the server on this port alone")
            self.is_stop_server = True
            return
        logger.info("Stopping  server service")
        if self.handle_process is not None:
            # signal only our own process group: whatever else is on the port may belong to another job
            self._stop_own_process()
        else:
            # instance created with is_start_server=True but without a handle: stop by port
            self._gfserver_cmd("stop")
            for proc in self._search_processing():
                with contextlib.suppress(psutil.NoSuchProcess, psutil.AccessDenied):
                    proc.kill()
        self.is_stop_server = True

    def _stop_own_process(self) -> None:
        """Terminate the gfServer we started (its own session/process group), then kill if it lingers."""
        if self.handle_process is None:
            return
        for sig in (signal.SIGTERM, signal.SIGKILL):
            with contextlib.suppress(ProcessLookupError, PermissionError):
                os.killpg(self.handle_process.pid, sig)
            with contextlib.suppress(subprocess.TimeoutExpired):
                self.handle_process.wait(timeout=10)
                return

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

        logger.trace(f"{self.ref_dir=}")
        cmd = (
            f"{self.gfclient} -minScore=20 -out=pslx -minIdentity={mini_identity} localhost {self.port} . "
            f"{in_fasta} {out_pslx}"
        )
        logger.trace(f"{cmd=}")
        # cwd= instead of os.chdir so a failed query cannot leave the process in ref_dir
        try:
            subprocess.check_call(
                cmd,
                stderr=subprocess.STDOUT,
                shell=True,
                stdout=subprocess.DEVNULL,
                cwd=self.ref_dir,
            )
        finally:
            self._remove(in_fasta)

        return out_pslx

    def _is_server_ready(self) -> bool:
        """Ready when our log says so or the server answers on our port."""
        if os.path.exists(self.log_file_path) and self.is_ready():
            return True
        return self._server_alive()

    def _check_if_self_ready(self, interval: int = 60 * 2) -> None:
        """Wait until the server this job started is ready.

        Gives up waiting (and drops ownership) when our server exited before becoming ready,
        so that ``query`` can reuse or restart a server instead of looping forever.

        :param interval: the interval time for checking the server service
        :raise RuntimeError: when the server is not ready within ``self.max_wait`` seconds
        """
        if not self.is_start_server or self._confirmed_ready:
            return
        deadline = time.monotonic() + self.max_wait
        while not self._is_server_ready():
            if self._own_server_gone():
                logger.warning(f"BLAT server on localhost:{self.port} exited before it was ready")
                self.is_start_server = False
                return
            if time.monotonic() > deadline:
                msg = f"BLAT server on localhost:{self.port} not ready after {self.max_wait} s"
                raise RuntimeError(msg)
            time.sleep(interval)
        self._confirmed_ready = True

    def query(self, in_seq: str, mini_identity: int = 90) -> str:
        """Function for querying the sequence to the server service.

        Retries while the server on ``self.port`` is starting; starts our own server when nothing
        serves the port any more (e.g. the job that owned it finished).

        :param in_seq: the sequence of input sequence
        :param mini_identity: the threshold of the identity for aligning
        :return: the path for PSLX file
        :raise RuntimeError: when no server answers on ``self.port`` within ``self.max_wait`` seconds
        """
        deadline = time.monotonic() + self.max_wait
        while True:
            try:
                self._check_if_self_ready(self.retry_interval)
                return self._query(in_seq, mini_identity)
            except subprocess.CalledProcessError:
                self._confirmed_ready = False
                if time.monotonic() > deadline:
                    msg = f"BLAT server on localhost:{self.port} did not answer within {self.max_wait} s; check --blat-port"
                    raise RuntimeError(msg) from None
                if not self.is_running():
                    logger.warning(f"no BLAT server on localhost:{self.port}; starting one")
                    this_lock = self.lock if self.lock is not None else contextlib.nullcontext()
                    with this_lock:  # one worker starts it, the others then find it on the port
                        self.start_server()
                time.sleep(self.retry_interval)

    @staticmethod
    def hsp_aligned_query_len(hsp: Any) -> int:
        """Number of query bases placed on the genome by an HSP (sum of its aligned query blocks)."""
        return sum(hsp.query_span_all)

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
        min_aligned_len: int = 0,
    ) -> Any:
        """Helper function for querying insertion sequence.

        An HSP is kept when its matched length covers more than ``threshold_identity`` of the query,
        it lies on a primary contig, and it aligns at least ``min_aligned_len`` query bases.
        Shorter alignments are treated as unaligned.
        """
        hsps = blat_result.hsps
        hsps.sort(key=lambda x: x.score, reverse=True)
        hsps = hsps[:top]
        keep_hsps = []
        for hsp in hsps:
            if (
                Blat.hsp_matched_len(hsp) / len(insert_seq) > threshold_identity
                and "_" not in hsp.hit_id
                and Blat.hsp_aligned_query_len(hsp) >= min_aligned_len
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

        # a query shorter than min_aligned_len can never yield an acceptable hit: skip BLAT
        if len(insert_seq) < max(align_len_threshold, self.min_aligned_len):
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
            min_aligned_len=self.min_aligned_len,
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
