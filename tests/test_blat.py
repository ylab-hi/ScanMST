"""Tests for `scanmst.blat` platform/download helpers and `scanmst.base.blat.Blat`."""

from __future__ import annotations

import contextlib
import os
import platform
import signal
import subprocess
from pathlib import Path
from types import SimpleNamespace

import psutil
import pytest
import requests
from Bio import SearchIO

from scanmst import Blat
from scanmst import blat as blat_module
from scanmst.base import Insertion, NovelInsertion


class TestBlatPureHelpers:
    """Platform- and binary-independent unit tests for BLAT static helpers and downloaders."""

    def test_get_platform_base_url(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr(platform, "system", lambda: "Linux")
        monkeypatch.setattr(platform, "machine", lambda: "x86_64")
        sys_name, url = blat_module.get_platform_base_url()
        assert sys_name == "linux"
        assert "linux.x86_64" in url

        monkeypatch.setattr(platform, "system", lambda: "Darwin")
        monkeypatch.setattr(platform, "machine", lambda: "arm64")
        assert blat_module.get_platform_base_url()[1].endswith("macOSX.arm64/")

        monkeypatch.setattr(platform, "machine", lambda: "x86_64")
        assert blat_module.get_platform_base_url()[1].endswith("macOSX.x86_64/")

        monkeypatch.setattr(platform, "system", lambda: "Windows")
        with pytest.raises(NotImplementedError, match="not supported"):
            blat_module.get_platform_base_url()

    def test_download_blat_tools_mocked(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, fake_logger) -> None:
        fake_pkg_file = tmp_path / "__init__.py"
        fake_pkg_file.write_text("")
        monkeypatch.setattr(blat_module, "__file__", str(fake_pkg_file))
        monkeypatch.setattr(blat_module, "get_platform_base_url", lambda: ("linux", "https://example.invalid/"))

        class _FakeResp:
            def __enter__(self) -> _FakeResp:
                return self

            def __exit__(self, *_args: object) -> None:
                return None

            def raise_for_status(self) -> None:
                return None

            def iter_content(self, chunk_size: int = 8192):
                yield b"#!/bin/sh\n"

        monkeypatch.setattr(requests, "get", lambda *_a, **_kw: _FakeResp())
        blat_module.download_blat_tools(fake_logger)
        for tool in ("blat", "gfServer", "gfClient", "faToTwoBit"):
            p = tmp_path / "linux" / tool
            assert p.exists()
            assert p.stat().st_mode & 0o100

        # Calling a second time hits the already-exists fast path without network calls
        monkeypatch.setattr(requests, "get", lambda *_a, **_kw: (_ for _ in ()).throw(AssertionError("should not download")))
        blat_module.download_blat_tools(fake_logger)

        # Download failure cleans up partial file and raises SystemExit
        (tmp_path / "linux" / "blat").unlink()
        monkeypatch.setattr(requests, "get", lambda *_a, **_kw: (_ for _ in ()).throw(OSError("network down")))
        with pytest.raises(SystemExit, match="Could not download blat"):
            blat_module.download_blat_tools(fake_logger)

    def test_hsp_len_and_mapq(self) -> None:
        hsp_del = SimpleNamespace(query_span_all=[30, 20], hit_span_all=[30, 20], mismatch_num=1, hit_gap_num=4, score=45)
        hsp_intron = SimpleNamespace(query_span_all=[30, 20], hit_span_all=[30, 20], mismatch_num=1, hit_gap_num=500, score=30)
        assert Blat.hsp_aligned_query_len(hsp_del) == 50
        assert Blat.hsp_matched_len(hsp_del) == 50 - 1 - 4
        assert Blat.hsp_matched_len(hsp_intron) == 50 - 1

        assert Blat._calculate_mapq([]) == 0
        assert Blat._calculate_mapq([hsp_del]) == 60
        assert Blat._calculate_mapq([hsp_del, hsp_intron]) == int(60 * (45 - 30) / 45)

    def test_psl2sam_forward_reverse_del_intron_and_ins(self) -> None:
        class _FakePSL:
            def __init__(
                self,
                *,
                strand: int,
                q_start: int,
                q_end: int,
                q_spans: list[int],
                q_ranges: list[tuple[int, int]],
                h_starts: list[int],
                h_range: tuple[int, int],
                mismatch: int = 1,
            ) -> None:
                self.query_start = q_start
                self.query_end = q_end
                self.query_strand_all = [strand]
                self.query_span_all = q_spans
                self.query_range_all = q_ranges
                self.hit_start_all = h_starts
                self.hit_range = h_range
                self.hit_range_all = [(s, s + sp) for s, sp in zip(h_starts, q_spans, strict=True)]
                self.hit_id = "chr1"
                self.mismatch_num = mismatch

            def __len__(self) -> int:
                return len(self.query_span_all)

        # Forward strand with 5bp 5' clip, 20M, 4D (<10 -> D), 15M, 200N (>=10 -> N), 10M, 5bp 3' clip
        fwd = _FakePSL(
            strand=1,
            q_start=5,
            q_end=50,
            q_spans=[20, 15, 10],
            q_ranges=[(5, 25), (25, 40), (40, 50)],
            h_starts=[100, 124, 339],
            h_range=(100, 349),
            mismatch=1,
        )
        chrom, pos, strand, cigar, nm = Blat.psl2sam(fwd, in_seq_len=55)
        assert (chrom, pos, strand, cigar, nm) == ("chr1", 100, "+", "5S20M4D15M200N10M5S", 5)

        subs_f, ins_f, del_f = Blat.obtain_variants_stats(fwd, in_seq_len=55, large_indel_len_threshold=4)
        # aligned hit blocks = 45 + 4bp del (<10) = 49 total_reference_length
        assert subs_f == pytest.approx(1 / 49)
        assert ins_f == 0.0
        assert del_f == pytest.approx(1 / 49)

        # Reverse strand with a 5bp query insertion between blocks
        rev = _FakePSL(
            strand=-1,
            q_start=0,
            q_end=45,
            q_spans=[20, 20],
            q_ranges=[(0, 20), (25, 45)],  # reversed on len=50: [(30,50), (5,25)] -> sorted block starts y=[30,5]?
            h_starts=[500, 520],
            h_range=(500, 540),
            mismatch=0,
        )
        # Note Biopython PSL on '-' orders blocks in target order, so reversed query ranges run [(5, 25), (30, 50)]
        rev.query_range_all = [(25, 45), (0, 20)]
        chrom_r, pos_r, strand_r, cigar_r, nm_r = Blat.psl2sam(rev, in_seq_len=50)
        assert (chrom_r, pos_r, strand_r, cigar_r, nm_r) == ("chr1", 500, "-", "5S20M5I20M", 5)
        subs_r, ins_r, del_r = Blat.obtain_variants_stats(rev, in_seq_len=50, large_indel_len_threshold=4)
        assert (subs_r, ins_r, del_r) == (0.0, pytest.approx(1 / 40), 0.0)

    def test_server_lifecycle_and_query_insertion_mocked(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr("scanmst.base.blat.load_gfserver", lambda: tmp_path / "gfServer")
        monkeypatch.setattr("scanmst.base.blat.load_gfclient", lambda: tmp_path / "gfClient")
        monkeypatch.setattr(os.path, "expanduser", lambda _p: str(tmp_path))
        monkeypatch.setattr("time.sleep", lambda _s: None)

        b_tilde = Blat(ref_2bit="~/sub/dummy.2bit", port=19999, output_dir=str(tmp_path), min_aligned_len=30)
        assert b_tilde.ref_dir == str(tmp_path / "sub")

        log_path = tmp_path / "gf.log"
        blat_inst = Blat(
            ref_2bit=str(tmp_path / "dummy.2bit"),
            port=19999,
            output_dir=str(tmp_path),
            fix_log_file=str(log_path),
            min_aligned_len=30,
        )

        # _remove existing and missing
        victim = tmp_path / "rm_me.txt"
        victim.write_text("x")
        Blat._remove(str(victim))
        assert not victim.exists()
        Blat._remove(str(victim))

        # query(): one CalledProcessError while a server is up on our port -> retry -> succeeds via _query
        blat_inst.is_start_server = False
        monkeypatch.setattr(blat_inst, "is_running", lambda: True)
        call_state = {"q": 0}

        def _fake_check_call(*_a: object, **_kw: object) -> int:
            call_state["q"] += 1
            if call_state["q"] == 1:
                raise subprocess.CalledProcessError(1, "gfClient")
            return 0

        monkeypatch.setattr(subprocess, "check_call", _fake_check_call)
        out_pslx = blat_inst.query("A" * 40)
        assert out_pslx.endswith(".pslx")
        assert call_state["q"] == 2

        # query_insertion: too short -> NovelInsertion
        flag_short, ins_short = blat_inst.query_insertion("ACGT", align_len_threshold=20)
        assert flag_short is False
        assert isinstance(ins_short, NovelInsertion)

        # query_insertion: SearchIO ValueError -> NovelInsertion
        monkeypatch.setattr(blat_inst, "query", lambda **_kw: str(tmp_path / "empty.pslx"))
        monkeypatch.setattr(SearchIO, "read", lambda *_a, **_kw: (_ for _ in ()).throw(ValueError("empty psl")))
        flag_err, ins_err = blat_inst.query_insertion("A" * 40)
        assert flag_err is False
        assert isinstance(ins_err, NovelInsertion)

        # query_insertion: 1 hit -> Insertion, and 0 hits -> NovelInsertion
        good_hsp = SimpleNamespace(
            score=40,
            hit_id="chr1",
            query_span_all=[40],
            hit_span_all=[40],
            mismatch_num=0,
            hit_gap_num=0,
            query_start=0,
            query_end=40,
            query_strand_all=[1],
            query_range_all=[(0, 40)],
            hit_start_all=[100],
            hit_range=(100, 140),
            hit_range_all=[(100, 140)],
            __len__=lambda self: 1,
        )

        # SimpleNamespace instance method __len__ needs class definition for len()
        class _OneBlockHSP:
            score = 40
            hit_id = "chr1"
            query_span_all = [40]
            hit_span_all = [40]
            mismatch_num = 0
            hit_gap_num = 0
            query_start = 0
            query_end = 40
            query_strand_all = [1]
            query_range_all = [(0, 40)]
            hit_start_all = [100]
            hit_range = (100, 140)
            hit_range_all = [(100, 140)]

            def __len__(self) -> int:
                return 1

        _ = good_hsp
        monkeypatch.setattr(SearchIO, "read", lambda *_a, **_kw: SimpleNamespace(hsps=[_OneBlockHSP()]))
        flag_ok, ins_ok = blat_inst.query_insertion("A" * 40)
        assert flag_ok is True
        assert isinstance(ins_ok, Insertion)
        assert (ins_ok.chrom, ins_ok.ref_start, ins_ok.cigarstring) == ("chr1", 100, "40M")

        monkeypatch.setattr(SearchIO, "read", lambda *_a, **_kw: SimpleNamespace(hsps=[]))
        flag_zero, ins_zero = blat_inst.query_insertion("A" * 40)
        assert flag_zero is False
        assert isinstance(ins_zero, NovelInsertion)


@pytest.mark.skipif(
    not blat_module.load_blat().exists(),
    reason="BLAT executables are downloaded from UCSC on demand and are not installed",
)
class TestBlatServerHelpers:
    """Tests that construct a real Blat instance using tmp_path so no files leak into the repo."""

    @pytest.fixture
    def blat(self, tmp_path: Path) -> Blat:
        two_bit = tmp_path / "dummy.2bit"
        two_bit.write_bytes(b"")
        return Blat(ref_2bit=str(two_bit), port=88888, output_dir=str(tmp_path))

    def test_paths_and_is_ready(self, blat: Blat) -> None:
        assert Path(blat.ref_dir).is_absolute()
        assert Path(blat.log_file_path).is_absolute()

        with pytest.raises(RuntimeError, match="log file does not exist"):
            blat.is_ready()

        Path(blat.log_file_path).write_text("Starting...\n")
        assert blat.is_ready() is False

        Path(blat.log_file_path).write_text("Server ready for queries!\n")
        assert blat.is_ready() is True

    def test_is_running(self, blat: Blat, mocker) -> None:
        killed: list[str] = []

        class _Proc:
            def __init__(self, name: str, status: str) -> None:
                self._name = name
                self._status = status

            def name(self) -> str:
                return self._name

            def cmdline(self) -> list[str]:
                return ["gfServer", "-canStop", "start", "localhost", str(blat.port), "dummy.2bit"]

            def status(self) -> str:
                return self._status

            def kill(self) -> None:
                killed.append(self._name)

        mocker.patch.object(blat, "_server_alive", return_value=False)
        mocker.patch("psutil.process_iter", return_value=[_Proc("gfServer", "running")])
        assert blat.is_running() is True

        mocker.patch("psutil.process_iter", return_value=[_Proc("gfServer", "stopped")])
        assert blat.is_running() is False
        assert killed == ["gfServer"]


def _start_server_in_worker(two_bit: str, port: int, out_dir: str) -> int:
    """Pool-worker entry point: start a (fake) gfServer and return its pid."""
    worker_blat = Blat(ref_2bit=two_bit, port=port, output_dir=out_dir)
    worker_blat._start_server()
    assert worker_blat.handle_process is not None
    return worker_blat.handle_process.pid


class _FakeGfProc:
    """psutil.Process stand-in for a gfServer started as ``gfServer <action> localhost <port> ...``."""

    def __init__(self, port: int, status: str = "sleeping", action: str = "start", name: str = "gfServer") -> None:
        self._port, self._status, self._action, self._name = port, status, action, name
        self.killed = False
        self.deny_kill = False

    def name(self) -> str:
        return self._name

    def cmdline(self) -> list[str]:
        return ["/x/gfServer", "-canStop", "-log=x.log", "-stepSize=5", self._action, "localhost", str(self._port), "hg38.2bit"]

    def status(self) -> str:
        return self._status

    def kill(self) -> None:
        if self.deny_kill:
            raise psutil.AccessDenied(4242)
        self.killed = True


class _FakePopen:
    """subprocess.Popen stand-in whose exit status is controlled by the test."""

    def __init__(self, returncode: int | None = None) -> None:
        self.returncode = returncode
        self.pid = 4242

    def poll(self) -> int | None:
        return self.returncode

    def wait(self, timeout: float | None = None) -> int | None:
        return self.returncode


class TestBlatPortScopedServer:
    """Server reuse, start, stop and retry are scoped to ``Blat.port`` (no binaries needed)."""

    PORT = 19999

    @pytest.fixture
    def blat_inst(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Blat:
        monkeypatch.setattr("scanmst.base.blat.load_gfserver", lambda: tmp_path / "gfServer")
        monkeypatch.setattr("scanmst.base.blat.load_gfclient", lambda: tmp_path / "gfClient")
        monkeypatch.setattr("time.sleep", lambda _s: None)
        return Blat(
            ref_2bit=str(tmp_path / "hg38.2bit"),
            port=self.PORT,
            output_dir=str(tmp_path),
            fix_log_file=str(tmp_path / "gf.log"),
        )

    def test_search_processing_ignores_other_ports_and_subcommands(self, blat_inst: Blat, monkeypatch: pytest.MonkeyPatch) -> None:
        ours = _FakeGfProc(self.PORT)
        procs = [ours, _FakeGfProc(55752), _FakeGfProc(self.PORT, action="status"), _FakeGfProc(self.PORT, name="python")]
        monkeypatch.setattr("psutil.process_iter", lambda *_a, **_kw: procs)
        assert blat_inst._search_processing() == [ours]

    def test_is_running_only_for_our_port(self, blat_inst: Blat, monkeypatch: pytest.MonkeyPatch) -> None:
        # the 22Rv1 hang: another job's server on a different port must not count as ours
        monkeypatch.setattr("psutil.process_iter", lambda *_a, **_kw: [_FakeGfProc(55752)])
        monkeypatch.setattr(blat_inst, "_server_alive", lambda: False)
        assert blat_inst.is_running() is False

        # a server psutil cannot see (other user) still counts when it answers on our port
        monkeypatch.setattr(blat_inst, "_server_alive", lambda: True)
        assert blat_inst.is_running() is True

        # stopped process on our port is killed and not counted; other ports are untouched
        stopped, other = _FakeGfProc(self.PORT, status="stopped"), _FakeGfProc(55752, status="stopped")
        monkeypatch.setattr("psutil.process_iter", lambda *_a, **_kw: [stopped, other])
        monkeypatch.setattr(blat_inst, "_server_alive", lambda: False)
        assert blat_inst.is_running() is False
        assert stopped.killed is True
        assert other.killed is False

    def test_server_alive_and_served_files(self, blat_inst: Blat, monkeypatch: pytest.MonkeyPatch) -> None:
        calls: list[list[str]] = []

        def _fake_run(cmd: list[str], **_kw: object) -> subprocess.CompletedProcess:
            calls.append(cmd)
            return subprocess.CompletedProcess(cmd, 0, stdout="hg38.2bit\n", stderr="")

        monkeypatch.setattr(subprocess, "run", _fake_run)
        assert blat_inst._server_alive() is True
        assert blat_inst._served_files() == ["hg38.2bit"]
        assert calls[0][1:] == ["status", "localhost", str(self.PORT)]
        assert calls[1][1:] == ["files", "localhost", str(self.PORT)]

        monkeypatch.setattr(subprocess, "run", lambda cmd, **_kw: subprocess.CompletedProcess(cmd, 255, stdout="", stderr="refused"))
        assert blat_inst._server_alive() is False
        assert blat_inst._served_files() == []

        def _timeout(cmd: list[str], **_kw: object) -> subprocess.CompletedProcess:
            raise subprocess.TimeoutExpired(cmd, 30)

        monkeypatch.setattr(subprocess, "run", _timeout)
        assert blat_inst._server_alive() is False

    def test_start_server_reuse_mismatch_wait_or_own(self, blat_inst: Blat, monkeypatch: pytest.MonkeyPatch) -> None:
        popen_calls: list[tuple[list[str], dict]] = []

        def _fake_popen(cmd: list[str], **kw: object) -> _FakePopen:
            popen_calls.append((cmd, kw))
            return _FakePopen()

        monkeypatch.setattr(subprocess, "Popen", _fake_popen)

        # alive on our port with our reference -> reuse, no new server
        monkeypatch.setattr(blat_inst, "_server_alive", lambda: True)
        monkeypatch.setattr(blat_inst, "_served_files", lambda: ["hg38.2bit"])
        blat_inst.start_server()
        assert blat_inst.is_start_server is False
        assert popen_calls == []

        # alive on our port with another reference -> refuse
        monkeypatch.setattr(blat_inst, "_served_files", lambda: ["mm10.2bit"])
        with pytest.raises(RuntimeError, match="different --blat-port"):
            blat_inst.start_server()

        # process on our port still building its index -> wait for it, no new server
        monkeypatch.setattr(blat_inst, "_server_alive", lambda: False)
        monkeypatch.setattr(blat_inst, "_port_process_alive", lambda: True)
        blat_inst.start_server()
        assert blat_inst.is_start_server is False
        assert popen_calls == []

        # nothing on our port -> start our own, detached, from ref_dir
        monkeypatch.setattr(blat_inst, "_port_process_alive", lambda: False)
        blat_inst.start_server()
        assert blat_inst.is_start_server is True
        (cmd, kw), = popen_calls
        assert cmd[-4:] == ["start", "localhost", str(self.PORT), "hg38.2bit"]
        assert kw["cwd"] == blat_inst.ref_dir
        assert kw["start_new_session"] is True

    def test_review_edge_cases(self, blat_inst: Blat, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        # a server on our port that is still indexing another reference must not be reused silently
        starting = _FakeGfProc(self.PORT)  # cmdline serves hg38.2bit
        monkeypatch.setattr("psutil.process_iter", lambda *_a, **_kw: [starting])
        monkeypatch.setattr(blat_inst, "_server_alive", lambda: False)
        mm10 = Blat(ref_2bit=str(tmp_path / "mm10.2bit"), port=self.PORT, output_dir=str(tmp_path))
        monkeypatch.setattr(mm10, "_server_alive", lambda: False)
        with pytest.raises(RuntimeError, match="different --blat-port"):
            mm10.start_server()
        blat_inst.start_server()  # same reference -> wait for it
        assert blat_inst.is_start_server is False

        # another user's stopped gfServer on our port: AccessDenied on kill must not crash is_running()
        foreign = _FakeGfProc(self.PORT, status="stopped")
        foreign.deny_kill = True
        monkeypatch.setattr("psutil.process_iter", lambda *_a, **_kw: [foreign])
        assert blat_inst.is_running() is False

        # restarting truncates the log so a stale "Server ready" from a dead server is not trusted
        Path(blat_inst.log_file_path).write_text("Server ready for queries!\n")
        monkeypatch.setattr(subprocess, "Popen", lambda *_a, **_kw: _FakePopen())
        blat_inst._start_server()
        assert Path(blat_inst.log_file_path).read_text() == ""

        # a failed query does not leave its temp FASTA behind
        def _refused(*_a: object, **_kw: object) -> int:
            raise subprocess.CalledProcessError(255, "gfClient")

        monkeypatch.setattr(subprocess, "check_call", _refused)
        with pytest.raises(subprocess.CalledProcessError):
            blat_inst._query("A" * 40)
        assert list(tmp_path.glob("*.fasta")) == []

    def test_stop_server_only_stops_our_port(self, blat_inst: Blat, monkeypatch: pytest.MonkeyPatch) -> None:
        ours, other = _FakeGfProc(self.PORT), _FakeGfProc(55752)
        monkeypatch.setattr("psutil.process_iter", lambda *_a, **_kw: [ours, other])
        run_calls: list[list[str]] = []
        monkeypatch.setattr(subprocess, "run", lambda cmd, **_kw: run_calls.append(cmd) or subprocess.CompletedProcess(cmd, 0))

        # reused server (not ours) -> untouched
        blat_inst.is_start_server = False
        blat_inst.stop_server()
        assert run_calls == []
        assert not ours.killed

        # our server exited early (lost a bind race) -> whatever is on the port is not ours
        blat_inst.is_start_server, blat_inst.handle_process = True, _FakePopen(returncode=1)
        blat_inst.stop_server()
        assert run_calls == []
        assert not ours.killed
        assert blat_inst.is_stop_server is True

        # our live server -> signal only our own process group; nothing else on the port is touched,
        # because a concurrent starter may own the process that actually listens there
        signals: list[tuple[int, int]] = []
        monkeypatch.setattr(os, "killpg", lambda pid, sig: signals.append((pid, sig)))
        blat_inst.is_start_server, blat_inst.handle_process = True, _FakePopen()
        blat_inst.stop_server()
        assert signals == [(4242, signal.SIGTERM)]
        assert run_calls == []
        assert not ours.killed
        assert not other.killed

        # flag-only owner (no handle) -> gfServer stop on our port, kill leftovers on our port only
        blat_inst.is_start_server, blat_inst.handle_process = True, None
        blat_inst.stop_server()
        assert run_calls == [[str(blat_inst.gfserver), "stop", "localhost", str(self.PORT)]]
        assert ours.killed is True
        assert other.killed is False

    def test_check_if_self_ready(self, blat_inst: Blat, monkeypatch: pytest.MonkeyPatch) -> None:
        blat_inst.is_start_server = True
        monkeypatch.setattr(blat_inst, "_is_server_ready", lambda: False)

        # our server died before it was ready -> stop waiting and drop ownership (no recursion, no hang)
        blat_inst.handle_process = _FakePopen(returncode=1)
        blat_inst._check_if_self_ready(interval=0)
        assert blat_inst.is_start_server is False

        # still alive but never ready -> bounded wait
        blat_inst.is_start_server, blat_inst.handle_process, blat_inst.max_wait = True, _FakePopen(), -1
        with pytest.raises(RuntimeError, match="not ready"):
            blat_inst._check_if_self_ready(interval=0)

        # ready -> returns and is remembered
        monkeypatch.setattr(blat_inst, "_is_server_ready", lambda: True)
        blat_inst._check_if_self_ready(interval=0)
        assert blat_inst._confirmed_ready is True

    def test_query_restarts_missing_server_then_succeeds(self, blat_inst: Blat, monkeypatch: pytest.MonkeyPatch) -> None:
        starts: list[str] = []
        monkeypatch.setattr(blat_inst, "is_running", lambda: False)
        monkeypatch.setattr(blat_inst, "start_server", lambda: starts.append("start"))
        attempts = {"n": 0}

        def _check_call(*_a: object, **_kw: object) -> int:
            attempts["n"] += 1
            if attempts["n"] == 1:
                raise subprocess.CalledProcessError(1, "gfClient")
            return 0

        monkeypatch.setattr(subprocess, "check_call", _check_call)
        assert blat_inst.query("A" * 40).endswith(".pslx")
        assert starts == ["start"]

    def test_query_gives_up_instead_of_looping_forever(self, blat_inst: Blat, monkeypatch: pytest.MonkeyPatch) -> None:
        def _refused(*_a: object, **_kw: object) -> int:
            raise subprocess.CalledProcessError(255, "gfClient")

        monkeypatch.setattr(subprocess, "check_call", _refused)
        monkeypatch.setattr(blat_inst, "is_running", lambda: True)
        blat_inst.max_wait = -1
        with pytest.raises(RuntimeError, match=f"localhost:{self.PORT} did not answer"):
            blat_inst.query("A" * 40)

    def test_pool_worker_that_starts_server_can_exit(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        # the DU145 hang: a ProcessPoolExecutor worker that starts gfServer must not block shutdown
        import multiprocessing
        import signal
        import threading
        from concurrent import futures

        fake_server = tmp_path / "gfServer"
        fake_server.write_text("#!/bin/sh\nexec sleep 600\n")
        fake_server.chmod(0o755)
        monkeypatch.setattr("scanmst.base.blat.load_gfserver", lambda: fake_server)
        monkeypatch.setattr("scanmst.base.blat.load_gfclient", lambda: tmp_path / "gfClient")

        executor = futures.ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("fork"))
        pid = executor.submit(_start_server_in_worker, str(tmp_path / "hg38.2bit"), self.PORT, str(tmp_path)).result(timeout=30)
        try:
            shutdown = threading.Thread(target=executor.shutdown, kwargs={"wait": True}, daemon=True)
            shutdown.start()
            shutdown.join(timeout=30)
            assert not shutdown.is_alive(), "pool worker hung waiting for the gfServer it started"
        finally:
            with contextlib.suppress(ProcessLookupError):
                os.kill(pid, signal.SIGKILL)
