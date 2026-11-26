"""Tests for `scanmst.blat` platform/download helpers and `scanmst.base.blat.Blat`."""

from __future__ import annotations

import os
import platform
import subprocess
from pathlib import Path
from types import SimpleNamespace

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

        # _run_cmd normal, CalledProcessError, and KeyboardInterrupt -> stop_server
        monkeypatch.setattr(subprocess, "check_call", lambda *_a, **_kw: 0)
        blat_inst._run_cmd("true")

        def _raise_cpe(*_a: object, **_kw: object) -> int:
            raise subprocess.CalledProcessError(1, "cmd")

        monkeypatch.setattr(subprocess, "check_call", _raise_cpe)
        blat_inst._run_cmd("false")

        stopped_procs: list[str] = []

        class _FakeProc:
            def start(self) -> None:
                stopped_procs.append("started")

            def kill(self) -> None:
                stopped_procs.append("killed")

        monkeypatch.setattr("scanmst.base.blat.Process", lambda **_kw: _FakeProc())
        monkeypatch.setattr(blat_inst, "_search_processing", lambda: [_FakeProc()])
        blat_inst.is_start_server = True

        def _raise_kbi(*_a: object, **_kw: object) -> int:
            raise KeyboardInterrupt

        monkeypatch.setattr(subprocess, "check_call", _raise_kbi)
        blat_inst._run_cmd("interrupt")
        assert blat_inst.is_stop_server is True
        assert "killed" in stopped_procs

        # start_server when already running vs not running
        monkeypatch.setattr(blat_inst, "is_running", lambda: True)
        blat_inst.start_server()
        assert blat_inst.is_start_server is False

        monkeypatch.setattr(blat_inst, "is_running", lambda: False)
        blat_inst.start_server()
        assert blat_inst.is_start_server is True
        assert "started" in stopped_procs

        # _check_if_self_ready: wait once when log missing, then succeed once Server ready is written
        ready_checks = {"n": 0}

        def _fake_ready() -> bool:
            ready_checks["n"] += 1
            return True

        def _write_log_on_sleep(_s: float) -> None:
            log_path.write_text("Server ready\n")

        monkeypatch.setattr("time.sleep", _write_log_on_sleep)
        monkeypatch.setattr(blat_inst, "is_ready", _fake_ready)
        blat_inst._check_if_self_ready(interval=0)
        assert ready_checks["n"] == 1

        # query(): first not running -> starts server -> CalledProcessError once -> succeeds via _query
        blat_inst.is_start_server = False
        call_state = {"q": 0}

        def _fake_check_call(*_a: object, **_kw: object) -> int:
            call_state["q"] += 1
            if call_state["q"] == 1:
                raise subprocess.CalledProcessError(1, "gfClient")
            return 0

        monkeypatch.setattr(subprocess, "check_call", _fake_check_call)
        out_pslx = blat_inst.query("A" * 40)
        assert out_pslx.endswith(".pslx")

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

            def cmdline(self) -> bool:
                return True

            def status(self) -> str:
                return self._status

            def kill(self) -> None:
                killed.append(self._name)

        mocker.patch("psutil.process_iter", return_value=[_Proc("gfServer", "running")])
        assert blat.is_running() is True

        mocker.patch("psutil.process_iter", return_value=[_Proc("gfServer", "stopped")])
        assert blat.is_running() is False
        assert killed == ["gfServer"]
