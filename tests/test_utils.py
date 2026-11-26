"""Tests for scanmst/utils.py."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from scanmst import ToolNotFoundError
from scanmst.base import BreakPoint, MappingMode
from scanmst.graph import Edge, EdgeData, VariationType
from scanmst.utils import (
    change_dir,
    change_dir_decorator,
    cigar_validity,
    cigarstring2cigartuples,
    determine_mst_link_type,
    external_tool_checking,
    find_2bit_file,
    get_current_time,
    get_longest_insertion_sequence,
    get_softclip_length,
    timeit,
    wait_for_aligner,
)


def _make_seg(cigar: str | None, seq: str | None, ref_start: int = 391287) -> pysam.AlignedSegment:
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr20", "LN": 1000000}]})
    seg = pysam.AlignedSegment(header)
    seg.query_name = "read_util"
    seg.reference_id = 0
    seg.reference_start = ref_start
    if cigar is not None:
        seg.cigarstring = cigar
    if seq is not None:
        seg.query_sequence = seq
    return seg


@pytest.mark.parametrize(
    ("mode", "cigar", "seq", "expected"),
    [
        (
            MappingMode.Type0,
            "20S20M10S",
            "TAAGACTAGGAATGGAGCAGTTCAGTCTAAAAAATATCACAGGTAACAGA",
            (20, "TAAGACTAGGAATGGAGCAG", 391287, MappingMode.SM),
        ),
        (
            MappingMode.Type0,
            "10S20M20S",
            "TAAGACTAGGAATGGAGCAGTTCAGTCTAAAAAATATCACAGGTAACAGA",
            (20, "AAAATATCACAGGTAACAGA", 391307, MappingMode.MS),
        ),
        (
            MappingMode.Type0,
            "15S20M15S",
            "A" * 50,
            None,
        ),
        (
            MappingMode.MS,
            "20S20M10S",
            "TAAGACTAGGAATGGAGCAGTTCAGTCTAAAAAATATCACAGGTAACAGA",
            (10, "AGGTAACAGA", 391307, MappingMode.MS),
        ),
        (
            MappingMode.SM,
            "20S20M10S",
            "TAAGACTAGGAATGGAGCAGTTCAGTCTAAAAAATATCACAGGTAACAGA",
            (20, "TAAGACTAGGAATGGAGCAG", 391287, MappingMode.SM),
        ),
    ],
)
def test_get_softclip_length(mode: MappingMode, cigar: str, seq: str, expected: object) -> None:
    assert get_softclip_length(_make_seg(cigar, seq), mode) == expected


def test_get_softclip_length_rejects_missing_seq_or_cigar() -> None:
    with pytest.raises(ValueError, match="query sequence or cigar is None"):
        get_softclip_length(_make_seg(None, None), MappingMode.SM)


def test_get_longest_insertion_sequence() -> None:
    with pytest.raises(ValueError, match="query sequence or cigar is None"):
        get_longest_insertion_sequence(_make_seg(None, None))

    # No insertion -> (0, "", 0)
    assert get_longest_insertion_sequence(_make_seg("50M", "A" * 50, ref_start=100)) == (0, "", 0)

    # Insertion below cutoff -> (0, "", 0)
    assert get_longest_insertion_sequence(_make_seg("10M20I10M", "A" * 40, ref_start=100), insertion_length_cutoff=50) == (
        0,
        "",
        0,
    )

    # Multiple insertions, picks longest; with adjacent 500N shifts ins_ref_pos by +500
    seq = ("A" * 5) + ("C" * 10) + ("G" * 20) + ("T" * 60) + ("A" * 10)
    seg = _make_seg("5S10M20I5D60I500N10M", seq, ref_start=100)
    # Before 60I: current_pos = 100 + 10 (M) + 5 (D) = 115; + 500 (adjacent N) = 615
    assert get_longest_insertion_sequence(seg, insertion_length_cutoff=50) == (615, "T" * 60, 60)


def test_cigar_conversions_and_validity() -> None:
    assert cigarstring2cigartuples("10S30M2I3D100N5H") == [(4, 10), (0, 30), (1, 2), (2, 3), (3, 100), (5, 5)]
    assert cigar_validity("1S2S5M0I5M3S0S2S") == "3S10M5S"


def test_determine_mst_link_type() -> None:
    def _edge(vtype: VariationType, m1: MappingMode, m2: MappingMode) -> Edge:
        return Edge(
            "n1",
            "n2",
            EdgeData(vtype, BreakPoint("chr1", 1), BreakPoint("chr1", 2), 1, ["r"], m1, m2, "G1", "G2", 0, 1),
        )

    assert determine_mst_link_type(_edge(VariationType.TDUP, MappingMode.SM, MappingMode.MS)) == "ICRL"
    assert determine_mst_link_type(_edge(VariationType.INV, MappingMode.MS, MappingMode.MS)) == "ICTL"
    assert determine_mst_link_type(_edge(VariationType.TRA, MappingMode.MS, MappingMode.MS)) == "ITTL"
    assert determine_mst_link_type(_edge(VariationType.TRA, MappingMode.MS, MappingMode.SM)) == "ITPL"
    assert determine_mst_link_type(_edge(VariationType.DEL, MappingMode.MS, MappingMode.SM)) == "NA"


def test_external_tool_checking(monkeypatch: pytest.MonkeyPatch, fake_logger) -> None:
    monkeypatch.setattr(shutil, "which", lambda _x: "/usr/bin/gfServer")
    external_tool_checking(["gfServer"], fake_logger)

    monkeypatch.setattr(shutil, "which", lambda _x: None)
    with pytest.raises(ToolNotFoundError):
        external_tool_checking(["gfServer"], fake_logger)


def test_time_dir_2bit_and_wait_for_aligner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    assert len(get_current_time()) == 15

    @timeit
    def _add(a: int, b: int) -> int:
        return a + b

    assert _add(2, 3) == 5

    orig_cwd = Path.cwd()

    def _enter_and_fail() -> None:
        with change_dir(str(tmp_path)):
            assert Path.cwd() == tmp_path.resolve()
            raise RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        _enter_and_fail()
    assert Path.cwd() == orig_cwd

    @change_dir_decorator(str(tmp_path))
    def _where() -> Path:
        return Path.cwd()

    assert _where() == tmp_path.resolve()
    assert Path.cwd() == orig_cwd

    # find_2bit_file: creates when absent, skips when present
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nACGT\n")
    calls: list[list[object]] = []

    def _fake_check_call(cmd: list[object]) -> int:
        calls.append(cmd)
        (tmp_path / "ref.2bit").write_bytes(b"2bit")
        return 0

    monkeypatch.setattr(subprocess, "check_call", _fake_check_call)
    monkeypatch.setattr("scanmst.utils.load_fa2bit", lambda: Path("/fake/faToTwoBit"))
    out_2bit = find_2bit_file(str(fa))
    assert out_2bit.endswith("ref.2bit")
    assert len(calls) == 1
    assert find_2bit_file(str(fa)) == out_2bit
    assert len(calls) == 1  # already exists -> not called again

    # wait_for_aligner: succeeds once running (note logger.trace also calls is_running()), raises TimeoutError when stuck
    states = iter([False, False, True])
    aligner_ok = SimpleNamespace(
        is_running=lambda: next(states),
        log_file_path=str(tmp_path / "blat.log"),
        is_start_server=True,
    )
    wait_for_aligner(aligner_ok, max_timeout=5, check_interval=0.001)

    aligner_stuck = SimpleNamespace(
        is_running=lambda: False,
        log_file_path=str(tmp_path / "blat.log"),
        is_start_server=True,
    )
    with pytest.raises(TimeoutError, match="did not start running"):
        wait_for_aligner(aligner_stuck, max_timeout=0.005, check_interval=0.002)
