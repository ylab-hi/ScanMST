"""Unit tests for base primitives, intervals, logging, parallel workers, and exceptions."""

from __future__ import annotations

import copy
from typing import Any

import pytest
from loguru import logger as loguru_logger

from scanmst.base import (
    AnnotationCode,
    BreakPoint,
    CigarCode,
    Event,
    Exon,
    Exons,
    Insertion,
    Interval,
    Intervals,
    Introns,
    MappingMode,
    MicroHomology,
    MyLogger,
    NovelInsertion,
    ParallelWorker,
    Read,
    Strand,
    reverse_complement,
)
from scanmst.exception import (
    AnnotationCodeNotFoundError,
    BreakpointNotFoundError,
    ExonsNotFoundError,
    GenesNotFoundError,
    ModesNotEqualError,
    ModesNotFoundError,
    ReadNotConnectedError,
    ReadNotFoundError,
    ScanmstExceptionError,
    SeqNotFoundError,
    SplicingCodeNotFoundError,
    ToolNotFoundError,
)
from scanmst.graph import Node, NodeIdentity


def _parallel_add_worker(key: str, _lock: Any, *, offset: int = 0) -> int:
    """Top-level picklable worker for ParallelWorker.run tests."""
    return int(key) + offset


def _parallel_square(val: int) -> int:
    """Top-level picklable worker for ParallelWorker.map tests."""
    return val * val


class TestEnums:
    """Test MappingMode, AnnotationCode, Strand, and CigarCode."""

    def test_mapping_mode(self) -> None:
        assert MappingMode.from_int(0) is MappingMode.Type0
        assert MappingMode.from_int(1) is MappingMode.MS
        assert MappingMode.from_int(2) is MappingMode.SM
        assert MappingMode.from_int(MappingMode.SM) is MappingMode.SM
        with pytest.raises(ValueError, match="Invalid mode"):
            MappingMode.from_int(99)

        assert MappingMode.MS.reversed() is MappingMode.SM
        assert MappingMode.SM.reversed() is MappingMode.MS
        assert MappingMode.SM.is_sm() is True
        assert MappingMode.SM.is_ms() is False
        assert MappingMode.MS.is_ms() is True
        assert MappingMode.SM.to_str() == "SM"
        assert MappingMode.MS.to_str() == "MS"

    def test_annotation_code_and_cigar_code(self) -> None:
        assert AnnotationCode.Type1.reversed() is AnnotationCode.Type2
        assert AnnotationCode.Type2.reversed() is AnnotationCode.Type1
        assert CigarCode.Match == 0
        assert CigarCode.Ref_skip == 3
        assert CigarCode.Soft_clip == 4
        assert CigarCode.Back == 9

    def test_strand(self) -> None:
        assert Strand.from_str("+") is Strand.Forward
        assert Strand.from_str("-") is Strand.Reverse
        assert Strand.from_str(Strand.Reverse) is Strand.Reverse
        with pytest.raises(ValueError, match="Invalid strand"):
            Strand.from_str("?")

        assert Strand.Forward.is_forward() is True
        assert Strand.Forward.is_reverse() is False
        assert Strand.Forward.reversed() is Strand.Reverse
        assert Strand.Reverse.reversed() is Strand.Forward
        assert str(Strand.Forward) == "+"
        assert str(Strand.Reverse) == "-"


class TestIntervalAndIntervals:
    """Test Interval and Intervals (Exons/Introns) operations."""

    def test_interval_basic_and_arithmetic(self) -> None:
        assert Exon is Interval
        assert Exons is Intervals
        assert Introns is Intervals

        with pytest.raises(ValueError, match="start: 10 > end: 5"):
            Interval(10, 5)

        iv = Interval(10, 20)
        assert len(iv) == 10
        assert repr(iv) == "[10, 20)"
        assert str(iv) == "[10, 20)"
        assert tuple(iv) == (10, 20)
        assert iv[0] == 10
        assert iv[1] == 20
        with pytest.raises(ValueError, match="is out"):
            _ = iv[2]
        with pytest.raises(IndexError, match="is out"):
            _ = iv[-1]

        assert iv == Interval.from_list([10, 20])
        assert iv != (10, 20)
        assert hash(iv) == hash(Interval(10, 20))

        assert iv + 5 == Interval(15, 25)
        assert iv + Interval(1, 2) == Interval(11, 22)
        with pytest.raises(TypeError, match="not Interval or int"):
            _ = iv + "x"  # type: ignore[operator]

        assert iv - 5 == Interval(5, 15)
        assert iv - Interval(1, 2) == Interval(9, 18)
        with pytest.raises(TypeError, match="not Interval or int"):
            _ = iv - "x"  # type: ignore[operator]

    def test_interval_contain_and_join(self) -> None:
        outer = Interval(10, 30)
        inner = Interval(15, 25)
        assert inner in outer
        assert outer not in inner
        assert outer.contain(inner) is True
        assert outer.contain(outer, same_left=True, same_right=True) is True
        assert outer.contain(inner, same_left=True, same_right=True) is False
        assert outer.contain(Interval(20, 30), same_left=True) is True
        assert outer.contain(Interval(10, 25), same_right=True) is True
        with pytest.raises(ValueError, match="not Interval"):
            outer.contain((15, 25))  # type: ignore[arg-type]

        # s |o o| s
        assert Interval(10, 30).join(Interval(15, 25)) == (Interval(15, 25), Interval(10, 30))
        # s |o s| o
        assert Interval(10, 20).join(Interval(15, 25)) == (Interval(15, 20), Interval(10, 25))
        # o |s s| o
        assert Interval(15, 25).join(Interval(10, 30)) == (Interval(15, 25), Interval(10, 30))
        # o |s o| s
        assert Interval(15, 25).join(Interval(10, 20)) == (Interval(15, 20), Interval(10, 25))
        # adjacent: self.start == other.end and self.end == other.start
        assert Interval(20, 30).join(Interval(10, 20)) == (Interval(20, 20), Interval(10, 30))
        assert Interval(10, 20).join(Interval(20, 30)) == (Interval(20, 20), Interval(10, 30))
        # disjoint
        assert Interval(10, 20).join(Interval(25, 30)) is None
        with pytest.raises(ValueError, match="not a Interval"):
            Interval(10, 20).join((15, 25))  # type: ignore[arg-type]

    def test_intervals_collection_methods(self) -> None:
        ivs = Intervals.from_list([[20, 30], [0, 10]])
        ivs.sort()
        assert ivs == Intervals([Interval(0, 10), Interval(20, 30)])
        assert ivs != [(0, 10), (20, 30)]
        assert hash(ivs) == hash((Interval(0, 10), Interval(20, 30)))
        assert len(ivs) == 2
        assert ivs.first == Interval(0, 10)
        assert ivs.last == Interval(20, 30)
        assert ivs[0] == Interval(0, 10)
        assert ivs[0:1] == [Interval(0, 10)]
        assert Interval(0, 10) in ivs
        assert Interval(0, 5) not in ivs
        assert str(ivs) == "[0-10,20-30]"
        assert "Exons(" in repr(ivs)
        assert ivs.reversed() == Intervals([Interval(20, 30), Interval(0, 10)])

        # Arithmetic and concatenation
        assert ivs + 5 == Intervals([Interval(5, 15), Interval(25, 35)])
        assert ivs + Intervals([Interval(1, 1), Interval(2, 2)]) == Intervals(
            [Interval(1, 11), Interval(22, 32)],
        )
        with pytest.raises(TypeError, match="not int or Intervals"):
            _ = ivs + "bad"  # type: ignore[operator]

        assert ivs - 2 == Intervals([Interval(-2, 8), Interval(18, 28)])
        assert ivs - Intervals([Interval(1, 1), Interval(2, 2)]) == Intervals(
            [Interval(-1, 9), Interval(18, 28)],
        )
        with pytest.raises(TypeError, match="not int or Intervals"):
            _ = ivs - "bad"  # type: ignore[operator]

        concat = ivs.__concat__(Intervals([Interval(40, 50)]))
        assert len(concat) == 3
        with pytest.raises(TypeError, match="Can only concatenate"):
            ivs.__concat__([Interval(40, 50)])  # type: ignore[arg-type]

        # Mutation: setitem, append, extend, introns, index
        mut = Intervals([])
        assert len(mut.introns()) == 0
        mut.append((0, 10))
        assert len(mut.introns()) == 0
        mut.append(Interval(20, 30))
        with pytest.raises(TypeError, match="not Interval or tuple"):
            mut.append([40, 50])  # type: ignore[arg-type]

        mut.extend(Intervals([Interval(40, 50), Interval(0, 10)]))
        with pytest.raises(TypeError, match="not Intervals"):
            mut.extend([Interval(60, 70)])  # type: ignore[arg-type]

        assert len(Intervals.from_list([[0, 10], [20, 30], [45, 55]]).introns()) == 2

        assert mut.index(Interval(0, 10)) == 0
        assert mut.index(Interval(0, 10), 1) == 3
        assert mut.index(Interval(20, 30), -3, -1) == 1
        with pytest.raises(ValueError, match="not in exon_list"):
            mut.index(Interval(99, 100))
        with pytest.raises(TypeError, match="not an Interval"):
            mut.index((0, 10))  # type: ignore[arg-type]

        mut[0] = Interval(2, 10)
        assert mut[0] == Interval(2, 10)


class TestBasicClasses:
    """Test BreakPoint, Insertion/NovelInsertion/MicroHomology, Event, and reverse_complement."""

    def test_breakpoint_methods(self) -> None:
        bp = BreakPoint.from_str("chr1:1000", depth=7)
        assert str(bp) == "chr1:1000"
        assert bp.to_tuple() == ("chr1", 1000)
        assert bp.depth == 7
        with pytest.raises(ValueError, match="empty string"):
            BreakPoint.from_str("")

        assert bp.equal(BreakPoint("chr1", 1005), threshold=5) is True
        assert bp.equal(BreakPoint("chr1", 1006), threshold=5) is False
        assert bp.equal(BreakPoint("chr2", 1000), threshold=5) is False

        fwd_node = Node("q1", "chr1", "+", 100, 200, NodeIdentity.HEAD, Intervals.from_list([[100, 200]]))
        rev_node = Node("q2", "chr1", "-", 100, 200, NodeIdentity.HEAD, Intervals.from_list([[100, 200]]))
        assert BreakPoint.from_node(fwd_node) == (BreakPoint("chr1", 100), BreakPoint("chr1", 200))
        assert BreakPoint.from_node(rev_node) == (BreakPoint("chr1", 200), BreakPoint("chr1", 100))

    def test_insertions_and_microhomology(self) -> None:
        assert reverse_complement("ACGT") == "ACGT"
        assert reverse_complement("AACG") == "CGTT"

        ni1 = NovelInsertion(hit_num=2, query_sequence="ACGT")
        ni2 = NovelInsertion(hit_num=0, query_sequence="ACGT")
        assert repr(ni1) == "NovelInsertion(ACGT:2)"
        assert ni1 == ni2
        assert ni1 != "ACGT"

        mh1 = MicroHomology("GGCC")
        mh2 = MicroHomology("GGCC")
        mh1.increment_ao(3)
        assert mh1.ao == 4
        assert repr(mh1) == "MicroHomology(GGCC)"
        assert mh1 == mh2
        assert mh1 != ni1

        ins = Insertion(
            hit_num=1,
            chrom="chr1",
            ref_start=100,
            strand="+",
            cigarstring="1S2M1S",
            mapq=60,
            nm=0,
            query_sequence="AGCT",
            query_qualities=[30, 30, 30, 30],
        )
        assert "Insertion(" in repr(ins)
        # source_s="right" and opposite strand swaps ls/rs
        ins.update_cigarstring_sms(sms=(5, 10, 8), source_s="right", source_strand=Strand.Reverse)
        assert ins.sms == (4, 4, 15)
        assert ins.cigarstring == "5S2M16S"

    def test_event_branches(self, event: Event, reads: list[Read]) -> None:
        ev = copy.deepcopy(event)
        ev.annotation_code = AnnotationCode.Type1
        ev.reverse()
        assert ev.annotation_code == AnnotationCode.Type2
        ev.reverse()
        assert ev.annotation_code == AnnotationCode.Type1
        ev.merge(event)
        assert "Event(" in repr(ev)
        assert ev.chrom2 == "chr17"
        assert ev.junction_pos1 == 7701655
        assert ev.junction_pos2 == 7708249
        assert ev.insertion_seq2 == "GG"
        assert ev.insertion_microhomology_len == 2

        ins_ev = copy.deepcopy(event)
        ins_ev.insertion_info = ("+ACGT", "+ACGT")
        assert ins_ev.has_insertion() is True
        assert ins_ev.has_microhomology() is False
        assert ins_ev.insertion_microhomology_len == 4

        blunt_ev = copy.deepcopy(event)
        blunt_ev.insertion_info = ("", "")
        assert blunt_ev.insertion_microhomology_len == 0

        na_ev = copy.deepcopy(event)
        na_ev.sv_type = "NA"
        assert na_ev.is_type_na() is True

        with pytest.raises(ReadNotFoundError):
            ev.read2([], 0)


class TestLoggerParallelAndExceptions:
    """Test MyLogger, ParallelWorker, and custom exception hierarchy."""

    def test_my_logger_prefixes_contig(self) -> None:
        records: list[str] = []
        sink_id = loguru_logger.add(records.append, level="TRACE", format="{level.name}:{message}")
        try:
            wrapped = MyLogger("chr17", loguru_logger)  # type: ignore[arg-type]
            wrapped.trace("t")
            wrapped.debug("d")
            wrapped.info("i")
            wrapped.warning("w")
            wrapped.error("e")
            wrapped.critical("c")
            wrapped.success("s")
            wrapped.complete()
        finally:
            loguru_logger.remove(sink_id)

        joined = "\n".join(records)
        for level, char in [
            ("TRACE", "t"),
            ("DEBUG", "d"),
            ("INFO", "i"),
            ("WARNING", "w"),
            ("ERROR", "e"),
            ("CRITICAL", "c"),
            ("SUCCESS", "s"),
        ]:
            assert f"{level}:chr17: {char}" in joined

    def test_parallel_worker_run_map_and_cap(self, fake_logger) -> None:
        worker = ParallelWorker(_parallel_add_worker, fake_logger, n_jobs=99999)
        assert worker.n_jobs >= 1

        worker2 = ParallelWorker(_parallel_add_worker, fake_logger, n_jobs=1)
        assert worker2.n_jobs == 1
        out = worker2.run("1", "2", "3", offset=10)
        assert out == {"1": 11, "2": 12, "3": 13}

        mapper = ParallelWorker(_parallel_square, fake_logger, n_jobs=1)
        assert list(mapper.map([2, 3, 5])) == [4, 9, 25]

    def test_all_custom_exceptions(self) -> None:
        instances = [
            ToolNotFoundError("samtools"),
            ReadNotFoundError(),
            ReadNotConnectedError(),
            SeqNotFoundError(),
            ExonsNotFoundError("n1"),
            ModesNotFoundError("n1"),
            ModesNotEqualError("r1/r2"),
            GenesNotFoundError("n1"),
            SplicingCodeNotFoundError("n1"),
            AnnotationCodeNotFoundError("n1"),
            BreakpointNotFoundError("n1"),
        ]
        for exc in instances:
            assert isinstance(exc, ScanmstExceptionError)
            assert len(str(exc)) > 5
