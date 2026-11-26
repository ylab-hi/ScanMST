"""Tests for the compiled C++ extension (`scanmst.cppext`) and `SRRescuer`."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pytest

from scanmst import cppext as cpp
from scanmst.base import BreakPoint as PyBreakPoint
from scanmst.base import Intervals, MappingMode
from scanmst.graph import EdgeData, NLGraph, Node, NodeIdentity, SRRescuer, VariationType
from scanmst.graph.sr_rescuer import is_middle_node, make_breakpoint

DATA_DIR = Path(__file__).parent / "data"
SR_BAM = (DATA_DIR / "test_calculate_sr.bam").as_posix()


class TestCigarAndStructs:
    """Test CigarResult, parseCigar, Region, Options, BreakPoint, and Seqs bindings."""

    def test_parse_cigar_softclip_flags_and_fields(self) -> None:
        both = cpp.parseCigar("10S50M5S")
        assert both.both_soft_clipped() is True
        assert both.is_left_soft_clipped() is True
        assert both.is_right_soft_clipped() is True
        assert both.lt_soft_len == 10
        assert both.rt_soft_len == 5
        assert both.read_match == 50
        assert both.ref_match == 50
        assert both.query_len == 65
        assert "lt_soft_len: 10" in both.to_string()
        assert "rt_soft_len: 5" in str(both)

        left_only = cpp.parseCigar("8S40M")
        assert left_only.both_soft_clipped() is False
        assert left_only.is_left_soft_clipped() is True
        assert left_only.is_right_soft_clipped() is False

        right_only = cpp.parseCigar("40M8S")
        assert right_only.both_soft_clipped() is False
        assert right_only.is_left_soft_clipped() is False
        assert right_only.is_right_soft_clipped() is True

        complex_cigar = cpp.parseCigar("10S20M5I3D100N30M")
        assert complex_cigar.read_match == 55  # 20M + 5I + 30M
        assert complex_cigar.ref_match == 153  # 20M + 3D + 100N + 30M
        assert complex_cigar.indel_len == 98  # 3 + 100 - 5
        assert complex_cigar.query_len == 65

        copied = cpp.CigarResult(both)
        assert copied.query_len == 65
        empty = cpp.CigarResult()
        empty.assign(both)
        assert empty.lt_soft_len == 10

    def test_region_and_options_and_breakpoint_and_seqs(self) -> None:
        reg_default = cpp.Region()
        assert reg_default.chrom == ""
        reg = cpp.Region("chr17", 100, 200)
        assert (reg.chrom, reg.start, reg.end) == ("chr17", 100, 200)
        assert "chr17:100-200" in reg.to_string()
        reg_copy = cpp.Region(reg)
        reg_copy.start = 150
        assert reg_copy.start == 150
        reg_default.assign(reg_copy)
        assert reg_default.start == 150

        opts = cpp.Options().file(SR_BAM).mapq(20).soft_len(6).mismatch(2).identity(0.85).min_seq_align_len(12).average_read_depth(30)
        assert opts.file_ == SR_BAM
        assert opts.mapq_ == 20
        assert opts.soft_len_ == 6
        assert opts.mismatch_ == 2
        assert opts.identity_ == pytest.approx(0.85)
        assert opts.min_seq_align_len_ == 12
        assert opts.average_read_depth_ == 30
        assert "mapq: 20" in opts.to_string()
        opts_copy = cpp.Options(opts)
        assert opts_copy.mapq_ == 20

        bp = cpp.BreakPoint(7705201, 7708250, 1, False, True)
        assert (bp.read_start, bp.read_end, bp.mode, bp.is_reverse, bp.is_middle) == (
            7705201,
            7708250,
            1,
            False,
            True,
        )
        assert "mode: 1" in bp.to_string()
        bp_default = cpp.BreakPoint()
        bp_default.assign(bp)
        assert bp_default.is_middle is True

        s1 = cpp.Seqs("ACGTACGTACGT")
        assert s1.seq1 == "ACGTACGTACGT"
        assert s1.seq2 is None
        s2 = cpp.Seqs("ACGTACGTACGT", "TGCATGCATGCA")
        assert s2.seq2 == "TGCATGCATGCA"
        assert "ACGTACGTACGT" in s2.to_string()
        assert "TGCATGCATGCA" in str(s2)


class TestBamReader:
    """Test BamReader and BamReader.Iterator against tests/data/test_calculate_sr.bam."""

    def test_count_and_query_overloads(self) -> None:
        reader = cpp.BamReader(SR_BAM)
        region = cpp.Region("chr17", 7708157, 7708250)
        assert reader.count("chr17", 7708157, 7708250) == 5
        assert reader.count(region) == 5
        assert reader.count("chr17", 1, 10) == 0

        records = []
        it = reader.query(cpp.Region("chr17", 7701655, 7702552))
        while not it.is_end():
            records.append(
                (
                    it.read_name(),
                    it.chrom(),
                    it.pos(),
                    it.end_pos(),
                    it.quality(),
                    len(it.sequence()),
                    it.cigar_string(),
                    it.cigar_length(),
                    it.is_reverse(),
                    it.same_strand_with(False),
                    it.quality_eq_than(60),
                    it.quality_eq_than(61),
                    bool(it.to_string()),
                ),
            )
            it.next()

        assert len(records) == 2
        names = [r[0] for r in records]
        assert names == ["one_hop_back", "CCS"]
        for rec in records:
            assert rec[1] == "chr17"
            assert rec[2] == 7701655
            assert rec[3] == 7702552
            assert rec[4] == 60
            assert rec[8] is False
            assert rec[9] is True
            assert rec[10] is True
            assert rec[11] is False
            assert rec[12] is True

        empty_it = reader.query("chr17", 1, 10)
        assert empty_it.is_end() is True


class TestRescuerAndSRRescuer:
    """Test C++ Rescuer alignment/rescue logic and Python SRRescuer integration."""

    @pytest.fixture
    def rescuer(self) -> cpp.Rescuer:
        opts = cpp.Options().file(SR_BAM).mapq(15).soft_len(5).mismatch(3).identity(0.8).min_seq_align_len(10)
        return cpp.Rescuer(opts)

    def test_static_helpers(self) -> None:
        assert cpp.Rescuer.get_align_seq_len(5, 10) is None
        assert cpp.Rescuer.get_align_seq_len(10, 10) == 10
        assert cpp.Rescuer.get_align_seq_len(42, 10) == 42
        assert cpp.Rescuer.get_align_seq_len(250, 10) == 100

        assert cpp.Rescuer.check_if_same_isform([0, 100], [0, 250]) is True
        assert cpp.Rescuer.check_if_same_isform([0, 134, 3, 554, 0, 209], [0, 50, 3, 554, 0, 80]) is True
        assert cpp.Rescuer.check_if_same_isform([0, 134, 3, 554, 0, 209], [0, 134, 3, 600, 0, 209]) is False
        assert cpp.Rescuer.check_if_same_isform([0, 100], [0, 50, 3, 200, 0, 50]) is False

    def test_identity_and_check_rescue(self, rescuer: cpp.Rescuer) -> None:
        seq_a = "ACGTACGTACGTACGTACGT"
        seq_b = "ACGTACGTACGTACGTACGT"
        seq_diff = "TTTTTTTTTTTTTTTTTTTT"

        ident = rescuer.calculate_identity(seq_a, seq_b)
        assert ident == pytest.approx(1.0)
        assert rescuer.check_identity(seq_a, seq_b, 0.9) is True
        assert rescuer.check_rescue_condition(seq_a, seq_b) is True
        assert rescuer.check_rescue_condition(seq_a, seq_diff) is False

        # 1-sequence end-node rescue check
        assert rescuer.check_rescue(cpp.Seqs(seq_a), cpp.Seqs(seq_b)) is True
        assert rescuer.check_rescue(cpp.Seqs(seq_a), cpp.Seqs(seq_diff)) is False

        # 2-sequence middle-node rescue check
        assert rescuer.check_rescue(cpp.Seqs(seq_a, seq_b), cpp.Seqs(seq_a, seq_b)) is True
        assert rescuer.check_rescue(cpp.Seqs(seq_a, seq_diff), cpp.Seqs(seq_a, seq_b)) is False

        assert rescuer.count_reads("chr17", 7708249, 7708250) == 5
        assert "min_mapq: 15" in rescuer.to_string()

    def test_calculate_sr_and_average_depth_cap(self, rescuer: cpp.Rescuer) -> None:
        # QQ, QQ1, QQ2 all align at chr17:7708157-7708250 with 93M390S (mode 1).
        # Using QQ as the node read rescues QQ1 and QQ2 (sr == 2).
        rescuer.reset_names_list(["QQ"])
        region = cpp.Region("chr17", 7708249, 7708250)
        bp = cpp.BreakPoint(7708157, 7708250, 1, False, False)
        assert rescuer.calculate_sr(region, bp, ["QQ"], [0, 93]) == 2

        # Once rescued, QQ1 and QQ2 are recorded in m_names_list so they are not double-counted
        assert rescuer.calculate_sr(region, bp, ["QQ"], [0, 93]) == 0

        # Capped rescuer stops early at average_read_depth=1
        capped_opts = cpp.Options().file(SR_BAM).mapq(15).soft_len(5).mismatch(3).identity(0.8).min_seq_align_len(10).average_read_depth(1)
        capped = cpp.Rescuer(capped_opts)
        capped.reset_names_list(["QQ"])
        assert capped.calculate_sr(region, bp, ["QQ"], [0, 93]) == 1

        # Mode 2 (left soft-clip) at chr17:7701655-7702552 (one_hop_back vs CCS)
        rescuer.reset_names_list(["one_hop_back"])
        reg_lt = cpp.Region("chr17", 7701655, 7701656)
        bp_lt = cpp.BreakPoint(7701655, 7702552, 2, False, False)
        sr_lt = rescuer.calculate_sr(reg_lt, bp_lt, ["one_hop_back"], [0, 134, 3, 554, 0, 209])
        assert sr_lt >= 0

    def test_python_sr_rescuer_on_nlgraph(self, fake_logger, tmp_path: Path) -> None:
        n_head = Node(
            query_name="QQ",
            chrom="chr17",
            strand="+",
            ref_start=7708157,
            ref_end=7708250,
            identity=NodeIdentity.HEAD,
            exons=Intervals.from_list([[7708157, 7708250]]),
            cigartuples_without_soft=[0, 93],
        )
        n_mid = Node(
            query_name="QQ",
            chrom="chr17",
            strand="-",
            ref_start=7701655,
            ref_end=7702552,
            identity=NodeIdentity.MID,
            exons=Intervals.from_list([[7701655, 7701789], [7702343, 7702552]]),
            cigartuples_without_soft=[0, 134, 3, 554, 0, 209],
        )
        n_tail = Node(
            query_name="one_hop_back",
            chrom="chr17",
            strand="+",
            ref_start=7701655,
            ref_end=7702552,
            identity=NodeIdentity.TAIL,
            exons=Intervals.from_list([[7701655, 7701789], [7702343, 7702552]]),
            cigartuples_without_soft=[0, 134, 3, 554, 0, 209],
        )

        assert is_middle_node(n_head) is False
        assert is_middle_node(n_mid) is True
        bp_cpp = make_breakpoint(n_mid, int(MappingMode.SM))
        assert bp_cpp.is_reverse is True
        assert bp_cpp.is_middle is True

        assert SRRescuer.obtain_region_for_rescue_sr(n_head, MappingMode.MS) == ("chr17", 7708250)
        assert SRRescuer.obtain_region_for_rescue_sr(n_head, MappingMode.SM) == ("chr17", 7708251)
        assert SRRescuer.obtain_region_for_rescue_sr(n_mid, MappingMode.MS) == ("chr17", 7701655)

        sr_rescuer = SRRescuer(
            input_bam_file=SR_BAM,
            mapq_cutoff=15,
            soft_len_cutoff=5,
            mismatch_cutoff=3,
            alignment_frac=0.8,
            node_rescued_sr_maximum=100,
            average_read_depth=50,
        )
        graph = NLGraph(
            logger=fake_logger,
            rescuer=sr_rescuer,
            merge_threshold=10,
            support_reads=1,
            input_bam_path=Path(SR_BAM),
            output_dir=tmp_path,
            rescue_sr=True,
        )
        graph.nodes = graph.dict_factory()
        graph.edges = defaultdict(list)
        n_head.is_in_graph = True
        n_tail.is_in_graph = True
        graph.add_node_with_similar_key(n_head)
        graph.add_node_with_similar_key(n_tail)

        edata = EdgeData(
            variantion_type=VariationType.TDUP,
            break_point1=PyBreakPoint("chr17", 7708250),
            break_point2=PyBreakPoint("chr17", 7701655),
            sr=1,
            read_ids=["QQ", "one_hop_back"],
            mode1=MappingMode.MS,
            mode2=MappingMode.SM,
            gene1="G1",
            gene2="G1",
            annotation_code=3,
            splicing_code=1,
        )
        n_head.add_successor(n_tail, graph, edata)
        edge = graph.find_edges(n_head, n_tail)[0]

        sr_rescuer.init(graph)
        assert edge.break_point1.depth == 5
        assert edge.break_point2.depth == 2

        sr_rescuer(graph)
        assert edge.original_sr == 1
        assert edge.sr >= 3  # rescued QQ1 + QQ2 at breakpoint 1

        # Missing cigartuples_without_soft raises ValueError
        n_tail.cigartuples_without_soft = None
        with pytest.raises(ValueError, match="None value"):
            sr_rescuer.update_breakpoint_depth(graph, n_head)
