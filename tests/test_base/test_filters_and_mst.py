"""Tests for ExonFilter, RTSwitchingFilter, CircRNAFilter, Read.splice_site_checker,
mst_inference, and helper alignment/chimeric functions."""

from __future__ import annotations

from pathlib import Path

import HTSeq
import pysam
import pytest
from pyfaidx import Fasta

from scanmst.base import (
    Blat,
    BreakPoint,
    CircRNAFilter,
    Event,
    ExonFilter,
    Exons,
    Insertion,
    Intervals,
    MappingMode,
    NovelInsertion,
    Read,
    RTSwitchingFilter,
    Strand,
    infer_mst_from_connected_reads,
)
from scanmst.base.filters import ExonInfo
from scanmst.base.helper import (
    blat2chimeric_alignment,
    diff_chrom_diff_strand_handler,
    diff_chrom_same_strand_handler,
    extract_splice_sites,
    find_match_length,
    gene_annotation,
    get_transcriptome_length,
    insertion2chimeric_alignment,
    obtain_bp_region_seq,
    obtain_insertion_surrouding_cigarstrings,
    obtain_read_segment_length_from_cigar_string,
    obtain_variants_stats,
    same_chrom_diff_strand_handler,
    same_chrom_same_strand_handler,
    softclipped_length_and_event_size_checker,
    strand_mode_checker,
)
from scanmst.exception import ModesNotEqualError
from scanmst.graph import EdgeData, NLPath, Node, NodeIdentity, VariationType
from tests import FakeBlat


def _make_event(
    sv_type: str = "TDUP",
    bp1: str = "chr1:1200",
    bp2: str = "chr1:1800",
    mode1: MappingMode = MappingMode.SM,
    mode2: MappingMode = MappingMode.MS,
    ins1: str = "",
    ins2: str = "",
    strands: tuple[str, str] = ("+", "+"),
) -> Event:
    return Event(
        (
            sv_type,
            0,
            1,
            (bp1, bp2, mode1, mode2),
            (1000, 1200, [[1000, 1200]]),
            (1800, 2000, [[1800, 2000]]),
            (ins1, ins2),
            strands,
            ["G1", "G1"],
            False,
        ),
    )


def _build_nlpath(nodes: list[Node], sv_types: list[VariationType]) -> NLPath:
    edges_data = [
        EdgeData(
            variantion_type=sv_t,
            break_point1=BreakPoint(nodes[i].chrom, nodes[i].ref_end),
            break_point2=BreakPoint(nodes[i + 1].chrom, nodes[i + 1].ref_start),
            sr=1,
            read_ids=["r1"],
            mode1=MappingMode.MS,
            mode2=MappingMode.SM,
            gene1="G1",
            gene2="G1",
            annotation_code=3,
            splicing_code=1,
        )
        for i, sv_t in enumerate(sv_types)
    ]
    return NLPath.from_nodes_and_edges_data(nodes, edges_data)


GTF_TEXT = (
    'chr1\tsrc\tgene\t1001\t9000\t.\t+\t.\tgene_id "G1"; gene_name "GENE1";\n'
    'chr1\tsrc\texon\t1001\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
    'chr1\tsrc\texon\t5001\t6000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
    'chr1\tsrc\texon\t8001\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "GENE1";\n'
    'chr1\tsrc\tgene\t20001\t29000\t.\t-\t.\tgene_id "G2"; gene_name "GENE2";\n'
    'chr1\tsrc\texon\t20001\t21000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENE2";\n'
    'chr1\tsrc\texon\t24001\t25000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENE2";\n'
    'chr1\tsrc\texon\t28001\t29000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "GENE2";\n'
    'chr1\tsrc\tgene\t40001\t40500\t.\t+\t.\tgene_id "G3"; gene_name "SINGLE";\n'
    'chr1\tsrc\texon\t40001\t40500\t.\t+\t.\tgene_id "G3"; transcript_id "T3"; gene_name "SINGLE";\n'
)


@pytest.fixture(scope="module")
def gtf_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    p = tmp_path_factory.mktemp("filters_gtf") / "annot.gtf"
    p.write_text(GTF_TEXT)
    return p


class TestExonRTAndCircRNAFilters:
    """Test ExonFilter, RTSwitchingFilter, and full CircRNAFilter.is_circrna paths."""

    def test_exon_info_and_exon_filter(self, gtf_path: Path) -> None:
        info = ExonInfo("chr1", Exons.from_list([[100, 200]])[0], Strand.Forward, "T1")
        info.start = 110
        info.end = 190
        assert (info.start, info.end) == (110, 190)
        assert info.obtain_trx_id() == "T1"
        assert "T1" in repr(info)
        assert isinstance(hash(info), int)

        ef = ExonFilter(str(gtf_path), boundary_size=10)
        # TRA and DEL always bypass ExonFilter
        assert ef.is_breakpoints_in_same_exon(_make_event("TRA", "chr1:1200", "chr1:1500")) is False
        assert ef.is_breakpoints_in_same_exon(_make_event("DEL", "chr1:1200", "chr1:1500")) is False
        # Both inside shrunk exon [1010, 1990)
        assert ef.is_breakpoints_in_same_exon(_make_event("TDUP", "chr1:1200", "chr1:1500")) is True
        # Breakpoint 1 at 1005 is outside shrunk exon [1010, 1990)
        assert ef.is_breakpoints_in_same_exon(_make_event("TDUP", "chr1:1005", "chr1:1500")) is False
        # Different exons of T1
        assert ef.is_breakpoints_in_same_exon(_make_event("INV", "chr1:1500", "chr1:5500")) is False

    def test_rt_switching_filter(self) -> None:
        rtf = RTSwitchingFilter(rt_switching_filter_len=5)
        assert rtf.is_from_rt_switching(_make_event(ins1="-ACGTAC", ins2="-ACGTAC")) is True
        assert rtf.is_from_rt_switching(_make_event(ins1="-ACGTA", ins2="-ACGTA")) is False
        assert rtf.is_from_rt_switching(_make_event(ins1="+ACGTACGT", ins2="+ACGTACGT")) is False
        assert rtf.is_from_rt_switching(_make_event(ins1="", ins2="")) is False

    def test_circrna_one_hop_conditions_and_multihop(self, gtf_path: Path) -> None:
        cf = CircRNAFilter(str(gtf_path), boundary_size=10, breakpoint_diff_threshold=10)

        # Non-TDUP one-hop -> False
        n1 = Node("q", "chr1", "+", 1000, 6000, NodeIdentity.HEAD, Intervals.from_list([[1000, 2000], [5000, 6000]]))
        n2 = Node("q", "chr1", "+", 1000, 2000, NodeIdentity.TAIL, Intervals.from_list([[1000, 2000]]))
        assert cf.is_circrna(_build_nlpath([n1, n2], [VariationType.TRA])) is False

        # Condition 1: TDUP with exon superset + superposed with annotated exons of T1
        assert cf.is_circrna(_build_nlpath([n1, n2], [VariationType.TDUP])) is True

        # Condition 2: partial loop within annotated transcript T1 (+ strand) and T2 (- strand)
        p_pos_1 = Node("q", "chr1", "+", 5200, 6000, NodeIdentity.HEAD, Intervals.from_list([[5200, 6000]]))
        p_pos_2 = Node("q", "chr1", "+", 1000, 1800, NodeIdentity.TAIL, Intervals.from_list([[1000, 1800]]))
        assert cf.is_circrna(_build_nlpath([p_pos_1, p_pos_2], [VariationType.TDUP])) is True

        p_neg_1 = Node("q", "chr1", "-", 20000, 20800, NodeIdentity.HEAD, Intervals.from_list([[20000, 20800]]))
        p_neg_2 = Node("q", "chr1", "-", 24200, 25000, NodeIdentity.TAIL, Intervals.from_list([[24200, 25000]]))
        assert cf.is_circrna(_build_nlpath([p_neg_1, p_neg_2], [VariationType.TDUP])) is True

        # Condition 3: matching introns and start within breakpoint_diff_threshold
        c3_a = Node("q", "chr1", "+", 70000, 70500, NodeIdentity.HEAD, Intervals.from_list([[70000, 70500]]))
        c3_b = Node("q", "chr1", "+", 70005, 70300, NodeIdentity.TAIL, Intervals.from_list([[70005, 70300]]))
        assert cf.is_circrna(_build_nlpath([c3_a, c3_b], [VariationType.TDUP])) is True

        # Multi-hop (4 nodes, 3 TDUP hops): [e2] -> [e1,e2] -> [e1,e2] -> [e1]
        m1 = Node("q", "chr1", "+", 5000, 6000, NodeIdentity.HEAD, Intervals.from_list([[5000, 6000]]))
        m2 = Node("q", "chr1", "+", 1000, 6000, NodeIdentity.MID, Intervals.from_list([[1000, 2000], [5000, 6000]]))
        m3 = Node("q", "chr1", "+", 1000, 6000, NodeIdentity.MID, Intervals.from_list([[1000, 2000], [5000, 6000]]))
        m4 = Node("q", "chr1", "+", 1000, 2000, NodeIdentity.TAIL, Intervals.from_list([[1000, 2000]]))
        assert cf.is_circrna(_build_nlpath([m1, m2, m3, m4], [VariationType.TDUP] * 3)) is True

        # Multi-hop ONT condition: one TDUP hop with near-identical ref_start + mixed TRA hop
        o1 = Node("q", "chr1", "+", 70000, 70500, NodeIdentity.HEAD, Intervals.from_list([[70000, 70500]]))
        o2 = Node("q", "chr1", "+", 70003, 70400, NodeIdentity.MID, Intervals.from_list([[70003, 70400]]))
        o3 = Node("q", "chr2", "+", 1000, 1500, NodeIdentity.TAIL, Intervals.from_list([[1000, 1500]]))
        assert cf.is_circrna(_build_nlpath([o1, o2, o3], [VariationType.TDUP, VariationType.TRA])) is True


class TestSpliceCheckerMSTAndHelpers:
    """Test Read.splice_site_checker, mst_inference routing, and helper functions."""

    @pytest.fixture
    def synthetic_ref(self, tmp_path: Path) -> Fasta:
        # Build a 3000bp synthetic chromosome with planted GT-AG and CT-AC splice sites:
        # Forward intron [200, 300): bases 200..202 = "GT", bases 298..300 = "AG"
        # Reverse intron [500, 600): bases 500..502 = "CT" (rc="AG"), bases 598..600 = "AC" (rc="GT")
        # Non-canonical intron [800, 900): all 'C'
        seq = list("C" * 3000)
        seq[200:202] = list("GT")
        seq[298:300] = list("AG")
        seq[500:502] = list("CT")
        seq[598:600] = list("AC")
        fa_path = tmp_path / "synth.fa"
        fa_path.write_text(">chr1\n" + "".join(seq) + "\n>chr2\n" + "".join(seq) + "\n")
        return Fasta(str(fa_path), sequence_always_upper=True)

    def test_read_splice_site_checker_and_hash_repr(self, synthetic_ref: Fasta) -> None:
        unspliced = Read.new("r0", "chr1", 100, "+", "100M", 60, 0, "C" * 100, [30] * 100)
        assert unspliced.splice_site_checker(synthetic_ref) is True
        assert isinstance(hash(unspliced), int)
        assert "Read(" in repr(unspliced)

        fwd_can = Read.new("r1", "chr1", 150, "+", "50M100N50M", 60, 0, "C" * 100, [30] * 100)
        assert fwd_can.splice_site_checker(synthetic_ref) is True

        rev_can = Read.new("r2", "chr1", 450, "-", "50M100N50M", 60, 0, "G" * 100, [30] * 100)
        assert rev_can.splice_site_checker(synthetic_ref) is True

        non_can = Read.new("r3", "chr1", 750, "+", "50M100N50M", 60, 0, "C" * 100, [30] * 100)
        assert non_can.splice_site_checker(synthetic_ref) is False

    def test_helper_small_functions_and_cs_stats(self, gtf_path: Path, synthetic_ref: Fasta) -> None:
        assert get_transcriptome_length("human") > 0
        assert get_transcriptome_length("mouse") > 0

        subs, ins_f, del_f = obtain_variants_stats(":50*ag+acgt-gggg:50", large_indel_len_threshold=4)
        # ref length = 50 + 1 + 4 + 50 = 105; 1 substitution, 1 ins >=4, 1 del >=4
        assert subs == pytest.approx(1 / 105)
        assert ins_f == pytest.approx(1 / 105)
        assert del_f == pytest.approx(1 / 105)

        assert strand_mode_checker("+", "+", MappingMode.MS, MappingMode.SM) is True
        assert strand_mode_checker("+", "+", MappingMode.MS, MappingMode.MS) is False
        assert strand_mode_checker("+", "-", MappingMode.MS, MappingMode.MS) is True

        r_sm = Read.new("q", "chr1", 100, "+", "20S80M", 60, 0, "A" * 20 + "C" * 80, [30] * 100)
        r_ms = Read.new("q", "chr1", 100, "+", "80M20S", 60, 0, "C" * 80 + "T" * 20, [30] * 100)
        assert softclipped_length_and_event_size_checker(r_sm, MappingMode.SM, 30) is True
        assert softclipped_length_and_event_size_checker(r_ms, MappingMode.MS, 10) is False

        assert obtain_bp_region_seq(r_sm, MappingMode.SM, 4, synthetic_ref) == "+AAAA"
        assert obtain_bp_region_seq(r_ms, MappingMode.MS, 4, synthetic_ref) == "+TTTT"
        assert obtain_bp_region_seq(r_sm, MappingMode.SM, -3, synthetic_ref) == "-CCC"
        assert obtain_bp_region_seq(r_ms, MappingMode.MS, -3, synthetic_ref) == "-CCC"

        assert find_match_length("", "ABC", "right") == 0
        assert find_match_length("XXATG", "YYATG", "right") == 3
        assert find_match_length("ATGXX", "ATGYY", "left") == 3
        assert find_match_length("AAA", "BBB", "right") == 0
        with pytest.raises(ValueError, match="Invalid value for left_or_right"):
            find_match_length("A", "A", "middle")

        assert obtain_insertion_surrouding_cigarstrings("20M100I30M", 100) == ("20M", "30M")
        assert obtain_insertion_surrouding_cigarstrings("50M", 100) == ("", "")
        assert obtain_read_segment_length_from_cigar_string("10S20M5I3D100N") == 35

        cvg, gene_iv = extract_splice_sites(str(gtf_path), 5)
        assert gene_annotation("chrUn", 10, "chr1", 40200, gene_iv) == ("INTERGENIC", "SINGLE")
        assert "XX" in cvg[HTSeq.GenomicPosition("chr1", 40000)]

    def test_mst_inference_del_idup_and_routing(self, synthetic_ref: Fasta, fake_logger) -> None:
        cvg = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        gene_iv = HTSeq.GenomicArrayOfSets("auto", stranded=False)

        # 1. Splice-site check rejection on same-chrom same-strand and diff-chrom
        bad_lt = Read.new("q", "chr1", 750, "+", "50M100N50M50S", 60, 0, "C" * 150, [30] * 150)
        bad_rt = Read.new("q", "chr1", 100, "+", "100S50M", 60, 0, "C" * 150, [30] * 150)
        assert (
            infer_mst_from_connected_reads(bad_lt, bad_rt, MappingMode.MS, MappingMode.SM, 5, synthetic_ref, cvg, gene_iv, motif_required=False)
            is None
        )
        bad_rt_chr2 = Read.new("q", "chr2", 100, "+", "100S50M", 60, 0, "C" * 150, [30] * 150)
        assert (
            infer_mst_from_connected_reads(bad_lt, bad_rt_chr2, MappingMode.MS, MappingMode.SM, 5, synthetic_ref, cvg, gene_iv, motif_required=False)
            is None
        )

        # 2. DEL branch (evt_size <= 0) via MS + SM (is_reverse=True in same_chrom_same_strand_handler)
        r_del_1 = Read.new("q", "chr1", 100, "+", "50M50S", 60, 0, "C" * 100, [30] * 100)
        r_del_2 = Read.new("q", "chr1", 500, "+", "50S50M", 60, 0, "C" * 100, [30] * 100)
        del_res = infer_mst_from_connected_reads(
            r_del_1, r_del_2, MappingMode.MS, MappingMode.SM, 5, synthetic_ref, cvg, gene_iv, motif_required=False
        )
        assert del_res is not None
        assert del_res[0] == "DEL"

        # 3. Large TDUP branch (evt_size >= query_offset) with motif_required=False
        r_tdup_1 = Read.new("q", "chr1", 100, "+", "50S50M", 60, 0, "C" * 100, [30] * 100)
        r_tdup_2 = Read.new("q", "chr1", 100, "+", "50M50S", 60, 0, "C" * 100, [30] * 100)
        tdup_res = same_chrom_same_strand_handler(
            r_tdup_1, r_tdup_2, MappingMode.SM, MappingMode.MS, 5, synthetic_ref, cvg, gene_iv, False, fake_logger
        )
        assert tdup_res is not None
        assert tdup_res[0] == "TDUP"

        # 4. IDUP branch (same_chrom_diff_strand_handler when ra_bp == sa_bp)
        r_idup_1 = Read.new("q", "chr1", 300, "+", "50S50M", 60, 0, "C" * 100, [30] * 100)
        r_idup_2 = Read.new("q", "chr1", 300, "-", "50S50M", 60, 0, "C" * 100, [30] * 100)
        idup_res = infer_mst_from_connected_reads(
            r_idup_1, r_idup_2, MappingMode.SM, MappingMode.SM, 5, synthetic_ref, cvg, gene_iv, motif_required=False
        )
        assert idup_res is not None
        assert idup_res[0] == "IDUP"

        # Mismatched modes in same_chrom_diff_strand returns None; diff_chrom_diff_strand raises ModesNotEqualError
        assert (
            same_chrom_diff_strand_handler(r_idup_1, r_idup_2, MappingMode.SM, MappingMode.MS, 5, synthetic_ref, cvg, gene_iv, False, fake_logger)
            is None
        )
        r_chr2_neg = Read.new("q", "chr2", 300, "-", "50S50M", 60, 0, "C" * 100, [30] * 100)
        with pytest.raises(ModesNotEqualError):
            diff_chrom_diff_strand_handler(r_idup_1, r_chr2_neg, MappingMode.SM, MappingMode.MS, 5, synthetic_ref, cvg, gene_iv, False, fake_logger)

        # diff_chrom_same_strand_handler routing (SM+MS and MS+SM)
        r_chr2_pos = Read.new("q", "chr2", 300, "+", "50M50S", 60, 0, "C" * 100, [30] * 100)
        tra_res = diff_chrom_same_strand_handler(
            r_idup_1, r_chr2_pos, MappingMode.SM, MappingMode.MS, 5, synthetic_ref, cvg, gene_iv, False, fake_logger, 20
        )
        assert tra_res is not None
        assert tra_res[0] == "TRA"

    def test_blat2chimeric_and_insertion2chimeric(self, synthetic_ref: Fasta) -> None:
        class HitBlat(FakeBlat):
            def __init__(self, flag: bool, ins_obj: Insertion | NovelInsertion) -> None:
                super().__init__()
                self._flag = flag
                self._ins = ins_obj

            def query_insertion(self, *_a: object, **_kw: object) -> tuple[bool, Insertion | NovelInsertion]:
                return self._flag, self._ins

        assert isinstance(HitBlat(False, NovelInsertion(0, "")), Blat)
        assert blat2chimeric_alignment("A" * 40, "T" * 60 + "A" * 40, 100, "+", MappingMode.MS, None, 20, 10, synthetic_ref, 10) == ""

        for r_strand, r_mode, ins_strand in [
            ("+", MappingMode.MS, "+"),
            ("+", MappingMode.SM, "+"),
            ("-", MappingMode.MS, "-"),
            ("-", MappingMode.SM, "-"),
            ("+", MappingMode.MS, "-"),
            ("+", MappingMode.SM, "-"),
        ]:
            ins_var = Insertion(1, "chr1", 500, ins_strand, "40M", 60, 0, "A" * 40, [30] * 40)
            res_sa = blat2chimeric_alignment(
                soft_seq_original="A" * 40,
                seq_original=("T" * 60 + "A" * 40) if r_mode == MappingMode.MS else ("A" * 40 + "T" * 60),
                read_length=100,
                read_strand=r_strand,
                read_mode=r_mode,
                aligner=HitBlat(True, ins_var),
                mapq_cutoff=20,
                max_allowed_nm=10,
                genome_fasta=synthetic_ref,
                rt_switching_filter_len=10,
            )
            assert res_sa.startswith(f"chr1,501,{ins_strand},")

        # insertion2chimeric_alignment with a real pysam.AlignedSegment (reference_end=1170 >= ref_end_aligner=1150)
        header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 3000}]})
        seg = pysam.AlignedSegment(header)
        seg.query_name = "ins_read"
        seg.reference_id = 0
        seg.reference_start = 1000
        seg.cigarstring = "50M100I120M"
        seg.query_sequence = "A" * 270
        seg.mapping_quality = 60
        seg.set_tag("NM", 100)

        assert insertion2chimeric_alignment(seg, 1050, "A" * 100, 270, "+", 10, None) == ("", "")

        dup_ins = Insertion(1, "chr1", 1050, "+", "100M", 60, 0, "A" * 100, [30] * 100)
        ra_cigar, sa_cigar = insertion2chimeric_alignment(seg, 1050, "A" * 100, 270, "+", 10, HitBlat(True, dup_ins), align_len_threshold=50)
        assert ra_cigar == "150S120M"
        assert sa_cigar == "chr1,1001,+,150M120S,60,0;"
