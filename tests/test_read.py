"""Tests for Read, ReadsConnector, and detect_read_read_connections_from_cigar."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest
from Bio import SearchIO
from pyfaidx import Fasta

from scanmst import Read, ReadsConnector, cppext
from scanmst.base import MappingMode, Strand, detect_read_read_connections_from_cigar
from tests import FakeBlat, FakeHsp


@pytest.fixture
def param_dict() -> dict[str, object]:
    """Create a dict of parameters for read init."""
    return {
        "query_name": "test_name",
        "chrom": "chr1",
        "strand": "+",
        "ref_start": 1,
        "cigar_str": "10S50M5S",
        "mapq": 10,
        "nm": 0,
        "query_seq": "A" * 65,
        "query_qualities": None,
    }


@pytest.fixture
def connector_ref(tmp_path: Path) -> Fasta:
    """Synthetic reference genome for ReadsConnector realignment tests."""
    fa = tmp_path / "conn_ref.fa"
    # Alternating 'GC' (self-reverse-complement) ensures 'A'/'T' read segments have 0 homology shift on both strands
    fa.write_text(">chr1\n" + ("GC" * 1000) + "\n>chrM\n" + ("GC" * 250) + "\n")
    return Fasta(str(fa), sequence_always_upper=True)


class TestRead:
    """Test the Read class and CIGAR feature extraction."""

    @pytest.mark.parametrize(
        ("cigar", "expected_result"),
        [
            ("10S50M5S", (10, 5, 50, 50, 0, [0, 50], 65, [4, 10, 0, 50, 4, 5])),
            (
                "15S5I1D80M",
                (
                    15,
                    0,
                    85,
                    81,
                    -4,
                    [1, 5, 2, 1, 0, 80],
                    100,
                    [4, 15, 1, 5, 2, 1, 0, 80],
                ),
            ),
        ],
    )
    def test_parse_cigar(self, cigar: str, expected_result: tuple[object, ...]) -> None:
        result = cppext.parseCigar(cigar)
        assert result.lt_soft_len == expected_result[0]
        assert result.rt_soft_len == expected_result[1]
        assert result.read_match == expected_result[2]
        assert result.ref_match == expected_result[3]
        assert result.indel_len == expected_result[4]
        assert result.cigartuples_without_soft == expected_result[5]
        assert result.query_len == expected_result[6]
        assert result.cigartuples == expected_result[7]

    def test_init(self, param_dict: dict[str, object]) -> None:
        read = Read.new(**param_dict)  # type: ignore[arg-type]
        assert read.query_name == param_dict["query_name"]
        assert read.lt_soft_len == 10
        assert read.rt_soft_len == 5
        assert read.query_length == 65

    def test_get_exons_and_introns(self, param_dict: dict[str, object]) -> None:
        read = Read.new(**param_dict)  # type: ignore[arg-type]
        assert [(e.start, e.end) for e in read.get_exons()] == [(1, 51)]
        param_dict["cigar_str"] = "15S2M2N3I50M2S"
        read = Read.new(**param_dict)  # type: ignore[arg-type]
        assert [(e.start, e.end) for e in read.get_exons()] == [(1, 3), (5, 55)]
        assert [(i.start, i.end) for i in read.get_exons().introns()] == [(3, 5)]


@pytest.fixture
def reads_connector(reads: list[Read], fake_blat: FakeBlat, fake_logger, connector_ref: Fasta) -> ReadsConnector:
    """Create a ReadsConnector backed by a real synthetic Fasta."""
    return ReadsConnector(
        reads,
        fake_blat,
        fake_logger,
        mapq_cutoff=0,
        genome_fasta=connector_ref,
        rt_switching_filter_len=10,
    )


class TestReadsConnector:
    """Test ReadsConnector connection, microhomology, and BLAT realignment logic."""

    def test_reset_and_increment_index(self, reads_connector: ReadsConnector) -> None:
        reads_connector.index = 1
        reads_connector.reset_index()
        assert reads_connector.index == 0
        reads_connector.increment_index()
        assert reads_connector.index == 1

    def test__get_mode_and_init_read_mode(self, reads_connector: ReadsConnector) -> None:
        assert reads_connector._get_mode((1, 2, 3)) == 1
        assert reads_connector._get_mode((1, 2, 1)) == 1
        assert reads_connector._get_mode((3, 1, 1)) == 2

        read1, read2 = reads_connector.aln_list
        read1.mode = read2.mode = None
        read1.adhocsms = read1.sms  # type: ignore[assignment]
        read2.adhocsms = read2.sms  # type: ignore[assignment]
        reads_connector.init_read_mode(read1, read2)
        assert read1.mode == MappingMode.SM
        assert read2.mode == MappingMode.MS

    @pytest.mark.parametrize(
        ("query_seq", "target_seq", "same_strand", "expected_result"),
        [
            ("ATCG", "ATCG", True, False),
            ("ATCGC", "TCATCGC", True, True),
            ("ATCGC", "GCGATGA", False, True),
        ],
    )
    def test_check_if_ms_match(
        self,
        reads_connector: ReadsConnector,
        query_seq: str,
        target_seq: str,
        same_strand: bool,
        expected_result: bool,
    ) -> None:
        assert (
            reads_connector.check_if_ms_match(
                query_seq,
                target_seq,
                same_strand=same_strand,
                minimum_s_length=4,
            )
            == expected_result
        )

    def test_microhomology_and_update_query_sequence(self, reads_connector: ReadsConnector) -> None:
        # bp_region_seq_len = 100 - 0 - 0 - 52 - 52 = -4 -> microhomology of 4 bp
        is_mh, mh_len = ReadsConnector._determine_microhomology_len(
            "A" * 52,
            100,
            (48, 52, 0),
            (0, 52, 48),
            MappingMode.SM,
            MappingMode.MS,
        )
        assert (is_mh, mh_len) == (True, 4)
        trimmed_sm = ReadsConnector.update_query_sequence(
            "CCCC" + "A" * 48,
            "X" * 100,
            (48, 52, 0),
            (0, 52, 48),
            MappingMode.SM,
            MappingMode.MS,
        )
        assert trimmed_sm == "A" * 48

        trimmed_ms = ReadsConnector.update_query_sequence(
            "A" * 48 + "CCCC",
            "X" * 100,
            (0, 52, 48),
            (48, 52, 0),
            MappingMode.MS,
            MappingMode.SM,
        )
        assert trimmed_ms == "A" * 48

    def test__check_if_strand_mode_for_compare_ms(self, reads_connector: ReadsConnector) -> None:
        read1, read2 = reads_connector.aln_list
        assert reads_connector._check_if_strand_mode_for_compare_ms(
            read1,
            read2,
            same_strand=True,
            first_is_matched=True,
            second_is_matched=False,
        )
        assert not reads_connector._check_if_strand_mode_for_compare_ms(
            read1,
            read2,
            same_strand=False,
            first_is_matched=False,
            second_is_matched=True,
        )

    def test_test_2case_and_determine_new_read_mode(self, reads_connector: ReadsConnector) -> None:
        read1, read2 = reads_connector.aln_list
        read1.adhocsms = read1.sms  # type: ignore[assignment]
        read2.adhocsms = read2.sms  # type: ignore[assignment]
        read1.adhocseq = read1.query_sequence  # type: ignore[assignment]
        read2.adhocseq = read2.query_sequence  # type: ignore[assignment]
        assert reads_connector.test_2case(read1, read2, is_compare_for_ms=False)
        assert reads_connector.test_2case(read2, read1, is_compare_for_ms=False)

        reads_connector.init_read_mode(read1, read2)
        assert read1.mode == MappingMode.SM
        new_mode = reads_connector._double_check_for_start_end_read_determine_new_read_mode(
            read1,
            Strand.Forward,
        )
        assert read1.mode == MappingMode.MS
        assert new_mode == MappingMode.SM

    def test__double_check_create_new_read_calculate_sms(
        self,
        connector_ref: Fasta,
        fake_logger,
        tmp_path: Path,
    ) -> None:
        assert ReadsConnector._find_match_length("", "A", "right") == 0
        with pytest.raises(ValueError, match="Invalid value"):
            ReadsConnector._find_match_length("AA", "AA", "middle")

        # Read has 40bp left soft-clip ('A'*40) and 60bp matched ('T'*60) in MS mode
        read = Read.new(
            query_name="realign_q",
            chrom="chr1",
            ref_start=1000,
            strand="+",
            cigar_str="40S60M",
            mapq=60,
            nm=0,
            query_seq=("A" * 40) + ("T" * 60),
            query_qualities=[30] * 100,
        )
        read.mode = MappingMode.MS
        fake_hsp = FakeHsp(query_start=0, query_end=40, query_seq="A" * 40, hit_range=(500, 540))
        blat_double = FakeBlat(psl2sam_return=("chr1", 500, "+", "40M", 0))
        rc = ReadsConnector([read], blat_double, fake_logger, 20, connector_ref, rt_switching_filter_len=10)

        new_read = rc._double_check_create_new_read_calculate_sms(fake_hsp, "A" * 40, read, 60)
        assert new_read is not None
        assert read.mode == MappingMode.SM
        assert new_read.mode == MappingMode.MS
        assert new_read.chrom == "chr1"
        assert new_read.ref_start == 500
        assert new_read.cigarstring == "40M60S"

        # Exceeding rt_switching_filter_len returns None (8 matching 'T's followed by 'G' mismatch at 540..550)
        poly_t_seq = list("G" * 2000)
        poly_t_seq[540:548] = list("T" * 8)
        poly_t_fa = tmp_path / "polyt.fa"
        poly_t_fa.write_text(">chr1\n" + "".join(poly_t_seq) + "\n")
        rc_rt = ReadsConnector(
            [read],
            blat_double,
            fake_logger,
            20,
            Fasta(str(poly_t_fa), sequence_always_upper=True),
            rt_switching_filter_len=5,
        )
        read.mode = MappingMode.MS
        assert rc_rt._double_check_create_new_read_calculate_sms(fake_hsp, "A" * 40, read, 60) is None

        # Exercise all (r_strand, init_mode, blat_strand) branches with both unclipped ('40M') and clipped ('5S30M5S') HSPs
        partial_hsp = FakeHsp(query_start=5, query_end=35, query_seq="A" * 30, hit_range=(500, 530))
        for hsp_obj, blat_cigar in [(fake_hsp, "40M"), (partial_hsp, "5S30M5S")]:
            for r_strand, init_mode, blat_strand, cigar_in in [
                ("+", MappingMode.MS, "+", "40S60M"),
                ("+", MappingMode.SM, "+", "60M40S"),
                ("-", MappingMode.MS, "-", "40S60M"),
                ("-", MappingMode.SM, "-", "60M40S"),
                ("+", MappingMode.SM, "-", "60M40S"),
                ("+", MappingMode.MS, "-", "40S60M"),
                ("-", MappingMode.SM, "+", "60M40S"),
                ("-", MappingMode.MS, "+", "40S60M"),
            ]:
                qseq = (("A" * 40) + ("T" * 60)) if cigar_in == "40S60M" else (("T" * 60) + ("A" * 40))
                r_var = Read.new("rq", "chr1", 1000, r_strand, cigar_in, 60, 0, qseq, [30] * 100)
                r_var.mode = init_mode
                bd = FakeBlat(psl2sam_return=("chr1", 500, blat_strand, blat_cigar, 0))
                rc_var = ReadsConnector([r_var], bd, fake_logger, 20, connector_ref, rt_switching_filter_len=10)
                realigned = rc_var._double_check_create_new_read_calculate_sms(hsp_obj, "A" * 40, r_var, 60)
                assert realigned is not None
                assert realigned.chrom == "chr1"

    def test__double_check_for_start_end_read(
        self,
        connector_ref: Fasta,
        fake_logger,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        monkeypatch.setattr(SearchIO, "read", lambda *_a, **_kw: SimpleNamespace(hsps=[]))
        read = Read.new(
            query_name="start_end_q",
            chrom="chr1",
            ref_start=1000,
            strand="+",
            cigar_str="40S60M",
            mapq=60,
            nm=0,
            query_seq=("A" * 40) + ("C" * 60),
            query_qualities=[30] * 100,
        )
        read.mode = MappingMode.MS

        fake_hsp = FakeHsp(query_start=0, query_end=40, query_seq="A" * 40, hit_range=(500, 540))
        blat_double = FakeBlat(
            query_return="dummy.pslx",
            query_insertion_return=(1, fake_hsp, 60),  # type: ignore[arg-type]
            psl2sam_return=("chr1", 500, "+", "40M", 0),
        )
        rc = ReadsConnector([read], blat_double, fake_logger, 20, connector_ref, rt_switching_filter_len=10)
        rc.reads_chain = [read]
        assert rc._double_check_for_start_end_read(read, "start") is True
        assert len(rc.reads_chain) == 2
        assert rc.reads_chain[0].cigarstring == "40M60S"
        assert (rc.reads_chain[0], read) in rc.read_pair_mode_dict

        # Long unaligned soft-clip (>= min_soft_seg_len) returns False
        long_read = Read.new(
            query_name="long_unaligned",
            chrom="chr1",
            ref_start=1000,
            strand="+",
            cigar_str="220S60M",
            mapq=60,
            nm=0,
            query_seq=("A" * 220) + ("C" * 60),
            query_qualities=[30] * 280,
        )
        long_read.mode = MappingMode.MS
        blat_miss = FakeBlat(query_insertion_return=(0, "", 0))  # type: ignore[arg-type]
        rc_miss = ReadsConnector([long_read], blat_miss, fake_logger, 20, connector_ref, rt_switching_filter_len=10)
        assert rc_miss._double_check_for_start_end_read(long_read, "start") is False

    def test_connect_two_and_three_segments(self, reads_connector: ReadsConnector, connector_ref: Fasta, fake_logger) -> None:
        assert reads_connector.connect() is True
        assert len(reads_connector.reads_chain) == 2

        # 3-segment multi-hop chain: seg1 (50M100S) -> seg2 (50S50M50S) -> seg3 (100S50M)
        s1 = "ACGT" * 12 + "AA"  # 50bp
        s2 = "TGCA" * 12 + "CC"  # 50bp
        s3 = "GGTT" * 12 + "GG"  # 50bp
        full_seq = s1 + s2 + s3
        r1 = Read.new("mhop", "chr1", 100, "+", "50M100S", 60, 0, full_seq, [30] * 150)
        r2 = Read.new("mhop", "chr1", 500, "+", "50S50M50S", 60, 0, full_seq, [30] * 150)
        r3 = Read.new("mhop", "chr1", 900, "+", "100S50M", 60, 0, full_seq, [30] * 150)
        rc3 = ReadsConnector([r1, r2, r3], None, fake_logger, 20, connector_ref, rt_switching_filter_len=10)
        assert rc3.connect() is True
        assert rc3.reads_chain == [r1, r2, r3]


class TestDetectReadReadConnectionsFromCigar:
    """Test detect_read_read_connections_from_cigar filtering and connection."""

    @staticmethod
    def _make_bam_seg(
        chrom: str = "chr1",
        pos: int = 100,
        cigar: str = "50M50S",
        seq: str | None = None,
        sa: str | None = "chr1,501,+,50S50M,60,0;",
        mapq: int = 60,
        nm: int = 0,
    ) -> pysam.AlignedSegment:
        header = pysam.AlignmentHeader.from_dict(
            {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 5000}, {"SN": "chrM", "LN": 1000}]},
        )
        seg = pysam.AlignedSegment(header)
        seg.query_name = "chim_q"
        seg.reference_id = 0 if chrom == "chr1" else 1
        seg.reference_start = pos
        seg.cigarstring = cigar
        seg.query_sequence = seq if seq is not None else ("A" * 50 + "C" * 50)
        seg.query_qualities = pysam.qualitystring_to_array("I" * len(seg.query_sequence))
        seg.mapping_quality = mapq
        seg.set_tag("NM", nm)
        if sa is not None:
            seg.set_tag("SA", sa)
        return seg

    def test_detect_connections_filters_and_success(self, connector_ref: Fasta, fake_logger) -> None:
        # Missing SA tag -> None
        no_sa = self._make_bam_seg(sa=None)
        assert detect_read_read_connections_from_cigar(no_sa, 20, 10, None, 0.9, connector_ref, 10, 200, fake_logger) is None

        # Mitochondrial primary -> None
        mt_seg = self._make_bam_seg(chrom="chrM")
        assert detect_read_read_connections_from_cigar(mt_seg, 20, 10, None, 0.9, connector_ref, 10, 200, fake_logger) is None

        # Mitochondrial SA -> None
        mt_sa = self._make_bam_seg(sa="chrM,101,+,50S50M,60,0;")
        assert detect_read_read_connections_from_cigar(mt_sa, 20, 10, None, 0.9, connector_ref, 10, 200, fake_logger) is None

        # Excessive NM -> None
        high_nm = self._make_bam_seg(nm=50)
        assert detect_read_read_connections_from_cigar(high_nm, 20, 10, None, 0.9, connector_ref, 10, 200, fake_logger) is None

        # Low MAPQ -> None
        low_mapq = self._make_bam_seg(mapq=5)
        assert detect_read_read_connections_from_cigar(low_mapq, 20, 10, None, 0.9, connector_ref, 10, 200, fake_logger) is None

        # Reverse transcription artifact (opposite strand on same chrom at identical start) -> None
        rt_art = self._make_bam_seg(pos=100, sa="chr1,105,-,50S50M,60,0;")
        assert detect_read_read_connections_from_cigar(rt_art, 20, 10, None, 0.9, connector_ref, 10, 200, fake_logger) is None

        # Valid 2-segment chimeric read -> returns (reads_chain, pair_mode_dict, num_added_reads)
        valid = self._make_bam_seg()
        ret = detect_read_read_connections_from_cigar(valid, 20, 10, None, 0.9, connector_ref, 10, 200, fake_logger)
        assert ret is not None
        chain, mode_dict, num_added = ret
        assert len(chain) == 2
        assert len(mode_dict) == 1
        assert num_added == 0
