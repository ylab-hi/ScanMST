"""Test the read module.

@Filename:    test_read.py
@license:     MIT Licence
@Time:        2/4/22 10:32 AM
"""
import pytest
from scannls import Read, ReadsConnector, cppext

from tests import FakeHsp


@pytest.fixture()
def param_dict():
    """Create a dict of parameters for read init."""
    return {
        "query_name": "test_name",
        "chrom": "chr1",
        "strand": "+",
        "ref_start": 1,
        "cigar_str": "10S50M5S",
        "mapq": 10,
        "nm": 0,
        "query_seq": "A",
    }


class TestRead:
    """Test the read module."""

    @pytest.mark.parametrize(
        "cigar, expected_result",
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
    def test_parse_cigar(self, cigar, expected_result):
        """Test the calculate_features function."""
        result = cppext.parseCigar(cigar)
        assert result.lt_soft_len == expected_result[0]
        assert result.rt_soft_len == expected_result[1]
        assert result.read_match == expected_result[2]
        assert result.ref_match == expected_result[3]
        assert result.indel_len == expected_result[4]
        assert result.cigartuples_without_soft == expected_result[5]
        assert result.query_len == expected_result[6]
        assert result.cigartuples == expected_result[7]

    def test_init(self, param_dict):
        """Test the init function."""
        read = Read.new(**param_dict)

        assert read.query_name == param_dict["query_name"]
        assert read.lt_soft_len == 10
        assert read.rt_soft_len == 5
        assert read.query_length == 65

    def test_get_exons_and_introns(self, param_dict):
        """Test the get_exons_and_introns function."""
        read = Read.new(**param_dict)
        assert read.get_exons_and_introns() == ([[1, 51]], [])
        param_dict["cigar_str"] = "15S2M2N3I50M2S"
        read = Read.new(**param_dict)
        assert read.get_exons_and_introns() == ([[1, 3], [5, 55]], [[3, 5]])


@pytest.fixture()
def reads_connector(reads, fake_blat, fake_logger):
    """Create a ReadsConnector object.

    reads from conftest.py
    """
    return ReadsConnector(reads, fake_blat, fake_logger)


class TestReadsConnector:
    """Test ReadsConnector."""

    def test_reset_index(self, reads_connector):
        """Test reset index."""
        reads_connector.index = 1
        reads_connector.reset_index()
        assert reads_connector.index == 0

    def test_increment_index(self, reads_connector):
        """Test increment index."""
        reads_connector.increment_index()
        assert reads_connector.index == 1

    def test__get_mode(self, reads_connector):
        """Test the _get_mode function."""
        assert reads_connector._get_mode((1, 2, 3)) == 1
        assert reads_connector._get_mode((1, 2, 1)) == 1
        assert reads_connector._get_mode((3, 1, 1)) == 2

    def test_init_mode_judge(self, reads_connector):
        """Test the init_mode_judge function."""
        read1, read2 = reads_connector.aln_list
        read1.mode = read2.mode = -1
        read1.adhocsms = read1.sms
        read2.adhocsms = read2.sms
        reads_connector.init_mode_judge(read1, read2)
        assert read1.mode == 2
        assert read2.mode == 1

    @pytest.mark.parametrize(
        "query_seq, target_seq, same_strand, expected_result",
        [
            ("ATCG", "ATCG", True, False),
            ("ATCGC", "TCATCGC", True, True),
        ],
    )
    def test_check_if_ms_match(
        self, reads_connector, query_seq, target_seq, same_strand, expected_result
    ):
        """Test the check_if_ms_match function."""
        assert (
            reads_connector.check_if_ms_match(
                query_seq, target_seq, same_strand, minimum_s_length=4
            )
            == expected_result
        )

    def test__determine_microhomology_len(self, reads_connector):
        """Test determine microhomology length."""
        read1, read2 = reads_connector.aln_list
        assert reads_connector._determine_microhomology_len(
            read1.read_match_size,
            read1.query_length,
            read1.sms,
            read2.sms,
            read1.mode,
            read2.mode,
        ) == (False, 0)

    def test_update_query_sequence(self, reads_connector):
        """Test update query sequence."""
        read1, read2 = reads_connector.aln_list
        assert (
            reads_connector.update_query_sequence(
                read1.query_sequence[read1.lt_soft_len :],
                read1.query_sequence,
                read1.sms,
                read2.sms,
                read1.mode,
                read2.mode,
            )
            == read1.query_sequence[read1.lt_soft_len :]
        )

    def test__check_if_strand_mode_for_compare_ms(self, reads_connector):
        """Test check if strand mode for compare ms."""
        read1, read2 = reads_connector.aln_list
        assert reads_connector._check_if_strand_mode_for_compare_ms(
            read1, read2, True, True, False
        )
        assert not reads_connector._check_if_strand_mode_for_compare_ms(
            read1, read2, False, False, True
        )

    def test_test_2case(self, reads_connector):
        """Test test 2case."""
        read1, read2 = reads_connector.aln_list
        read1.adhocsms = read1.sms
        read2.adhocsms = read2.sms
        read1.adhocseq = read1.query_sequence
        read2.adhocseq = read2.query_sequence
        assert reads_connector.test_2case(read1, read2, False)
        assert reads_connector.test_2case(read2, read1, False)

    def test__double_check_for_start_end_read_determine_new_read_mode(
        self, reads_connector
    ):
        """Test double check for start end read."""
        read1, read2 = reads_connector.aln_list
        read1.adhocsms = read1.sms
        read2.adhocsms = read2.sms
        reads_connector.init_mode_judge(read1, read2)
        assert read1.mode == 2
        assert read2.mode == 1
        read1.mode = 1
        reads_connector._double_check_for_start_end_read_determine_new_read_mode(
            read1, read2
        )
        assert read2.mode == 1

    def test__double_check_create_new_read_calculate_sms(self, reads_connector):
        """Test double check create new read calculate sms."""
        read1, _ = reads_connector.aln_list
        read1.mode = 2
        fake_hsp = FakeHsp(query_start=1, query_end=3, query_seq="ATC")
        reads_connector.blat.psl2sam_return = ("chr1", 1, "+", "2S1M1S", 2)

        new_read = reads_connector._double_check_create_new_read_calculate_sms(
            fake_hsp, "ATCGC", read1
        )

        assert read1.mode == 1
        assert new_read.mode == 2
        assert new_read.cigarstring == "1448S1M3S"
        assert new_read.query_sequence == read1.query_sequence
        assert new_read.chrom == "chr1"

    @pytest.mark.skip(reason="not implemented")
    def test__double_check_for_start_end_read(
        self,
        reads_connector,
        monkeypatch,
        align_len_threshold=1,
        threshold_identity=0.9,
        top=1,
    ):
        """Test double check for start end read."""
        from Bio import SearchIO  # type: ignore

        monkeypatch.setattr(SearchIO, "read", lambda x, y: None)
        read1, read2 = reads_connector.aln_list
        read1.mode = 2
        fake_hsp = FakeHsp(query_start=1, query_end=3, query_seq="ATC")
        reads_connector.blat.query_insertion_return = (1, [fake_hsp])
        reads_connector.align_len_threshold = align_len_threshold
        reads_connector.threshold_identity = threshold_identity
        reads_connector.top = top

    def test_connect(self, reads_connector):
        """Test connect."""
        assert reads_connector.connect()
        assert len(reads_connector.reads_chain) == len(reads_connector.aln_list)
