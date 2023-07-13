"""Connecter Reads.
@Time:        12/15/21 2:14 PM.
"""
import re
from itertools import combinations
from typing import Any

from Bio import SearchIO
from loguru import logger

from scannls.cli.helper import cigar_validity
from scannls.type import LoggerType

from .basic_class import reverse_complement
from .basic_read import Read
from .blat import Blat


class ReadsConnector:
    """ReadsConnector class is used to connect the reads and identify the mode of the reads.

    :param read_list: the list of the alignment
    :param blat: :class: `class.Blat` for the BLAT search
    :param logger: :class: `loguru.logger` for logging

    :Example:

    >>> from loguru import  logger
    >>> read_list = []
    >>> blat = Blat(ref_2bit='reference.2bit', logger= logger, port=88888, output_dir='/tmp')
    >>> readconnector = ReadsConnector(read_list=read_list, blat=blat, logger=logger)
    >>> readconnector.connect()
    >>> readconnector.reads_chain
    [Read(chr1, 6524193, 6524850, +, 60, 8), Read(chr1, 6522473, 6522883, +, 60, 4)]
    >>> readconnector.read_pair_mode_dict
    {(Read(chr1, 6524193, 6524850, +, 60, 8), Read(chr1, 6522473, 6522883, +, 60, 4)): (1, 2)}
    """

    def __init__(
        self,
        read_list: list[Read],
        blat: Blat,
        logger: LoggerType,
        align_len_threshold: int = 20,
        threshold_identity: float = 0.99,
        top: int = 3,
    ) -> None:
        """Initialize the ReadsConnector class."""
        self.candidate_nodes: list[Read] = []
        self.reads_chain: list[Read] = []
        self.read_pair_mode_dict, self.insertion_dict = {}, {}  # type: ignore
        self.aln_list = read_list
        self.logger = logger
        self.blat = blat
        self.index = 0
        self.num_added_reads = 0

        self.align_len_threshold: int = align_len_threshold
        self.threshold_identity: float = threshold_identity
        self.top: int = top

    def reset_index(self) -> None:
        """Reset the index in order to fetch read in  candidate reads in new iteration."""
        self.index = 0

    def increment_index(self) -> None:
        """Increment the index in order to fetch read in  candidate nodes."""
        self.index += 1

    @staticmethod
    def _get_mode(sms: Any) -> int:
        """Get the mode of the reads.

        :param sms:
        :return: mode

        .. note::
            the length of left s more than right s, the mode is 2
            the length of left s less than right s, the mode is 1
        """
        _lt, _, _rt = sms
        # SM
        if _lt > _rt:
            return 2
        # MS
        return 1

    @staticmethod
    def init_mode_judge(read1: Read, read2: Read) -> None:
        """Initialize the mode of the reads."""
        sorted_by_s_length = sorted(
            [read1, read2],
            key=lambda x: min(x.lt_soft_len, x.rt_soft_len),
        )

        if read1 == sorted_by_s_length[0]:
            read1.mode = ReadsConnector._get_mode(read1.adhocsms)
            if read1.strand == read2.strand:
                read2.mode = 1 if read1.mode == 2 else 2
            else:
                read2.mode = read1.mode

        else:
            read2.mode = ReadsConnector._get_mode(read2.sms)
            if read1.strand == read2.strand:
                read1.mode = 1 if read2.mode == 2 else 2
            else:
                read1.mode = read2.mode

    def check_if_ms_match(
        self,
        query_seq: str,
        target_seq: str,
        *,
        same_strand: bool,
        minimum_s_length: int = 30,
        minimum_terminal_length: int = 5,
    ) -> bool:
        """Compare m of start read with s of read.

        query_seq: M  target_seq: S

        :param minimum_terminal_length: the minimum length away from the terminal
        :param minimum_s_length: the minimum length of the s
        :param same_strand: whether the start read and read are on the same strand
        :param target_seq: the s of the read
        :param query_seq: the m of the start read
        """
        match_flag = False

        self.logger.trace(
            f"query length ={len(query_seq)} target length ={len(target_seq)}",
        )
        if len(target_seq) <= minimum_s_length:
            return match_flag

        if not same_strand:
            target_seq = reverse_complement(target_seq)

        pattern = re.compile(f"({query_seq})")
        temp_indices = [item.span() for item in re.finditer(pattern, target_seq)]
        if temp_indices:
            min_indices = [
                min(temp_index[0], len(target_seq) - temp_index[1])
                for temp_index in temp_indices
            ]
            index = min(min_indices)
            if index <= minimum_terminal_length:
                match_flag = True

        return match_flag

    @staticmethod
    def _determine_microhomology_len(
        read_match_sequence: str,
        read_query_length: int,
        prev_sms: tuple[int, int, int],
        next_sms: tuple[int, int, int],
        prev_read_mode: int,
        next_read_mode: int,
    ) -> Any:
        """Determine the length of the microhomology.

        :param read_match_sequence: the match sequence of the read
        :param read_query_length: whole reads length
        :param prev_sms: previous read sms
        :param next_sms: next read sms
        :param prev_read_mode: previous predicted connected read mode
        :param next_read_mode: next predicted connected read mode
        """
        bp_region_seq_len = 0
        _lt_len_r1, _read_match_r1, _rt_len_r1 = prev_sms
        _lt_len_r2, _read_match_r2, _rt_len_r2 = next_sms
        if prev_read_mode == 2:
            if next_read_mode == 2:
                bp_region_seq_len = (
                    read_query_length
                    - _rt_len_r1
                    - _rt_len_r2
                    - _read_match_r1
                    - _read_match_r2
                )
            elif next_read_mode == 1:
                bp_region_seq_len = (
                    read_query_length
                    - _rt_len_r1
                    - _lt_len_r2
                    - _read_match_r1
                    - _read_match_r2
                )
        elif next_read_mode == 2:
            bp_region_seq_len = (
                read_query_length
                - _lt_len_r1
                - _rt_len_r2
                - _read_match_r1
                - _read_match_r2
            )
        elif next_read_mode == 1:
            bp_region_seq_len = (
                read_query_length
                - _lt_len_r1
                - _lt_len_r2
                - _read_match_r1
                - _read_match_r2
            )
        is_microhomology = False
        microhomology_length = 0

        if bp_region_seq_len < 0:
            microhomology_length = -bp_region_seq_len
            if microhomology_length < len(read_match_sequence):
                is_microhomology = True
        return is_microhomology, microhomology_length

    @staticmethod
    def update_query_sequence(
        read_match_sequence: str,
        read_query_sequence: str,
        prev_sms: tuple[int, int, int],
        next_sms: tuple[int, int, int],
        prev_read_mode: int,
        next_read_mode: int,
    ) -> Any:
        """Update the query sequence based on the length of the microhomology."""
        (
            is_microhomology,
            microhomology_length,
        ) = ReadsConnector._determine_microhomology_len(
            read_match_sequence,
            len(read_query_sequence),
            prev_sms,
            next_sms,
            prev_read_mode,
            next_read_mode,
        )
        if is_microhomology:
            if prev_read_mode == 2:
                return read_match_sequence[microhomology_length:]

            if prev_read_mode == 1:
                return read_match_sequence[:-microhomology_length]
            return None

        return read_match_sequence

    @staticmethod
    def _check_if_strand_mode_for_compare_ms(
        start_read: Read,
        read: Read,
        *,
        same_strand: bool,
        first_is_matched: bool = False,
        second_is_matched: bool = False,
    ) -> bool:
        """Check the strand mode for compare ms."""
        flag = True
        mode1 = start_read.mode
        mode2 = read.mode
        if first_is_matched:
            mode2 = 2
        if second_is_matched:
            mode2 = 1
        if same_strand:
            if mode1 == mode2:
                flag = False
        elif mode1 != mode2:
            flag = False
        return flag

    def test_2case(
        self,
        start_read: Read,
        read: Read,
        *,
        is_compare_for_ms: bool,
    ) -> Any:
        """Test 2 case for two reads to check if they are connected.

        start read -> read

        :param start_read: start read
        :param read: next read
        :param is_compare_for_ms: is compare for ms
        """
        condition1, condition2 = False, False
        self.logger.debug(f"{start_read.mode=}, {read.mode=}, {start_read.query_name=}")
        self.logger.debug(f"{start_read.adhocsms=}, {read.sms=}")
        if not is_compare_for_ms:  # one hop
            self.read_pair_mode_dict[(start_read, read)] = (start_read.mode, read.mode)
            self.reads_chain.append(read)
            return True, start_read

        _lt_len_r1, _read_match_r1, _ = start_read.adhocsms
        _lt_len_r2, _read_match_r2, _rt_len_r2 = read.sms
        read_query_sequence = read.query_sequence

        same_strand = start_read.adhocseq == read.query_sequence

        # first case
        self.logger.debug("testing first case M vs LS")

        next_read_mode = 2
        read_match_sequence = start_read.adhocseq[
            _lt_len_r1 : _lt_len_r1 + _read_match_r1
        ]
        read_match_sequence = ReadsConnector.update_query_sequence(
            read_match_sequence,
            read_query_sequence,
            start_read.adhocsms,
            read.sms,
            start_read.mode,
            next_read_mode,
        )

        match_flag1 = self.check_if_ms_match(
            read_match_sequence,
            read.query_sequence[:_lt_len_r2],
            same_strand=same_strand,
        )
        if match_flag1:
            condition1 = self._check_if_strand_mode_for_compare_ms(
                start_read,
                read,
                same_strand=same_strand,
                first_is_matched=match_flag1,
            )
        if match_flag1 and condition1:  # may same
            read.mode = 2
            self.logger.debug(f"{start_read.mode=}, {read.mode=}")
            self.read_pair_mode_dict[(start_read, read)] = (start_read.mode, read.mode)
            self.reads_chain.append(read)

            if read in self.candidate_nodes:
                self.candidate_nodes.remove(read)
                self.reset_index()

            start_read = read
            start_read.adhocsms = 0, _lt_len_r2 + _read_match_r2, _rt_len_r2
            start_read.adhocseq = read.query_sequence

            return True, start_read

        self.logger.debug("testing second case M vs RS")
        # second case
        next_read_mode = 1
        read_match_sequence = start_read.adhocseq[
            _lt_len_r1 : _lt_len_r1 + _read_match_r1
        ]
        read_match_sequence = ReadsConnector.update_query_sequence(
            read_match_sequence,
            read_query_sequence,
            start_read.adhocsms,
            read.sms,
            start_read.mode,
            next_read_mode,
        )

        match_flag2 = self.check_if_ms_match(
            read_match_sequence,
            read.query_sequence[-_rt_len_r2:],
            same_strand=same_strand,
        )

        if match_flag2:
            condition2 = self._check_if_strand_mode_for_compare_ms(
                start_read,
                read,
                same_strand=same_strand,
                first_is_matched=False,
                second_is_matched=match_flag2,
            )
        if match_flag2 and condition2:  # may same
            read.mode = 1

            self.logger.debug(f"{start_read.mode=}, {read.mode=}")

            self.read_pair_mode_dict[(start_read, read)] = (start_read.mode, read.mode)

            self.reads_chain.append(read)

            if read in self.candidate_nodes:
                self.candidate_nodes.remove(read)
                self.reset_index()

            start_read = read
            start_read.adhocsms = (
                _lt_len_r2,
                _read_match_r2 + _rt_len_r2,
                0,
            )
            start_read.adhocseq = read.query_sequence

            return True, start_read
        self.logger.debug(
            "start read cannot connect with read and try to connect other reads",
        )
        return False, start_read  # not match

    @staticmethod
    def _double_check_for_start_end_read_determine_new_read_mode(
        read: Read,
        new_read_strand: str,
    ) -> int:
        """Double check for start and end read determine new read mode."""
        read.mode = 1 if read.mode == 2 else 2
        if new_read_strand == read.strand:
            return 1 if read.mode == 2 else 2
        return 1 if read.mode == 1 else 2

    def _double_check_create_new_read_calculate_sms(
        self,
        hsp: Any,
        query_seq: str,
        read: Read,
    ) -> Read:
        """Double check creat new read and calculate sms."""
        mapq = 60
        chrom, position, strand, cigar_str, num_of_mismatch = self.blat.psl2sam(
            hsp,
            len(query_seq),
        )

        lt_s_len = hsp.query_start
        rt_s_len = len(query_seq) - hsp.query_end
        new_read_mode = (
            ReadsConnector._double_check_for_start_end_read_determine_new_read_mode(
                read,
                strand,
            )
        )

        if new_read_mode == 1:
            cigar_str = (
                f"{lt_s_len}S"
                + cigar_str
                + f"{rt_s_len}S"
                + f"{read.read_match_size + read.rt_soft_len}S"
            )
            cigar_str = cigar_str[2:] if cigar_str.startswith("0S") else cigar_str
        else:
            cigar_str = (
                f"{read.lt_soft_len + read.read_match_size}S"
                f"{lt_s_len}S"
                f"{cigar_str}"
                f"{rt_s_len}S"
            )
            cigar_str = cigar_str[:-2] if cigar_str.endswith("0S") else cigar_str

        assert read.query_qualities is not None
        new_read = Read.new(
            read.query_name,
            chrom,
            position,
            strand,
            cigar_validity(cigar_str),
            mapq,
            num_of_mismatch,
            read.query_sequence,
            read.query_qualities,
        )

        new_read.mode = new_read_mode

        return new_read

    def __double_check_blat_query(
        self,
        query_sequence,
        align_len_threshold,
        threshold_identity,
        top,
    ):
        """Double check blat query."""
        flag = False

        if len(query_sequence) < align_len_threshold:
            return flag, None, None

        out_blat = self.blat.query(in_seq=query_sequence)
        try:
            blat_result = SearchIO.read(out_blat, "blat-psl")
        except ValueError:
            return flag, None, None

        hit, keep_hsp = self.blat._query_insertion(
            blat_result,
            query_sequence,
            threshold_identity,
            top=top,
        )
        return True, hit, keep_hsp

    def _double_check_for_start_end_read(
        self,
        read: Read,
        read_type: str,
    ) -> None:
        """Double check for start and end read.

        To see if there are True first read or True end read.
        """
        query_sequence = (
            read.query_sequence[: read.lt_soft_len]
            if read.mode == 1
            else read.query_sequence[len(read.query_sequence) - read.rt_soft_len :]
        )

        flag, hit, keep_hsp = self.__double_check_blat_query(
            query_sequence,
            self.align_len_threshold,
            self.threshold_identity,
            self.top,
        )

        assert keep_hsp is not None

        if flag and hit == 1:
            self.num_added_reads += 1
            hsp = keep_hsp[0]
            new_read = self._double_check_create_new_read_calculate_sms(
                hsp,
                query_sequence,
                read,
            )

            if read_type == "start":
                self.reads_chain.insert(0, new_read)
                self.read_pair_mode_dict[(new_read, read)] = (new_read.mode, read.mode)
            else:
                self.reads_chain.append(new_read)
                self.read_pair_mode_dict[(read, new_read)] = (read.mode, new_read.mode)

    @staticmethod
    def __sort_candidate_reads_key(read: Read, start_read: Read) -> int:
        """Sort candidate reads."""
        start_read_match_sequence = start_read.adhocsms[1]

        return min(
            abs(read.lt_soft_len - start_read_match_sequence),
            abs(read.rt_soft_len - start_read_match_sequence),
        )

    def connect(self) -> bool:
        """Find the best connected paths for a list of chimeric alignments.

        .. note::
            Read-to-Read chain scenarios
            * [[Read1, Read2, Read3]]
            * [[Read1, Read2, Read3],[Read4,Read5]]

            Dictionary of Read-pair scenarios
            * (Read1, Read2) => mode-of-Read1, mode-of-Read2
            * (Read2, Read1) => mode-of-Read2, mode-of-Read1
        """
        # find start node and end node
        temp_list = sorted(
            self.aln_list,
            key=lambda x: min(x.lt_soft_len, x.rt_soft_len),
        )
        start_nodes = temp_list[:2]
        self.candidate_nodes = temp_list[2:]

        start_read = start_nodes[0]
        end_read = start_nodes[1]

        if start_read.lt_soft_len > start_read.rt_soft_len:
            start_read.adhocsms = (
                start_read.lt_soft_len,
                start_read.rt_soft_len + start_read.read_match_size,
                0,
            )
        else:
            start_read.adhocsms = (
                0,
                start_read.lt_soft_len + start_read.read_match_size,
                start_read.rt_soft_len,
            )

        self.candidate_nodes.sort(
            key=lambda x: self.__sort_candidate_reads_key(x, start_read),
        )
        start_read.adhocseq = start_read.query_sequence

        self.reads_chain.append(start_read)
        flag = True
        if not self.candidate_nodes:  # []
            self.logger.debug("ReadsConnector: candidate_nodes is []")
            ReadsConnector.init_mode_judge(start_read, end_read)
            _, _ = self.test_2case(start_read, end_read, is_compare_for_ms=False)
            self._double_check_for_start_end_read(start_read, "start")
            self._double_check_for_start_end_read(end_read, "end")

        else:
            candidate_read_len = len(self.candidate_nodes)
            while self.candidate_nodes:
                if self.index == len(self.candidate_nodes):
                    logger.warning(
                        f"ReadsConnector: cannot connect all reads in candidate_nodes "
                        f"{start_read.query_name}",
                    )
                    return False
                read = self.candidate_nodes[self.index]
                ReadsConnector.init_mode_judge(start_read, read)
                flag, start_read = self.test_2case(
                    start_read,
                    read,
                    is_compare_for_ms=True,
                )

                if not flag:  # False
                    self.increment_index()

                if flag and len(self.candidate_nodes) + 1 == candidate_read_len:
                    self._double_check_for_start_end_read(start_nodes[0], "start")

            ReadsConnector.init_mode_judge(start_read, end_read)
            _, start_read = self.test_2case(
                start_read,
                end_read,
                is_compare_for_ms=True,
            )
            self._double_check_for_start_end_read(end_read, "end")

        return flag


def detect_read_read_connections_from_cigar(
    read: Any,
    mapq_cutoff: int,
    max_allowed_nm: int,
    blat: Blat,
    logger: LoggerType,
) -> Any:
    """Detecting read-read connections with chimeric alignments CIGAR string.

    :param logger:
    :param blat:
    :param mapq_cutoff: MAPQ cutoff
    :param max_allowed_nm: NM cutoff
    :type read: pysam.AlignedSegment object
    :type mapq_cutoff: int
    :type max_allowed_nm: int
    :return: Read-to-Read chain (a list of lists), a dictionary of Read-pair(Read1, Read2) =>
        mode-of-Read1, mode-of-Read2
    :rtype: tuple

    .. note::
        Read-to-Read chain scenarios
        * [[Read1, Read2, Read3]]
        * [[Read1, Read2, Read3],[Read4,Read5]]

        Dictionary of Read-pair scenarios
        * (Read1, Read2) => mode-of-Read1, mode-of-Read2
        * (Read2, Read1) => mode-of-Read2, mode-of-Read1

    .. important::
        If no 'SA' tag is found in this read, read-to-read chain and the read-pair =>
        mode dictionary will become empty.

    """

    def format_sa_tag(in_str: str) -> Any:
        """To keep read.reference_start and start position of SA alignment consistent.

        start position of SA alignment need to subtract 1

        :param in_str: string of supplementary read item in the SA tag
        :type in_str: str
        :return: chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa
        :rtype: tuple

        .. note::
             pos_sa, mapq_sa and nm_sa are integral variables now.
        """
        chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = in_str.split(",")
        pos_sa = int(pos_sa) - 1  # type: ignore
        mapq_sa = int(mapq_sa)  # type: ignore
        nm_sa = int(nm_sa)  # type: ignore
        return chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa

    def obtain_sa_query_seq_from_ra(
        query_seq_ra: str,
        strand_ra: str,
        strand_sa: str,
    ) -> str:
        """Helper function to define query_seq for the supplementary alignment.

        :param query_seq_ra: query sequence of representative alignment
        :param strand_ra: direction of representative read (-|+)
        :type strand_ra: str
        :param strand_sa: direction of supplementary read (-|+)
        :type strand_sa: str
        :return: query sequence of supplementary alignment
        :rtype: str
        """
        return (
            query_seq_ra if strand_ra == strand_sa else reverse_complement(query_seq_ra)
        )

    def mean(in_list: list[int]) -> float:
        """Helper function to calculate mean value of a list."""
        return sum(in_list) / len(in_list)

    def is_reverse_transcription_artifacts(
        read_list: list[Read],
        minimum_cutoff: int = 20,
        maximum_cutoff: int = 200,
        base_quality_cutoff: int = 10,
    ) -> bool:
        """Helper function to determine the reverse transcription (RT) artifacts.

        :param read_list: a list of chimeric read, the first one should be the representative alignment
        :param minimum_cutoff: minimum allowed distance between transcript and the corresponding RT artifact
        :param maximum_cutoff: maximum allowed distance between transcript and the corresponding RT artifact
        :param base_quality_cutoff: the difference in averaged base quality between transcript and RT artifact
        :return: the input read list contains reverse transcription artifacts
        """
        is_artifact = False
        for read1, read2 in combinations(read_list, 2):
            if read1.query_qualities is None or read2.query_qualities is None:
                mean_qualities_read1_match = 40.0
                mean_qualities_read2_match = 40.0
            else:
                mean_qualities_read1_match = mean(
                    read1.query_qualities[
                        read1.lt_soft_len : (read1.query_length - read1.rt_soft_len)
                    ],
                )
                mean_qualities_read2_match = mean(
                    read2.query_qualities[
                        read2.lt_soft_len : (read2.query_length - read2.rt_soft_len)
                    ],
                )
            if (
                read1.strand != read2.strand
                and read1.chrom == read2.chrom
                and (
                    (
                        abs(read1.ref_start - read2.ref_start) <= minimum_cutoff
                        or abs(read1.ref_end - read2.ref_end) <= minimum_cutoff
                    )
                    or (
                        (
                            minimum_cutoff
                            < abs(read1.ref_start - read2.ref_start)
                            < maximum_cutoff
                            or minimum_cutoff
                            < abs(read1.ref_end - read2.ref_end)
                            < maximum_cutoff
                        )
                        and (
                            abs(mean_qualities_read1_match - mean_qualities_read2_match)
                            > base_quality_cutoff
                        )
                    )
                )
            ):
                return True

        return is_artifact

    noreturn = [], {}, 0  # type: ignore

    # if no 'SA' tag was found, read-to-read chain will be empty
    try:
        chimeric_aln = read.get_tag("SA")[:-1].split(";")  # type: ignore
    except KeyError:
        return noreturn

    # chimeric alignments for a chimeric read
    # a chimeric read can have multiple chimeric alignments
    chimeric_aln_list = []
    mapq_list = []

    chrm_ra = read.reference_name
    pos_ra = read.reference_start
    strand_ra = "-" if read.is_reverse else "+"

    cigar_ra = read.cigarstring
    mapq_ra = read.mapping_quality
    nm_ra = read.get_tag("NM")
    seq_ra = read.query_sequence
    query_qualities_ra = read.query_qualities

    if (
        chrm_ra is None
        or read.query_name is None
        or seq_ra is None
        or nm_ra is None
        or cigar_ra is None
    ):
        msg = "None value found in read"
        raise ValueError(msg)

    # filter reads in uncommon chromosome and mitochondrion
    if "_" in chrm_ra or chrm_ra in {"chrM", "MT"}:
        return noreturn

    if nm_ra < max_allowed_nm:  # type: ignore
        mapq_list.append(mapq_ra)
        chimeric_aln_list.append(
            Read.new(
                read.query_name,
                chrm_ra,
                pos_ra,
                strand_ra,
                cigar_ra,
                mapq_ra,
                nm_ra,  # type: ignore
                seq_ra,
                query_qualities_ra,
            ),
        )

    for sa_string in chimeric_aln:
        chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = format_sa_tag(sa_string)
        seq_sa = obtain_sa_query_seq_from_ra(seq_ra, strand_ra, strand_sa)
        if strand_sa == strand_ra:
            query_qualities_sa = query_qualities_ra
        elif query_qualities_ra is None:
            query_qualities_sa = None
        else:
            query_qualities_sa = query_qualities_ra[::-1]

        # filter reads in uncommon chromosome and mitochondrion
        if "_" in chrm_sa or chrm_sa in {"chrM", "MT"}:
            return noreturn

        if nm_sa < max_allowed_nm:
            mapq_list.append(mapq_sa)
            assert query_qualities_sa is not None
            chimeric_aln_list.append(
                Read.new(
                    read.query_name,
                    chrm_sa,
                    pos_sa,
                    strand_sa,
                    cigar_sa,
                    mapq_sa,
                    nm_sa,
                    seq_sa,
                    query_qualities_sa,
                ),
            )

    if (len(chimeric_aln_list) < 1 + len(chimeric_aln)) or (
        min(mapq_list) < mapq_cutoff
    ):
        return noreturn

    if is_reverse_transcription_artifacts(chimeric_aln_list):
        logger.debug(f"{chimeric_aln_list=} has reverse transcription artifacts")
        return noreturn

    read_connector = ReadsConnector(
        read_list=chimeric_aln_list,
        blat=blat,
        logger=logger,
    )

    flag = read_connector.connect()

    if flag:
        logger.debug(
            f"reads chain: {read_connector.reads_chain};"
            f" reads pair mode: {read_connector.read_pair_mode_dict}",
        )
        return (
            read_connector.reads_chain,
            read_connector.read_pair_mode_dict,
            read_connector.num_added_reads,
        )

    return noreturn
