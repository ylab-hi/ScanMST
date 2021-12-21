# !/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Filename:    readConnector.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        12/15/21 2:14 PM
"""
from typing import Any
from typing import List

from Bio import SearchIO  # type: ignore
from loguru import logger
from loguru._logger import Logger
from pysam import AlignedSegment  # type: ignore

from ..draft.helper import cigar_validity  # type: ignore
from ..utils import reverse_complement
from .basicClass import Read
from .blat import Blat
from .exception import ReadNotConnectedError


class ReadsConnector(object):
    """the ReadsConnector class is used to connect the reads and identify the mode of the reads

    :param aln_list: the list of the alignment
    :param blat: :class: `class.Blat` for the BLAT search
    :param logger: :class: `loguru.logger` for logging

    :Example:

    >>> from loguru import  logger
    >>> aln_list = []
    >>> blat = Blat(ref_2bit='reference.2bit', logger= logger, port=88888, output_dir='/tmp')
    >>> readconnector = ReadsConnector(aln_list=aln_list, blat=blat, logger=logger)
    >>> readconnector.run()
    >>> readconnector.reads_chain
    [Read(chr1, 6524193, 6524850, +, 60, 8), Read(chr1, 6522473, 6522883, +, 60, 4)]
    >>> readconnector.read_pair_mode_dict
    {(Read(chr1, 6524193, 6524850, +, 60, 8), Read(chr1, 6522473, 6522883, +, 60, 4)): (1, 2)}
    """

    def __init__(
        self,
        aln_list: List[Read],
        blat: Blat,
        logger: Logger,
        align_len_threshold: int = 20,
        threshold_identity: float = 0.99,
        top: int = 3,
    ) -> None:

        self.candidate_nodes: List = []
        self.reads_chain: List = []
        self.read_pair_mode_dict, self.insertion_dict = {}, {}  # type: ignore
        self.aln_list = aln_list
        self.logger = logger
        self.blat = blat
        self.index = 0
        self.num_added_reads = 0

        self.align_len_threshold: int = align_len_threshold
        self.threshold_identity: float = threshold_identity
        self.top: int = top

    def reset_index(self):
        """ reset the index in order to fetch read in  candidate reads in new iteration """
        self.index = 0

    def increment_index(self):
        """ increment the index in order to fetch read in  candidate nodes """
        self.index += 1

    @staticmethod
    def _get_mode(sms: Any) -> int:
        """get the mode of the reads

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
        else:
            return 1

    @staticmethod
    def init_mode_judge(read1: Read, read2: Read) -> None:
        """initialize the mode of the reads"""
        sorted_by_s_length = sorted(
            [read1, read2], key=lambda x: min(x.lt_soft_len, x.rt_soft_len)
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

    def compare_ms(
        self,
        query_seq: str,
        target_seq: str,
        same_strand: bool,
        minimum_s_length: int = 30,
        minimum_terminal_length: int = 5,
    ) -> bool:
        """compare m of start read with s of read

        query_seq: M  target_seq: S

        :param minimum_terminal_length: the minimum length away from the terminal
        :param minimum_s_length: the minimum length of the s
        :param same_strand: whether the start read and read are on the same strand
        :param target_seq: the s of the read
        :param query_seq: the m of the start read

        """
        match_flag = False

        self.logger.trace(
            f"query length ={len(query_seq)} target length ={len(target_seq)}"
        )  # type: ignore

        if len(target_seq) <= minimum_s_length:
            return match_flag

        if not same_strand:
            target_seq = reverse_complement(target_seq)

        if query_seq in target_seq:
            temp_index = target_seq.index(query_seq)
            index = min(temp_index, len(target_seq) - len(query_seq) - temp_index)
            if index <= minimum_terminal_length:
                match_flag = True

        return match_flag

    @staticmethod
    def _determine_microhomology_len(
        read_match_sequence,
        read_query_length,
        prev_sms,
        next_sms,
        prev_read_mode,
        next_read_mode,
    ):
        """
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
        else:
            if next_read_mode == 2:
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
        read_match_sequence,
        read_query_sequence,
        prev_sms,
        next_sms,
        prev_read_mode,
        next_read_mode,
    ):
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
            elif prev_read_mode == 1:
                return read_match_sequence[:-microhomology_length]
        else:
            return read_match_sequence

    def test_2case(self, start_read: Read, read: Read, is_compare_for_ms: bool) -> Any:

        """

        :param start_read:
        :param read:
        :param is_compare_for_ms:
        :return:
        """
        compare_mode = None

        self.logger.debug(f"{start_read.mode=}, {read.mode=}")  # type: ignore
        self.logger.debug(f"{start_read.adhocsms=}, {read.sms=}")  # type: ignore

        if not is_compare_for_ms:  # one hop
            self.read_pair_mode_dict[(start_read, read)] = (start_read.mode, read.mode)
            self.reads_chain.append(read)
            return True, start_read, compare_mode

        _lt_len_r1, _read_match_r1, _rt_len_r1 = start_read.adhocsms  # type: ignore
        _lt_len_r2, _read_match_r2, _rt_len_r2 = read.sms
        read_query_sequence = read.query_sequence

        same_strand = True if start_read.adhocseq == read.query_sequence else False

        # first case
        self.logger.debug("testing first case M vs LS")

        next_read_mode = 2
        read_match_sequence = start_read.adhocseq[
            _lt_len_r1 : _lt_len_r1 + _read_match_r1
        ]  # type: ignore
        read_match_sequence = ReadsConnector.update_query_sequence(
            read_match_sequence,
            read_query_sequence,
            start_read.adhocsms,
            read.sms,
            start_read.mode,
            next_read_mode,
        )

        match_flag = self.compare_ms(
            read_match_sequence, read.query_sequence[:_lt_len_r2], same_strand
        )

        if match_flag:  # may same
            compare_mode = "ls"
            read.mode = 2
            self.logger.debug(f"{start_read.mode=}, {read.mode=}")
            self.read_pair_mode_dict[(start_read, read)] = (start_read.mode, read.mode)
            self.reads_chain.append(read)

            if read in self.candidate_nodes:
                self.candidate_nodes.remove(read)
                self.reset_index()

            start_read = read
            start_read.adhocsms = 0, _lt_len_r2 + _read_match_r2, _rt_len_r2  # type: ignore
            start_read.adhocseq = read.query_sequence

            return True, start_read, compare_mode

        self.logger.debug("testing second case M vs RS")  # type: ignore
        # second case
        next_read_mode = 1
        read_match_sequence = start_read.adhocseq[
            _lt_len_r1 : _lt_len_r1 + _read_match_r1
        ]  # type: ignore
        read_match_sequence = ReadsConnector.update_query_sequence(
            read_match_sequence,
            read_query_sequence,
            start_read.adhocsms,
            read.sms,
            start_read.mode,
            next_read_mode,
        )

        match_flag = self.compare_ms(
            read_match_sequence, read.query_sequence[-_rt_len_r2:], same_strand
        )

        if match_flag:
            compare_mode = "rs"

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

            return True, start_read, compare_mode
        self.logger.debug("start read cannot connect with read and try to connect other reads")  # type: ignore
        return False, start_read, compare_mode  # not match

    def _double_check_for_one_hop_for_end_read_add_new_read(
        self, read: Read, new_read: Read
    ) -> None:
        self.reads_chain.append(new_read)
        if new_read.strand == read.strand:
            new_read.mode = 1 if read.mode == 2 else 2
        else:
            new_read.mode = 1 if read.mode == 1 else 2
        self.read_pair_mode_dict[(read, new_read)] = (read.mode, new_read.mode)

    def _double_check_for_one_hop_for_start_read_add_new_read(
        self, read: Read, new_read: Read
    ) -> None:
        self.reads_chain.insert(0, new_read)
        if new_read.strand == read.strand:
            new_read.mode = 1 if read.mode == 2 else 2
        else:
            new_read.mode = 1 if read.mode == 1 else 2
        self.read_pair_mode_dict[(read, new_read)] = (read.mode, new_read.mode)

    def _double_check_for_one_hop_creat_new_read_and_calculate_sms(
        self, hsp: Any, query_seq: str, read: Read
    ) -> Read:

        mapq = 60
        chrom, position, strand, cigar_str, num_of_mismatch = self.blat.psl2sam(
            hsp, len(query_seq)
        )

        lt_s_len = hsp.query_start
        rt_s_len = len(query_seq) - hsp.query_end - 1

        if read.mode == 1:
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
                + f"{lt_s_len}S"
                + cigar_str
                + f"{rt_s_len}S"
            )
            cigar_str = cigar_str[:-2] if cigar_str.endswith("0S") else cigar_str

        new_read = Read.init(
            chrom,
            position,
            strand,
            cigar_validity(cigar_str),
            mapq,
            num_of_mismatch,
            read.query_sequence,
        )

        return new_read

    def __double_check_for_one_hop_blat_query(
        self, query_sequence, align_len_threshold, threshold_identity, top
    ):
        flag = False

        if len(query_sequence) < align_len_threshold:
            return flag, None, None

        out_blat = self.blat.query(in_seq=query_sequence)
        try:
            blat_result = SearchIO.read(out_blat, "blat-psl")
        except ValueError:
            return flag, None, None

        hit, keep_hsp = self.blat._query_insertion(
            blat_result, query_sequence, threshold_identity, top=top
        )
        return True, hit, keep_hsp

    def _double_check_for_one_hop_for_end_read(
        self,
        read: Read,
    ) -> None:

        query_sequence = (
            read.query_sequence[: read.lt_soft_len]
            if read.mode == 1
            else read.query_sequence[-read.rt_soft_len :]
        )

        flag, hit, keep_hsp = self.__double_check_for_one_hop_blat_query(
            query_sequence, self.align_len_threshold, self.threshold_identity, self.top
        )
        if flag and hit == 1:
            self.num_added_reads += 1
            hsp = keep_hsp[0]
            new_read = self._double_check_for_one_hop_creat_new_read_and_calculate_sms(
                hsp, query_sequence, read
            )
            self._double_check_for_one_hop_for_end_read_add_new_read(read, new_read)

    @staticmethod
    def _double_check_for_start_read_determine_s_source_for_blat(
        start_read, read, compare_mode: str
    ) -> str:

        if compare_mode == "ls":
            return "ls" if start_read.strand == read.strand else "rs"

        elif compare_mode == "rs":
            return "rs" if start_read.strand == read.strand else "ls"

    def _double_check_for_start_read(
        self, start_read: Read, read: Read, compare_mode: str
    ) -> None:
        s_source_blat = self._double_check_for_start_read_determine_s_source_for_blat(
            start_read, read, compare_mode
        )
        query_sequence = (
            start_read.query_sequence[: start_read.lt_soft_len]
            if s_source_blat == "ls"
            else start_read.query_sequence[-start_read.rt_soft_len :]
        )

        flag, hit, keep_hsp = self.__double_check_for_one_hop_blat_query(
            query_sequence, self.align_len_threshold, self.threshold_identity, self.top
        )

        if flag and hit == 1:
            self.num_added_reads += 1
            hsp = keep_hsp[0]
            new_read = self._double_check_for_one_hop_creat_new_read_and_calculate_sms(
                hsp, query_sequence, start_read
            )
            self._double_check_for_one_hop_for_start_read_add_new_read(
                start_read, new_read
            )

    def run(self) -> bool:
        """Find the best connected paths for a list of chimeric alignments
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
            self.aln_list, key=lambda x: min(x.lt_soft_len, x.rt_soft_len)
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

        start_read.adhocseq = start_read.query_sequence

        self.reads_chain.append(start_read)
        flag = True
        if not self.candidate_nodes:  # []

            self.logger.debug("ReadsConnector: candidate_nodes is []")
            ReadsConnector.init_mode_judge(start_read, end_read)
            _, _ = self.test_2case(start_read, end_read, is_compare_for_ms=False)
            self._double_check_for_one_hop_for_end_read(end_read)

        else:
            candidate_read_len = len(self.candidate_nodes)
            while self.candidate_nodes:
                if self.index == candidate_read_len:
                    logger.error(
                        "ReadsConnector: cannot connect all reads in candidate_nodes"
                    )
                    raise ReadNotConnectedError
                read = self.candidate_nodes[self.index]
                ReadsConnector.init_mode_judge(start_read, read)
                flag, start_read, compare_mode = self.test_2case(
                    start_read, read, is_compare_for_ms=True
                )

                if not flag:  # False
                    self.increment_index()

                if flag and len(self.candidate_nodes) + 1 == candidate_read_len:
                    self._double_check_for_start_read(
                        start_nodes[0], start_read, compare_mode
                    )

            ReadsConnector.init_mode_judge(start_read, end_read)
            _, start_read, compare_mode = self.test_2case(
                start_read, end_read, is_compare_for_ms=True
            )
            self._double_check_for_one_hop_for_end_read(end_read)

        return flag


def detect_read_read_connections_from_cigar(
    read: AlignedSegment,
    mapq_cutoff: int,
    max_allowed_nm: int,
    blat: Blat,
    logger: Logger,
) -> Any:
    """Detecting read-read connections with chimeric alignments CIGAR string

    :param logger:
    :param blat:
    :param mapq_cutoff: MAPQ cutoff
    :param max_allowed_nm: NM cutoff
    :type read: pysam.AlignedSegment object
    :type mapq_cutoff: int
    :type max_allowed_nm: int
    :return: Read-to-Read chain (a list of lists), a dictionary of Read-pair(Read1, Read2) => mode-of-Read1, mode-of-Read2
    :rtype: tuple
    .. note::
        Read-to-Read chain scenarios
        * [[Read1, Read2, Read3]]
        * [[Read1, Read2, Read3],[Read4,Read5]]

        Dictionary of Read-pair scenarios
        * (Read1, Read2) => mode-of-Read1, mode-of-Read2
        * (Read2, Read1) => mode-of-Read2, mode-of-Read1

    .. important::
        If no 'SA' tag is found in this read, read-to-read chain and the read-pair => mode dictionary will become empty.

    #return: NLS_type(TDUP/INV), exon_boundary(0/1/2/3), canonical_or_not (1/0), [position, size, rep_aln_mode, sup_aln_mode], [++]
    #        TRA, canonical_or_not (1/0), [position, sup_position, rep_aln_mode, sup_aln_mode], [+-]
    #        e.g., INV,1,43947377,181934993,1,1,++
    #              TRA,1,160289623,chr17:17189212,1,1,+-
    """

    def format_sa_tag(in_str: str) -> Any:
        """
        To keep read.reference_start and start position of SA alignment consistent, start position of SA alignment need to substract 1
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
        query_seq_ra: str, strand_ra: str, strand_sa: str
    ) -> str:
        """a helper function to define query_seq for the supplementary alignment
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

    if read.has_tag("SV"):
        return [], {}

    if read.is_supplementary:
        return [], {}

    # if no 'SA' tag was found, read-to-read chain will be empty
    try:
        chimeric_aln = read.get_tag("SA")[:-1].split(";")
    except KeyError:
        return [], {}

    # chimeric alignments for a chimeric read
    # a chimeric read can have multiple chimeric alignments
    chimeric_aln_list = []

    chrm_ra = read.reference_name
    pos_ra = read.reference_start
    strand_ra = "-" if read.is_reverse else "+"

    cigar_ra = read.cigarstring
    mapq_ra = read.mapping_quality
    nm_ra = read.get_tag("NM")
    seq_ra = read.query_sequence

    if mapq_ra > mapq_cutoff and nm_ra < max_allowed_nm:
        chimeric_aln_list.append(
            Read.init(chrm_ra, pos_ra, strand_ra, cigar_ra, mapq_ra, nm_ra, seq_ra)
        )

    for sa_string in chimeric_aln:
        chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = format_sa_tag(sa_string)
        seq_sa = obtain_sa_query_seq_from_ra(seq_ra, strand_ra, strand_sa)
        if mapq_sa > mapq_cutoff and nm_sa < max_allowed_nm:
            chimeric_aln_list.append(
                Read.init(chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa, seq_sa)
            )

    if len(chimeric_aln_list) < 1 + len(chimeric_aln):
        return [], {}
    else:
        read_connector = ReadsConnector(
            aln_list=chimeric_aln_list, blat=blat, logger=logger
        )
        flag = read_connector.run()
        if flag:
            logger.debug(
                f"reads chain: {read_connector.reads_chain};"
                f" reads pair mode: {read_connector.read_pair_mode_dict}"
            )
            return (
                read_connector.reads_chain,
                read_connector.read_pair_mode_dict,
                read_connector.num_added_reads,
            )
        else:
            return [], {}
