# !/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Filename:    readConnecter.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        12/15/21 2:14 PM
"""
from typing import Any
from typing import List

from loguru import logger
from loguru._logger import Logger
from pysam import AlignedSegment  # type: ignore

from ..utils import reverse_complement
from .basicClass import Read
from .blat import Blat
from .exception import ReadNotConnectedError


class ReadsConnecter(object):
    """the ReadsConnecter class is used to connect the reads and identify the mode of the reads

    :param aln_list: the list of the alignment
    :param blat: `class.Blat` for the BLAT search
    :param logger: `loguru.logger` for logging

    :Example:

    >>> from loguru import  logger
    >>> aln_list = []
    >>> blat = Blat(ref_2bit='reference.2bit', logger= logger, port=88888, output_dir='/tmp')
    >>> readconnecter = ReadsConnecter(aln_list=aln_list, blat=blat, logger=logger)
    >>> readconnecter.run()
    >>> readconnecter.reads_chain
    [Read(chr1, 6524193, 6524850, +, 60, 8), Read(chr1, 6522473, 6522883, +, 60, 4)]
    >>> readconnecter.read_pair_mode_dict
    {(Read(chr1, 6524193, 6524850, +, 60, 8), Read(chr1, 6522473, 6522883, +, 60, 4)): (1, 2)}
    """

    def __init__(
        self,
        aln_list: List[Read],
        blat: Blat,
        logger: Logger,
    ) -> None:

        self.candidate_nodes: List = []
        self.reads_chain: List = []
        self.read_pair_mode_dict, self.insertion_dict = {}, {}  # type: ignore
        self.aln_list = aln_list
        self.logger = logger
        self.blat = blat
        self.index = 0

    def reset_index(self):
        """ reset the index in order to fetch read in  candidate reads in new iteration """
        self.index = 0

    def increment_index(self):
        """ increment the index in order to fetch read in  candidate nodes """
        self.index += 1

    @staticmethod
    def init_mode_judge(sms: Any) -> int:
        """initialize the mode of the reads

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

    # @staticmethod
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
        :param is_compare_by_in: whether to compare m and s by in
        :param same_strand: whether the start read and read are on the same strand
        :param target_seq: the s of the read
        :param query_seq: the m of the start read

        .. note::

            `is_compare_by_in=False` means candidate reads is empty -> two hop

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
        ) = ReadsConnecter._determine_microhomology_len(
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

    def test_2case(self, start_read: Read, read: Read, is_align_for_ms: bool) -> Any:

        """

        :param start_read:
        :param read:
        :param is_align_for_ms:
        :return:
        """
        _lt_len_r1, _read_match_r1, _rt_len_r1 = start_read.adhocsms  # type: ignore
        _lt_len_r2, _read_match_r2, _rt_len_r2 = read.sms
        read_query_sequence = read.query_sequence

        self.logger.debug(f"{start_read.mode=}, {read.mode=}")  # type: ignore
        self.logger.debug(f"{start_read.adhocsms=}, {read.sms=}")  # type: ignore

        if not is_align_for_ms:
            self.read_pair_mode_dict[(start_read, read)] = (start_read.mode, read.mode)
            self.reads_chain.append(read)
            return True, start_read

        same_strand = True if start_read.adhocseq == read.query_sequence else False

        # first case
        self.logger.debug("testing first case M vs LS")

        next_read_mode = 2
        read_match_sequence = start_read.adhocseq[
            _lt_len_r1 : _lt_len_r1 + _read_match_r1
        ]  # type: ignore
        read_match_sequence = ReadsConnecter.update_query_sequence(
            read_match_sequence,
            read_query_sequence,
            start_read.adhocsms,
            read.sms,
            start_read.mode,
            next_read_mode,
        )

        match_flag = self.compare_ms(
            read_match_sequence,
            read.query_sequence[:_lt_len_r2],
            same_strand,
            is_align_for_ms,
        )

        if match_flag:  # may same

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

            return True, start_read

        self.logger.debug("testing second case M vs RS")  # type: ignore
        # second case
        next_read_mode = 1
        read_match_sequence = start_read.adhocseq[
            _lt_len_r1 : _lt_len_r1 + _read_match_r1
        ]  # type: ignore
        read_match_sequence = ReadsConnecter.update_query_sequence(
            read_match_sequence,
            read_query_sequence,
            start_read.adhocsms,
            read.sms,
            start_read.mode,
            next_read_mode,
        )

        match_flag = self.compare_ms(
            read_match_sequence,
            read.query_sequence[-_rt_len_r2:],
            same_strand,
            is_align_for_ms,
        )

        if match_flag:

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
        self.logger.debug("start read cannot connect with read")  # type: ignore
        return False, start_read  # not match

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

            self.logger.debug("ReadsConnecter: candidate_nodes is []")
            start_read.mode, end_read.mode = (
                ReadsConnecter.init_mode_judge(start_read.adhocsms),
                ReadsConnecter.init_mode_judge(end_read.sms),
            )

            flag, _ = self.test_2case(start_read, end_read, is_align_for_ms=False)

        else:
            candidate_read_len = len(self.candidate_nodes)
            while self.candidate_nodes:
                if self.index == candidate_read_len:
                    logger.error(
                        "ReadsConnecter: cannot connect all reads in candidate_nodes"
                    )
                    raise ReadNotConnectedError
                read = self.candidate_nodes[self.index]
                start_read.mode, read.mode = (
                    ReadsConnecter.init_mode_judge(start_read.adhocsms),
                    ReadsConnecter.init_mode_judge(read.sms),
                )
                flag, start_read = self.test_2case(
                    start_read, read, is_align_for_ms=True
                )

                if not flag:  # False
                    self.increment_index()

            start_read.mode, end_read.mode = (
                ReadsConnecter.init_mode_judge(start_read.adhocsms),
                ReadsConnecter.init_mode_judge(end_read.sms),
            )
            _, start_read = self.test_2case(start_read, end_read, is_align_for_ms=True)

        return flag


def detect_read_read_connections_from_cigar(
    read: AlignedSegment, mapq_cutoff: int, blat: Blat, logger: Logger
) -> Any:
    """Detecting read-read connections with chimeric alignments CIGAR string

    :param logger:
    :param blat:
    :param mapq_cutoff: MAPQ cutoff
    :type read: pysam.AlignedSegment object
    :type mapq_cutoff: int
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
    if read.is_reverse:
        strand_ra = "-"
    else:
        strand_ra = "+"
    cigar_ra = read.cigarstring
    mapq_ra = read.mapping_quality
    nm_ra = read.get_tag("NM")
    seq_ra = read.query_sequence

    if mapq_ra > mapq_cutoff:
        chimeric_aln_list.append(
            Read.init(chrm_ra, pos_ra, strand_ra, cigar_ra, mapq_ra, nm_ra, seq_ra)
        )

    for sa_string in chimeric_aln:
        chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = format_sa_tag(sa_string)
        seq_sa = obtain_sa_query_seq_from_ra(seq_ra, strand_ra, strand_sa)
        if mapq_sa > mapq_cutoff:
            chimeric_aln_list.append(
                Read.init(chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa, seq_sa)
            )

    if len(chimeric_aln_list) < 1 + len(chimeric_aln):
        return [], {}
    else:

        read_connecter = ReadsConnecter(
            aln_list=chimeric_aln_list, blat=blat, logger=logger
        )
        flag = read_connecter.run()
        if flag:
            logger.debug(
                f"reads chain: {read_connecter.reads_chain};"
                f" reads pair mode: {read_connecter.read_pair_mode_dict}"
            )
            return (
                read_connecter.reads_chain,
                read_connecter.read_pair_mode_dict,
            )
        else:
            return [], {}
