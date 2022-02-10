# !/usr/bin/env python
"""Basic Read Class.

@Filename:    basicRead.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        1/9/22 12:13 PM
"""
import re
from typing import Any
from typing import List

from scannls import cppext


class Read:
    """Build a read class for storing information of every junction read.

    :param chrom: chromosome of genome
    :param ref_start: start position of chimeric read
    :param strand: direction of chimeric read (-|+)
    :param cigarstring: cigar string of chimeric read (-|+)
    :param mapq: MAPQ of chimeric read
    :param nm: number of mismatches of chimeric read
    :param query_sequence: read sequence in the BAM file
    :param linked_paths: linked paths for the read
    :param lt_soft_len: softclipped segment length on the left side
    :param rt_soft_len: softclipped segment length on the right side
    :param read_match_size: M+I
    :param reference_match_size: M+D+N
    :param indel_size: D+N-I
    :param cigartuples: cigarstring tuple version: [ (operation code, length) ];
        operation code: {'M':0,'I':1,'D':2,'N':3,'S':4,'H':5}
    :param cigartuples_without_soft: cigarstring tuple verion [exclude softclipping]:
        [(operation code, length)]; operation code: {'M':0,'I':1,'D':2,'N':3}
    :param query_length: length of the chimeric read

    :Example:

    >>> chrm_ra, pos_ra, strand_ra, cigar_ra, mapq_ra, nm_ra, seq_ra = ('chr1', 6524193,
    ...     '+', '5S10M2I5M10N10M15S', 60, 0, 'ATCGAAATTAGCTGGGTGTAGTGGCAGGTACCTATGGTCCTGGCTAC')
    >>> read = Read.init(chrm_ra, pos_ra, strand_ra, cigar_ra, mapq_ra, nm_ra, seq_ra)
    >>> read
    Read(chr1, 6524193, 6524213, +, 60, 0)
    >>> read.read_match_size
    27
    >>> read.reference_match_size
    35
    >>> read.sms
    5, 27, 15
    """

    __slots__ = (
        "chrom",
        "ref_start",
        "strand",
        "cigarstring",
        "mapq",
        "nm",
        "query_sequence",
        "query_name",
        "lt_soft_len",
        "rt_soft_len",
        "read_match_size",
        "reference_match_size",
        "indel_size",
        "cigartuples_without_soft",
        "cigartuples",
        "query_length",
        "adhocsms",
        "adhocseq",
        "mode",
        "sms",
        "ref_end",
    )

    def __init__(
        self,
        query_name: str,
        chrom: str,
        ref_start: int,
        strand: str,
        cigarstring: str,
        mapq: int,
        nm: int,
        query_sequence: str,
        lt_soft_len: int,
        rt_soft_len: int,
        read_match_size: int,
        reference_match_size: int,
        indel_size: int,
        cigartuples_without_soft: List[int],
        query_length: int,
    ) -> None:
        """Initialize a read class."""
        self.query_name = query_name
        self.chrom = chrom
        self.ref_start = ref_start
        self.strand = strand
        self.cigarstring = cigarstring
        self.mapq = mapq
        self.nm = nm
        self.query_sequence = query_sequence
        self.lt_soft_len = lt_soft_len
        self.rt_soft_len = rt_soft_len
        self.read_match_size = read_match_size
        self.reference_match_size = reference_match_size
        self.indel_size = indel_size
        self.cigartuples_without_soft = cigartuples_without_soft
        self.query_length = query_length
        self.ref_end = self.ref_start + self.reference_match_size

        self.sms = self.lt_soft_len, self.read_match_size, self.rt_soft_len
        self.adhocsms: Any = None
        self.adhocseq: Any = None
        self.mode: Any = None

    def __hash__(self) -> int:
        """Get the hash value of the read.

        :return: hash value of the read
        """
        return (
            hash(self.chrom)
            ^ hash(self.ref_start)
            ^ hash(self.ref_end)
            ^ hash(self.strand)
            ^ hash(self.mapq)
            ^ hash(self.nm)
        )

    def __repr__(self) -> str:
        """Get the representation of the read.

        :return: representation of the read
        """
        return (
            f"Read({self.chrom}, {self.ref_start}, {self.ref_end}, "
            f"{self.strand}, {self.mapq}, {self.nm})"
        )

    @staticmethod
    def _calculate_features(cigar_str: str) -> Any:
        """Calculate the features of the read.

        :param cigar_str: cigar string of the read
        """
        cigar_char_dict = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5}
        # 'length', 'operation char'
        len_type_tuple = re.findall(r"(\d+)(\w)", cigar_str)
        # (operation code, length)
        cigartuples = [(cigar_char_dict[j], int(i)) for i, j in len_type_tuple]

        query_length = 0
        indel_size = 0
        reference_match_size = 0
        read_match_size = 0
        cigartuples_without_soft = []

        for op_code, _len_ in cigartuples:
            if op_code == 0:  # M
                reference_match_size += _len_
                read_match_size += _len_
                query_length += _len_
                cigartuples_without_soft.append([0, _len_])
            elif op_code == 1:  # I
                indel_size += -_len_
                read_match_size += _len_
                query_length += _len_
                cigartuples_without_soft.append([1, _len_])
            elif op_code == 2:  # D
                indel_size += _len_
                reference_match_size += _len_
                cigartuples_without_soft.append([2, _len_])
            elif op_code == 3:  # N
                indel_size += _len_
                reference_match_size += _len_
                cigartuples_without_soft.append([3, _len_])
            elif op_code == 4:  # S
                query_length += _len_

        lt_soft_len = 0
        rt_soft_len = 0
        lt_op, lt_len = cigartuples[0]
        rt_op, rt_len = cigartuples[-1]
        if lt_op == 4:
            lt_soft_len = lt_len
        if rt_op == 4:
            rt_soft_len = rt_len
        return (
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigartuples_without_soft,
            query_length,
            cigartuples,
        )

    @classmethod
    def init(
        cls,
        query_name: str,
        chrom: str,
        ref_start: int,
        strand: str,
        cigar_str: str,
        mapq: int,
        nm: int,
        query_seq: str,
    ) -> "Read":
        """Calculate the features of the read and initialize the read."""
        parse_cigar_result = cppext.parse_cigar(cigar_str)

        return cls(
            query_name,
            chrom,
            ref_start,
            strand,
            cigar_str,
            mapq,
            nm,
            query_seq,
            parse_cigar_result.lt_soft_len,
            parse_cigar_result.rt_soft_len,
            parse_cigar_result.read_match,
            parse_cigar_result.ref_match,
            parse_cigar_result.indel_len,
            parse_cigar_result.cigartuples_without_soft,
            parse_cigar_result.query_len,
        )

    def get_exons_and_introns(self) -> Any:
        """Get the coordinates for reads matched part (without softclipping).

        :return: exons coordinates and introns coordinates
        :rtype: tuple
        """
        exons = []
        current_pos = self.ref_start
        start_pos = self.ref_start

        for ind in range(len(self.cigartuples_without_soft), 2):
            op_code = self.cigartuples_without_soft[ind]
            _len = self.cigartuples_without_soft[ind + 1]

            if op_code in {0, 2}:  # M, D
                current_pos = current_pos + _len
            elif op_code == 3:  # N
                exons.append([start_pos, current_pos])
                current_pos = current_pos + _len
                start_pos = current_pos

        exons.append([start_pos, current_pos])

        introns = []
        # No 'N' in the cigar
        if len(exons) > 1:
            _positions = []
            for i, j in exons:
                _positions.extend([i, j])
            _positions.sort()
            _positions.pop(0)
            _positions.pop(-1)
            introns = [[x, y] for x, y in zip(_positions[::2], _positions[1::2])]
        return exons, introns
