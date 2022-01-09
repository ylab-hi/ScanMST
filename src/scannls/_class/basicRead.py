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
from typing import Tuple

import pyfaidx  # type: ignore


class Read:
    """Build a read class for storing information of every junction read.

    :param chrom: chromosome of genome
    :type chrom: str
    :param ref_start: start position of chimeric read
    :type ref_start: int
    :param strand: direction of chimeric read (-|+)
    :type strand: str
    :param cigarstring: cigar string of chimeric read (-|+)
    :type cigarstring: str
    :param mapq: MAPQ of chimeric read
    :type mapq: int
    :param nm: number of mismatches of chimeric read
    :type nm: int
    :param query_sequence: read sequence in the BAM file
    :type query_sequence: str
    :param linked_paths: linked paths for the read
    :type linked_paths: list
    :param lt_soft_len: softclipped segment length on the left side
    :type lt_soft_len: int
    :param rt_soft_len: softclipped segment length on the right side
    :type rt_soft_len: int
    :param read_match_size: M+I
    :type read_match_size: int
    :param reference_match_size: M+D+N
    :type reference_match_size: int
    :param indel_size: D+N-I
    :type indel_size: int
    :param cigartuples: cigarstring tuple version: [ (operation code, length) ];
        operation code: {'M':0,'I':1,'D':2,'N':3,'S':4,'H':5}
    :type cigartuples: list
    :param cigartuples_without_soft: cigarstring tuple verion [exclude softclipping]:
        [(operation code, length)]; operation code: {'M':0,'I':1,'D':2,'N':3}
    :type cigartuples_without_soft: list
    :param query_length: length of the chimeric read
    :type query_length: int

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
        "linked_paths",
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
        cigartuples_without_soft: List[Tuple[str, int]],
        query_length: int,
        cigartuples: Any,
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
        self.linked_paths = []  # type: Any
        self.lt_soft_len = lt_soft_len
        self.rt_soft_len = rt_soft_len
        self.read_match_size = read_match_size
        self.reference_match_size = reference_match_size
        self.indel_size = indel_size
        self.cigartuples_without_soft = cigartuples_without_soft
        self.query_length = query_length
        self.cigartuples = cigartuples
        self.ref_end = self.ref_start + self.reference_match_size

        self.sms = self.lt_soft_len, self.read_match_size, self.rt_soft_len
        self.adhocsms: Any = None
        self.adhocseq: Any = None
        self.mode: Any = None

    def __eq__(self, other: Any) -> bool:
        """Compare two reads.

        :param other: other read to compare
        :return: True if two reads are equal, False otherwise
        """
        if isinstance(other, Read) and (
            self.chrom == other.chrom
            and self.ref_start == other.ref_start
            and self.ref_end == other.ref_end
            and self.strand == other.strand
            and self.mapq == other.mapq
            and self.nm == other.nm
        ):
            return True
        return False

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

    def __lt__(self, other) -> bool:
        """Compare two reads.

        :param other:  other read to compare
        :return: True if self is smaller than other, False otherwise
        """
        if self.chrom == other.chorm:
            return bool(self.ref_start < other.ref_start)
        else:
            chrm_dict = {"chrM": 0, "MT": 0, "chrX": 23, "chrY": 24, "X": 23, "Y": 24}
            for i in range(1, 23):
                chrm_dict.update({f"chr{i}": i})
                chrm_dict.update({f"{i}": i})
            return bool(chrm_dict[self.chrom] < chrm_dict[other.chrom])

    def __repr__(self) -> str:
        """Get the representation of the read.

        :return: representation of the read
        """
        return (
            f"Read({self.chrom}, {self.ref_start}, {self.ref_end}, "
            f"{self.strand}, {self.mapq}, {self.nm})"
        )

    def __str__(self) -> str:
        """Get the string representation of the read.

        :return: string representation of the read
        """
        return ",".join(
            map(
                str,
                [
                    self.chrom,
                    self.ref_start,
                    self.ref_end,
                    self.strand,
                    self.mapq,
                    self.nm,
                    f"{self.lt_soft_len}:{self.read_match_size}:{self.rt_soft_len}",
                ],
            ),
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
        (
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigartuples_without_soft,
            query_length,
            cigartuples,
        ) = Read._calculate_features(cigar_str)

        return cls(
            query_name,
            chrom,
            ref_start,
            strand,
            cigar_str,
            mapq,
            nm,
            query_seq,
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigartuples_without_soft,
            query_length,
            cigartuples,
        )

    @property
    def reference_span(self) -> int:
        """M+N+D."""
        return self.reference_match_size

    def add_path(self, path) -> None:
        """Path is an instance of Path class."""
        self.linked_paths.append(path)

    def get_exons_and_introns(self) -> Any:
        """Get the coordinates for reads matched part (without softclipping).

        :return: exons coordinates and introns coordinates
        :rtype: tuple
        """
        cigartuples_without_soft = self.cigartuples_without_soft
        exons = []
        current_pos = self.ref_start
        start_pos = self.ref_start

        for op_code, _len_ in cigartuples_without_soft:

            if op_code in {0, 2}:  # M, D
                current_pos = current_pos + _len_
            elif op_code == 3:  # N
                exons.append([start_pos, current_pos])
                current_pos = current_pos + _len_
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

    def splice_site_checker(
        self, genome_fasta: pyfaidx.Fasta, fraction_cutoff=0
    ) -> bool:
        """Check whether the fraction of canonical splice site usage in read reference.

        Which matched part is bigger than 'fraction_cutoff' or not

        :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
        :param fraction_cutoff: fraction of canonical splice sites used in the putative
            introns inferred from the CIGAR
        :type genome_fasta: pyfaidx.Fasta
        :type fraction_cutoff: float
        :return: using canonical splice sites OR not
        :rtype: bool
        """
        exons, introns = self.get_exons_and_introns()

        if len(introns) == 0:
            return True

        intron_count = 0
        can_count = 0
        can_sites = {"GT-AG", "GC-AG", "AT-AC"}
        for start, end in introns:
            if end - start >= 10:
                intron_count += 1
                if self.strand == "-":
                    left_site = genome_fasta[self.chrom][
                        end - 2 : end
                    ].reverse.complement.seq
                    right_site = genome_fasta[self.chrom][
                        start : start + 2
                    ].reverse.complement.seq
                else:
                    left_site = genome_fasta[self.chrom][start : start + 2].seq
                    right_site = genome_fasta[self.chrom][end - 2 : end].seq
                if f"{left_site}-{right_site}" in can_sites:
                    can_count += 1
        return bool(can_count / intron_count >= fraction_cutoff)
