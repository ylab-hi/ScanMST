"""Basic Read Class.

@Filename:    basicRead.py
@Author:      YangyangLi
@contact:     yangyang.li@northwestern.edu
@Time:        1/9/22 12:13 PM
"""

from __future__ import annotations

from scannls import cppext

from .basic import CigarCode, Intervals, MappingMode, Strand


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
    :param cigar_tuples: cigarstring tuple version: [ (operation code, length) ];
        operation code: {'M':0,'I':1,'D':2,'N':3,'S':4,'H':5}
    :param cigartuples_without_soft: cigarstring tuple verion [exclude softclipping]:
        [(operation code, length)]; operation code: {'M':0,'I':1,'D':2,'N':3}
    :param query_length: length of the chimeric read

    :Example:

    >>> chrm_ra, pos_ra, strand_ra, cigar_ra, mapq_ra, nm_ra, seq_ra = ('chr1', 6524193,
    ...     '+', '5S10M2I5M10N10M15S', 60, 0, 'ATCGAAATTAGCTGGGTGTAGTGGCAGGTACCTATGGTCCTGGCTAC')
    >>> from array import array
    >>> query_qualities = array('B', [24, 23, 24, 24, 25, 27, 29, 31, 34, 33, 32, 34, 19,
            19, 19, 19, 33, 30, 29, 30, 30, 31, 28, 25, 26, 27, 14, 11, 29, 31, 34, 33, 32, 34, 19,
            29, 31, 34, 33, 32, 34, 19, 29, 31, 34, 33, 32])
    >>> read = Read.init(chrm_ra, pos_ra, strand_ra, cigar_ra, mapq_ra, nm_ra, seq_ra, query_qualities)
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
        "adhocseq",
        "adhocsms",
        "chrom",
        "cigarstring",
        "cigartuples_without_soft",
        "indel_size",
        "lt_soft_len",
        "mapq",
        "mode",
        "nm",
        "query_length",
        "query_name",
        "query_qualities",
        "query_sequence",
        "read_match_size",
        "ref_end",
        "ref_start",
        "reference_match_size",
        "rt_soft_len",
        "sms",
        "strand",
    )

    def __init__(
        self,
        query_name: str,
        chrom: str,
        ref_start: int,
        strand: str | Strand,
        cigarstring: str,
        mapq: int,
        nm: int,
        query_sequence: str,
        lt_soft_len: int,
        rt_soft_len: int,
        read_match_size: int,
        reference_match_size: int,
        indel_size: int,
        cigartuples_without_soft: list[int],
        query_length: int,
        query_qualities: list[int] | None = None,
    ) -> None:
        """Initialize a read class."""
        self.query_name = query_name
        self.chrom = chrom
        self.ref_start = ref_start
        self.strand = Strand.from_str(strand)
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
        self.query_qualities = query_qualities
        self.ref_end = self.ref_start + self.reference_match_size

        self.sms = self.lt_soft_len, self.read_match_size, self.rt_soft_len
        self.adhocsms = None
        self.adhocseq = None
        self.mode: MappingMode | None = None

    def __hash__(self) -> int:
        """Get the hash value of the read.

        :return: hash value of the read
        """
        return hash(self.chrom) ^ hash(self.ref_start) ^ hash(self.ref_end) ^ hash(self.strand) ^ hash(self.mapq) ^ hash(self.nm)

    def __repr__(self) -> str:
        """Get the representation of the read.

        :return: representation of the read
        """
        return (
            f"Read({self.query_name=}, {self.chrom=}, {self.ref_start=}, {self.ref_end=}, {self.sms=} "
            f"{self.mode=}, {self.strand=}, {self.mapq=}, {self.nm=})"
        )

    @classmethod
    def new(
        cls,
        query_name: str,
        chrom: str,
        ref_start: int,
        strand: str | Strand,
        cigar_str: str,
        mapq: int,
        nm: int,
        query_seq: str,
        query_qualities: list[int] | None,
    ) -> Read:
        """Calculate the features of the read and initialize the read."""
        parse_cigar_result = cppext.parseCigar(cigar_str)

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
            query_qualities,
        )

    def get_exons(self) -> Intervals:
        """Get the coordinates for reads matched part (without softclipping).

        :return: exons coordinates and introns coordinates
        :rtype: tuple
        """
        exons = []
        current_pos = self.ref_start
        start_pos = self.ref_start

        for ind in range(0, len(self.cigartuples_without_soft), 2):
            op_code = self.cigartuples_without_soft[ind]
            _len = self.cigartuples_without_soft[ind + 1]

            if op_code in {CigarCode.Match, CigarCode.Del}:  # M, D
                current_pos += _len
            elif op_code == CigarCode.Ref_skip:  # N
                exons.append((start_pos, current_pos))
                current_pos += _len
                start_pos = current_pos

        exons.append((start_pos, current_pos))
        return Intervals.from_list(exons)

    def splice_site_checker(self, genome_fasta, fraction_cutoff=0.6) -> bool:
        """Check whether the fraction of canonical splice site usage.

        in the read is bigger than 'fraction_cutoff' or not.
        :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
        :param fraction_cutoff: fraction of canonical splice sites used in the putative
            introns inferred from the CIGAR
        :type genome_fasta: pyfaidx.Fasta
        :type fraction_cutoff: float
        :return: using canonical splice sites OR not
        :rtype: bool
        """
        exons = self.get_exons()
        introns = exons.introns()

        if introns is None:
            return True

        intron_count = len(introns)
        can_count = 0
        can_sites = {"GT-AG", "GC-AG", "AT-AC"}
        for start, end in introns:
            donor_site = (
                genome_fasta[self.chrom][end - 2 : end].reverse.complement.seq
                if self.strand.is_reverse()
                else genome_fasta[self.chrom][start : start + 2].seq
            )
            acceptor_site = (
                genome_fasta[self.chrom][start : start + 2].reverse.complement.seq
                if self.strand.is_reverse()
                else genome_fasta[self.chrom][end - 2 : end].seq
            )

            if f"{donor_site}-{acceptor_site}" in can_sites:
                can_count += 1
        return can_count / intron_count >= fraction_cutoff
