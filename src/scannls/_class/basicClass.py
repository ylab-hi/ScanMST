#!/usr/bin/env python3
import re
from typing import Any
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

from loguru._logger import Logger  # type: ignore

from ..draft.helper import cigar_validity  # type: ignore
from ..draft.nls_inference import infer_nls_from_connected_reads  # type: ignore
from ..utils import reverse_complement  # type: ignore
from .exception import ReadNotFoundError  # type: ignore


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
        chrom,
        ref_start,
        strand,
        cigarstring,
        mapq,
        nm,
        query_sequence,
        lt_soft_len,
        rt_soft_len,
        read_match_size,
        reference_match_size,
        indel_size,
        cigartuples_without_soft,
        query_length,
        cigartuples,
    ) -> None:
        """Initialize a read class."""
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

    def __eq__(self, other) -> bool:
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
            return self.ref_start < other.ref_start
        else:
            chrm_dict = {"chrM": 0, "MT": 0, "chrX": 23, "chrY": 24, "X": 23, "Y": 24}
            for i in range(1, 23):
                chrm_dict.update({f"chr{i}": i})
                chrm_dict.update({f"{i}": i})
            return chrm_dict[self.chrom] < chrm_dict[other.chrom]

    def __repr__(self) -> str:
        """Get the representation of the read.

        :return: representation of the read
        """
        return (
            f"Read({self.chrom}, {self.ref_start}, {self.ref_end}"
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
    def _calculate_features(cigar_str):
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
    def init(cls, chrom, ref_start, strand, cigar_str, mapq, nm, query_seq):
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
        exons = []
        current_pos = self.ref_start
        start_pos = self.ref_start
        for op_code, _len_ in self.cigartuples_without_soft:
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

    def splice_site_checker(self, genome_fasta, fraction_cutoff=0) -> bool:
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
        if can_count / intron_count >= fraction_cutoff:
            return True
        else:
            return False


class NovelInsertion(Read):
    """NovelInsertion is used to represent reads insertion whose hit is 0 or >1.

    :Example:

    >>> novel_insertion = NovelInsertion(hit_num=0, query_sequence='ATCA')
    >>> novel_insertion
    NovelInsertion(ATCA:0)
    >>> novel_insertion.query_sequence
    ATCA
    >>> novel_insertion.hit_num
    0

    .. note::
        `NovelInsertion` is a subclass of :class:`Read`, and siblings of :class:`Insertion`

    .. seealso::
        :class:`Insertion`
    """

    def __init__(self, hit_num: int, query_sequence: str):
        """Initialize NovelInsertion."""
        self.query_sequence = query_sequence
        self.hit_num = hit_num
        self.insertion_info = None

    def __repr__(self):
        """Represent NovelInsertion object."""
        return f"NovelInsertion({self.query_sequence}:{self.hit_num})"

    def reverse_completement_query(self):
        """Reverse complement query sequence."""
        self.query_sequence = reverse_complement(self.query_sequence)


class MicroHomology:
    """MicroHomology is used to represent microhomology.

    :Example:

    >>> microhomology = MicroHomology(query_sequence="ATCA")
    >>> microhomology
    MicroHomology(ATCA)
    >>> microhomology.query_sequence
    ATCA

    .. seealso::
        :class:`Insertion` and :class:`NovelInsertion`
    """

    def __init__(self, query_sequence: str):
        """Initialize MicroHomology."""
        self.query_sequence = query_sequence

    def __repr__(self):
        """Represent MicroHomology object."""
        return f"MicroHomology({self.query_sequence})"

    def reverse_completement_query(self):
        """Reverse complement query sequence."""
        self.query_sequence = reverse_complement(self.query_sequence)


class Insertion(Read):
    """Insertion is used to represent reads insertion whose hit is 1.

    :param chrom: chromosome of genome
    :param ref_start: start position of chimeric read
    :param strand: direction of chimeric read (-|+)
    :param cigarstring: cigar string of chimeric read (-|+)
    :param mapq: MAPQ of chimeric read
    :param nm: number of mismatches of chimeric read
    :param query_sequence: read sequence in the BAM file

    :Example:

    >>> insertion = Insertion(hit_num=1, chrom= '1', ref_start=1, strand='+',
    ...                 cigarstring='1S1M1S',mapq=60, nm=0, query_sequence='ATCA')
    >>> insertion
    Insertion(1:1-4:+, 1-2|2-3, TPA, 1, 4)

    .. note::
        `Insertion` is a subclass of :class:`Read`, and siblings of :class:`NovelInsertion`
        `Insertion` includes the attributes of :class:`Node` in order to enable us to
        manipulate the attributes of `Insertion` same as :class:`Node` in :class:`Series`

    .. seealso:: :class:`NovelInsertion`, :class:`Node` and :class:`Read`
    """

    def __init__(
        self,
        hit_num: int,
        chrom: str,
        ref_start: int,
        strand: str,
        cigarstring: str,
        mapq: int,
        nm: int,
        query_sequence: str,
    ):
        """Initialize Insertion."""
        (
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigartuples_without_soft,
            query_length,
            cigartuples,
        ) = Read._calculate_features(cigarstring)
        super().__init__(
            chrom,
            ref_start,
            strand,
            cigarstring,
            mapq,
            nm,
            query_sequence,
            lt_soft_len,
            rt_soft_len,
            read_match_size,
            reference_match_size,
            indel_size,
            cigartuples_without_soft,
            query_length,
            cigartuples,
        )
        self.hit_num = hit_num
        self.sv_type = None

        # add attributes for insertion in order to be compatible with the class Node
        self.prev_breakpoint: Optional[str] = None
        self.next_breakpoint: Optional[str] = None
        self.modes = None
        self.genes = None
        self.annotation_code = None
        self.splicing_code = None
        self.sr = 1
        self.insertion_info = None

        self.exons, self.introns = self.get_exons_and_introns()
        self.successors: Any = []
        self.predecessors: Any = []
        self.unique_key = None
        self.merged_child_nodes: List = []
        self.merged_parent_nodes: List = []

        self.next_node_in_series = None
        self.previous_node_in_series = None
        self.is_merged, self.is_in_graph = False, False

    def __repr__(self):
        """Represent Insertion object."""
        exons_repr = "|".join([f"{i}-{j}" for i, j in self.exons])
        return (
            f"Insertion({self.chrom}:{self.ref_start}-{self.ref_end}:{self.strand}"
            f"{exons_repr}, {self.sv_type}, {self.prev_breakpoint}, "
            f"{self.next_breakpoint}, SR={self.sr})"
        )

    def __hash__(self) -> int:
        """Hash Insertion object."""
        return (
            hash(self.chrom)
            ^ hash(self.ref_start)
            ^ hash(self.ref_end)
            ^ hash(self.sv_type)
            ^ hash(self.prev_breakpoint)
            ^ hash(self.next_breakpoint)
            ^ hash(self.strand)
        )

    def update_cigarstring_sms(self, sms, source_s, source_strand):
        """Update cigarstring and sms of Insertion object."""
        _ls, _m, _rs = sms
        if source_s == "left":
            ls = _ls - self.query_length
            rs = _rs + _m
        else:
            ls = _ls + _m
            rs = _rs - self.query_length

        if source_strand != self.strand:
            rs, ls = ls, rs

        self.lt_soft_len = ls
        self.rt_soft_len = rs
        self.sms = (ls, self.query_length, rs)
        self.query_length = self.query_length + ls + rs
        self.cigarstring = cigar_validity(f"{ls}S{self.cigarstring}{rs}S")

    def reverse_completement_query(self):
        """Reverse complement query sequence of Insertion object."""
        self.query_sequence = reverse_complement(self.query_sequence)

    def reverse_strand(self):
        """Reverse strand of Insertion object."""
        self.strand = "-" if self.strand == "+" else "+"

    @property
    def similar_key(self):
        """Return similar key of Insertion object.

        .. note::
            similar key is consisted of chrom and introns info.
        """
        introns = self.introns

        key = "-".join([f"{i - j}" for i, j in introns]) if introns else "None"
        key = f"{self.chrom}-{key}"

        return key

    def get_unique_key(self):
        """Return unique key of Insertion object.

        .. note::
            unique key is consisted of chrom, introns, and breakpoints info
        """
        introns = self.introns

        key = "-".join([f"{i - j}" for i, j in introns]) if introns else "None"

        key = f"{self.chrom}-{key}-{self.prev_breakpoint}-{self.next_breakpoint}"

        self.unique_key = key

        return key

    def is_start_node(self):
        """Return True if Insertion object is start node."""
        return True if not self.has_predecessor() else False

    def is_end_node(self):
        """Return True if Insertion object is end node."""
        return True if not self.has_successor() else False

    def has_predecessor(self):
        """Return True if Insertion object has predecessor."""
        return True if self.predecessors else False

    def has_successor(self):
        """Return True if Insertion object has successor."""
        return True if self.successors else False

    def add_successor_from_list(self, successors):
        """Add successor from list of Insertion object."""
        for successor in successors:
            self.add_successor(successor)

    def add_predecessor_from_list(self, predecessors):
        """Add predecessor from list of Insertion object."""
        for predecessor in predecessors:
            self.add_predecessor(predecessor)

    def _add_successor(self, successor):
        """Helper function to add successor to Insertion object."""
        self.successors.append(successor)
        successor.add_predecessor(self)

    def _add_predecessor(self, predecessor):
        """Helper function to add predecessor to Insertion object."""
        self.predecessors.append(predecessor)
        predecessor.add_successor(self)

    def add_successor(self, successor):
        """Add successor to Insertion object."""
        if successor is not None and successor not in self.successors:
            if successor.is_in_graph:
                self._add_successor(successor)
            else:
                self.add_successor_from_list(successor.merged_parent_nodes)

    def add_predecessor(self, predecessor):
        """Node must be in the graph if the function is called.

        :param predecessor: predecessor of Insertion object
        """
        if predecessor is not None and predecessor not in self.predecessors:
            if predecessor.is_in_graph:
                self._add_predecessor(predecessor)
            else:
                self.add_predecessor_from_list(predecessor.merged_parent_nodes)

    def update_sr(self, key=1):
        """Update sr."""
        self.sr += key

    def update_next_and_previous_node_in_series(self, index, series):
        """Update next and previous node in series."""
        if index == 0:
            self.next_node_in_series = series[index + 1]
        elif index == len(series) - 1:
            self.previous_node_in_series = series[index - 1]
        else:
            self.next_node_in_series = series[index + 1]
            self.previous_node_in_series = series[index - 1]


class Node:
    """Build a breakpoint node class for storing information of every breakpoint.

    :param prev_breakpoint: breakpoint for the previous breakpoints connections
    :param next_bp: breakpoint for the next breakpoints connections
    :param strand: direction of chimeric read (-|+)
    :param chrom: chromosome
    :param ref_start: reference start position
    :param ref_end: reference end position
    :param exons: CIGAR inferred exons in the read. e.g., [(100, 200), (300, 500)]
    :param sv_type: one of the SV types (TDUP/INV/TRA)
    :param annot: gene annotation code
    :param canonical: canonical splice site code {1: canonical, 0: noncanonical}
    :param modes: read modes of connected breakpoints
    :param genes: overlapped genes of connected breakpoints
    :param insertion_info: insertion information, (True, Insertion) or
        (False, NovelInsertion) or (False, MicroHomology)

    .. note::
        connection-level fields:
        * sv_type
        * modes
        * genes
        * annotation_code
        * splicing_code

    :Example:

    >>> node1 = Node(
                prev_bp=None,
                next_bp='chr10:93636994',
                strand='+',
                chrom='chr10',
                ref_start=93636994
                ref_end=93637094,
                exons=[(93636994, 93637094)],
            )
    >>> node1
    Node(chr10:93636994-93637094:-, 93636994-93637094, TRA, None, chr10:93636994)
    """

    __slots__ = (
        "next_breakpoint",
        "prev_breakpoint",
        "strand",
        "chrom",
        "ref_start",
        "ref_end",
        "exons",
        "sv_type",
        "modes",
        "genes",
        "annotation_code",
        "splicing_code",
        "sr",
        "insertion_info",
        "predecessors",
        "successors",
        "unique_key",
        "merged_child_nodes",
        "merged_parent_nodes",
        "next_node_in_series",
        "previous_node_in_series",
        "is_in_graph",
        "is_merged",
    )

    def __init__(
        self,
        prev_bp: Optional[str] = None,
        next_bp: Optional[str] = None,
        strand: Optional[str] = None,
        chrom: Optional[str] = None,
        ref_start: Optional[int] = None,
        ref_end: Optional[int] = None,
        exons: List[Any] = None,
        sv_type: Optional[str] = None,
        annot: Optional[int] = None,
        canonical: Optional[int] = None,
        modes: Optional[Tuple[int]] = None,
        genes: Optional[Tuple[str]] = None,
        sr: Optional[int] = 1,
        insertion_info: Optional[Tuple[bool, Union[Insertion, NovelInsertion]]] = None,
        # type: ignore
    ) -> None:
        """Initialize a Node object."""
        self.chrom = chrom
        self.prev_breakpoint = prev_bp
        self.next_breakpoint = next_bp
        self.strand = strand
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.exons = exons
        self.sv_type = sv_type
        self.modes = modes
        self.genes = genes
        self.annotation_code = annot
        self.splicing_code = canonical
        self.sr = sr
        self.insertion_info = insertion_info
        self.successors: Any = []
        self.predecessors: Any = []

        self.unique_key = None
        self.merged_child_nodes: List = []
        self.merged_parent_nodes: List = []

        self.next_node_in_series = None
        self.previous_node_in_series = None
        self.is_merged, self.is_in_graph = False, False

    def __eq__(self, other) -> bool:
        """Compare two nodes."""
        if isinstance(other, Node) and (
            self.chrom == other.chrom
            and self.ref_start == other.ref_start
            and self.ref_end == other.ref_end
            and self.sv_type == other.sv_type
            and self.prev_breakpoint == other.prev_breakpoint
            and self.next_breakpoint == other.next_breakpoint
            and self.strand == other.strand
        ):
            return True
        return False

    def __hash__(self) -> int:
        """Hash a node."""
        return (
            hash(self.chrom)
            ^ hash(self.ref_start)
            ^ hash(self.ref_end)
            ^ hash(self.sv_type)
            ^ hash(self.prev_breakpoint)
            ^ hash(self.next_breakpoint)
            ^ hash(self.strand)
        )

    def __repr__(self) -> str:
        """Get a string representation of a node."""
        exons_repr = "|".join([f"{i}-{j}" for i, j in self.exons])  # type: ignore
        return (
            f"Node({self.chrom}:{self.ref_start}-{self.ref_end}:{self.strand}, "
            f"{exons_repr}, {self.sv_type}, {self.prev_breakpoint}, "
            f"{self.next_breakpoint}, SR={self.sr}) "
        )

    __str__ = __repr__

    @classmethod
    def create_nodes(cls, number):
        """Create a list of nodes.

        :param number: number of nodes to create
        :return: list of nodes
        """
        return [cls() for _ in range(number)]

    @property
    def introns(self):
        """Get introns of a node."""
        if len(self.exons) <= 1:
            return []
        else:
            _positions = []
            for i, j in self.exons:
                _positions.extend([i, j])
            _positions.pop(0)
            _positions.pop(-1)
            _introns = list(zip(_positions[::2], _positions[1::2]))
            return _introns

    @property
    def similar_key(self):
        """Get similar key of a node."""
        introns = self.introns

        key = "-".join([f"{i - j}" for i, j in introns]) if introns else "None"
        key = f"{self.chrom}-{key}"

        return key

    def get_unique_key(self):
        """Get unique key of a node."""
        introns = self.introns

        key = "-".join([f"{i - j}" for i, j in introns]) if introns else "None"

        key = f"{self.chrom}-{key}-{self.sv_type}-{self.prev_breakpoint}-{self.next_breakpoint}"

        if self.insertion_info is not None:
            _, insertion_type = self.insertion_info
            if insertion_type.__class__.__name__ in ["NovelInsertion", "MicroHomology"]:
                key = f"{insertion_type.query_sequence}-{key}"

        self.unique_key = key
        return key

    def is_start_node(self):
        """Check if a node is a start node."""
        return False if self.has_predecessor() else True

    def is_end_node(self):
        """Check if a node is an end node."""
        return False if self.has_successor() else True

    def has_predecessor(self):
        """Check if a node has a predecessor."""
        return True if self.predecessors else False

    def has_successor(self):
        """Check if a node has a successor."""
        return True if self.successors else False

    def add_successor_from_list(self, successors):
        """Add successors from a list."""
        for successor in successors:
            self.add_successor(successor)

    def add_predecessor_from_list(self, predecessors):
        """Add predecessors from a list."""
        for predecessor in predecessors:
            self.add_predecessor(predecessor)

    def _add_successor(self, successor):
        """Helper function to add a successor."""
        self.successors.append(successor)
        successor.add_predecessor(self)

    def _add_predecessor(self, predecessor):
        """Helper function to add a predecessor."""
        self.predecessors.append(predecessor)
        predecessor.add_successor(self)

    def add_successor(self, successor):
        """Add a successor."""
        if successor is not None and successor not in self.successors:
            if successor.is_in_graph:
                self._add_successor(successor)
            else:
                self.add_successor_from_list(successor.merged_parent_nodes)

    def add_predecessor(self, predecessor):
        """Node must be in the graph if the function is called.

        :param predecessor: predecessor node
        """
        if predecessor is not None and predecessor not in self.predecessors:
            if predecessor.is_in_graph:
                self._add_predecessor(predecessor)
            else:
                self.add_predecessor_from_list(predecessor.merged_parent_nodes)

    def update_sr(self, key=1):
        """Update the sr of a node."""
        self.sr += key

    def update_next_and_previous_node_in_series(self, index, series):
        """Update the next and previous node in series."""
        if index == 0:
            self.next_node_in_series = series[index + 1]
        elif index == len(series) - 1:
            self.previous_node_in_series = series[index - 1]
        else:
            self.next_node_in_series = series[index + 1]
            self.previous_node_in_series = series[index - 1]


NodeType = Union[Node, Insertion]


class Series:
    """Construct a sequence of Nodes for storing information of connected breakpoints.

    :param nodes: sequence of Nodes
    :type nodes: list
    :param assemblied: The series is from assembly of reads (True) or a single read (False)
    :type assemblied: bool or None

    .. note::
        [('TDUP', 0, 1, ('chr17:7708250', 'chr17:7701656', 1, 2), ('+', '+'),
         ['INTERGENIC', 'INTERGENIC']),

        ('TRA', 0, 1, ('chr17:7702552', 'chr1:15872815', 1, 2), ('+', '+'),
        ['INTERGENIC', 'INTERGENIC']),

        ('TDUP', 0, 1, ('chr1:15876678', 'chr1:15777169', 1, 2), ('+', '+'),
        ['INTERGENIC', 'INTERGENIC'])]

        Node(TDUP, None, chr17:7708250, +);Node(TRA, chr17:7701656, chr17:7702552, +)
        Node(TDUP, chr1:15872815, chr1:15876678, +);Node(None, chr1:15777169, None, +)


                         bp1               bp2  bp3                bp4
                ---------|------    -------|----|------    --------|---------
                       Node1                 Node2                Node3
    prev_breakpoint:   None                 bp2                   bp4
    next_breakpoint:    bp1                 bp3                   None
    sv_type:        TDUP/INV/TRA        TDUP/INV/TRA              None

    :Example:

    >>> series = Series(blat=None, logger=logger)
    >>> series.add_node(Node(prev_bp=None,next_bp='chr17:7708250',strand='+',
    ... chrom='chr17',ref_start=7706250,ref_end=7708250,exons=[[7706250,7708250]],sv_type='TDUP'))
    >>> series.add_node(Node(prev_bp='chr17:7701656',next_bp='chr17:7702552',strand='+',
    ... chrom='chr17',ref_start=7701656,ref_end=7702552,exons=[[7701656, 7702552]], sv_type='TRA'))
    >>> series.add_node(Node(prev_bp='chr1:15872815',next_bp='chr1:15876678',strand='+',
    ... chrom='chr1',ref_start=15872815,ref_end=15876678,exons=[[15872815,15876678]],
    ... sv_type='TDUP'))
    >>> series.add_node(Node(prev_bp='chr1:15777169',next_bp=None,strand='+',
    ... chrom='chr1',ref_start=15777169,ref_end=15777589,exons=[[15777169,15777589]],
    ... sv_type=None))
    >>> series
    Series(
        Node(chr17:7706250-7708250:+, 7706250-7708250, TDUP, None, chr17:7708250)
        Node(chr17:7701656-7702552:+, 7701656-7702552, TRA, chr17:7701656, chr17:7702552)
        Node(chr1:15872815-15876678:+, 15872815-15876678, TDUP, chr1:15872815, chr1:15876678)
        Node(chr1:15777169-15777589:+, 15777169-15777589, None, chr1:15777169, None) )

    >>> series_with_novel_insertion = Series(blat=None, logger=logger)
    >>> series_with_novel_insertion.nodes = [ Node(prev_bp=None,next_bp='chr17:7702552',
    ... strand='+',chrom='chr17',ref_start=7701656,ref_end=7702552,exons=[[7701656, 7702552]],
    ... sv_type='TRA', insertion_info=(False, NovelInsertion(hit_num=1,
    ... query_sequence='ATCGATCG'))), Node(prev_bp='chr1:15872815',next_bp=None,strand='+',
    ... chrom='chr1',ref_start=15872815,ref_end=15876678,exons=[[15872815,15876678]],
    ... sv_type=None)]
    >>> series_with_novel_insertion
    Series(
        Node(chr17:7701656-7702552:+, 7701656-7702552, TRA, None, chr17:7702552)
        Node(chr1:15872815-15876678:+, 15872815-15876678, None, chr1:15872815, None) )
    """

    def __init__(self, blat: Any, logger: Logger) -> None:
        """Initialize a Series object."""
        self.nodes: List[Union[Node, Insertion]] = []
        self.is_in_graph = False
        self.blat = blat
        self.logger = logger

    def add_node(self, node: Union[Node, Insertion]) -> None:
        """Add a node to the series."""
        self.nodes.append(node)

    @classmethod
    def create_series_from_node_list(
        cls, node_list: List[NodeType], logger: Logger
    ) -> "Series":
        """Create a series from a list of nodes."""
        series_instance = cls(None, logger)
        for node in node_list:
            series_instance.add_node(node)
        return series_instance

    def init(
        self,
        event_list,
        read_chains,
        splice_bin,
        genome_fasta,
        cvg,
        gene_iv,
        motif_required,
    ) -> None:
        """Add event list as Node to self.nodes."""
        event_list = [
            Event(event)
            for event in self.order_events_by_trancription_direction(event_list)
            if event[0] != "NA"
        ]
        event_list_len = len(event_list)
        previous_breakpoint = None
        for index, event in enumerate(event_list):

            read1_node = Node(
                prev_bp=previous_breakpoint,
                next_bp=event.bp1,
                strand=event.strand1,
                chrom=event.chrom1,
                ref_start=event.read1_ref_start,
                ref_end=event.read1_ref_end,
                exons=event.read1_exons,
            )
            previous_breakpoint = event.bp2

            # is insertions
            if event.has_insertion():

                insertion_seq = event.insertion_seq1  # pick from the first read
                insertion_seq = (
                    reverse_complement(insertion_seq)
                    if event.strand1 == "-"
                    else insertion_seq
                )
                flag, insertion = self.blat.query_insertion(insertion_seq)  # type: ignore
                if flag:  # only one hit
                    # add first node and insertion node
                    source_s = event.source_s1

                    # get type of insertion between first node and insertion node
                    read1 = event.read1(read_chains)
                    insertion.update_cigarstring_sms(
                        read1.sms, source_s=source_s, source_strand=event.strand1
                    )
                    self.logger.trace(f"{insertion.strand=}, {insertion.cigarstring}")
                    if event.strand1 == insertion.strand:
                        insertion_mode = 2 if event.mode1 == 1 else 1
                    else:
                        insertion_mode = event.mode1
                    self.logger.trace("nls reference for read1 and insertion")
                    read1_insertion_event = Event(
                        infer_nls_from_connected_reads(
                            read_lt=read1,
                            read_rt=insertion,
                            lt_mode=event.mode1,
                            rt_mode=insertion_mode,
                            splice_bin=splice_bin,
                            genome_fasta=genome_fasta,
                            cvg=cvg,
                            gene_iv=gene_iv,
                            motif_required=motif_required,
                            logger=self.logger,
                        )
                    )

                    # get type of insertion between insertion node and second node
                    read2 = event.read2(read_chains)

                    if insertion.strand == read2.strand:
                        insertion_mode = 2 if event.mode2 == 1 else 1
                    else:
                        insertion_mode = event.mode2

                    self.logger.trace("nls reference for read1 and insertion")
                    insertion_read2_event = Event(
                        infer_nls_from_connected_reads(
                            read_lt=insertion,
                            read_rt=read2,
                            lt_mode=insertion_mode,
                            rt_mode=event.mode2,
                            splice_bin=splice_bin,
                            genome_fasta=genome_fasta,
                            cvg=cvg,
                            gene_iv=gene_iv,
                            motif_required=motif_required,
                            logger=self.logger,
                        )
                    )

                    if (
                        read1_insertion_event.is_type_na()
                        or insertion_read2_event.is_type_na()
                    ):
                        # only add read1, False means that the insertion type (hit 1 insertion)
                        # are not added in series
                        read1_node = event.update_node_info(
                            False, read1_node, insertion
                        )
                        self.add_node(read1_node)
                    else:
                        # add read1 and insertion
                        # True means that the insertion type(hit 1 insertion) are added in series
                        read1_node = read1_insertion_event.update_node_info(
                            True, read1_node, insertion
                        )
                        self.add_node(read1_node)

                        insertion = insertion_read2_event.update_insertion_info(
                            insertion
                        )

                        self.logger.trace(f"Add {insertion=} to Series")

                        self.add_node(insertion)

                else:  # no hits or multiple hits

                    self.logger.trace(f"Add Novel Insertion {insertion=} to read1")
                    # only add read1 with insertion info
                    # False means that the insertion type (hit more insertion) are
                    # not added in series
                    read1_node = event.update_node_info(False, read1_node, insertion)
                    self.add_node(read1_node)
            # no insertion
            elif event.has_microhomology():
                # add read 1 with on insertion
                microhomology = MicroHomology(event.insertion_seq1)

                self.logger.trace(f"Add MicroHomology {microhomology=} to read1")
                if event.strand1 == "-":
                    microhomology.reverse_completement_query()

                read1_node = event.update_node_info(False, read1_node, microhomology)
                self.add_node(read1_node)

            else:
                read1_node = event.update_node_info(False, read1_node, None, False)
                self.add_node(read1_node)

            # add final node
            if index == event_list_len - 1:
                final_node = Node(
                    prev_bp=previous_breakpoint,
                    strand=event.strand2,
                    chrom=event.chrom2,
                    ref_start=event.read2_ref_start,
                    ref_end=event.read2_ref_end,
                    exons=event.read2_exons,
                )

                self.add_node(final_node)

    @staticmethod
    def reorder_event(event):
        """Order breakpoints pairs following the transcription direction using.

        information of reads 'mode' and 'strand'
        +1;-1 => up;down
        +2;-2 => down;up
        """
        (
            sv_type,
            annot,
            canonical,
            _positions,
            read1_info,
            read2_info,
            insertion_info,
            strands,
            genes,
        ) = event
        bp1 = _positions[0]
        bp2 = _positions[1]
        mode1 = _positions[2]
        mode2 = _positions[3]
        strand1 = strands[0]
        strand2 = strands[1]

        is_bp1_upstream = None
        if strand1 == "+" and strand2 == "-":
            if mode1 == 1 and mode2 == 1:
                is_bp1_upstream = True
            elif mode1 == 2 and mode2 == 2:
                is_bp1_upstream = False
        elif strand1 == "-" and strand2 == "+":
            if mode1 == 1 and mode2 == 1:
                is_bp1_upstream = False
            elif mode1 == 2 and mode2 == 2:
                is_bp1_upstream = True
        elif strand1 == "+" and strand2 == "+":
            if mode1 == 1 and mode2 == 2:
                is_bp1_upstream = True
            elif mode1 == 2 and mode2 == 1:
                is_bp1_upstream = False
        elif strand1 == "-" and strand2 == "-":
            if mode1 == 1 and mode2 == 2:
                is_bp1_upstream = False
            elif mode1 == 2 and mode2 == 1:
                is_bp1_upstream = True

        if not is_bp1_upstream:
            if annot == 1:
                annot = 2
            elif annot == 2:
                annot = 1
            _positions = (bp2, bp1, mode2, mode1)
            strands = (strand2, strand1)
            genes = list(reversed(genes))
            return (
                sv_type,
                annot,
                canonical,
                _positions,
                read2_info,
                read1_info,
                insertion_info[::-1],
                strands,
                genes,
            )
        else:
            return event

    @staticmethod
    def order_events_by_trancription_direction(event_list):
        """Construct breakpoints order following transcription direction.

         for multiple-hop events or one-hop events
                bp1                bp2   bp3               bp4
        ---------|------    -------|----|------    --------|---------
              Node1                 Node2                Node3
        ..note ::
               requirements
               * bp2 and bp3 at the same chromosome
               * if strand(+): bp3 > bp2
                 if strand(-): bp3 < bp2
        """
        kept_right_pos = None
        kept_right_chrm = None
        kept_right_strand = None

        keep_event_list_order = True
        output_event_list = []
        for evt in event_list:
            ordered_evt = Series.reorder_event(evt)
            output_event_list.append(ordered_evt)
            (
                sv_type,
                annot,
                canonical,
                _positions,
                read1_info,
                read2_info,
                insertion_info,
                strands,
                genes,
            ) = ordered_evt
            chrm1, _pos1 = _positions[0].split(":")
            chrm2, _pos2 = _positions[1].split(":")
            pos1 = int(_pos1)
            pos2 = int(_pos2)
            strand1 = strands[0]
            strand2 = strands[1]
            if not kept_right_pos:
                kept_right_pos = pos2
                kept_right_chrm = chrm2
                kept_right_strand = strand2
            else:
                if chrm1 == kept_right_chrm and strand1 == kept_right_strand:
                    if (strand1 == "+" and pos1 > kept_right_pos) or (
                        strand1 == "-" and pos1 < kept_right_pos
                    ):
                        kept_right_pos = pos2
                        kept_right_chrm = chrm2
                        kept_right_strand = strand2
                    else:
                        keep_event_list_order = False
                else:
                    keep_event_list_order = False

        if not keep_event_list_order:
            output_event_list = list(reversed(output_event_list))

        return output_event_list

    def __getitem__(self, index):
        """Return the event at the given index."""
        return self.nodes[index]

    def __hash__(self) -> int:
        """Return the hash of the event."""
        return hash(";".join(map(str, self.nodes)))

    def __eq__(self, other) -> bool:
        """Return True if the events are equal."""
        return ";".join(map(str, self.nodes)) == ";".join(map(str, other.nodes))

    def __len__(self) -> int:
        """Return the number of events."""
        return len(self.nodes)

    def __lt__(self, other) -> bool:
        """Return True if the event is less than the other event."""
        return len(self.nodes) < len(other.nodes)

    def __repr__(self) -> str:
        """Return the string representation of the event."""
        _repr = "\nSeries("
        space = " " * 4
        for n in self.nodes:
            _repr += f"\n{space}{n!r}"

        _repr += ")"
        return _repr

    def __iter__(self):
        """Return an iterator over the events."""
        yield from self.nodes

    __str__ = __repr__

    def disable_blat_logger(self):
        """Disable blat logger."""
        self.blat, self.logger = None, None

    @property
    def start_node(self):
        """Return the start node of the event."""
        return self.nodes[0]

    @start_node.setter
    def start_node(self, node):
        """Set the start node of the event."""
        self.nodes[0] = node

    @property
    def end_node(self):
        """Return the end node of the event."""
        return self.nodes[-1]

    @end_node.setter
    def end_node(self, node):
        """Set the end node of the event."""
        self.nodes[-1] = node

    @property
    def unique_key(self):
        """Return the unique key of the event."""
        return "".join([node.get_unique_key() for node in self.nodes])


class Event:
    """Event class is used to parse the return value from the function nls_inference.

    :Example:

    >>> args, kwargs = [], {}
    >>> event = Event(infer_nls_from_connected_reads(*args, **kwargs))
    >>> event.sv_type
    TRA
    >>> event
    Event(TRA, )

    .. todo::
        add more examples
    """

    def __init__(self, event):
        """Initialize the event."""
        (
            sv_type,
            annot,
            canonical,
            _positions,
            read1_info,
            read2_info,
            insertion_info,
            strands,
            genes,
        ) = event

        self.sv_type = sv_type
        if self.sv_type != "NA":
            self.annotation_code = annot
            self.splicing_code = canonical
            self.genes = genes
            self.insertion_info = insertion_info
            self.bp1, self.bp2 = _positions[:2]
            self.mode1, self.mode2 = _positions[2:]
            self.strand1, self.strand2 = strands
            self.read1_ref_start, self.read1_ref_end, self.read1_exons = read1_info
            self.read2_ref_start, self.read2_ref_end, self.read2_exons = read2_info

    def __repr__(self):
        """Return the string representation of the event."""
        if self.sv_type != "NA":
            return (
                f"Event({self.sv_type}, {self.annotation_code}, {self.splicing_code})"
            )

    @property
    def modes(self) -> List[int]:
        """Return the modes of the event.

        :return: the mode of read1 and read2 in the event
        """
        return [self.mode1, self.mode2]

    @property
    def chrom1(self) -> str:
        """Return chromosome of read1."""
        return self.bp1.split(":")[0]

    @property
    def chrom2(self) -> str:
        """Return chromosome of read2."""
        return self.bp2.split(":")[0]

    @property
    def insertion_seq1(self) -> str:
        """Return the insertion sequence of read1."""
        return self.insertion_info[0][1:]

    @property
    def insertion_seq2(self) -> str:
        """Return the insertion sequence of read2."""
        return self.insertion_info[1][1:]

    @property
    def source_s1(self) -> str:
        """Source of insertion of read1."""
        return "left" if self.mode1 == 2 else "right"

    @property
    def source_s2(self) -> str:
        """Source of insertion of read2."""
        return "left" if self.mode2 == 2 else "right"

    def is_type_na(self) -> bool:
        """Return True if the event is NA."""
        return True if self.sv_type == "NA" else False

    def has_insertion(self) -> bool:
        """Return True if the event has insertion."""
        return True if self.insertion_info[0].startswith("+") else False

    def has_microhomology(self) -> bool:
        """Return True if the event has microhomology."""
        return True if self.insertion_info[0].startswith("-") else False

    def is_same_strand(self) -> bool:
        """Return True if the event is same strand."""
        return self.strand1 == self.strand2

    def read1(self, read_chains: List[Read]) -> Read:
        """Return the read1 of the event."""
        for read in read_chains:
            if read.ref_start == self.read1_ref_start:
                return read
        else:
            raise ReadNotFoundError

    def read2(self, read_chains: List[Read]) -> Read:
        """Return the read2 of the event."""
        for read in read_chains:
            if read.ref_start == self.read2_ref_start:
                return read
        else:
            raise ReadNotFoundError

    def update_specific_info_within_event(
        self, node: Union[Node, Insertion], info_key_list: List[str]
    ) -> Union[Node, Insertion]:
        """Update node info from the event by the info_key_list.

        :param node:  Node
        :param info_key_list: [key1, key2, ...]
        :return: Node with updated info
        """
        for key in info_key_list:
            setattr(node, key, getattr(self, key))

        return node

    def update_node_info(
        self,
        flag: bool,
        new_node: Node,
        insertion: Union[Insertion, None, MicroHomology],
        is_update_insertion_info: bool = True,
    ) -> Node:
        """Update the common info the node in the front, and the common info includes.

        sv_type, annot, canonical, genes, insertion_info, and the breakpoints, mode

        :param flag: the flag indicates whether there is a insertion
        :param new_node: the new node to be updated
        :param insertion: the insertion to be updated
        :param is_update_insertion_info: whether to update the insertion info
        :return: the updated node
        """
        new_node = self.update_specific_info_within_event(
            new_node, ["sv_type", "annotation_code", "splicing_code", "modes", "genes"]
        )  # type: ignore
        if is_update_insertion_info:
            new_node.insertion_info = (flag, insertion)  # type: ignore
        return new_node

    def update_insertion_info(self, insertion: Insertion) -> Union[Node, Insertion]:
        """Update the information of insertion.

        :param insertion: the insertion to be updated
        :return: the updated insertion
        """
        insertion.prev_breakpoint = f"{insertion.chrom}:{insertion.ref_start}"
        insertion.next_breakpoint = f"{insertion.chrom}:{insertion.ref_end}"

        insertion.exons, _ = insertion.get_exons_and_introns()
        return self.update_specific_info_within_event(
            insertion, ["sv_type", "annotation_code", "splicing_code", "modes", "genes"]
        )
