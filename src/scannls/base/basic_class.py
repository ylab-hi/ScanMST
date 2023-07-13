"""Type of the scannls.

@Author:      YangyangLi
@Time:        12/30/21 2:20 PM
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from scannls import cppext
from scannls.cli.helper import cigar_validity
from scannls.exception import ReadNotFoundError

from .basic import AnnotationCode, Mode, Strand
from .basic_read import Read

if TYPE_CHECKING:
    from collections.abc import Iterable

    from scannls.type import EventType


class NovelInsertion:
    """NovelInsertion is used to represent reads insertion whose hit is 0 or >1.

    :Example:

    >>> novel_insertion = NovelInsertion(hit_num=0, query_sequence='ATCA')
    >>> novel_insertion
    NovelInsertion(ATCA:0)
    >>> novel_insertion.query_sequence
    ATCA
    >>> novel_insertion.hit_num
    0
    >>> novel_insertion.ao
    1

    .. note::
        `NovelInsertion` is a subclass of :class:`Read`, and siblings of :class:`Insertion`

    .. seealso::
        :class: `Insertion`
    """

    def __init__(self, hit_num: int, query_sequence: str) -> None:
        """Initialize NovelInsertion."""
        self.query_sequence = query_sequence
        self.hit_num = hit_num
        self.insertion_info = None
        self.ao = 1

    def __repr__(self) -> str:
        """Represent NovelInsertion object."""
        return f"NovelInsertion({self.query_sequence}:{self.hit_num})"

    def reverse_completement_query(self) -> None:
        """Reverse complement query sequence."""
        self.query_sequence = reverse_complement(self.query_sequence)

    def increment_ao(self, num=1) -> None:
        """Increment ao."""
        self.ao += num


class MicroHomology:
    """MicroHomology is used to represent microhomology.

    :Example:

    >>> microhomology = MicroHomology(query_sequence="ATCA")
    >>> microhomology
    MicroHomology(ATCA)
    >>> microhomology.query_sequence
    ATCA
    >>> microhomology.ao
    1

    .. seealso::
        :class:`Insertion` and :class:`NovelInsertion`
    """

    def __init__(self, query_sequence: str) -> None:
        """Initialize MicroHomology."""
        self.query_sequence = query_sequence
        self.ao = 1

    def __repr__(self) -> str:
        """Represent MicroHomology object."""
        return f"MicroHomology({self.query_sequence})"

    def reverse_completement_query(self):
        """Reverse complement query sequence."""
        self.query_sequence = reverse_complement(self.query_sequence)

    def increment_ao(self, num=1) -> None:
        """Increment ao."""
        self.ao += num


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

    >>> import array
    >>> insertion = Insertion(hit_num=1, chrom='1', ref_start=1, strand='+',
    ...                 cigarstring='1S1M1S',mapq=60, nm=0, query_sequence='ATCA',
    ...                 query_qualities=array.array('B', [10,20,10,9]))
    >>> insertion
    Insertion(1:1-4:+, 1-2|2-3, TPA, 1, 4)

    .. note::
        `Insertion` is a subclass of :class:`Read`, and siblings of :class:`NovelInsertion`
        `Insertion` includes the attributes of :class:`Node` in order to enable us to
        manipulate the attributes of `Insertion` same as :class:`Node` in :class:`Series`

    .. seealso:: :class:`NovelInsertion`, :class:`Node` and :class:`Read`
    """

    __slots__ = ("hit_num", *Read.__slots__)

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
        query_qualities: list[int],
    ) -> None:
        """Initialize Insertion."""
        parse_cigar_result = cppext.parseCigar(cigarstring)
        super().__init__(
            "",  # query_name
            chrom,
            ref_start,
            strand,
            cigarstring,
            mapq,
            nm,
            query_sequence,
            parse_cigar_result.lt_soft_len,
            parse_cigar_result.rt_soft_len,
            parse_cigar_result.read_match,
            parse_cigar_result.ref_match,
            parse_cigar_result.indel_len,
            parse_cigar_result.cigartuples_without_soft,
            parse_cigar_result.query_len,
            query_qualities,
        )

        self.hit_num = hit_num

    def __repr__(self) -> str:
        """Represent Insertion object."""
        return (
            f"{self.__class__.__name__}({self.hit_num=},"
            f"{self.chrom}:{self.ref_start}-{self.ref_end}:{self.strand})"
        )

    def update_cigarstring_sms(
        self,
        sms: Iterable[int],
        source_s: str,
        source_strand: str | Strand,
    ) -> None:
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

    def reverse_completement_query(self) -> None:
        """Reverse complement query sequence of Insertion object."""
        self.query_sequence = reverse_complement(self.query_sequence)

    def reverse_strand(self) -> None:
        """Reverse strand of Insertion object."""
        self.strand.reverse()


@dataclass(unsafe_hash=True)
class BreakPoint:
    """BreakPoint is used to represent breakpoints."""

    chrom: str
    pos: int
    depth: int = 1

    def __str__(self) -> str:
        """Return string representation of BreakPoint object."""
        return f"{self.chrom}:{self.pos}"

    def to_tuple(self) -> tuple[str, int]:
        """Return tuple representation of BreakPoint object."""
        return self.chrom, self.pos

    @classmethod
    def from_str(cls, breakpoint_str: str, depth: int | None = None) -> BreakPoint:
        """Create BreakPoint object from string."""
        if breakpoint_str == "":
            msg = "breakpoint_str can not be empty string"
            raise ValueError(msg)
        chrom, pos = breakpoint_str.split(":")
        depth = int(depth) if depth is not None else cls.depth
        return cls(chrom, int(pos), depth)

    @classmethod
    def from_node(cls, node) -> tuple[BreakPoint, BreakPoint]:
        if node.strand.is_forward():
            prev_breakpoint = BreakPoint.from_str(f"{node.chrom}:{node.ref_start}")
            next_breakpoint = BreakPoint.from_str(f"{node.chrom}:{node.ref_end}")
        else:
            prev_breakpoint = BreakPoint.from_str(f"{node.chrom}:{node.ref_end}")
            next_breakpoint = BreakPoint.from_str(f"{node.chrom}:{node.ref_start}")

        return prev_breakpoint, next_breakpoint

    def equal(
        self,
        other: BreakPoint,
        threshold: int,
    ) -> bool:
        """Check if two breakpoints are different.
        :param breakpoint1:  breakpoint1
        :param breakpoint2:  breakpoint2
        :param threshold:  threshold for checking if two breakpoints are different
        :return:  True if two breakpoints are different, otherwise False.
        """
        if self.chrom != other.chrom:
            return False

        return abs(self.pos - other.pos) <= threshold


class Event:
    """Event class is used to parse the return value from the function nls_inference.

    :Example:

    >>> args, kwargs = [], {}
    >>> event = Event(infer_nls_from_connected_reads(*args, **kwargs))
    >>> event.sv_type
    TRA
    >>> event
    Event(TRA, )
    """

    def __init__(self, event: EventType) -> None:
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
            self.positions = _positions
            self.bp1, self.bp2 = _positions[:2]
            self.mode1, self.mode2 = _positions[2:]
            self.strand1 = Strand.from_str(strands[0])
            self.strand2 = Strand.from_str(strands[1])
            self.read1_ref_start, self.read1_ref_end, self.read1_exons = read1_info
            self.read2_ref_start, self.read2_ref_end, self.read2_exons = read2_info

    def __repr__(self) -> str:
        """Return the string representation of the event."""
        if self.sv_type != "NA":
            return (
                f"Event({self.sv_type}, {self.annotation_code}, {self.splicing_code} ({self.bp1} "
                f"{self.bp2} {self.mode1} {self.mode2}) "
                f"{self.strand1} {self.read1_ref_start} {self.read1_ref_end} {self.read1_exons} "
                f"{self.strand2} {self.read2_ref_start} {self.read2_ref_end} {self.read2_exons} "
                f"{self.insertion_info})"
            )
        return f"Event({self.sv_type})"

    def reverse(self) -> None:
        """Reverse breakpoint1 and breakpoint2."""
        if self.annotation_code == AnnotationCode.Type1:
            self.annotation_code = AnnotationCode.Type2
        elif self.annotation_code == AnnotationCode.Type2:
            self.annotation_code = AnnotationCode.Type1

        self.bp1, self.bp2 = self.bp2, self.bp1
        self.mode1, self.mode2 = self.mode2, self.mode1
        self.strand1, self.strand2 = self.strand2, self.strand1
        self.genes = self.genes[::-1]
        self.read1_ref_start, self.read2_ref_start = (
            self.read2_ref_start,
            self.read1_ref_start,
        )
        self.read1_ref_end, self.read2_ref_end = self.read2_ref_end, self.read1_ref_end
        self.read1_exons, self.read2_exons = self.read2_exons, self.read1_exons
        self.insertion_info = self.insertion_info[::-1]

    @property
    def modes(self) -> list[int]:
        """Return the modes of the event.

        :return: the mode of read1 and read2 in the event
        """
        return [self.mode1, self.mode2]

    @property
    def chrom1(self) -> str:
        """Return chromosome of read1."""
        return str(self.bp1.split(":")[0])

    @property
    def chrom2(self) -> str:
        """Return chromosome of read2."""
        return str(self.bp2.split(":")[0])

    @property
    def junction_pos1(self) -> int:
        """Return junction position of read1."""
        return int(self.bp1.split(":")[1])

    @property
    def junction_pos2(self) -> int:
        """Return junction position of read2."""
        return int(self.bp2.split(":")[1])

    @property
    def insertion_seq1(self) -> str:
        """Return the insertion sequence of read1."""
        return str(self.insertion_info[0][1:])

    @property
    def insertion_seq2(self) -> str:
        """Return the insertion sequence of read2."""
        return str(self.insertion_info[1][1:])

    @property
    def source_s1(self) -> str:
        """Source of insertion of read1."""
        return "left" if self.mode1 == Mode.Type2 else "right"

    @property
    def source_s2(self) -> str:
        """Source of insertion of read2."""
        return "left" if self.mode2 == Mode.Type2 else "right"

    def is_type_na(self) -> bool:
        """Return True if the event is NA."""
        return self.sv_type == "NA"

    def has_insertion(self) -> bool:
        """Return True if the event has insertion."""
        return bool(self.insertion_info[0].startswith("+"))

    def has_microhomology(self) -> bool:
        """Return True if the event has microhomology."""
        return bool(self.insertion_info[0].startswith("-"))

    @property
    def insertion_microhomology_len(self) -> int:
        """Return the length of insertion or microhomology."""
        if self.has_insertion() or self.has_microhomology():
            return len(self.insertion_info[0]) - 1
        return 0

    def is_same_strand(self) -> bool:
        """Return True if the event is same strand."""
        return self.strand1 == self.strand2

    def read1(self, read_chains: list[Read]) -> Read:
        """Return the read1 of the event."""
        for read in read_chains:
            if (
                read.ref_start == self.read1_ref_start
                and read.ref_end == self.read1_ref_end
                and read.strand == self.strand1
            ):
                return read
        raise ReadNotFoundError

    def read2(self, read_chains: list[Read]) -> Read:
        """Return the read2 of the event."""
        for read in read_chains:
            if (
                read.ref_start == self.read2_ref_start
                and read.ref_end == self.read2_ref_end
                and read.strand == self.strand2
            ):
                return read
        raise ReadNotFoundError

    def update_specific_info_within_event(self, node, info_key_list: list[str]):
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
        new_node,
    ):
        """Update the common info the node in the front, and the common info includes.

        sv_type, annot, canonical, genes, insertion_info, and the breakpoints, mode

        :param flag: the flag indicates whether there is a insertion
        :param new_node: the new node to be updated
        :param insertion: the insertion to be updated
        :param is_update_insertion_info: whether to update the insertion info
        :return: the updated node
        """
        new_node = self.update_specific_info_within_event(
            new_node,
            ["annotation_code", "splicing_code", "modes", "genes"],
        )

    def update_insertion_node_info(self, insertion_node):
        """Update the information of insertion.

        :param insertion_node: the insertion to be updated
        :return: the updated insertion
        """
        self.update_specific_info_within_event(
            insertion_node,
            ["annotation_code", "splicing_code", "modes", "genes"],
        )


def reverse_complement(seq: str) -> str:
    """Obtain reverse complement sequence."""
    rctrans = str.maketrans("ACGT", "TGCA")
    return str.translate(seq, rctrans)[::-1]
