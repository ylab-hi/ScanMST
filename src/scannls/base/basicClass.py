"""Type of the scannls.

@Author:      YangyangLi
@license:     MIT Licence
@Time:        12/30/21 2:20 PM
"""
from collections import Counter
from collections.abc import Iterable
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Any
from typing import Optional

import pyfaidx
from scannls import cppext

from ..cli.helper import cigar_validity
from ..cli.nls_inference import infer_nls_from_connected_reads
from .basicRead import Read
from .exception import ReadNotFoundError
from .type import EventType
from .type import LoggerType


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

    def __repr__(self):
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
        ),

        self.hit_num = hit_num

    def __repr__(self):
        """Represent Insertion object."""
        return (
            f"{self.__class__.__name__}({self.hit_num=},"
            f"{self.chrom}:{self.ref_start}-{self.ref_end}:{self.strand})"
        )

    def update_cigarstring_sms(
        self, sms: Iterable[int], source_s: str, source_strand: str
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
        self.strand = "-" if self.strand == "+" else "+"


class BasicNode:
    """BasicNode is used to represent nodes in the splice graph."""

    __slots__ = (
        # parent class: basic node fields
        "successors",
        "predecessors",
        "merged_child_nodes",
        "merged_parent_nodes",
        "next_node_in_series",
        "previous_node_in_series",
        "is_merged",
        "is_in_graph",
        "is_traced",
        "trace_id",
        "sr",
        "harmonic_mean_sr",
        "original_sr",
    )

    def __init__(self) -> None:
        """Initialize BasicNode object."""
        self.successors: list[Node] = []
        self.predecessors: list[Node] = []

        self.merged_child_nodes: list[Node] = []
        self.merged_parent_nodes: list[Node] = []

        self.next_node_in_series: Node | None = None
        self.previous_node_in_series: Node | None = None
        self.is_merged, self.is_in_graph, self.is_traced = False, False, False
        self.trace_id: int = -1
        self.sr: int = 1
        self.original_sr: int = -1
        self.harmonic_mean_sr: float = 1

    def __eq__(self, other) -> bool:
        """Compare two nodes in strict mode same memory address."""
        return id(self) == id(other)

    def update_sr(self, key=1) -> None:
        """Update the sr of a node."""
        self.sr += key

    def set_harmoic_mean_sr(self, sr: int) -> None:
        """Get harmonic mean of current sr and predecessor.sr."""
        harmoic_mean_sr = 1 / sr
        _node_numbers = 1

        for _node_numbers, node in enumerate(self.predecessors, 2):
            harmoic_mean_sr += 1 / node.sr

        self.harmonic_mean_sr = _node_numbers / harmoic_mean_sr

    def set_trace_id(self, trace_id: int) -> None:
        """Set trace_id."""
        if not self.is_traced:
            self.trace_id = trace_id
            self.is_traced = True

    def reset_trace_id(self) -> None:
        """Reset trace_id."""
        self.trace_id = -1
        self.is_traced = False

    def is_start_node(self) -> bool:
        """Return True if Insertion object is start node."""
        return bool(not self.has_predecessor())

    def is_end_node(self) -> bool:
        """Return True if Insertion object is end node."""
        return bool(not self.has_successor())

    def has_predecessor(self) -> bool:
        """Return True if Insertion object has predecessor."""
        return bool(self.predecessors)

    def has_successor(self) -> bool:
        """Return True if Insertion object has successor."""
        return bool(self.successors)

    def add_successor_from_list(self, successors) -> None:
        """Add successor from list of Insertion object."""
        for successor in successors:
            self.add_successor(successor)

    def add_predecessor_from_list(self, predecessors) -> None:
        """Add predecessor from list of Insertion object."""
        for predecessor in predecessors:
            self.add_predecessor(predecessor)

    def _add_successor(self, successor) -> None:
        """Helper function to add successor to Insertion object."""
        self.successors.append(successor)
        successor.add_predecessor(self)

    def _add_predecessor(self, predecessor) -> None:
        """Helper function to add predecessor to Insertion object."""
        self.predecessors.append(predecessor)
        predecessor.add_successor(self)

    def add_successor(self, successor) -> None:
        """Add successor to Insertion object."""
        if successor is not None and successor not in self.successors:
            if successor.is_in_graph:
                self._add_successor(successor)
            else:
                self.add_successor_from_list(successor.merged_parent_nodes)

    def add_predecessor(self, predecessor) -> None:
        """Node must be in the graph if the function is called.

        :param predecessor: predecessor of Insertion object
        """
        if predecessor is not None and predecessor not in self.predecessors:
            if predecessor.is_in_graph:
                self._add_predecessor(predecessor)
            else:
                self.add_predecessor_from_list(predecessor.merged_parent_nodes)

    def update_next_and_previous_node_in_series(
        self, index: int, series: "Series"
    ) -> None:
        """Update next and previous node in series."""
        if index == 0:
            self.next_node_in_series = series[index + 1]
        elif index == len(series) - 1:
            self.previous_node_in_series = series[index - 1]
        else:
            self.next_node_in_series = series[index + 1]
            self.previous_node_in_series = series[index - 1]

    def clear_next_and_previous_node_in_series(self) -> None:
        """Clear next and previous node in series."""
        self.next_node_in_series = None
        self.previous_node_in_series = None


@dataclass(unsafe_hash=True)
class BreakPoint:
    """BreakPoint is used to represent breakpoints."""

    chrom: str
    pos: int

    def __str__(self):
        """Return string representation of BreakPoint object."""
        return f"{self.chrom}:{self.pos}"

    def to_tuple(self) -> tuple[str, int]:
        """Return tuple representation of BreakPoint object."""
        return self.chrom, self.pos

    @classmethod
    def from_str(cls, breakpoint_str: str | None) -> Optional["BreakPoint"]:
        """Create BreakPoint object from string."""
        if breakpoint_str is None:
            return None
        chrom, pos = breakpoint_str.split(":")
        return cls(chrom, int(pos))


class Node(BasicNode):
    """Build a breakpoint node class for storing information of every breakpoint.

    :param prev_breakpoint: breakpoint for the previous breakpoints connections
    :param next_bp: breakpoint for the next breakpoints connections
    :param strand: direction of chimeric read (-|+)
    :param chrom: chromosome
    :param ref_start: reference start position
    :param ref_end: reference end position
    :param exons: CIGAR inferred exons in the read. e.g., [(100, 200), (300, 500)]
    :param sv_type: one of the NLS types (TDUP/INV/TRA)
    :param annot: gene annotation code
    :param canonical: canonical splice site code {1: canonical, 0: noncanonical}
    :param modes: read modes of connected breakpoints
    :param genes: overlapped genes of connected breakpoints

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
        "prev_breakpoint_depth",
        "next_breakpoint_depth",
        "strand",
        "chrom",
        "ref_start",
        "ref_end",
        "exons",
        "_introns",
        "sv_type",
        "prev_sv_type",
        "modes",
        "genes",
        "query_name",
        "annotation_code",
        "splicing_code",
        "insertion_info",
        "unique_key",
        "is_polya",
        "_exon_repr",
        "cigartuples_without_soft",
        *BasicNode.__slots__,
    )

    def __init__(
        self,
        prev_bp: str | None = None,
        next_bp: str | None = None,
        strand: str | None = None,
        chrom: str | None = None,
        ref_start: int | None = None,
        ref_end: int | None = None,
        exons: list[Any] | None = None,
        sv_type: str | None = None,
        annot: int | None = None,
        canonical: int | None = None,
        modes: list[int] | None = None,
        genes: tuple[str, str] | None = None,
        query_name: str = "",
    ) -> None:
        """Initialize a Node object."""
        super().__init__()  # initialize BasicNode object
        self._introns = None
        self.chrom = chrom
        self.query_name = query_name
        self.prev_breakpoint = BreakPoint.from_str(prev_bp)
        self.next_breakpoint = BreakPoint.from_str(next_bp)
        self.prev_breakpoint_depth: int | None = None
        self.next_breakpoint_depth: int | None = None
        self.strand = strand
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.exons = exons
        self._exon_repr = ""
        self.sv_type = sv_type
        self.prev_sv_type = None
        self.modes = modes
        self.genes = genes
        self.annotation_code = annot
        self.splicing_code = canonical
        self.insertion_info = None
        self.unique_key = None
        self.is_polya = False
        self.cigartuples_without_soft: list[int] | None = None

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
        return (
            f"{self.__class__.__name__}({self.chrom}:{self.ref_start}-{self.ref_end}:{self.strand}, "
            f"{self.exons_repr}, {self.prev_sv_type}, {self.sv_type}, "
            f"{self.prev_breakpoint}|DP:{self.prev_breakpoint_depth}, "
            f"{self.next_breakpoint}|DP:{self.next_breakpoint_depth}, modes={self.modes}, "
            f"SR={self.sr}, query_name={self.query_name.split(',')[:3]}, trace_id={self.trace_id})"
        )

    @classmethod
    def create_nodes(cls, number):
        """Create a list of nodes.

        :param number: number of nodes to create
        :return: list of nodes
        """
        return [cls() for _ in range(number)]

    @property
    def exons_repr(self) -> str:
        """Get a string representation of exons."""
        if self._exon_repr == "":
            self._exon_repr = "|".join([f"{i}-{j}" for i, j in self.exons])  # type: ignore
        return self._exon_repr

    @property
    def introns(self):
        """Get introns of a node."""
        if self._introns is not None:
            return self._introns

        if len(self.exons) <= 1:
            return []

        positions = []
        for i, j in self.exons:
            positions.extend([i, j])
        positions.pop(0)
        positions.pop(-1)
        self._introns = list(zip(positions[::2], positions[1::2]))
        return self._introns

    @property
    def similar_key(self) -> str:
        """Get similar key of a node."""
        introns = self.introns
        chosen_intron = None
        if introns:
            chosen_intron = introns[-1] if self.strand == "+" else introns[0]

        key = f"{chosen_intron[0]}-{chosen_intron[1]}" if chosen_intron else "None"
        return f"{self.chrom}_{key}"

    def get_unique_key(self):
        """Get unique key of a node."""
        introns = self.introns

        key = "-".join([f"{i}-{j}" for i, j in introns]) if introns else "None"

        key = f"{self.chrom}-{key}-{self.sv_type}-{self.prev_breakpoint}-{self.next_breakpoint}"

        if self.insertion_info is not None:
            _, insertion_type = self.insertion_info
            if insertion_type.__class__.__name__ in ("NovelInsertion", "MicroHomology"):
                key = f"{insertion_type.query_sequence}-{key}"

        self.unique_key = key
        return key

    def get_breakpoint_depth_pos(self, mode: int, direc: str) -> tuple[str, Any]:
        """Get update breakpoint depth and position of a node."""
        break_point = self.prev_breakpoint if direc == "prev" else self.next_breakpoint
        if break_point is not None:
            chrom, pos = break_point.to_tuple()
            if mode == 1:
                pos -= 1
            return chrom, pos
        return " ", 1

    def is_reverse(self) -> bool:
        """Check if a node is reverse."""
        return self.strand == "-"


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

    >>> from loguru import logger
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
    Series(
            >>> series_with_novel_insertion
        Node(chr17:7701656-7702552:+, 7701656-7702552, TRA, None, chr17:7702552)
        Node(chr1:15872815-15876678:+, 15872815-15876678, None, chr1:15872815, None) )
    """

    reorder_conditions_dict = {
        "+-11": True,
        "+-22": False,
        "-+11": False,
        "-+22": True,
        "++12": True,
        "++21": False,
        "--12": False,
        "--21": True,
    }

    def __init__(self, blat: Any, logger: LoggerType) -> None:
        """Initialize a Series object."""
        self.nodes: list[Node] = []
        self.is_in_graph = False
        self.blat = blat
        self.logger = logger
        self.id = -1

    def add_node(self, node: Node) -> None:
        """Add a node to the series."""
        self.nodes.append(node)

    def is_all_type_del(self) -> bool:
        """Check if sv_type of all nodes in the series are DEL."""
        return all(
            node.sv_type == "DEL" for node in self.nodes if node.sv_type is not None
        )

    def is_all_node_sr_higher_than_threshold(self, threshold: int) -> bool:
        """Check if all nodes in the series have sr > threshold."""
        return all(node.sr >= threshold for node in self.nodes[:-1])

    def get_sr_sum_for_all_node(self) -> int:
        """Get sum of sr for all nodes in the series."""
        return sum(node.sr for node in self.nodes)

    @classmethod
    def create_series_from_node_list(
        cls,
        node_list: list[Node],
        logger: LoggerType,
        nodes_keys: set[str],
        is_add_key: bool = True,
    ) -> "Series":
        """Create a series from a list of nodes."""
        series_instance = cls(None, logger)
        for node in node_list:
            if is_add_key and (key := node.unique_key) is not None:
                nodes_keys.add(key)
            series_instance.add_node(node)
        series_instance.disable_blat_logger()  # support parallel processing
        return series_instance

    def __getitem__(self, index: int) -> Node:
        """Return the event at the given index."""
        return self.nodes[index]

    def __hash__(self) -> int:
        """Return the hash of the event."""
        return hash(";".join(map(str, self.nodes)))

    def __len__(self) -> int:
        """Return the number of events."""
        return len(self.nodes)

    def __lt__(self, other: Any) -> bool:
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

    def __iter__(self) -> Iterator[Node]:
        """Return an iterator over the events."""
        yield from self.nodes

    def disable_blat_logger(self) -> None:
        """Disable blat logger."""
        self.blat, self.logger = None, None  # type: ignore

    @property
    def unique_key(self) -> str:
        """Return the unique key of the event."""
        return "".join([node.get_unique_key() for node in self.nodes])

    @staticmethod
    def reorder_event(evt: "Event"):
        """Order breakpoint pairs following the transcription direction using.

        information of reads 'mode' and 'strand'
        +1;-1 => up;down
        +2;-2 => down;up
        """
        is_bp1_upstream = Series.reorder_conditions_dict.get(
            f"{evt.strand1}{evt.strand2}{evt.mode1}{evt.mode2}", None
        )

        if not is_bp1_upstream:
            evt.reverse()

        return evt

    @staticmethod
    def order_events_by_trancription_direction(event_list: list["Event"]):
        """Construct breakpoints order following transcription direction.

         for multiple-hop events or one-hop events
                bp1                bp2   bp3               bp4
        ---------|------    -------|----|------    --------|---------
              Node1      |         Node2        |        Node3
                       event1                 event2
        ..note ::
               requirements
               * the read where `breakpoint2` of event1 habors and
               the read where `breakpoint1` of event2 habors should be the identical
        """
        is_reversed = False
        output_event_list: list["Event"] = []
        for index, evt in enumerate(event_list):
            parsed_evt = Series.reorder_event(evt)

            if index == 1:
                last_event = output_event_list[-1]
                if not (
                    last_event.chrom2 == parsed_evt.chrom1
                    and last_event.strand2 == parsed_evt.strand1
                    and last_event.read2_ref_start == parsed_evt.read1_ref_start
                    and last_event.read2_ref_end == parsed_evt.read1_ref_end
                ):
                    is_reversed = True

            output_event_list.append(parsed_evt)
        if is_reversed:
            output_event_list = output_event_list[::-1]

        return output_event_list

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
        event_list = self.order_events_by_trancription_direction(event_list)

        event_list_len = len(event_list)
        previous_breakpoint = None
        prev_sv_type = None
        for index, event in enumerate(event_list):
            read1: Read = event.read1(read_chains)
            read2: Read = event.read2(read_chains)

            read1_node: Node = Node(
                prev_bp=previous_breakpoint,
                next_bp=event.bp1,
                strand=event.strand1,
                chrom=event.chrom1,
                ref_start=event.read1_ref_start,
                ref_end=event.read1_ref_end,
                exons=event.read1_exons,  # type: ignore
                query_name=read1.query_name,
            )

            read1_node.prev_sv_type = prev_sv_type
            read1_node.cigartuples_without_soft = read1.cigartuples_without_soft

            previous_breakpoint = event.bp2
            prev_sv_type = event.sv_type
            self.logger.trace(f"{read1=} {read2=}")
            # is insertions
            if event.has_insertion():
                insertion_seq = event.insertion_seq1  # pick from the first read
                insertion_seq = (
                    reverse_complement(insertion_seq)
                    if event.strand1 == "-"
                    else insertion_seq
                )

                flag, insertion = self.blat.query_insertion(insertion_seq)

                insertion.query_name = read1.query_name
                if flag:  # only one hit
                    # add first node and insertion node
                    source_s = event.source_s1

                    # get type of insertion between first node and insertion node
                    insertion.update_cigarstring_sms(
                        read1.sms, source_s=source_s, source_strand=event.strand1
                    )
                    self.logger.trace(f"{insertion.strand=}, {insertion.cigarstring}")
                    insertion_mode = (
                        (2 if event.mode1 == 1 else 1)
                        if event.strand1 == insertion.strand
                        else event.mode1
                    )
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
                    insertion_mode = (
                        (2 if event.mode2 == 1 else 1)
                        if insertion.strand == read2.strand
                        else event.mode2
                    )

                    self.logger.trace("nls reference for read2 and insertion")
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
                        # change prev sv type for next node or Insertion
                        prev_sv_type = read1_node.sv_type
                        self.add_node(read1_node)
                        #  creat node for insertion
                        insertion_node = Node(
                            chrom=insertion.chrom,
                            strand=insertion.strand,
                            ref_start=insertion.ref_start,
                            ref_end=insertion.ref_end,
                            query_name=insertion.query_name,
                        )
                        (
                            insertion_node.exons,
                            insertion_node._introns,
                        ) = insertion.get_exons_and_introns()

                        insertion_node = (
                            insertion_read2_event.update_insertion_node_info(
                                insertion_node
                            )
                        )
                        # change prev sv type for next node or Insertion
                        insertion_node.prev_sv_type, prev_sv_type = (
                            prev_sv_type,
                            insertion_node.sv_type,
                        )

                        # add cigartuples_without_soft for insertion node
                        insertion_node.cigartuples_without_soft = (
                            insertion.cigartuples_without_soft
                        )

                        self.logger.trace(f"Add Insertion {insertion_node=} to Series")

                        self.add_node(insertion_node)

                else:  # no hits or multiple hits
                    self.logger.trace(f"Add Novel Insertion {insertion=} to read1")
                    # only add read1 with insertion info
                    # False means that the insertion type (hit more insertion) are
                    # not added in series
                    read1_node = event.update_node_info(False, read1_node, insertion)
                    self.add_node(read1_node)
            # no insertion and has microhomology
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
                    exons=event.read2_exons,  # type: ignore
                )
                final_node.query_name = read2.query_name
                final_node.prev_sv_type = prev_sv_type
                final_node.cigartuples_without_soft = read2.cigartuples_without_soft

                check_end_node_is_ploya(final_node, genome_fasta)

                self.add_node(final_node)


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
            self.strand1, self.strand2 = strands
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
        if self.annotation_code == 1:
            self.annotation_code = 2
        elif self.annotation_code == 2:
            self.annotation_code = 1
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
        return "left" if self.mode1 == 2 else "right"

    @property
    def source_s2(self) -> str:
        """Source of insertion of read2."""
        return "left" if self.mode2 == 2 else "right"

    def is_type_na(self) -> bool:
        """Return True if the event is NA."""
        return self.sv_type == "NA"

    def has_insertion(self) -> bool:
        """Return True if the event has insertion."""
        return bool(self.insertion_info[0].startswith("+"))

    def has_microhomology(self) -> bool:
        """Return True if the event has microhomology."""
        return bool(self.insertion_info[0].startswith("-"))

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

    def update_specific_info_within_event(
        self, node: Node, info_key_list: list[str]
    ) -> Node:
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
        insertion: Insertion | None | MicroHomology,
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
        )
        if is_update_insertion_info:
            new_node.insertion_info = (flag, insertion)  # type: ignore
        return new_node

    def update_insertion_node_info(self, insertion_node: Node) -> Node:
        """Update the information of insertion.

        :param insertion_node: the insertion to be updated
        :return: the updated insertion
        """
        if insertion_node.strand == "+":
            insertion_node.prev_breakpoint = BreakPoint.from_str(
                f"{insertion_node.chrom}:{insertion_node.ref_start}"
            )
            insertion_node.next_breakpoint = BreakPoint.from_str(
                f"{insertion_node.chrom}:{insertion_node.ref_end}"
            )
        else:
            insertion_node.prev_breakpoint = BreakPoint.from_str(
                f"{insertion_node.chrom}:{insertion_node.ref_end}"
            )
            insertion_node.next_breakpoint = BreakPoint.from_str(
                f"{insertion_node.chrom}:{insertion_node.ref_start}"
            )

        return self.update_specific_info_within_event(
            insertion_node,
            ["sv_type", "annotation_code", "splicing_code", "modes", "genes"],
        )


def reverse_complement(seq: str) -> str:
    """Obtain reverse complement sequence."""
    rctrans = str.maketrans("ACGT", "TGCA")
    return str.translate(seq, rctrans)[::-1]


def check_end_node_is_ploya(
    node: Node, genome_fasta: pyfaidx.Fasta, ratio: float = 0.7, length: int = 20
) -> None:
    """Check whether the node is bona fide polyA or internal priming events."""
    if node.ref_end is None or node.ref_start is None:
        raise SystemExit(f"{node} has no start or end position")

    if node.strand == "+":
        seq = genome_fasta[node.chrom][node.ref_end : node.ref_end + length].seq
    else:
        seq = genome_fasta[node.chrom][
            node.ref_start - length : node.ref_start
        ].reverse.complement.seq

    counter: dict[str, int] = Counter(seq)
    if counter["A"] <= ratio * len(seq):
        node.is_polya = True
