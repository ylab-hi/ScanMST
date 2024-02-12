from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import Enum, auto
from itertools import combinations
from typing import TYPE_CHECKING, Any, ClassVar

from loguru import logger

from scannls.base import (
    BreakPoint,
    Event,
    Exons,
    Introns,
    MappingMode,
    MicroHomology,
    NovelInsertion,
    Strand,
    reverse_complement,
)
from scannls.cli import infer_nls_from_connected_reads

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    import pyfaidx

    from scannls.base.basic_read import Read


class BasicNode:
    """BasicNode is used to represent nodes in the nlgraph."""

    __slots__ = (
        "successors",
        "predecessors",
        "merged_child_nodes",
        "merged_parent_nodes",
        "next_node_in_nlpath",
        "previous_node_in_nlpath",
        "previous_edge_in_nlapth",
        "is_merged",
        "is_in_graph",
        "is_traced",
        "trace_id",
    )

    def __init__(self) -> None:
        """Initialize BasicNode object."""
        self.successors: list[Node] = []
        self.predecessors: list[Node] = []

        self.merged_child_nodes: list[Node] = []
        self.merged_parent_nodes: list[Node] = []

        self.next_node_in_nlpath: Node | None = None
        self.previous_node_in_nlpath: Node | None = None
        self.previous_edge_in_nlapth: Edge | None = None
        self.is_merged, self.is_in_graph, self.is_traced = False, False, False
        self.trace_id: int = -1

    def __eq__(self, other) -> bool:
        """Compare two nodes in strict mode same memory address."""
        return id(self) == id(other)

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

    def add_successor_from_list(self, successors, graph, edge_data) -> None:
        """Add successor from list of Insertion object."""
        for successor in successors:
            self.add_successor(successor, graph, edge_data)

    def add_predecessor_from_list(self, predecessors, graph, edge_data) -> None:
        """Add predecessor from list of Insertion object."""
        for predecessor in predecessors:
            self.add_predecessor(predecessor, graph, edge_data)

    def _add_successor(self, successor, graph, edge_data) -> None:
        """Helper function to add successor to Insertion object."""
        self.successors.append(successor)

        if graph is not None and edge_data is not None:
            graph.add_edge(self, successor, edge_data)

        successor.add_predecessor(self, graph, edge_data)

    def _add_predecessor(self, predecessor, graph, edge_data) -> None:
        """Helper function to add predecessor to Insertion object."""
        self.predecessors.append(predecessor)

        predecessor.add_successor(self, graph, edge_data)

    def add_successor(self, successor, graph=None, edge_data=None) -> None:
        """Add successor to Insertion object."""
        if successor is not None and successor not in self.successors:
            if successor.is_in_graph:
                self._add_successor(successor, graph, edge_data)
            else:
                self.add_successor_from_list(
                    successor.merged_parent_nodes,
                    graph,
                    edge_data,
                )

            # WARN: miss only add edge <Yangyang Li>

    def add_predecessor(self, predecessor, graph=None, edge_data=None) -> None:
        """Node must be in the graph if the function is called.

        :param predecessor: predecessor of Insertion object
        """
        if predecessor is not None:
            if predecessor not in self.predecessors:
                if predecessor.is_in_graph:
                    self._add_predecessor(predecessor, graph, edge_data)
                else:
                    self.add_predecessor_from_list(
                        predecessor.merged_parent_nodes,
                        graph,
                        edge_data,
                    )
            elif graph is not None and edge_data is not None:
                # only merge nodes in the same graph
                graph.add_edge(predecessor, self, edge_data)

    def update_next_and_previous_node_in_nlpath(self, index: int, nlpath) -> None:
        """Update next and previous node in series."""
        if index == 0:
            self.next_node_in_nlpath = nlpath[index + 1]
        elif index == len(nlpath) - 1:
            self.previous_node_in_nlpath = nlpath[index - 1]
            self.previous_edge_in_nlapth = nlpath.edges[
                Edge.create_key_from_node(nlpath[index - 1], self)
            ]
        else:
            self.next_node_in_nlpath = nlpath[index + 1]
            self.previous_node_in_nlpath = nlpath[index - 1]
            self.previous_edge_in_nlapth = nlpath.edges[
                Edge.create_key_from_node(nlpath[index - 1], self)
            ]

    def clear_next_and_previous_node_in_series(self) -> None:
        """Clear next and previous node in series."""
        self.next_node_in_nlpath = None
        self.previous_node_in_nlpath = None


class NodeIdentity(Enum):
    HEAD = auto()
    TAIL = auto()
    MID = auto()

    @classmethod
    def from_str(cls, s) -> NodeIdentity:
        if s == "HEAD":
            return cls.HEAD
        if s == "TAIL":
            return cls.TAIL
        if s == "MID":
            return cls.MID

        msg = f"Invalid value for NodeIdentity: {s}"
        raise ValueError(msg)

    # fmt: off
    def is_head(self) -> bool: return self == NodeIdentity.HEAD
    def is_tail(self) -> bool: return self == NodeIdentity.TAIL
    def is_mid(self) -> bool: return self == NodeIdentity.MID
    def __hash__(self): return hash(self.name)
    # fmt: on


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

    .. example::

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
        "strand",
        "chrom",
        "_ref_start",
        "_ref_end",
        "exons",
        "_introns",
        "gene_names",
        "query_name",
        "_unique_key",
        "is_polya",
        "cigartuples_without_soft",
        "identities",
        *BasicNode.__slots__,
    )

    def __init__(
        self,
        query_name: str,
        chrom: str,
        strand: Strand | str,
        ref_start: int,
        ref_end: int,
        identity: NodeIdentity,
        exons: Exons,
        cigartuples_without_soft: list[int] | None = None,
    ) -> None:
        """Initialize a Node object."""
        super().__init__()
        self.query_name = query_name
        self.chrom = chrom
        self.strand = Strand.from_str(strand)
        self._ref_start = ref_start
        self._ref_end = ref_end

        self.exons = exons
        self._introns = exons.introns()

        self.gene_names: list[str] = []

        self.is_polya = False
        self.cigartuples_without_soft = cigartuples_without_soft
        self.identities: dict[str, NodeIdentity] = {self.query_name: identity}

        self._unique_key = f"{self.chrom}-{self.introns}-{self.ref_start}-{self.ref_end}-{self.strand}-{self.query_name}"

    @property
    def ref_start(self) -> int:
        return self._ref_start

    @ref_start.setter
    def ref_start(self, value: int) -> None:
        self._ref_start = value
        self.exons.first.start = value

    @property
    def ref_end(self) -> int:
        return self._ref_end

    @ref_end.setter
    def ref_end(self, value: int) -> None:
        self._ref_end = value
        self.exons.last.end = value

    def __repr__(self) -> str:
        """Get a string representation of a node."""
        return (
            f"Node({self.chrom}:{self.ref_start}-{self.ref_end}:{self.strand}, {self.trace_id=} {self.self_identity} "
            f"{self.exons!s}, read_ids={self.read_ids})"
        )

    __str__ = __repr__

    def __hash__(self) -> int:
        return (
            hash(self.chrom)
            ^ hash(self.ref_start)
            ^ hash(self.ref_end)
            ^ hash(self.strand)
            ^ hash(self.exons)
        )

    @property
    def read_ids(self) -> list[str]:
        return list(self.identities.keys())

    @property
    def self_identity(self) -> NodeIdentity:
        return self.identities[self.query_name]

    @self_identity.setter
    def self_identity(self, identity: NodeIdentity) -> None:
        self.identities[self.query_name] = identity

    def identity(self, query_name: str) -> NodeIdentity | None:
        ret = self.identities.get(query_name)
        if ret is None:
            logger.warning(f"Identity for {query_name} not found")

        return ret

    @property
    def introns(self) -> Introns | None:
        """Get introns of a node."""
        if self._introns is not None:
            return self._introns

        self._introns = self.exons.introns()
        return self._introns

    @property
    def similar_key(self) -> str:
        """Get similar key of a node."""
        introns = self.introns
        chosen_intron = None
        if introns:
            chosen_intron = introns.last if self.strand.is_forward() else introns.first

        key = f"{chosen_intron!s}"
        return f"{self.chrom}_{key}"

    @property
    def exons_length(self) -> int:
        """Get total length of exon of a node."""
        return sum(len(exon) for exon in self.exons)

    @property
    def unique_key(self) -> str:
        return self._unique_key

    def merge(
        self,
        other: Node,
    ) -> None:
        """Merge two nodes.

        Merge two nodes for ref_start, ref_end, and identities.
        """
        if isinstance(other, Node):
            # update exon coordinates
            self.ref_start = min(
                self.exons.first.start,
                other.exons.first.start,
            )
            # WARN: ref_end may be not consistent with prev_breakpoint of next edge <Yangyang Li>
            self.ref_end = max(
                self.exons.last.end,
                other.exons.last.end,
            )

            if other.query_name in self.read_ids:
                logger.warning(
                    f"A circle in a path is detectd {other.query_name}",
                )

            # WARN: do not check if they have same key <Yangyang Li>
            self.identities.update(other.identities)

            return

        msg = f"Cannot merge {self!r} and {other!r}"
        raise TypeError(msg)

    def contains(
        self,
        other: Node,
        *,
        same_left=False,
        same_right=False,
        check_introns=False,
        threshold=0,
    ):
        if isinstance(other, Node):
            if check_introns and self.introns != other.introns:
                return False

            if not same_left and not same_right:
                return (
                    self.exons.first.start
                    <= other.exons.first.start
                    <= other.exons.last.end
                    <= self.exons.last.end
                )

            if same_left and same_right:
                return (
                    abs(self.exons.first.start - other.exons.first.start) <= threshold
                    and self.exons.last.end == other.exons.last.end
                )

            if same_left:
                return (
                    abs(self.exons.first.start - other.exons.first.start) <= threshold
                    and self.exons.last.end - other.exons.last.end >= -threshold
                )

            if same_right:
                return (
                    -threshold <= self.exons.first.start - other.exons.first.start
                    and abs(self.exons.last.end - other.exons.last.end) <= threshold
                )

        msg = f"{other} is not Node"
        raise ValueError(msg)


class VariationType(Enum):
    TRA = auto()
    DEL = auto()
    TDUP = auto()
    INV = auto()
    IDUP = auto()

    @classmethod
    def from_str(cls, s):
        if s == "TRA":
            return cls.TRA
        if s == "DEL":
            return cls.DEL
        if s == "TDUP":
            return cls.TDUP
        if s == "INV":
            return cls.INV
        if s == "IDUP":
            return cls.IDUP

        msg = f"Invalid SV type: {s}"
        raise ValueError(msg)

    # fmt: off
    def __str__(self) -> str: return self.name
    def __repr__(self) -> str: return self.name
    def is_tra(self) -> bool: return self == self.TRA
    def is_del(self) -> bool: return self == self.DEL
    def is_tdup(self) -> bool: return self == self.TDUP
    def is_inv(self) -> bool: return self == self.INV
    def is_idup(self) -> bool: return self == self.IDUP
    # fmt: on


@dataclass
class EdgeData:
    variantion_type: VariationType
    break_point1: BreakPoint
    break_point2: BreakPoint
    sr: int
    read_ids: list[str]
    mode1: MappingMode
    mode2: MappingMode

    gene1: str | None  # find by breakpoint1
    gene2: str | None  # find by breakpoint2
    annotation_code: int
    splicing_code: int

    insertion_info: Any | None = None
    original_sr: int = 1

    @classmethod
    def from_event(cls, event: Event, read_id: str) -> EdgeData:
        return cls(
            variantion_type=VariationType.from_str(event.sv_type),
            break_point1=BreakPoint.from_str(event.bp1),
            break_point2=BreakPoint.from_str(event.bp2),
            sr=1,
            read_ids=[read_id],
            mode1=MappingMode.from_int(event.mode1),
            mode2=MappingMode.from_int(event.mode2),
            gene1=event.genes[0],
            gene2=event.genes[1],
            annotation_code=event.annotation_code,
            splicing_code=event.splicing_code,
        )

    def equal(
        self,
        other: EdgeData,
        *,
        compared_break_point: bool,
        break_point_threshold: int,
    ):
        """Check if two edges are equal."""
        flag = self.variantion_type == self.variantion_type

        if compared_break_point:
            return (
                flag
                and self.break_point1.equal(
                    other.break_point1,
                    break_point_threshold,
                )
                and self.break_point2.equal(other.break_point2, break_point_threshold)
            )

        return flag


class Edge:
    def __init__(self, node1_key: str, node2_key: str, edge_data: EdgeData) -> None:
        """Initializes a new instance of the Edge class.

        Args:
            node1_key (str): The key of the first node connected by the edge.
            node2_key (str): The key of the second node connected by the edge.
            edge_data (EdgeData): The data of the edge.
        """
        self.node1_key = node1_key
        self.node2_key = node2_key
        self.edge_data = edge_data

    def __repr__(self) -> str:
        return f"Edge(data={self.edge_data})"

    # fmt: off
    @property
    def key(self): return f"{self.node1_key}-{self.node2_key}"
    @property
    def break_point1(self): return self.edge_data.break_point1
    @break_point1.setter
    def break_point1(self, value: BreakPoint): self.edge_data.break_point1 = value
    @property
    def break_point2(self): return self.edge_data.break_point2
    @break_point2.setter
    def break_point2(self, value: BreakPoint): self.edge_data.break_point2 = value
    @property
    def insertion_info(self): return self.edge_data.insertion_info
    @insertion_info.setter
    def insertion_info(self, value): self.edge_data.insertion_info = value
    @property
    def variation_type(self): return self.edge_data.variantion_type
    @property
    def sr(self): return self.edge_data.sr
    @sr.setter
    def sr(self, value): self.edge_data.sr = value
    @property
    def original_sr(self): return self.edge_data.original_sr
    @original_sr.setter
    def original_sr(self, value): self.edge_data.original_sr = value
    @property
    def read_ids(self): return self.edge_data.read_ids
    @property
    def mode1(self): return self.edge_data.mode1
    @property
    def mode2(self): return self.edge_data.mode2
    @property
    def modes(self): return self.mode1, self.mode2
    @property
    def gene1(self): return self.edge_data.gene1
    @gene1.setter
    def gene1(self, value): self.edge_data.gene1 = value
    @property
    def gene2(self): return self.edge_data.gene2
    @gene2.setter
    def gene2(self, value): self.edge_data.gene2 = value
    @property
    def annotation_code(self): return self.edge_data.annotation_code
    @annotation_code.setter
    def annotation_code(self, value): self.edge_data.annotation_code = value
    @property
    def splicing_code(self): return self.edge_data.splicing_code
    @splicing_code.setter
    def splicing_code(self, value): self.edge_data.splicing_code = value
    # fmt: on

    @staticmethod
    def create_key_from_node(node1: Node, node2: Node) -> str:
        return f"{node1.unique_key}-{node2.unique_key}"

    def add_read_id(self, read_id: str):
        if read_id not in self.edge_data.read_ids:
            self.edge_data.read_ids.append(read_id)

    def merge(self, other: Edge, pnode_strand: Strand, nnode_strand: Strand):
        # WARN: update breakpoint in cmparing way <07-03-23, Yangyang Li>

        self.beak_point1 = (
            min(
                [self.break_point1, other.break_point1],
                key=lambda x: x.pos,
            )
            if pnode_strand.is_forward()
            else max(
                [self.break_point1, other.break_point1],
                key=lambda x: x.pos,
            )
        )

        self.break_point2 = (
            min(
                [self.break_point2, other.break_point2],
                key=lambda x: x.pos,
            )
            if nnode_strand.is_forward()
            else max(
                [self.break_point2, other.break_point2],
                key=lambda x: x.pos,
            )
        )

        merge_insertion(self, other)

        self.sr += other.sr
        self.edge_data.read_ids.extend(other.read_ids)

        self.gene1 = other.gene1
        self.gene2 = other.gene2
        self.annotation_code = other.annotation_code
        self.splicing_code = other.splicing_code

        logger.trace(
            f"Merge edge {self.key} {self.read_ids=} with {other.key} {other.read_ids=}.",
        )

    def get_nodes(self, graph):
        node1 = graph.get_node_with_unique_key(self.node1_key)
        node2 = graph.get_node_with_unique_key(self.node2_key)

        if node1 is None or node2 is None:
            msg = "Node not found due to edge is invalidated"
            raise ValueError(msg)

        return node1, node2

    @classmethod
    def from_nodes(
        cls,
        node1: Node,
        node2: Node,
        edge_data: EdgeData,
    ):
        return cls(node1.unique_key, node2.unique_key, edge_data)

    def merged(
        self,
        other_edge: Edge,
        *,
        compared_break_point: bool = True,
        break_point_threshold: int = 10,
    ) -> bool:
        """Check if two edges are merged."""
        return self.edge_data.equal(
            other_edge.edge_data,
            compared_break_point=compared_break_point,
            break_point_threshold=break_point_threshold,
        )


class NLPath:
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
    ... chrom='chr17',ref_start=7706250,ref_end=7708250,exons=[(7706250,7708250)],sv_type='TDUP'))
    >>> series.add_node(Node(prev_bp='chr17:7701656',next_bp='chr17:7702552',strand='+',
    ... chrom='chr17',ref_start=7701656,ref_end=7702552,exons=[(7701656, 7702552)], sv_type='TRA'))
    >>> series.add_node(Node(prev_bp='chr1:15872815',next_bp='chr1:15876678',strand='+',
    ... chrom='chr1',ref_start=15872815,ref_end=15876678,exons=[(15872815,15876678)],
    ... sv_type='TDUP'))
    >>> series.add_node(Node(prev_bp='chr1:15777169',next_bp=None,strand='+',
    ... chrom='chr1',ref_start=15777169,ref_end=15777589,exons=[(15777169,15777589)],
    ... sv_type=None))
    >>> series
    Series(
        Node(chr17:7706250-7708250:+, 7706250-7708250, TDUP, None, chr17:7708250)
        Node(chr17:7701656-7702552:+, 7701656-7702552, TRA, chr17:7701656, chr17:7702552)
        Node(chr1:15872815-15876678:+, 15872815-15876678, TDUP, chr1:15872815, chr1:15876678)
        Node(chr1:15777169-15777589:+, 15777169-15777589, None, chr1:15777169, None) )

    >>> series_with_novel_insertion = Series(blat=None, logger=logger)
    >>> series_with_novel_insertion.nodes = [ Node(prev_bp=None,next_bp='chr17:7702552',
    ... strand='+',chrom='chr17',ref_start=7701656,ref_end=7702552,exons=[(7701656, 7702552)],
    ... sv_type='TRA', insertion_info=(False, NovelInsertion(hit_num=1,
    ... query_sequence='ATCGATCG'))), Node(prev_bp='chr1:15872815',next_bp=None,strand='+',
    ... chrom='chr1',ref_start=15872815,ref_end=15876678,exons=[(15872815,15876678)],
    ... sv_type=None)]
    >>> series_with_novel_insertion
    Series(
        Node(chr17:7701656-7702552:+, 7701656-7702552, TRA, None, chr17:7702552)
        Node(chr1:15872815-15876678:+, 15872815-15876678, None, chr1:15872815, None) )
    """

    reorder_conditions_dict: ClassVar[dict[str, bool]] = {
        "+-11": True,
        "+-22": False,
        "-+11": False,
        "-+22": True,
        "++12": True,
        "++21": False,
        "--12": False,
        "--21": True,
    }

    def __init__(
        self,
        nodes: list[Node],
    ) -> None:
        """Initialize a nlpath object."""
        self.nodes = nodes
        self.edges: dict[str, Edge] = {}
        self.is_in_graph = False
        self.id = -1
        self.merge_factor = 1
        self.extension = False

    def add_edge(self, nodes: Node, nodet: Node, edge: Edge | None = None) -> None:
        """Add edge to the path."""
        if nodes not in self.nodes:
            self.nodes.append(nodes)

        if nodet not in self.nodes:
            self.nodes.append(nodet)

        if edge is not None:
            self.edges[edge.key] = edge

    def only_add_edge(self, nodes: Node, nodet: Node, edge: Edge) -> None:
        key = Edge.create_key_from_node(nodes, nodet)
        if key not in self.edges:
            self.edges[key] = edge

    def next_edge(self, nodes: Node, nodes_idx: int | None = None) -> Edge | None:
        if nodes_idx is None:
            nodes_idx = self.nodes.index(nodes)

        if nodes_idx < len(self.nodes) - 1:
            return self.edges.get(
                Edge.create_key_from_node(nodes, self.nodes[nodes_idx + 1]),
            )

        return None

    def remove_edge(
        self,
        nodes: Node,
        nodet: Node | None = None,
        *,
        nodes_idx: int | None = None,
    ):
        if nodet is not None:
            key = Edge.create_key_from_node(nodes, nodet)
            self.edges.pop(key)
        else:
            if nodes_idx is None:
                nodes_idx = self.nodes.index(nodes)

            if nodes_idx < len(self.nodes) - 1:
                key = Edge.create_key_from_node(nodes, self.nodes[nodes_idx + 1])
                self.edges.pop(key)

    def get_edge(
        self,
        nodes: Node,
        nodet: Node | None = None,
        *,
        nodes_idx: int | None = None,
    ) -> Edge | None:
        """Get edge from the path."""
        if nodet is not None:
            key = Edge.create_key_from_node(nodes, nodet)
            return self.edges.get(key)

        return self.next_edge(nodes, nodes_idx)

    def is_all_type_del(self) -> bool:
        """Check if sv_type of all nodes in the series are DEL."""
        return all(
            edge.variation_type == VariationType.DEL for edge in self.edges.values()
        )

    def squeeze(self) -> None:
        """Squeeze nodes whose edge is del in the path."""
        logger.trace(f"Squeeze {self!r}")

        if not self.nodes:
            return

        if not any(
            edge.variation_type == VariationType.DEL for edge in self.edges.values()
        ):
            return

        edges = []
        for idx, node in enumerate(self.nodes):
            edge = self.get_edge(node, nodes_idx=idx)
            if edge is not None:
                edges.append(edge)

        new_edges = []
        new_nodes = []
        new_nodes.append(self.nodes[0])

        for idx, edge in enumerate(edges):
            prev_node = new_nodes[-1]
            next_node = self.nodes[idx + 1]

            if edge.variation_type.is_del():
                if prev_node.strand.is_forward() and next_node.strand.is_forward():
                    prev_node.exons.extend(next_node.exons)
                    prev_node.ref_end = next_node.ref_end
                    prev_node._introns = None
                else:
                    next_node.exons.extend(prev_node.exons)
                    next_node.ref_end = prev_node.ref_end
                    next_node._introns = None
                    new_nodes[-1] = next_node

                # WARN: cigartuples_without_soft is not update <Yangyang Li>
            else:
                new_edges.append(edge)
                new_nodes.append(next_node)

        for idx, node in enumerate(new_nodes):
            if idx == 0:
                node.self_identity = NodeIdentity.HEAD
            elif idx == len(new_nodes) - 1:
                node.self_identity = NodeIdentity.TAIL
            else:
                node.self_identity = NodeIdentity.MID

        self.nodes = new_nodes
        self.edges.clear()

        for idx, node in enumerate(self.nodes[:-1]):
            next_node = self.nodes[idx + 1]
            edge_key = Edge.create_key_from_node(node, next_node)
            edge = new_edges[idx]
            self.edges[edge_key] = edge

    def is_forming_circle(self, threshold: int = 20) -> bool:
        """Check if this nlpath itself can form a circle."""
        pair_indices = combinations(range(len(self.nodes)), 2)
        for _a, _b in pair_indices:
            node_a = self.nodes[_a]
            node_b = self.nodes[_b]
            intron_condition = (
                node_a.introns
                and node_b.introns
                and len(node_a.introns) > 0
                and len(node_b.introns) > 0
            ) or (not node_a.introns and not node_b.introns)

            if (
                node_a.strand == node_b.strand
                and node_a.chrom == node_b.chrom
                and abs(node_a.ref_start - node_b.ref_start) <= threshold
                and abs(node_a.ref_end - node_b.ref_end) <= threshold
            ) and intron_condition:
                return True
            # middle node vs. tail node
            if (
                _a > 0
                and _b == len(self.nodes) - 1
                and node_a.strand == node_b.strand
                and node_a.chrom == node_b.chrom
                and abs(node_a.ref_start - node_b.ref_start) <= threshold
                or abs(node_a.ref_end - node_b.ref_end) <= threshold
            ) and intron_condition:
                return True

        return False

    def is_maximum_novel_insertion_length_valid(self, threshold: int = 50) -> bool:
        """Check if nlpath with maximum insertion length > threshold, which indicates sequencing artifacts."""
        maximum_insertion_length = 0
        for event_id, _node in enumerate(self.nodes[:-1], 1):
            _edge = self.next_edge(_node, event_id - 1)
            if _edge.insertion_info and isinstance(
                _edge.insertion_info[1], NovelInsertion
            ):
                insertion = _edge.insertion_info[1]
                _insertion_length = len(insertion.query_sequence)
                if _insertion_length > maximum_insertion_length:
                    maximum_insertion_length = _insertion_length

        return maximum_insertion_length <= threshold

    def is_minimum_node_length_larger_than_threshold(self, threshold: int = 10) -> bool:
        """Check if minimum length of all nodes in the series > threshold."""
        return min(_node.exons_length for _node in self.nodes) > threshold

    def sum_sr(self) -> int:
        """Get sum of sr for all nodes in the series."""
        return sum(edge.sr for edge in self.edges.values())

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
        string = "\nNLPath("
        space = " " * 4
        for idx, n in enumerate(self.nodes):
            string += f"\n{space}\u001b[36m{n!r}"
            if idx < len(self.nodes) - 1:
                edge = self.get_edge(n, nodes_idx=idx)
                string += f"\n{space}\u001b[31m{edge!r}\u001b[36m"
        string += ")"
        return string

    def __iter__(self) -> Iterator[Node]:
        """Return an iterator over the events."""
        yield from self.nodes

    @property
    def unique_key(self) -> str:
        """Return the unique key of the event."""
        return "".join([node.unique_key for node in self.nodes])

    @staticmethod
    def reorder_event(evt: Event):
        """Order breakpoint pairs following the transcription direction using.

        information of reads 'mode' and 'strand'
        +1;-1 => up;down
        +2;-2 => down;up
        """
        is_bp1_upstream = NLPath.reorder_conditions_dict.get(
            f"{evt.strand1}{evt.strand2}{evt.mode1}{evt.mode2}",
            None,
        )

        if not is_bp1_upstream:
            logger.trace(f"reorder event {evt=} {is_bp1_upstream=}")
            evt.reverse()

        return evt

    @staticmethod
    def order_events_by_trancription_direction(event_list: list[Event]):
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
        output_event_list: list[Event] = []
        for index, evt in enumerate(event_list):
            parsed_evt = NLPath.reorder_event(evt)

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

    def anno_extension(self):
        read_ids_set = []
        for edge in self.edges.values():
            read_ids_set.append(set(edge.read_ids))

        if not set.intersection(*read_ids_set):
            self.extension = True

    @classmethod
    def create_path_from_node_edge_list(cls, node_edges):
        """Create a path from a list of nodes and edges.

        Args:
            node_edges: list of nodes and edges
        """
        instance = cls(nodes=[])

        for idx in range(0, len(node_edges), 2):
            current_node = node_edges[idx]
            if not isinstance(current_node, Node):
                msg = f"Expected Node, got {type(current_node)} from {current_node}"
                raise TypeError(msg)

            instance.nodes.append(current_node)

            if idx < len(node_edges) - 1:
                current_edge = node_edges[idx + 1]
                instance.edges[current_edge.key] = current_edge

        instance.anno_extension()
        return instance

    @classmethod
    def from_nodes_and_edges_data(cls, nodes, edges_data):
        instance = cls(nodes=nodes)

        for idx in range(len(nodes) - 1):
            edge = Edge.from_nodes(nodes[idx], nodes[idx + 1], edges_data[idx])
            instance.edges[edge.key] = edge

        return instance

    @classmethod
    def new(
        cls,
        events: list[Event],
        read_chains,
        splice_bin,
        genome_fasta,
        cvg,
        gene_iv,
        motif_required,
        aligner,
    ) -> NLPath:
        """Create a nlpath from a list of events."""

        events = NLPath.order_events_by_trancription_direction(events)

        logger.trace(f"reordered event:{events=}")
        edges_data: list[EdgeData] = []
        nodes = []

        events_len = len(events)

        for index, event in enumerate(events):
            read1: Read = event.read1(read_chains)
            read2: Read = event.read2(read_chains)

            read1_node = Node(
                query_name=read1.query_name,
                chrom=event.chrom1,
                strand=event.strand1,
                ref_start=event.read1_ref_start,
                ref_end=event.read1_ref_end,
                identity=NodeIdentity.HEAD if index == 0 else NodeIdentity.MID,
                exons=Exons.from_list(event.read1_exons),
                cigartuples_without_soft=read1.cigartuples_without_soft,
            )

            edge_data = EdgeData.from_event(event, read_id=read1.query_name)

            logger.trace(f"{read1=} {read2=}")

            # is insertions
            if event.has_insertion():
                insertion_seq = event.insertion_seq1  # pick from the first read
                insertion_seq = (
                    reverse_complement(insertion_seq)
                    if event.strand1.is_reverse()
                    else insertion_seq
                )

                if aligner is None:
                    flag, insertion = False, NovelInsertion(
                        hit_num=0, query_sequence=insertion_seq
                    )
                else:
                    flag, insertion = aligner.query_insertion(insertion_seq)

                insertion.query_name = read1.query_name

                if flag:  # only one hit
                    # add first node and insertion node
                    source_s = event.source_s1

                    # get type of insertion between first node and insertion node
                    insertion.update_cigarstring_sms(
                        read1.sms,
                        source_s=source_s,
                        source_strand=event.strand1,
                    )

                    logger.trace(f"{insertion.strand=}, {insertion.cigarstring}")
                    insertion_mode = (
                        (2 if event.mode1 == 1 else 1)
                        if event.strand1 == insertion.strand
                        else event.mode1
                    )
                    logger.trace("nls reference for read1 and insertion")

                    read1_insertion_event = infer_nls_from_connected_reads(
                        read_lt=read1,
                        read_rt=insertion,
                        lt_mode=event.mode1,
                        rt_mode=insertion_mode,
                        splice_bin=splice_bin,
                        genome_fasta=genome_fasta,
                        cvg=cvg,
                        gene_iv=gene_iv,
                        motif_required=motif_required,
                    )

                    # get type of insertion between insertion node and second node
                    insertion_mode = (
                        (2 if event.mode2 == 1 else 1)
                        if insertion.strand == read2.strand
                        else event.mode2
                    )

                    logger.trace("nls reference for read2 and insertion")
                    insertion_read2_event = infer_nls_from_connected_reads(
                        read_lt=insertion,
                        read_rt=read2,
                        lt_mode=insertion_mode,
                        rt_mode=event.mode2,
                        splice_bin=splice_bin,
                        genome_fasta=genome_fasta,
                        cvg=cvg,
                        gene_iv=gene_iv,
                        motif_required=motif_required,
                    )

                    if read1_insertion_event is None or insertion_read2_event is None:
                        # only add read1, False means that the insertion type (hit 1 insertion)
                        # are not added in series

                        nodes.append(read1_node)

                        edge_data.insertion_info = (False, insertion)
                        edges_data.append(edge_data)

                    else:
                        # add read1 and insertion
                        # True means that the insertion type(hit 1 insertion) are added in series
                        read1_insertion_event = Event(read1_insertion_event)
                        insertion_read2_event = Event(insertion_read2_event)

                        edge_data = EdgeData.from_event(
                            read1_insertion_event,
                            read1.query_name,
                        )
                        edge_data.insertion_info = (True, insertion)

                        #  creat node for insertion
                        insertion_node = Node(
                            query_name=insertion.query_name,
                            chrom=insertion.chrom,
                            strand=insertion.strand,
                            ref_start=insertion.ref_start,
                            ref_end=insertion.ref_end,
                            identity=NodeIdentity.MID,
                            exons=insertion.get_exons(),
                            cigartuples_without_soft=insertion.cigartuples_without_soft,
                        )

                        insertion_edge_data = EdgeData.from_event(
                            insertion_read2_event,
                            read1.query_name,
                        )

                        nodes.append(read1_node)
                        logger.trace(
                            f"auxiliary alignment[4] is effective here. reads_name:{read1.query_name} query_sequence:{insertion_seq}"
                        )
                        logger.trace(f"Add Insertion {insertion_node=} to path")
                        nodes.append(insertion_node)

                        edges_data.append(edge_data)
                        edges_data.append(insertion_edge_data)

                else:  # no hits or multiple hits
                    logger.trace(f"Add Novel Insertion {insertion=} to read1")
                    # only add read1 with insertion info
                    # False means that the insertion type (hit more insertion) are
                    # not added in series
                    nodes.append(read1_node)
                    edge_data.insertion_info = (False, insertion)
                    edges_data.append(edge_data)

            # no insertion and has microhomology
            elif event.has_microhomology():
                # add read 1 with on insertion
                microhomology = MicroHomology(event.insertion_seq1)

                logger.trace(f"Add MicroHomology {microhomology=} to read1")
                if event.strand1.is_reverse():
                    microhomology.reverse_completement_query()

                nodes.append(read1_node)
                edge_data.insertion_info = (False, microhomology)
                edges_data.append(edge_data)

            else:
                nodes.append(read1_node)
                edges_data.append(edge_data)

            # add final node
            if index == events_len - 1:
                final_node = Node(
                    query_name=read2.query_name,
                    chrom=event.chrom2,
                    strand=event.strand2,
                    ref_start=event.read2_ref_start,
                    ref_end=event.read2_ref_end,
                    identity=NodeIdentity.TAIL,
                    exons=Exons.from_list(event.read2_exons),
                    cigartuples_without_soft=read2.cigartuples_without_soft,
                )
                check_end_node_is_ploya(final_node, genome_fasta)

                nodes.append(final_node)

        return cls.from_nodes_and_edges_data(nodes, edges_data)


def update_node_with_other_node(
    node: Node,
    other_node: Node,
    features: Iterable[str],
) -> None:
    """Update node with another node.

    if current feature of node is None, then use another node's feature.
    """
    for feature in features:
        if getattr(node, feature) is None:
            setattr(node, feature, getattr(other_node, feature))


def check_end_node_is_ploya(
    node: Node,
    genome_fasta: pyfaidx.Fasta,
    ratio: float = 0.7,
    length: int = 20,
) -> None:
    """Check whether the node is bona fide polyA or internal priming events."""
    if node.ref_end is None or node.ref_start is None:
        msg = f"{node} has no start or end position"
        raise SystemExit(msg)

    if node.strand.is_forward():
        seq = genome_fasta[node.chrom][node.ref_end : node.ref_end + length].seq
    else:
        seq = genome_fasta[node.chrom][
            node.ref_start - length : node.ref_start
        ].reverse.complement.seq

    counter: dict[str, int] = Counter(seq)
    if counter["A"] <= ratio * len(seq):
        node.is_polya = True


def _check_insertion_conditions_for_compare_insertion(
    insertion_info1,
    insertion_info2,
) -> bool:
    if insertion_info1 is None and insertion_info2 is None:
        return True

    if insertion_info1 is not None and insertion_info2 is not None:
        if insertion_info1[0] and insertion_info2[0]:
            # 1 hit insertion that is added in the series
            return True

        if not insertion_info1[0] and not insertion_info2[0]:
            if isinstance(insertion_info1[1], MicroHomology) and isinstance(
                insertion_info2[1],
                MicroHomology,
            ):
                return True

            if (
                isinstance(insertion_info1[1], NovelInsertion)
                and isinstance(insertion_info2[1], NovelInsertion)
                and (
                    insertion_info1[1].query_sequence
                    == insertion_info2[1].query_sequence
                )
            ):
                return True

    return False


def merge_insertion(edge1: Edge, edge2: Edge):
    """edge1 merge edge2."""
    if edge1.insertion_info is not None and edge2.insertion_info is not None:
        if edge1.sr > edge2.sr:
            return

        if (edge1.sr < edge2.sr) or (
            isinstance(edge1.insertion_info[1], MicroHomology)
            and isinstance(
                edge2.insertion_info[1],
                NovelInsertion,
            )
        ):
            edge1.insertion_info = edge2.insertion_info
