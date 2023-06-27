from typing import Optional, Any
from enum import auto
from enum import Enum
from dataclasses import dataclass
from ..base.basicClass import NovelInsertion, MicroHomology


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

        self.next_node_in_series: Optional[Node] = None
        self.previous_node_in_series: Optional[Node] = None
        self.is_merged, self.is_in_graph, self.is_traced = False, False, False
        self.trace_id: int = -1
        self.sr: int = 1
        self.original_sr: int = 1
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
                    successor.merged_parent_nodes, graph, edge_data
                )

    def add_predecessor(self, predecessor, graph=None, edge_data=None) -> None:
        """Node must be in the graph if the function is called.

        :param predecessor: predecessor of Insertion object
        """
        if predecessor is not None and predecessor not in self.predecessors:
            if predecessor.is_in_graph:
                self._add_predecessor(predecessor, graph, edge_data)
            else:
                self.add_predecessor_from_list(
                    predecessor.merged_parent_nodes, graph, edge_data
                )

    def update_next_and_previous_node_in_series(self, index: int, series) -> None:
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


class NodeIdentity(Enum):
    HEAD = auto()
    TAIL = auto()
    MID = auto()

    @classmethod
    def from_str(cls, s) -> "NodeIdentity":
        if s == "HEAD":
            return cls.HEAD
        elif s == "TAIL":
            return cls.TAIL
        elif s == "MID":
            return cls.MID
        else:
            raise ValueError("Invalid value for NodeIdentity: {}".format(s))

    @classmethod
    def from_node(cls, node: "Node") -> "NodeIdentity":
        if node.prev_breakpoint is not None and node.next_breakpoint is not None:
            return cls.MID
        elif node.prev_breakpoint is None:
            return cls.HEAD
        elif node.next_breakpoint is None:
            return cls.TAIL
        else:
            raise ValueError("Invalid node identity: {}".format(node))

    def is_head(self) -> bool:
        return self == NodeIdentity.HEAD

    def is_tail(self) -> bool:
        return self == NodeIdentity.TAIL

    def is_mid(self) -> bool:
        return self == NodeIdentity.MID


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
        "identity",
        "read_names",
        *BasicNode.__slots__,
    )

    def __init__(
        self,
        prev_bp: Optional[str] = None,
        next_bp: Optional[str] = None,
        strand: Optional[str] = None,
        chrom: Optional[str] = None,
        ref_start: Optional[int] = None,
        ref_end: Optional[int] = None,
        exons: Optional[list[Any]] = None,
        sv_type: Optional[str] = None,
        annot: Optional[int] = None,
        canonical: Optional[int] = None,
        modes: Optional[list[int]] = None,
        genes: Optional[tuple[str, str]] = None,
        query_name: str = "",
    ) -> None:
        """Initialize a Node object."""
        super().__init__()  # initialize BasicNode object
        self._introns = None
        self.chrom = chrom
        self.query_name = query_name
        self.prev_breakpoint = BreakPoint.from_str(prev_bp)
        self.next_breakpoint = BreakPoint.from_str(next_bp)
        self.prev_breakpoint_depth: Optional[int] = None
        self.next_breakpoint_depth: Optional[int] = None
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
        self.unique_key: Optional[str] = None
        self.is_polya = False
        self.cigartuples_without_soft: Optional[list[int]] = None
        self.identity: dict[str, NodeIdentity] = {}
        self.read_names = [self.query_name]

        if (
            self.query_name != ""
            and self.prev_breakpoint is not None
            and self.next_breakpoint is not None
        ):
            self.identity[self.query_name] = NodeIdentity.from_node(self)

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

    @property
    def self_identity(self) -> NodeIdentity:
        assert self.query_name != ""

        if self.query_name not in self.identity:
            self.identity[self.query_name] = NodeIdentity.from_node(self)

        return self.identity[self.query_name]

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

        assert self.exons is not None

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

    @property
    def length(self) -> int:
        """Get total length of exon of a node."""
        assert self.exons is not None
        total_len = 0
        for i, j in self.exons:
            total_len += j - i
        return total_len

    def get_unique_key(self):
        """Get unique key of a node."""
        introns = self.introns

        key = "-".join([f"{i}-{j}" for i, j in introns]) if introns else "None"

        key = f"{self.chrom}-{key}-{self.sv_type}-{self.prev_breakpoint}-{self.next_breakpoint}"

        if self.insertion_info is not None:
            _, insertion_type = self.insertion_info  # type: ignore
            if insertion_type.__class__.__name__ in ("NovelInsertion", "MicroHomology"):
                key = f"{insertion_type.query_sequence}-{key}"

        self.unique_key = key
        return key

    def update_identity(self):
        self.identity[self.query_name] = NodeIdentity.from_node(self)

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

    # NOTE: may be removed in the future  <04-17-23, Yangyang Li>


class SpliceType(Enum):
    """Splice Type.

    used in prune
    """

    forward = auto()
    backward = auto()


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
        elif s == "DEL":
            return cls.DEL
        elif s == "TDUP":
            return cls.TDUP
        elif s == "INV":
            return cls.INV
        elif s == "IDUP":
            return cls.IDUP
        else:
            raise ValueError("Invalid SV type: {}".format(s))

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name


class Variation:
    def __init__(
        self, types: VariationType, break_point: BreakPoint, break_point_depth: int
    ):
        self.types = types
        self.break_point = break_point
        self.break_point_depth = break_point_depth

    def __repr__(self):
        return f"Variation({self.types=} {self.break_point=} {self.break_point_depth=})"

    @classmethod
    def from_node(cls, node: Node):
        assert node.next_breakpoint is not None
        next_breakpoint_depth = (
            0 if node.next_breakpoint_depth is None else node.next_breakpoint_depth
        )

        return cls(
            VariationType.from_str(node.sv_type),
            node.next_breakpoint,
            next_breakpoint_depth,
        )

    @staticmethod
    def is_merged(
        variation1: "Variation", variation2: "Variation", threshold: int
    ) -> bool:
        return variation1.types == variation2.types and BreakPoint.equal(
            variation1.break_point, variation2.break_point, threshold
        )


@dataclass
class EdgeData:
    variation: Variation
    sr: int
    insertion: Any
    read_ids: list[str]

    @classmethod
    def from_node(cls, node: Node):
        return cls(
            Variation.from_node(node), node.sr, node.insertion_info, [node.query_name]
        )


class Edge:
    def __init__(self, node1_key: str, node2_key: str, edge_data: EdgeData) -> None:
        """
        Initializes a new instance of the Edge class.

        Args:
            node1_key (str): The key of the first node connected by the edge.
            node2_key (str): The key of the second node connected by the edge.
            sv (Variation): The data associated with the edge.
        """
        self.node1_key = node1_key
        self.node2_key = node2_key
        self.edge_data = edge_data

    def __repr__(self) -> str:
        return (
            f"Edge(variation={self.variation}, sr={self.sr}, read_ids={self.read_ids})"
        )

    @property
    def key(self):
        return f"{self.node1_key}-{self.node2_key}"

    @property
    def insertion(self):
        return self.edge_data.insertion

    @property
    def variation(self):
        return self.edge_data.variation

    @property
    def sr(self):
        return self.edge_data.sr

    @property
    def read_ids(self):
        return self.edge_data.read_ids

    @staticmethod
    def create_key_from_node(node1: Node, node2: Node) -> str:
        if node1.unique_key is None and node2.unique_key is None:
            raise ValueError("Both nodes have no unique key")
        return f"{node1.unique_key}-{node2.unique_key}"

    def add_read_id(self, read_id: str):
        if read_id not in self.edge_data.read_ids:
            self.edge_data.read_ids.append(read_id)

    def updated(self, other: "Edge"):
        # WARN:  Do not update variation with break point <06-12-23>
        self.edge_data.sr += other.sr
        self.edge_data.read_ids.extend(other.read_ids)
        if self.edge_data.insertion and isinstance(
            self.edge_data.insertion[1], (NovelInsertion, MicroHomology)
        ):
            self.edge_data.insertion[1].increment_ao()

    def get_nodes(self, graph):
        node1 = graph.get_node_with_unique_key(self.node1_key)
        node2 = graph.get_node_with_unique_key(self.node2_key)

        if node1 is None or node2 is None:
            raise ValueError("Node not found due to edge is invalidated")

        return node1, node2

    @classmethod
    def from_nodes(
        cls,
        node1: Node,
        node2: Node,
        edge_data: Optional[EdgeData],
    ):
        assert node1.unique_key is not None
        assert node2.unique_key is not None

        if edge_data is None:
            edge_data = EdgeData.from_node(node1)

        return cls(node1.unique_key, node2.unique_key, edge_data)

    @classmethod
    def from_node_key(
        cls,
        node1_key: str,
        node2_key: str,
        edge_data: EdgeData,
    ):
        return cls(node1_key, node2_key, edge_data)

    @staticmethod
    def is_merged(edge1: "Edge", edge2: "Edge", break_point_threshold: int) -> bool:
        return (
            edge1.key == edge2.key
            and Variation.is_merged(
                edge1.variation, edge2.variation, break_point_threshold
            )
            and _check_insertion_conditions_for_compare_insertion(
                edge1.insertion, edge2.insertion
            )
        )
