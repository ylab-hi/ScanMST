"""
@Author:      YangyangLi
@Filename:    spliceGraph.py
@Time:        12/15/21 10:42 AM
"""
import copy
import types
from collections import defaultdict
from collections.abc import Iterable
from collections.abc import Iterator
from dataclasses import dataclass
from enum import auto
from enum import Enum
from typing import Any
from typing import Optional
from typing import Union

from .basicClass import BreakPoint
from .basicClass import MicroHomology
from .basicClass import Node
from .basicClass import NodeIdentity
from .basicClass import NovelInsertion
from .basicClass import Series
from .mergeCondition import (
    _compare_is_merged_helper_check_condition_for_head_and_middle_nodes_mode,
)
from .mergeCondition import (
    _compare_is_merged_helper_check_condition_for_head_and_tail_nodes_mode,
)
from .mergeCondition import (
    _compare_is_merged_helper_check_condition_for_tail_and_middle_nodes_mode,
)
from .mergeCondition import (
    _compare_is_merged_helper_check_condition_for_two_heads_nodes_mode,
)
from .mergeCondition import (
    _compare_is_merged_helper_check_condition_for_two_tail_nodes_mode,
)
from .plotGraph import plot_graph
from .srRescuer import SRRescuer
from .type import LoggerType


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
        return variation1.types == variation2.types and BreakPoint.equals(
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
        self.edge_data.sr += other.sr
        self.edge_data.read_ids.extend(other.read_ids)
        if self.edge_data.insertion and isinstance(
            self.edge_data.insertion[1], (NovelInsertion, MicroHomology)
        ):
            self.edge_data.insertion[1].increment_ao()

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


class SpliceGraph:
    """SpliceGraph class is used to trace the path of the splice graph."""

    dict_factory = dict
    list_factory = list

    def __init__(
        self, logger: LoggerType, rescuer: Any, prune_threshold: int = 10
    ) -> None:
        """Initialize SpliceGraph."""
        self.logger = logger
        self.prune_threshold = prune_threshold
        self.dict_factory = SpliceGraph.dict_factory  # type: ignore
        self.list_factory = SpliceGraph.list_factory  # type: ignore
        self.rescuer = rescuer

    def __call__(
        self,
        series_list: Iterable[Series],
        clique_ind: int,
        is_plot: bool = True,
        is_check_circle: bool = False,
    ) -> Iterable[Series]:
        """Find a specific path based on splice graph.

        :param series_list: series list

        :Example:

        >>> from loguru import logger
        >>> splice_graph = SpliceGraph(logger)
        >>> splice_graph(series_list)
        """
        if isinstance(series_list, types.GeneratorType):
            series_list = list(series_list)

        self.series_list = copy.deepcopy(series_list)

        del series_list  # remove reference to series_list
        self.nodes: dict[str, list[Node]] = self.dict_factory()

        self.edges: dict[str, list[Edge]] = defaultdict(list)

        # construct splice graph
        self.construct()
        # sr rescuer
        # self.rescuer(self)

        self.logger.trace(f"Splice Graph Node: {sum(1 for _ in self)}")

        # self.prune()

        # export_graph(self, Path(f"graph_{clique_ind}.adj"))

        if is_plot:
            plot_graph(self, f"_{clique_ind}", False)

        # trace path
        current_nodes_keys: set[str] = set()
        for node_list in self.trace():
            yield Series.create_path_from_node_edge_list(
                node_list,
            )

        if is_check_circle:
            # check circle in graph
            for node_list in self.check_circle_in_graph(current_nodes_keys):
                yield Series.create_series_from_node_list(
                    node_list, self.logger, set(), is_add_key=False
                )

    @classmethod
    def create_graph(
        cls,
        input_bam: str,
        mapq: int,
        soft_len: int,
        mismatch: int,
        alignment_fraction: float,
        logger: LoggerType,
        prune_threshold: int,
        node_rescued_sr_maximum: int,
        average_read_depth: Optional[int],
    ) -> "SpliceGraph":
        """Create splice graph."""
        rescuer = SRRescuer(
            input_bam,
            mapq,
            soft_len,
            mismatch,
            alignment_fraction,
            node_rescued_sr_maximum,
            average_read_depth,
        )

        return cls(logger, rescuer, prune_threshold)

    def add_edge(self, node1: Node, node2: Node, edge_data=None):
        """Add edge from node1 -> node2"""
        edge = Edge.from_nodes(node1, node2, edge_data)

        if self.edges.get(edge.key) is None:
            self.edges[edge.key].append(edge)
            return

        is_merged = False
        for current_edge in self.edges[edge.key]:
            if Edge.is_merged(edge, current_edge, self.prune_threshold):
                current_edge.updated(edge)
                is_merged = True
                break

        if not is_merged:
            self.edges[edge.key].append(edge)

    def find_edges(self, node1: Node, node2: Node):
        edge_key = Edge.create_key_from_node(node1, node2)
        if self.edges.get(edge_key) is None:
            raise KeyError(f"Edge {edge_key} not found in splice graph.")
        else:
            return self.edges[edge_key]

    @staticmethod
    def get_node_identity_base_edge(
        edge: Edge, node: Node
    ) -> dict[NodeIdentity, list[str]]:
        result = defaultdict(list)
        for read_id in edge.read_ids:
            result[node.identity[read_id]].append(read_id)
        return result

    @staticmethod
    def determine_edge(previous_edge_node_identity, next_edge_node_identity) -> bool:
        previous_edge_read_id_mid = set(previous_edge_node_identity[NodeIdentity.MID])
        next_edge_read_id_mid = set(next_edge_node_identity[NodeIdentity.MID])

        if not previous_edge_read_id_mid and not next_edge_read_id_mid:
            # two mid from different reads
            return len(previous_edge_node_identity & next_edge_node_identity) > 0

        return True

    def get_possible_edges(
        self, current_path: list[Union[Node, Edge]], current_node: Node, successor: Node
    ) -> list[Edge]:
        if not current_path:
            return self.find_edges(current_node, successor)

        previous_edge = current_path[-1]
        assert isinstance(previous_edge, Edge)

        previous_edge_node_identity = self.get_node_identity_base_edge(
            previous_edge, current_node
        )

        possible_edges = []

        for edge in self.find_edges(current_node, successor):
            edge_node_identity = self.get_node_identity_base_edge(edge, current_node)
            if self.determine_edge(previous_edge_node_identity, edge_node_identity):
                possible_edges.append(edge)

        return possible_edges

    def __contains__(self, node: Node) -> bool:
        """Check if node is in a graph.

        :param node: node to be checked
        :return: True if node is in graph, otherwise False

        .. note::
            please ensure that node has run :function: `Node.get_unique_key` before. So
            node.unique_key is not None.

        """
        return any(
            other_node.unique_key == node.unique_key
            for other_node in self.get_nodes_with_similar_key(node.similar_key)
        )

    def __iter__(self) -> Iterator[Node]:
        """Iterate over all nodes in graph."""
        for nodes in self.nodes.values():
            yield from nodes

    def print_path(self) -> None:
        """Print path based on splice graph."""
        for _node_list in self.trace():
            pass

    def get_start_nodes(self) -> Iterable[Node]:
        """Get start nodes based if node has predecessors."""
        return (
            node
            for nodes in self.nodes.values()
            for node in nodes
            if node.is_start_node()
        )

    def get_end_nodes(self) -> Iterable[Node]:
        """Get end nodes based if node has successors."""
        return (
            node
            for nodes in self.nodes.values()
            for node in nodes
            if node.is_end_node()
        )

    def remove_node(self, node: Node) -> None:
        """Remove node from graph.

        :param node: node to be removed
        """
        if node.similar_key is None:
            raise ValueError(f"node.similar_key is None, {node.query_name}")
        self.get_nodes_with_similar_key(node.similar_key).remove(node)

    def reset_trace_id(self) -> None:
        """Reset trace id for all nodes."""
        for node in self:
            node.reset_trace_id()

    def get_node_with_unique_key(self, unique_key: str) -> Optional[Node]:
        """Get node with unique key.

        :param unique_key: unique key
        :return: node with unique key
        """
        for node in self:
            if node.unique_key == unique_key:
                return node
        return None

    def get_nodes_with_similar_key(self, similar_key: str) -> list[Node]:
        """Get nodes in graph with similar key."""
        return self.nodes.get(similar_key, [])

    def add_node_with_similar_key(self, node: Node) -> None:
        """Add node to the splice graph.

        :param node: node to be added
        """
        if similar_nodes := self.get_nodes_with_similar_key(node.similar_key):
            similar_nodes.append(node)
        else:
            self.nodes[node.similar_key] = [node]

    @staticmethod
    def _compare_is_merged_helper(node1: Node, node2: Node, threshold: int) -> bool:
        """Check if node1 and node2 can be merged.

        .. note::
            End note will not merge with start/middle node,
            since every end node has polyA tail in library preparation.
        """
        if node1.strand != node2.strand:
            return False

        node1_self_identity: NodeIdentity = node1.self_identity
        node2_self_identity: NodeIdentity = node2.self_identity
        if (
            node1_self_identity.is_tail() and node2_self_identity.is_head()
        ):  # node1 is end node, node2 is start node
            return (
                _compare_is_merged_helper_check_condition_for_head_and_tail_nodes_mode(
                    node1, node2
                )
            )
        elif (
            node1_self_identity.is_tail() and node2_self_identity.is_mid()
        ):  # node1 is end node, node2 is middle node
            return _compare_is_merged_helper_check_condition_for_tail_and_middle_nodes_mode(
                node1, node2
            )

        elif (
            node1_self_identity.is_head() and node2_self_identity.is_head()
        ):  # both are start nodes
            return _compare_is_merged_helper_check_condition_for_two_heads_nodes_mode(
                node1, node2, threshold
            )

        elif (
            node1_self_identity.is_tail() and node2_self_identity.is_tail()
        ):  # both are end nodes  # check first exon start
            return _compare_is_merged_helper_check_condition_for_two_tail_nodes_mode(
                node1, node2
            )

        elif (
            node1_self_identity.is_head() and node2_self_identity.is_mid()
        ):  # node1 is start node, node2 is middle node
            return _compare_is_merged_helper_check_condition_for_head_and_middle_nodes_mode(
                node1, node2
            )

        elif (
            node1_self_identity.is_mid() and node2_self_identity.is_mid()
        ):  # both are middle nodes
            return (
                node1.exons[0][0] == node2.exons[0][0]  # type: ignore
                and node1.exons[-1][1] == node2.exons[-1][1]  # type: ignore
            )

        # swap node1 and node2 to check if they can be merged again
        return False

    @staticmethod
    def _compare_is_merged(node1: Node, node2: Node, threshold: int) -> bool:
        """Node1 is similar as node2 is precommit of the function.

         compare if node1 can merge node2
        :param node1: node1
        :param node2: node2
        :return:
        """

        # NOTE: may not swap order <04-24-23, Yangyang Li>
        if SpliceGraph._compare_is_merged_helper(node1, node2, threshold):
            return True

        if SpliceGraph._compare_is_merged_helper(node2, node1, threshold):
            return True

        return False

    def _check_if_current_node_is_merged_in_similar_nodes_in_graph(
        self, current_node: Node, similar_key: str, merged_nodes_pool: set[Node]
    ) -> None:
        """Check if current node is merged in similar nodes in graph."""
        # get similar nodes in the graph
        similar_nodes_in_graph = self.get_nodes_with_similar_key(similar_key)

        # iterate all similar nodes in the graph
        for similar_node_in_graph in similar_nodes_in_graph:
            # check if the current node is merged into a similar node in the graph
            if SpliceGraph._compare_is_merged(
                similar_node_in_graph, current_node, self.prune_threshold
            ):
                current_node.is_merged = True

                update_exon_coord_name_mode(similar_node_in_graph, current_node)

                # nodes in merged_parent_nodes are all in the graph
                current_node.merged_parent_nodes.append(similar_node_in_graph)
                merged_nodes_pool.add(current_node)

                # only consider nodes that have been processed: previous node in current series
                # keeps in mind the next node in current series is not processed yet!!!!
                # a -> b and b <- a
                edge_data = (
                    EdgeData.from_node(current_node.previous_node_in_series)
                    if current_node.previous_node_in_series is not None
                    else None
                )

                similar_node_in_graph.add_predecessor(
                    current_node.previous_node_in_series, self, edge_data
                )

    def _check_if_current_node_added_in_graph_and_update_predecessor_successor(
        self,
        current_node: Node,
        similar_key: str,
        merged_nodes_pool: set[Node],
    ) -> None:
        """Check if the current node is added in graph and update a predecessor and successor."""
        if not current_node.is_merged:
            current_node.is_in_graph = True  # check if node is in graph
            self.add_node_with_similar_key(current_node)

            # only consider nodes that have been processed: previous node in series
            # keeps in mind the next node in series is not processed yet!!!!
            edge_data = (
                EdgeData.from_node(current_node.previous_node_in_series)
                if current_node.previous_node_in_series is not None
                else None
            )

            current_node.add_predecessor(
                current_node.previous_node_in_series,
                self,
                edge_data,
            )

    def construct(self) -> None:
        """Main function to construct graph."""
        # iterate all series
        merged_nodes_pool: set[Node] = set()
        self.logger.debug(f"Input Clique {self.series_list=}")
        for series in self.series_list:
            # iterate all nodes in series
            for index, current_node in enumerate(series):
                self.logger.trace(f"{current_node=}")

                # add information about next and previous node in series to current node
                current_node.update_next_and_previous_node_in_series(index, series)

                # initialize and get unique key of current node and set node.unique_key
                # if not set when you reach node.unique_key, will return None
                _ = current_node.get_unique_key()

                # get node identity mid, tail, head
                current_node.update_identity()

                # get a similar key(chrom and first intron) of current node
                similar_key = current_node.similar_key

                self._check_if_current_node_is_merged_in_similar_nodes_in_graph(
                    current_node, similar_key, merged_nodes_pool
                )

                self._check_if_current_node_added_in_graph_and_update_predecessor_successor(
                    current_node, similar_key, merged_nodes_pool
                )

                current_node.clear_next_and_previous_node_in_series()

    def _trace_forward(
        self,
        start_node: Node,
        trace_id: int,
        path: list[Union[Node, Edge]],
        group_paths: list[list[Union[Node, Edge]]],
    ) -> None:
        """Helper function to trace through graph and find all paths.

        .. seealso::
            :func:`SpliceGraph.trace`
        """

        if start_node in path:
            self.logger.warning(
                f"A circle is found in the graph {start_node} in {path}"
            )

        if not start_node or start_node in path:
            # successor be [] or None
            group_paths.append(path)

        else:
            if successors := start_node.successors:
                for successor in successors:
                    successor.set_trace_id(trace_id)

                    # successor.set_harmoic_mean_sr(successor.sr)
                    for edge in self.get_possible_edges(path, start_node, successor):
                        self._trace_forward(
                            successor,
                            trace_id + 1,
                            [*path, start_node, edge],
                            group_paths,
                        )
            else:
                # successor be [] or None
                self._trace_forward(
                    successors,  # type: ignore
                    trace_id + 1,
                    [*path, start_node],
                    group_paths,
                )

    def _trace_backward(
        self,
        end_node: Node,
        trace_id: int,
        path: list[Union[Node, Edge]],
        group_paths: list[list[Union[Node, Edge]]],
    ) -> None:
        """Helper function to trace through graph and find all paths.

        .. seealso::
            :func:`SpliceGraph.trace`
        """
        if not end_node or end_node in path:
            # successor be [] or None
            group_paths.append(path)
        else:
            if predecessors := end_node.successors:
                for predecessor in predecessors:
                    predecessor.set_trace_id(trace_id)
                    # successor.set_harmoic_mean_sr(successor.sr)
                    for edge in self.get_possible_edges(path, end_node, predecessor):
                        self._trace_backward(
                            predecessor,
                            trace_id + 1,
                            [*path, end_node, edge],
                            group_paths,
                        )
            else:
                # successor be [] or None
                self._trace_backward(
                    predecessors,  # type: ignore
                    trace_id + 1,
                    [*path, end_node],
                    group_paths,
                )

    def _trace(self, direction: Enum) -> None:
        """Trace splice graph but only mark node with trace_id.

        :param direction: direction of trace, forward or backward

        .. note::
            before trace, you should call :func:`SpliceGraph.reset_trace_id` if
            you have already traced splice graph.
        """
        if direction == SpliceType.forward:
            for start_node in self.get_start_nodes():
                start_node.set_trace_id(1)
                # start_node.set_harmoic_mean_sr(start_node.sr)
                self._trace_forward(start_node, 2, [], [])
            return
        elif direction == SpliceType.backward:
            for end_node in self.get_end_nodes():
                end_node.set_trace_id(1)
                end_node.set_harmoic_mean_sr(end_node.sr)
                self._trace_backward(end_node, 2, [], [])
            return

        raise ValueError(f"{direction=} is not a valid direction[forward, backward]")

    def trace(self) -> Any:
        """Trace forward through graph and find all paths."""
        result_series_list = []

        if not self.get_start_nodes() and len(self.nodes.values()) > 0:
            self.logger.warning(f"A circle may exist in graph {self.nodes.values()}")

        for start_node in self.get_start_nodes():
            start_node.set_trace_id(1)
            group_paths: Any = []
            self._trace_forward(start_node, 2, [], group_paths)
            result_series_list.extend(group_paths)

        return result_series_list

    def create_same_level_node_list(self) -> list[list[Node]]:
        """Create node trace id dict.

        :return: trace_id: List[node] dict
        """
        node_trace_id_dict: dict[int, list[Node]] = defaultdict(list)
        for node in self:
            node_trace_id_dict[node.trace_id].append(node)

        return [i for i in node_trace_id_dict.values() if len(i) > 1]

    def _rule_out(self, winner: Node, loser: Node) -> None:
        """Rule out the loser and add sr to the winner.

        Loser is out, and its sr, successors, predecessors are added to the winner.
        """
        winner.update_sr(loser.sr)
        winner.add_successor_from_list(loser.successors, None, None)
        winner.add_predecessor_from_list(loser.predecessors, None, None)
        update_node_with_other_node(
            winner,
            loser,
            (
                "splicing_code",
                "annotation_code",
                "genes",
            ),
        )

        for loser_predecessor in loser.predecessors:
            loser_predecessor.successors.remove(loser)

        for loser_successor in loser.successors:
            loser_successor.predecessors.remove(loser)
        self.remove_node(loser)

    @staticmethod
    def _check_can_battle_condition(
        breakpoint1: Optional[BreakPoint],
        breakpoint2: Optional[BreakPoint],
        threshold: int,
    ) -> bool:
        """Check if two breakpoints are in the threshold.

        :param breakpoint1
        :param breakpoint2
        :return: True or False
        """
        if breakpoint1 is None or breakpoint2 is None:
            return False
        return abs(breakpoint1.pos - breakpoint2.pos) < threshold

    def check_can_battle(self, node_a: Node, node_b: Node) -> bool:
        """Check if two nodes can battle."""
        self.logger.trace(f"{node_a=}\n{node_b=}")
        if node_a.prev_sv_type != node_b.prev_sv_type:
            return False

        if node_a.introns != node_b.introns:
            return False

        if node_a.prev_breakpoint is None and node_b.prev_breakpoint is None:
            return SpliceGraph._check_can_battle_condition(
                node_a.next_breakpoint, node_b.next_breakpoint, self.prune_threshold
            )

        if node_a.next_breakpoint is None and node_b.next_breakpoint is None:
            return SpliceGraph._check_can_battle_condition(
                node_a.prev_breakpoint, node_b.prev_breakpoint, self.prune_threshold
            )

        return SpliceGraph._check_can_battle_condition(
            node_a.prev_breakpoint, node_b.prev_breakpoint, self.prune_threshold
        ) and SpliceGraph._check_can_battle_condition(
            node_a.next_breakpoint, node_b.next_breakpoint, self.prune_threshold
        )

    def _begin_battle(self, node_a: Node, node_b: Node) -> tuple[bool, ...]:
        """Begin battle between two nodes.

        :return: Two bool values:
                value1: True if node_a and node_b can battle.
                value2: if node_a is winner, return True, else return False
        """
        self.logger.trace(
            f"{node_a.harmonic_mean_sr=:.2f}\t{node_b.harmonic_mean_sr=:.2f}"
        )
        if (
            node_a.harmonic_mean_sr == node_b.harmonic_mean_sr
            or not self.check_can_battle(node_a, node_b)
        ):
            return False, False

        if node_a.harmonic_mean_sr > node_b.harmonic_mean_sr:
            self._rule_out(node_a, node_b)
            return True, True

        # node_a.original_sr < node_b.original_sr
        self._rule_out(node_b, node_a)
        return True, False

    def battle(self, same_level_node_list: list[list[Node]]) -> None:
        """Nodes with same trace id battle each other.

        Node with larger number of sr wins, otherwise lose.

        Loser will be rule out.
        """
        for node_list in same_level_node_list:
            while node_list:
                current_node = node_list.pop()
                for other_node in node_list:
                    can_battle, is_winner = self._begin_battle(current_node, other_node)
                    if can_battle:
                        if not is_winner:
                            # other_node is winner
                            break
                        # winner is last_node
                        node_list.remove(other_node)

    def _prune(self, direction: Enum) -> None:
        """Implement function to prune graph.

        :param direction: direction of prune, forward or backward
        """
        # 1. trace and mark node with trace_id
        self._trace(direction)
        # 2. save every trace_id and its corresponding node to be Dict
        # 3. check the length of a node list, if the number of nodes is less than 2, remove it
        same_level_node_list = self.create_same_level_node_list()
        # 4. compare them and rule out loser
        self.battle(same_level_node_list)

    def prune(self) -> None:
        """Prune graph.

        Algorithm:
        1. trace graph and mark every node with trace id
        2. save every node with trace id as a Dict[int, List[Node]]
        3. check nodes of Dict in terms of trace id if len(values)>1,
        4. then compare them and rule out loser in terms of sr number
        """
        self._prune(SpliceType.forward)
        self.reset_trace_id()
        self._prune(SpliceType.backward)

    def check_circle_in_graph(self, nodes_keys: set[str]):
        """Check if there is a circle in graph."""
        all_nodes_keys: set[str] = set()
        result_paths: list[list[Node]] = []

        for node in self:
            if (key := node.unique_key) is not None:
                all_nodes_keys.add(key)

        self.check_circle_in_graph_helper(all_nodes_keys - nodes_keys, result_paths)
        return result_paths

    def check_circle_in_graph_helper(
        self, nodes_keys: set[str], result_paths: list[list[Node]]
    ) -> None:
        """Check if there is a circle in graph."""
        if (
            nodes_keys
            and (start_node := self.get_node_with_unique_key(nodes_keys.pop()))
            is not None
        ):
            current_nodes_keys: set[str] = set()
            self._trace_forward_record_node_unique_keys(
                start_node, [], result_paths, current_nodes_keys
            )
            self.check_circle_in_graph_helper(
                nodes_keys - current_nodes_keys, result_paths
            )

    def _trace_forward_record_node_unique_keys(
        self,
        start_node: Node,
        path: list[Node],
        group_paths: list[list[Node]],
        nodes_keys: set[str],
    ) -> None:
        """Helper function to trace through graph and find all paths.

        .. seealso::
            :func:`SpliceGraph.trace`
        """
        if not start_node or start_node in path:
            group_paths.append(path)
        else:
            if successors := start_node.successors:
                for successor in successors:
                    if (key := successor.unique_key) is not None:
                        nodes_keys.add(key)
                    self._trace_forward_record_node_unique_keys(
                        successor, [*path, start_node], group_paths, nodes_keys
                    )
            else:
                self._trace_forward_record_node_unique_keys(
                    successors,  # type: ignore
                    [*path, start_node],
                    group_paths,
                    nodes_keys,
                )


def update_exon_coord_name_mode(updated_node: Node, current_node: Node) -> None:
    """Update exon coordinates of the updated node based on the current node.

    :param updated_node: node has been inserted into graph
    :param current_node: node has not been inserted into graph
    :return: None
    """
    if updated_node.exons is None or current_node.exons is None:
        raise ValueError(f"{updated_node} or {current_node} has no exons")

    _update_exon_coord_sr_svtype_breakpoints_name_mode_in_same_exons(
        updated_node, current_node
    )


def _update_exon_coord_sr_svtype_breakpoints_name_mode_in_same_exons(
    updated_node: Node, current_node: Node
) -> None:
    """Update exon coordinates of the updated node based on the current node.

    Two nodes with the same number of exons

    :param updated_node:  node has been inserted into graph
    :param current_node: node has not been inserted into graph
    :return: None
    """

    # update exon coordinates
    updated_node.ref_start = min(  # type: ignore
        updated_node.exons[0][0], current_node.exons[0][0]  # type: ignore
    )

    updated_node.exons[0] = updated_node.ref_start, updated_node.exons[0][1]  # type: ignore

    updated_node.ref_end = max(  # type: ignore
        updated_node.exons[-1][1], current_node.exons[-1][1]  # type: ignore
    )

    updated_node.exons[-1] = updated_node.exons[-1][0], updated_node.ref_end  # type: ignore

    update_node_with_other_node(
        updated_node,
        current_node,
        (
            "splicing_code",
            "annotation_code",
            "genes",
            "modes",
        ),
    )

    updated_node.read_names.append(current_node.query_name)
    updated_node.identity[current_node.query_name] = current_node.self_identity


def _check_insertion_conditions_for_compare_insertion(
    insertion_info1, insertion_info2
) -> bool:
    if insertion_info1 is None and insertion_info2 is None:
        return True

    elif insertion_info1 is not None and insertion_info2 is not None:
        if insertion_info1[0] and insertion_info2[0]:
            # 1 hit insertion that is added in the series
            return True

        if not insertion_info1[0] and not insertion_info2[0]:
            if (
                isinstance(insertion_info1[1], NovelInsertion)
                and isinstance(insertion_info2[1], NovelInsertion)
                and (
                    insertion_info1[1].query_sequence
                    == insertion_info2[1].query_sequence
                )
            ):
                return True

            elif isinstance(insertion_info1[1], MicroHomology) and isinstance(
                insertion_info2[1], MicroHomology
            ):
                return True

    return False


def _check_insertion_conditions_for_compare(node1: Node, node2: Node) -> bool:
    """Check if node1 and node2 can be merged based on insertion info."""
    insertion_info1 = node1.insertion_info
    insertion_info2 = node2.insertion_info

    return _check_insertion_conditions_for_compare_insertion(
        insertion_info1, insertion_info2
    )


def update_node_with_other_node(
    node: Node, other_node: Node, features: Iterable[str]
) -> None:
    """Update node with another node.

    if current feature of node is None, then use another node's feature.
    """
    for feature in features:
        if getattr(node, feature) is None:
            setattr(node, feature, getattr(other_node, feature))
