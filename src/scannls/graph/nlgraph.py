"""
@Author:      YangyangLi
@Filename:    nlgraph.py
@Time:        12/15/21 10:42 AM.
"""
from __future__ import annotations

import copy
import types
from collections import defaultdict
from typing import TYPE_CHECKING, Any

from loguru import logger

from .basic_graph import (
    Edge,
    NLPath,
    Node,
    NodeIdentity,
    update_node_with_other_node,
)
from .graphvis import plot_graph
from .merge_condition import MergeCondition
from .sr_rescuer import SRRescuer

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from scannls.type import LoggerType


class NLGraph:
    """SpliceGraph class is used to trace the path of the splice graph."""

    dict_factory = dict
    list_factory = list

    def __init__(
        self,
        logger: LoggerType,
        rescuer: Any,
        prune_threshold: int = 10,
        support_reads: int = 1,
    ) -> None:
        """Initialize SpliceGraph."""
        self.logger = logger
        self.prune_threshold = prune_threshold
        self.support_reads = support_reads
        self.dict_factory = NLGraph.dict_factory  # type: ignore
        self.list_factory = NLGraph.list_factory  # type: ignore
        self.rescuer = rescuer

    def __call__(
        self,
        nlpath_list: Iterable[NLPath],
        clique_ind: int,
        *,
        is_plot: bool,
    ) -> Iterable[NLPath]:
        """Find a specific path based on splice graph.

        :param series_list: series list

        .. example:

        >>> from loguru import logger
        >>> splice_graph = SpliceGraph(logger)
        >>> splice_graph(series_list)
        """
        if isinstance(nlpath_list, types.GeneratorType):
            nlpath_list = list(nlpath_list)

        self.nlpaths = copy.deepcopy(nlpath_list)

        del nlpath_list  # remove reference to series_list
        self.nodes: dict[str, list[Node]] = self.dict_factory()

        self.edges: dict[str, list[Edge]] = defaultdict(list)

        # construct splice graph
        self.construct()

        # sr rescuer
        self.logger.trace(f"NLGraph Node: {sum(1 for _ in self)}")

        if is_plot:
            plot_graph(self, f"_{clique_ind}", is_matplotlib=False)

        # trace path
        for node_list in self.trace():
            yield NLPath.create_path_from_node_edge_list(
                node_list,
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
        support_reads: int,
        node_rescued_sr_maximum: int,
        average_read_depth: int | None,
    ) -> NLGraph:
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

        return cls(logger, rescuer, prune_threshold, support_reads)

    def add_edge(self, node1: Node, node2: Node, edge_data):
        """Add edge from node1 -> node2."""
        edge = Edge.from_nodes(node1, node2, edge_data)

        if self.edges.get(edge.key) is None:
            self.edges[edge.key].append(edge)
            return

        is_merged = False
        for current_edge in self.edges[edge.key]:
            if current_edge.merged(
                edge,
                compared_break_point=True,
                break_point_threshold=self.prune_threshold,
            ):
                current_edge.merge(edge)
                is_merged = True
                break

        if not is_merged:
            self.edges[edge.key].append(edge)

    def find_edges(self, node1: Node, node2: Node):
        edge_key = Edge.create_key_from_node(node1, node2)
        if self.edges.get(edge_key) is None:
            msg = f"Edge {edge_key} not found in splice graph."
            raise KeyError(msg)
        return self.edges[edge_key]

    @staticmethod
    def get_node_identity_base_edge(
        edge: Edge,
        node: Node,
    ) -> dict[NodeIdentity, list[str]]:
        result = defaultdict(list)
        logger.debug(f"{node=}")
        for read_id in edge.read_ids:
            node_identity = node.identities[read_id]
            logger.debug(f"{node_identity=}")
            result[node_identity].append(read_id)
        return result

    @staticmethod
    def determine_edge(previous_edge_node_identity, next_edge_node_identity) -> bool:
        previous_edge_read_id_mid = set(previous_edge_node_identity[NodeIdentity.MID])
        next_edge_read_id_mid = set(next_edge_node_identity[NodeIdentity.MID])

        if not previous_edge_read_id_mid and not next_edge_read_id_mid:
            # check two mid from different reads
            return len(previous_edge_read_id_mid & next_edge_read_id_mid) > 0

        return True

    def get_possible_edges(
        self,
        current_path: list[Node | Edge],
        current_node: Node,
        successor: Node,
    ) -> Iterable[Edge]:
        if len(current_path) == 1:
            previous_edge_node_identity = None
        else:
            previous_edge = current_path[-2]
            if not isinstance(previous_edge, Edge):
                msg = f"previous_edge is not instance of Edge: {previous_edge}"
                raise TypeError(msg)

            previous_edge_node_identity = self.get_node_identity_base_edge(
                previous_edge,
                current_node,
            )

        edges = []
        for edge in self.find_edges(current_node, successor):
            if edge.sr >= self.support_reads:
                edge_node_identity = self.get_node_identity_base_edge(
                    edge,
                    current_node,
                )

                if previous_edge_node_identity is not None:
                    if self.determine_edge(
                        previous_edge_node_identity,
                        edge_node_identity,
                    ):
                        edges.append(edge)
                else:
                    edges.append(edge)

        return edges

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
            msg = f"node.similar_key is None, {node.query_name}"
            raise ValueError(msg)
        self.get_nodes_with_similar_key(node.similar_key).remove(node)

    def get_node_with_unique_key(self, unique_key: str) -> Node | None:
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
    def _compare_is_merged(node1: Node, node2: Node, threshold: int) -> bool:
        """Check if node1 and node2 can be merged.

        .. note::
            End note will not merge with start/middle node,
            since every end node has polyA tail in library preparation.
        """
        if node1.strand != node2.strand:
            return False
        merge_condition = MergeCondition(threshold)
        return merge_condition.merged(node1, node2)

    def _check_if_current_node_is_merged_in_similar_nodes_in_graph(
        self,
        current_node: Node,
        similar_key: str,
    ) -> None:
        """Check if current node is merged in similar nodes in graph."""

        # iterate all similar nodes in the graph
        for similar_node_in_graph in self.get_nodes_with_similar_key(similar_key):
            # check if the current node is merged into a similar node in the graph
            if NLGraph._compare_is_merged(
                similar_node_in_graph,
                current_node,
                self.prune_threshold,
            ):
                current_node.is_merged = True

                self.logger.warning(
                    f"merging node {similar_node_in_graph} and {current_node})",
                )

                update_exon_coord_name_mode(similar_node_in_graph, current_node)

                # nodes in merged_parent_nodes are all in the graph
                current_node.merged_parent_nodes.append(similar_node_in_graph)

                # only consider nodes that have been processed: previous node in current series
                # keeps in mind the next node in current series is not processed yet!!!!
                # a -> b and b <- a

                edge_data = (
                    current_node.previous_edge_in_nlapth.edge_data
                    if current_node.previous_edge_in_nlapth is not None
                    else None
                )

                similar_node_in_graph.add_predecessor(
                    current_node.previous_node_in_nlpath,
                    self,
                    edge_data,
                )

    def _check_if_current_node_added_in_graph_and_update_predecessor_successor(
        self,
        current_node: Node,
    ) -> None:
        """Check if the current node is added in graph and update a predecessor and successor."""
        if not current_node.is_merged:
            current_node.is_in_graph = True  # check if node is in graph
            self.add_node_with_similar_key(current_node)

            # only consider nodes that have been processed: previous node in series
            # keeps in mind the next node in series is not processed yet!!!!

            edge_data = (
                current_node.previous_edge_in_nlapth.edge_data
                if current_node.previous_edge_in_nlapth is not None
                else None
            )

            current_node.add_predecessor(
                current_node.previous_node_in_nlpath,
                self,
                edge_data,
            )

    def construct(self) -> None:
        """Main function to construct graph."""
        # iterate all series
        self.logger.debug(f"Input Cluster {self.nlpaths=}")
        for nlpath in self.nlpaths:
            # iterate all nodes in series
            for index, current_node in enumerate(nlpath):
                self.logger.trace(f"{current_node=}")

                # add information about next and previous node in series to current node
                current_node.update_next_and_previous_node_in_nlpath(index, nlpath)

                # get a similar key(chrom and first intron) of current node
                similar_key = current_node.similar_key

                self._check_if_current_node_is_merged_in_similar_nodes_in_graph(
                    current_node,
                    similar_key,
                )
                self._check_if_current_node_added_in_graph_and_update_predecessor_successor(
                    current_node,
                )

                current_node.clear_next_and_previous_node_in_series()

    def _trace_forward(
        self,
        start_node: Node,
        trace_id: int,
        path: list[Node | Edge],
        group_paths: list[list[Node | Edge]],
    ) -> None:
        """Helper function to trace through graph and find all paths.

        .. seealso::
            :func:`SpliceGraph.trace`
        """

        if not start_node:
            # successor be [] or None
            if not isinstance(path[-1], Node):
                msg = f"Last element in path is not a node {path[-1]}"
                raise TypeError(msg)
            group_paths.append(path)

        elif successors := start_node.successors:
            for successor in successors:
                successor.set_trace_id(trace_id)
                for edge in self.get_possible_edges(path, start_node, successor):
                    self._trace_forward(
                        successor,
                        trace_id + 1,
                        [*path, edge, successor],
                        group_paths,
                    )
        else:
            # successor be [] or None
            self._trace_forward(
                successors,  # type: ignore
                trace_id + 1,
                [*path],
                group_paths,
            )

    def trace(self) -> Any:
        """Trace forward through graph and find all paths."""
        result_series_list = []

        if not self.get_start_nodes() and len(self.nodes.values()) > 0:
            self.logger.warning(f"A circle may exist in graph {self.nodes.values()}")

        for start_node in self.get_start_nodes():
            start_node.set_trace_id(1)
            group_paths = []
            self._trace_forward(start_node, 2, [start_node], group_paths)
            result_series_list.extend(group_paths)

        if not result_series_list:
            self.logger.info(
                f"No path is found in graph {self.nodes.values()}",
            )

        return result_series_list

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
        self,
        nodes_keys: set[str],
        result_paths: list[list[Node]],
    ) -> None:
        """Check if there is a circle in graph."""
        if (
            nodes_keys
            and (start_node := self.get_node_with_unique_key(nodes_keys.pop()))
            is not None
        ):
            current_nodes_keys: set[str] = set()
            self._trace_forward_record_node_unique_keys(
                start_node,
                [],
                result_paths,
                current_nodes_keys,
            )
            self.check_circle_in_graph_helper(
                nodes_keys - current_nodes_keys,
                result_paths,
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
        elif successors := start_node.successors:
            for successor in successors:
                if (key := successor.unique_key) is not None:
                    nodes_keys.add(key)
                self._trace_forward_record_node_unique_keys(
                    successor,
                    [*path, start_node],
                    group_paths,
                    nodes_keys,
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
        msg = f"{updated_node} or {current_node} has no exons"
        raise ValueError(msg)

    _update_exon_coord_sr_svtype_breakpoints_name_mode_in_same_exons(
        updated_node,
        current_node,
    )


def _update_exon_coord_sr_svtype_breakpoints_name_mode_in_same_exons(
    updated_node: Node,
    current_node: Node,
) -> None:
    """Update exon coordinates of the updated node based on the current node.

    Two nodes with the same number of exons

    :param updated_node:  node has been inserted into graph
    :param current_node: node has not been inserted into graph
    :return: None
    """
    if updated_node.exons is None or current_node.exons is None:
        msg = f"{current_node} or {updated_node} exons is None"
        raise ValueError(msg)

    # update exon coordinates
    updated_node.ref_start = min(
        updated_node.exons.first.start,
        current_node.exons.first.start,
    )

    # WARN: ref_end may be not consistent with prev_breakpoint of next edge <Yangyang Li>
    updated_node.ref_end = max(
        updated_node.exons.last.end,
        current_node.exons.last.end,
    )

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

    updated_node.read_ids.append(current_node.query_name)

    if current_node.self_identity is None:
        msg = f"{current_node}'s self_identity is None"
        raise ValueError(msg)

    updated_node.identities[current_node.query_name] = current_node.self_identity
