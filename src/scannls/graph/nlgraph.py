"""NLGraph operations."""

from __future__ import annotations

import copy
from itertools import combinations
import types
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any

from loguru import logger

from .annotate import is_weakly_connected
from .basic_graph import (
    Edge,
    NLPath,
    Node,
    NodeIdentity,
    to_hash_identifier,
)
from .graphvis import default_visitors
from .merge_condition import MergeCondition
from .sr_rescuer import SRRescuer

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from scannls.mtype import LoggerType


class NLGraph:
    """SpliceGraph class is used to trace the path of the splice graph."""

    dict_factory = dict
    list_factory = list

    def __init__(
        self,
        logger: LoggerType,
        rescuer: Any,
        merge_threshold,
        support_reads,
        input_bam_path: Path,
        output_dir: Path,
        rescue_sr: bool,
        ignore_circle: bool = False,
        if_refine: bool = False,
    ) -> None:
        """Initialize SpliceGraph."""
        self.logger = logger
        self.merge_threshold = merge_threshold  # 10 is tolerance compared cluster phase
        self.support_reads = support_reads
        self.dict_factory = NLGraph.dict_factory  # type: ignore
        self.list_factory = NLGraph.list_factory  # type: ignore
        self.rescuer = rescuer
        self.rescue_sr = rescue_sr
        self.input_bam_path = input_bam_path
        self.output_dir = output_dir
        self.has_circle = False
        self.ignore_circle = ignore_circle
        self.if_refine = if_refine

    def __call__(
        self,
        nlpath_list: Iterable[NLPath],
        cluster_ind: int | str,
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
        self.logger.trace(f"NLGraph Node: {len(self)}")

        # caluclate the depth on breakpoints only
        self.rescuer.init(self)

        if self.rescue_sr:
            self.rescuer(self)

        if self.if_refine:
            self.refine()

        if not is_weakly_connected(self):
            logger.warning(f"Graph {self.nodes=} is not weakly connected")

        self.polish_edges()

        node_list: list[Node | Edge] = []

        all_paths = []
        # trace path
        for _idx, node_list in enumerate(self.trace(), 1):
            current_path = NLPath.create_path_from_node_edge_list(
                node_list,
            )
            # current_path.polish_edges()
            all_paths.append(current_path)

        if is_plot and not self.has_circle and node_list:
            plot_result = self.output_dir / Path(f"graph_{self.input_bam_path.stem}")
            plot_result.mkdir(exist_ok=True)
            cluster_name = f"{self.input_bam_path.stem}_{cluster_ind}" if self.input_bam_path is not None else f"{cluster_ind}"
            default_visitors(self, (plot_result / cluster_name).as_posix(), support_reads=1).visualize()

        update_node_ptf_in_path(all_paths, len(all_paths))
        return all_paths

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
        output_dir: Path,
        *,
        ignore_circle: bool,
        rescue_sr: bool,
        refine: bool = False,
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

        return cls(logger, rescuer, prune_threshold, support_reads, Path(input_bam), output_dir, rescue_sr, ignore_circle, refine)

    def add_edge(self, node1: Node, node2: Node, edge_data):
        """Add edge from node1 -> node2."""
        edge = Edge.from_nodes(node1, node2, edge_data)

        if self.edges.get(edge.key) is None:
            self.logger.info(f"add edge with {edge_data=}")
            self.edges[edge.key].append(edge)
            return

        is_merged = False

        for current_edge in self.edges[edge.key]:
            if current_edge.read_ids == edge.read_ids:
                # same edge
                return

            if current_edge.merged(
                edge,
                compared_break_point=False,
                break_point_threshold=self.merge_threshold,
            ):
                logger.info(f"merging {current_edge} and {edge}")
                current_edge.merge(edge, node1.strand, node2.strand)
                is_merged = True
                break

        if not is_merged:
            self.logger.info(f"add edge with {edge_data=}")

            if len(self.edges[edge.key]) > 1:
                self.logger.warning(f"Multiple edges {self.edges[edge.key]} found between {node1} and {node2}")

            self.edges[edge.key].append(edge)

    def find_edges(self, node1: Node, node2: Node) -> list[Edge]:
        edge_key = Edge.create_key_from_node(node1, node2)
        if self.edges.get(edge_key) is None:
            msg = f"Edge {edge_key} not found in graph."
            raise KeyError(msg)
        return self.edges[edge_key]

    @staticmethod
    def get_node_identity_base_edge(
        node: Node,
        edge: Edge,
    ) -> dict[NodeIdentity, list[str]]:
        result = defaultdict(list)
        for read_id in edge.read_ids:
            node_identity = node.identity(read_id)
            if node_identity is None:
                msg = f"node={node!r} cannot find identiy for read_id={read_id!r}"
                raise ValueError(msg)
            result[node_identity].append(read_id)
        return result

    @staticmethod
    def determine_edge(previous_edge_node_identity, next_edge_node_identity) -> bool:
        previous_edge_read_id_mid = set(previous_edge_node_identity[NodeIdentity.MID])
        next_edge_read_id_mid = set(next_edge_node_identity[NodeIdentity.MID])

        if previous_edge_read_id_mid and next_edge_read_id_mid:
            # check two mid from different reads
            return len(previous_edge_read_id_mid & next_edge_read_id_mid) > 0

        return True

    def get_possible_edges(
        self,
        current_path: list[Node | Edge],
        current_node: Node,
        successor: Node,
        support_reads: int,
        *,
        filter_edges: bool = True,
    ) -> Iterable[Edge]:
        if len(current_path) == 1:
            previous_edge_node_identity = None
        else:
            previous_edge = current_path[-2]
            if not isinstance(previous_edge, Edge):
                msg = f"previous_edge is not instance of Edge: {previous_edge}"
                raise TypeError(msg)

            previous_edge_node_identity = self.get_node_identity_base_edge(
                current_node,
                previous_edge,
            )

        edges = []
        for edge in self.find_edges(current_node, successor):
            if not filter_edges:
                edges.append(edge)
                continue

            if edge.sr >= support_reads:
                edge_node_identity = self.get_node_identity_base_edge(
                    current_node,
                    edge,
                )

                if previous_edge_node_identity is not None:
                    if self.determine_edge(
                        previous_edge_node_identity,
                        edge_node_identity,
                    ):
                        edges.append(edge)
                    else:
                        logger.warning(f"Read guided: {current_node} and {successor}")
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
        return any(other_node.unique_key == node.unique_key for other_node in self.get_nodes_with_similar_key(node.similar_key))

    def __iter__(self) -> Iterator[Node]:
        """Iterate over all nodes in graph."""
        for nodes in self.nodes.values():
            yield from nodes

    def __len__(self) -> int:
        """Get number of nodes in graph."""
        return sum(1 for _ in self)

    def get_start_nodes(self) -> Iterable[Node]:
        """Get start nodes based if node has predecessors."""
        return (node for nodes in self.nodes.values() for node in nodes if node.is_start_node())

    def get_end_nodes(self) -> Iterable[Node]:
        """Get end nodes based if node has successors."""
        return (node for nodes in self.nodes.values() for node in nodes if node.is_end_node())

    def remove_edge(self, edge: Edge) -> None:
        """Remove edge from graph.

        :param edge: edge to be removed
        """
        if edge.key is None:
            msg = f"edge.key is None, {edge}"
            raise ValueError(msg)
        self.edges[edge.key].remove(edge)

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
        logger.trace(f"compare merge: {node1=}, {node2=}")

        if node1.strand != node2.strand or node1.chrom != node2.chrom:
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
                self.merge_threshold,
            ):
                current_node.is_merged = True

                self.logger.info(
                    f"merging node {similar_node_in_graph} and {current_node})",
                )

                merge_nodes(similar_node_in_graph, current_node)

                # nodes in merged_parent_nodes are all in the graph
                current_node.merged_parent_nodes.append(similar_node_in_graph)

                # only consider nodes that have been processed: previous node in current series
                # keeps in mind the next node in current series is not processed yet!!!!
                # a -> b and b <- a

                edge_data = current_node.previous_edge_in_nlapth.edge_data if current_node.previous_edge_in_nlapth is not None else None

                similar_node_in_graph.add_predecessor(
                    current_node.previous_node_in_nlpath,
                    self,
                    edge_data,
                )
                break

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

            edge_data = current_node.previous_edge_in_nlapth.edge_data if current_node.previous_edge_in_nlapth is not None else None

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
        path: list[Node | Edge],
        group_paths: list[list[Node | Edge]],
    ):
        """Helper function to trace through graph and find all paths.

        .. seealso::
            :func:`SpliceGraph.trace`
        """
        if not start_node or self.has_circle:
            # successor be [] or None
            update_node_ptc_in_path(path)
            group_paths.append(path)

        elif successors := start_node.successors:
            for successor in successors:
                if successor in path:
                    self.logger.warning(
                        f"A circle exist in graph with nodes {path}",
                    )
                    self.has_circle = True
                    continue

                for edge_ind, edge in enumerate(
                    self.get_possible_edges(
                        path,
                        start_node,
                        successor,
                        self.support_reads,
                    )
                ):
                    if edge_ind > 0:
                        self.logger.warning(f"Multiple edges {edge} found between {start_node} and {successor}")

                    if start_node.is_start_node():
                        new_break_point = start_node.update_breakpoint()
                        path[-1].update_breakpoint()
                        edge.break_point1.pos = new_break_point

                    self._trace_forward(
                        successor,
                        [*path, edge, successor],
                        group_paths,
                    )

        else:
            if start_node.is_end_node():
                new_break_point = start_node.update_breakpoint()
                path[-1].update_breakpoint()
                path[-2].break_point2.pos = new_break_point

            # successor be [] or None
            self._trace_forward(
                successors,  # type: ignore
                [*path],
                group_paths,
            )

    def trace(self) -> Any:
        """Trace forward through graph and find all paths."""
        self.has_circle = False

        result_paths_list = []

        if not self.get_start_nodes() and len(self.nodes.values()) > 0:
            self.logger.warning(f"A circle may exist in graph {self.nodes.values()}")

        for start_node in self.get_start_nodes():
            group_paths = []
            self._trace_forward(
                start_node,
                [start_node],
                group_paths,
            )

            if not self.ignore_circle and self.has_circle:
                result_paths_list.clear()
                break

            result_paths_list.extend(group_paths)

        if not result_paths_list:
            self.logger.info(
                f"No path is found in graph {self.nodes.values()}",
            )

        return result_paths_list

    @property
    def id(self) -> str:
        """Get id of graph."""
        identifier = to_hash_identifier("-".join(node.id for node in self))
        return f"TSN{identifier}"

    def refine(self):
        """Refine the graph by merging nodes and edges.

        1. check every pair of nodes in the graph
        2. if the nodes are same information including chr, strand, start, end, exons.
        3. use one node with higher sr to merge another node, and update the edges of the nodes (by merging edges)
        4. remove the merged nodes from the graph, and remove the edges of the merged nodes
        5. update the graph structure to reflect the changes
        6. Only iterate the nodes once, so the time complexity is O(n^2)
        """
        logger.trace(f"Refine graph: number of nodes before refine: {len(self)}")

        # Find similar nodes to merge
        nodes_to_merge = []

        # Create a dictionary for faster node lookup
        node_dict = {node.unique_key: node for node in self}

        # Find pairs of nodes to merge - O(n²) operation
        for node1_key, node2_key in combinations(node_dict.keys(), 2):
            if node1_key == node2_key:
                continue

            node1 = node_dict[node1_key]
            node2 = node_dict[node2_key]

            if compare_node_when_refine(node1, node2):
                # Select the node with more read IDs as the primary node
                if len(node1.read_ids) >= len(node2.read_ids):
                    nodes_to_merge.append((node1, node2))
                else:
                    nodes_to_merge.append((node2, node1))

        # Merge the nodes
        for node1, node2 in nodes_to_merge:
            logger.trace(f"Refine graph: merging {node1} and {node2}")

            # Process predecessors
            self._merge_node_predecessors(node1, node2)

            # Process successors
            self._merge_node_successors(node1, node2)

            # Merge node data
            node1.merge(node2)

            # Remove the merged node from the graph
            self.remove_node(node2)

        self.logger.info(f"Refinement process completed: number of nodes after refine: {len(self)}")

    def _merge_node_predecessors(self, primary_node, secondary_node):
        """Handle predecessor merging during node refinement.

        Args:
            primary_node: The node that will remain in the graph
            secondary_node: The node that will be removed
        """
        # Process common predecessors
        common_predecessors = set(primary_node.predecessors) & set(secondary_node.predecessors)
        for predecessor in common_predecessors:
            try:
                # Merge edges for common predecessors
                primary_edges = self.find_edges(predecessor, primary_node)
                secondary_edges = self.find_edges(predecessor, secondary_node)

                if len(primary_edges) > 1:
                    logger.error(f"Multiple edges {primary_edges} found between {predecessor} and {primary_node}")
                if len(secondary_edges) > 1:
                    logger.error(f"Multiple edges {secondary_edges} found between {predecessor} and {secondary_node}")

                primary_edge = primary_edges[0]
                secondary_edge = secondary_edges[0]
                primary_edge.merge(secondary_edge, predecessor.strand, primary_node.strand)

                # Remove the secondary connection
                predecessor.successors.remove(secondary_node)
                self.remove_edge(secondary_edge)
            except (IndexError, KeyError) as e:
                logger.error(f"Error merging edges from {predecessor}: {e}")

        # Process unique predecessors
        unique_predecessors = set(secondary_node.predecessors) - set(primary_node.predecessors)
        for predecessor in unique_predecessors:
            try:
                # Get the edge data and create a new edge
                edge_data = self.find_edges(predecessor, secondary_node)[0].edge_data

                # Add to primary node's predecessors
                primary_node.predecessors.append(predecessor)

                # Create edge between predecessor and primary node
                self.add_edge(predecessor, primary_node, edge_data)

                # Update predecessor's successors list
                predecessor.successors.remove(secondary_node)
                predecessor.successors.append(primary_node)
            except (IndexError, KeyError) as e:
                logger.error(f"Error creating edge from {predecessor} to {primary_node}: {e}")

    def _merge_node_successors(self, primary_node, secondary_node):
        """Handle successor merging during node refinement.

        Args:
            primary_node: The node that will remain in the graph
            secondary_node: The node that will be removed
        """
        # Process common successors
        common_successors = set(primary_node.successors) & set(secondary_node.successors)
        for successor in common_successors:
            try:
                # Merge edges for common successors
                primary_edges = self.find_edges(primary_node, successor)
                secondary_edges = self.find_edges(secondary_node, successor)

                if len(primary_edges) > 1:
                    logger.error(f"Multiple edges {primary_edges} found between {primary_node} and {successor}")
                if len(secondary_edges) > 1:
                    logger.error(f"Multiple edges {secondary_edges} found between {secondary_node} and {successor}")

                primary_edge = primary_edges[0]
                secondary_edge = secondary_edges[0]
                primary_edge.merge(secondary_edge, primary_node.strand, successor.strand)

                # Remove the secondary connection
                successor.predecessors.remove(secondary_node)
                self.remove_edge(secondary_edge)
            except (IndexError, KeyError) as e:
                logger.error(f"Error merging edges to {successor}: {e}")

        # Process unique successors
        unique_successors = set(secondary_node.successors) - set(primary_node.successors)
        for successor in unique_successors:
            try:
                # Get the edge data and create a new edge
                edge_data = self.find_edges(secondary_node, successor)[0].edge_data

                # Add to primary node's successors
                primary_node.successors.append(successor)

                # Create edge between primary node and successor
                self.add_edge(primary_node, successor, edge_data)

                # Update successor's predecessors list
                successor.predecessors.remove(secondary_node)
                successor.predecessors.append(primary_node)
            except (IndexError, KeyError) as e:
                logger.error(f"Error creating edge from {primary_node} to {successor}: {e}")

    def polish_edges(self) -> None:
        """Polish edges in the graph."""

        # iterate all nodes in the graph
        # get every predecessor and successor of the node
        # update the break point of the edge
        for node in self:
            for successor in node.successors:
                # update edge from node to successor
                edges = self.find_edges(node, successor)
                if len(edges) > 1:
                    logger.error(f"Multiple edges {edges} found between {node} and {successor}")

                for edge in edges:
                    edge.break_point1.pos = node.ref_end if node.strand.is_forward() else node.ref_start
                    edge.break_point2.pos = successor.ref_start if successor.strand.is_forward() else successor.ref_end


def compare_node_when_refine(node1: Node, node2: Node, threshold=0) -> bool:
    """Compare two nodes.

    :param node1: node1
    :param node2: node2
    :return: True if the nodes are equal, otherwise False
    """
    return (
        node1.chrom == node2.chrom
        and node1.strand == node2.strand
        and node1.introns == node2.introns
        and abs(node1.ref_start - node2.ref_start) <= threshold
        and abs(node1.ref_end - node2.ref_end) <= threshold
    )


def merge_nodes(
    updated_node: Node,
    current_node: Node,
) -> None:
    """Update exon coordinates of the updated node based on the current node.

    Two nodes with the same number of exons

    :param updated_node:  node has been inserted into graph
    :param current_node: node has not been inserted into graph
    :return: None
    """
    updated_node.merge(current_node)


def update_node_ptc_in_path(path: list[Node | Edge]) -> None:
    """Update PTC information in path."""
    for node in path:
        if isinstance(node, Node):
            node.ptc += 1


def update_node_ptf_in_path(paths: list[NLPath], total_path: int) -> None:
    """Update PTF information in path."""
    for path in paths:
        for node in path:
            if isinstance(node, Node):
                node.ptf = node.ptc / total_path
