"""NLGraph."""

from __future__ import annotations

import copy
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
)
from .graphvis import default_visitors
from .merge_condition import MergeCondition

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
        junction_support_reads,
        input_bam_path: Path,
        output_dir: Path,
        *,
        rescue_sr: bool,
        refine_threshold: int,
        ignore_circle: bool = False,
        if_refine: bool = False,
        cluster_ind: int | str | None = None,
        is_plot: bool = False,
    ) -> None:
        """Initialize SpliceGraph."""
        self.logger = logger
        self.merge_threshold = merge_threshold  # 10 is tolerance compared cluster phase
        self.refine_threshold = refine_threshold

        self.support_reads = support_reads
        self.junction_support_reads = junction_support_reads

        self.use_junction_support = self.junction_support_reads > 0
        self.edge_support_threshold = self.junction_support_reads if self.use_junction_support else self.support_reads

        self.dict_factory = NLGraph.dict_factory  # type: ignore
        self.list_factory = NLGraph.list_factory  # type: ignore
        self.rescuer = rescuer
        self.rescue_sr = rescue_sr
        self.output_dir = output_dir
        self.has_circle = False
        self.ignore_circle = ignore_circle
        self.if_refine = if_refine
        self.possible_paths = None

        self.input_bam_path = input_bam_path
        self.cluster_ind: int | str | None = cluster_ind
        self.is_plot = is_plot

    def generate_paths(self) -> Iterable[NLPath]:
        """Generate paths from the nlgraph."""
        node_list: list[Node | Edge] = []
        all_paths = []
        # trace path
        for _idx, node_list in enumerate(self.trace(), 1):
            current_path = NLPath.create_path_from_node_edge_list(
                node_list,
            )
            all_paths.append(current_path)

        self.possible_paths = gather_possible_paths(all_paths)

        if self.is_plot and not self.has_circle and node_list:
            plot_result = self.output_dir / Path(f"graph_{self.input_bam_path.stem}")
            plot_result.mkdir(exist_ok=True)
            cluster_name = f"{self.input_bam_path.stem}_{self.cluster_ind}" if self.input_bam_path is not None else f"{cluster_ind}"
            default_visitors(self, (plot_result / cluster_name).as_posix(), support_reads=1, possible_paths=self.possible_paths).visualize()

        update_node_ptf_in_path(all_paths, len(all_paths))
        return all_paths

    def __call__(
        self,
        nlpath_list: Iterable[NLPath],
    ):
        """Find a specific path based on splice graph."""
        if isinstance(nlpath_list, types.GeneratorType):
            nlpath_list = list(nlpath_list)

        self.nlpaths = copy.deepcopy(nlpath_list)

        del nlpath_list  # remove reference to series_list
        self.nodes: dict[str, list[Node]] = self.dict_factory()
        self.edges: dict[str, list[Edge]] = defaultdict(list)

        # construct splice graph
        self.construct()
        # sr rescuer
        self.logger.trace(f"TSGraph Node: {len(self)}")
        # caluclate the depth on breakpoints only
        self.rescuer.init(self)

        if self.rescue_sr:
            self.rescuer(self)

        if self.if_refine:
            self.refine(self.refine_threshold)

        if not is_weakly_connected(self):
            logger.warning(f"Graph {self.nodes=} is not weakly connected")

        self.polish_edges()

    @classmethod
    def create_graph(
        cls,
        rescuer: Any,
        input_bam: str,
        logger: LoggerType,
        prune_threshold: int,
        refine_threshold: int,
        support_reads: int,
        junction_support_reads: int,
        output_dir: Path,
        *,
        ignore_circle: bool,
        rescue_sr: bool,
        refine: bool = False,
        cluster_ind: int | str | None = None,
        is_plot: bool = False,
    ) -> NLGraph:
        """Create  nlgraph."""
        return cls(
            logger,
            rescuer,
            prune_threshold,
            support_reads,
            junction_support_reads,
            Path(input_bam),
            output_dir,
            rescue_sr=rescue_sr,
            ignore_circle=ignore_circle,
            if_refine=refine,
            refine_threshold=refine_threshold,
            cluster_ind=cluster_ind,
            is_plot=is_plot,
        )

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

            edge_sr = edge.junction_sr if self.use_junction_support else edge.sr

            if edge_sr >= self.edge_support_threshold:
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

        # Get the list of nodes with the same similar_key
        nodes_list = self.get_nodes_with_similar_key(node.similar_key)

        # Check if the node is actually in the list
        if node in nodes_list:
            nodes_list.remove(node)
        else:
            # If the node is not in the list, try to find it by unique_key
            logger.warning(f"Node {node} not found in similar_key list. Attempting to find by unique_key.")
            for n in list(nodes_list):
                if n.unique_key == node.unique_key:
                    nodes_list.remove(n)
                    logger.info(f"Removed node with matching unique_key instead: {n}")
                    return

            # If we get here, the node wasn't found by either method
            logger.warning(f"Failed to remove node {node} (similar_key={node.similar_key}): not found in graph")

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

        merge_condition = MergeCondition(threshold, None)
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
            :func:`NLGraph.trace`
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

    def refine(self, threshold):
        """Refine the graph by merging nodes and edges.

        1. Group nodes by merge signature (chrom, strand, introns, start/end within threshold)
        2. For each group, compare pairs and merge as needed, skipping already merged/removed nodes
        3. For each merge:
            - Merge common predecessors/successors and their edges
            - For unique predecessors/successors, rewire edges
            - Remove the merged node from the graph
        4. Avoid duplicate predecessors/successors
        5. Add robust error handling and logging
        """
        logger.trace(f"Refine graph: number of nodes before refine: {len(self)} with {threshold=}")

        # Helper: create two signatures for grouping nodes that could be merged.
        # Using two overlapping bin offsets avoids boundary issues where nodes
        # within threshold distance land in different bins due to integer division.
        def merge_signatures(node, threshold):
            bin_size = threshold + 1
            offset = bin_size // 2
            sig1 = (
                node.chrom,
                node.strand,
                tuple(node.introns),
                int(node.ref_start // bin_size),
                int(node.ref_end // bin_size),
            )
            sig2 = (
                node.chrom,
                node.strand,
                tuple(node.introns),
                int((node.ref_start + offset) // bin_size),
                int((node.ref_end + offset) // bin_size),
            )
            return sig1, sig2

        # Group nodes by signature (each node may appear in up to two groups)
        from collections import defaultdict

        signature_to_nodes = defaultdict(list)
        for node in self:
            sig1, sig2 = merge_signatures(node, threshold=threshold)
            signature_to_nodes[sig1].append(node)
            if sig2 != sig1:
                signature_to_nodes[sig2].append(node)

        removed_nodes = set()  # Track nodes that have been merged/removed
        seen_pairs = set()  # Track compared pairs to avoid duplicates across overlapping bins

        # For each group, compare pairs
        for group in signature_to_nodes.values():
            n = len(group)
            for i in range(n):
                node1 = group[i]
                if node1 in removed_nodes:
                    continue
                for j in range(i + 1, n):
                    node2 = group[j]
                    if node2 in removed_nodes:
                        continue
                    pair_key = (id(node1), id(node2)) if id(node1) < id(node2) else (id(node2), id(node1))
                    if pair_key in seen_pairs:
                        continue
                    seen_pairs.add(pair_key)
                    if compare_node_when_refine(node1, node2, threshold):
                        # Select the node with more read IDs as the primary node
                        if len(node1.read_ids) >= len(node2.read_ids):
                            primary, secondary = node1, node2
                        else:
                            primary, secondary = node2, node1
                        logger.trace(f"Refine graph: merging {primary} and {secondary}")
                        try:
                            self._merge_node_predecessors(primary, secondary, removed_nodes)
                            self._merge_node_successors(primary, secondary, removed_nodes)
                            primary.merge(secondary)
                            self.remove_node(secondary)
                            removed_nodes.add(secondary)
                        except (IndexError, KeyError, ValueError) as e:
                            logger.error(f"Error merging nodes {primary} and {secondary}: {e}")

        self.logger.info(f"Refinement process completed: number of nodes after refine: {len(self)}")

    def _merge_node_predecessors(self, primary_node, secondary_node, removed_nodes=None):
        """Handle predecessor merging during node refinement.

        Args:
            primary_node: The node that will remain in the graph
            secondary_node: The node that will be removed
            removed_nodes: Set of nodes already removed (optional)
        """
        if removed_nodes is None:
            removed_nodes = set()
        # Process common predecessors
        common_predecessors = set(primary_node.predecessors) & set(secondary_node.predecessors)
        for predecessor in list(common_predecessors):
            if predecessor in removed_nodes:
                continue
            try:
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
                if secondary_node in predecessor.successors:
                    predecessor.successors.remove(secondary_node)
                self.remove_edge(secondary_edge)
            except (IndexError, KeyError) as e:
                logger.error(f"Error merging edges from {predecessor}: {e}")
        # Process unique predecessors
        unique_predecessors = set(secondary_node.predecessors) - set(primary_node.predecessors)
        for predecessor in list(unique_predecessors):
            if predecessor in removed_nodes:
                continue
            try:
                edge_data = self.find_edges(predecessor, secondary_node)[0].edge_data
                # Add to primary node's predecessors if not already present
                if predecessor not in primary_node.predecessors:
                    primary_node.predecessors.append(predecessor)
                # Create edge between predecessor and primary node
                self.add_edge(predecessor, primary_node, edge_data)
                # Update predecessor's successors list
                if secondary_node in predecessor.successors:
                    predecessor.successors.remove(secondary_node)
                if primary_node not in predecessor.successors:
                    predecessor.successors.append(primary_node)
            except (IndexError, KeyError) as e:
                logger.error(f"Error creating edge from {predecessor} to {primary_node}: {e}")

    def _merge_node_successors(self, primary_node, secondary_node, removed_nodes=None):
        """Handle successor merging during node refinement.

        Args:
            primary_node: The node that will remain in the graph
            secondary_node: The node that will be removed
            removed_nodes: Set of nodes already removed (optional)
        """
        if removed_nodes is None:
            removed_nodes = set()
        # Process common successors
        common_successors = set(primary_node.successors) & set(secondary_node.successors)
        for successor in list(common_successors):
            if successor in removed_nodes:
                continue
            try:
                primary_edges = self.find_edges(primary_node, successor)
                secondary_edges = self.find_edges(secondary_node, successor)
                if len(primary_edges) > 1:
                    logger.error(f"Multiple edges {primary_edges} found between {primary_node} and {successor}")
                if len(secondary_edges) > 1:
                    logger.error(f"Multiple edges {secondary_edges} found between {secondary_node} and {successor}")
                primary_edge = primary_edges[0]
                secondary_edge = secondary_edges[0]
                primary_edge.merge(secondary_edge, primary_node.strand, successor.strand)
                # Remove the secondary connection - use safe removal
                if secondary_node in successor.predecessors:
                    successor.predecessors.remove(secondary_node)
                self.remove_edge(secondary_edge)
            except (IndexError, KeyError) as e:
                logger.error(f"Error merging edges to {successor}: {e}")
        # Process unique successors
        unique_successors = set(secondary_node.successors) - set(primary_node.successors)
        for successor in list(unique_successors):
            if successor in removed_nodes:
                continue
            try:
                edge_data = self.find_edges(secondary_node, successor)[0].edge_data
                # Add to primary node's successors if not already present
                if successor not in primary_node.successors:
                    primary_node.successors.append(successor)
                # Create edge between primary node and successor
                self.add_edge(primary_node, successor, edge_data)
                # Update successor's predecessors list - use safe removal
                if secondary_node in successor.predecessors:
                    successor.predecessors.remove(secondary_node)
                if primary_node not in successor.predecessors:
                    successor.predecessors.append(primary_node)
            except (IndexError, KeyError) as e:
                logger.error(f"Error creating edge from {primary_node} to {successor}: {e}")

    def polish_edges(self) -> None:
        """Polish edges in the graph.

        For each node and its successors, update the breakpoints of the connecting edge(s).
        If multiple edges exist between a node pair, log a warning and only update the first edge.
        """
        logger.trace("Polish edges in the graph")
        for node in self:
            for successor in node.successors:
                try:
                    edges = self.find_edges(node, successor)
                except (KeyError, IndexError) as e:
                    logger.error(f"Cannot find edge(s) between {node} and {successor}: {e}")
                    continue
                if not edges:
                    logger.warning(f"No edge found between {node} and {successor}")
                    continue
                if len(edges) > 1:
                    logger.error(f"Multiple edges {edges} found between {node} and {successor}")
                # Only update the first edge if multiple exist
                edge = edges[0]
                try:
                    edge.break_point1.pos = node.ref_end if node.strand.is_forward() else node.ref_start
                    edge.break_point2.pos = successor.ref_start if successor.strand.is_forward() else successor.ref_end
                except Exception as e:  # noqa: BLE001
                    logger.error(f"Error updating breakpoints for edge {edge} between {node} and {successor}: {e}")


def compare_node_when_refine(node1: Node, node2: Node, threshold) -> bool:
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


def gather_possible_paths(nlpaths: list[NLPath]) -> dict[str, list[str]]:
    """Gather possible paths from NLPaths."""
    possible_paths = defaultdict(list)

    for path in nlpaths:
        for idx, node in enumerate(path.nodes[:-1]):
            possible_paths[path.id].append(node.id)

            next_node = path.nodes[idx + 1]
            edge = path.get_edge(node, next_node)

            possible_paths[path.id].append(edge.id)

        possible_paths[path.id].append(next_node.id)

    return possible_paths


def update_junction_support(graphs: Iterable[NLGraph]) -> None:
    """Update junction support for all graphs.

    Args:
        graphs: Iterable of NLGraph objects
        threshold: Position threshold for binning breakpoints (default: 10).
                   Positions within this distance will be grouped together.
    """
    edge_dict: dict[str, int] = defaultdict(int)

    # summary junction support
    for graph in graphs:
        for edge_list in graph.edges.values():
            for edge in edge_list:
                edge_key = f"{edge.break_point1.chrom}-{edge.break_point1.pos}-{edge.break_point2.chrom}-{edge.break_point2.pos}"
                edge_dict[edge_key] += edge.sr

    # update junction support
    for graph in graphs:
        for edge_list in graph.edges.values():
            for edge in edge_list:
                edge.edge_data.junction_sr = edge_dict[
                    f"{edge.break_point1.chrom}-{edge.break_point1.pos}-{edge.break_point2.chrom}-{edge.break_point2.pos}"
                ]
