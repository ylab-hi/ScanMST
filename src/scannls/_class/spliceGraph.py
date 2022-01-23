# !/usr/bin/env python
"""Splice Graph.

@Filename:    spliceGraph.py
@license:     MIT Licence
@Time:        12/15/21 10:42 AM
"""
import copy
import types
from typing import Any
from typing import Dict
from typing import Iterable
from typing import List
from typing import Set
from typing import Union

from ..type import LoggerType
from .basicClass import Insertion
from .basicClass import MicroHomology
from .basicClass import Node
from .basicClass import NovelInsertion
from .basicClass import Series
from .exception import ExonsNotFoundError
from .srRescuer import SRRescuer

NodeType = Union[Node, Insertion]


class SpliceGraph:
    """SpliceGraph class is used to trace the path of splice graph."""

    dict_factory = dict
    list_factory = list

    def __init__(self, logger: LoggerType):
        """Initialize SpliceGraph."""
        self.logger = logger
        self.dict_factory = SpliceGraph.dict_factory  # type: ignore
        self.list_factory = SpliceGraph.list_factory  # type: ignore

    def __call__(
        self, series_list: Iterable[Series], rescuer: SRRescuer
    ) -> Iterable[Series]:
        """Find specific path based on splice graph.

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
        self.nodes: Dict[str, List[NodeType]] = self.dict_factory()
        # construct splice graph
        self.construct()
        # sr rescuer
        rescuer(self)
        # prun the graph

        # trace path
        for node_list in self.trace():
            yield Series.create_series_from_node_list(node_list, self.logger)

    def get_start_nodes(self):
        """Get start nodes based if node has predecessors."""
        return [
            node
            for nodes in self.nodes.values()
            for node in nodes
            if node.is_start_node()
        ]

    def get_nodes_with_similar_key(self, similar_key: str) -> List[NodeType]:
        """Get nodes in graph with similar key."""
        return self.nodes.get(similar_key, [])

    def add_node_with_similar_key(self, node: NodeType) -> None:
        """Add node to the splice graph.

        :param node: node to be added
        """
        if similar_nodes := self.get_nodes_with_similar_key(node.similar_key):
            similar_nodes.append(node)
        else:
            self.nodes[node.similar_key] = [node]

    def __contains__(self, node: NodeType) -> bool:
        """Check if node is in graph.

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

    def __iter__(self):
        """Iterate over all nodes in graph."""
        for nodes in self.nodes.values():
            yield from nodes

    @staticmethod
    def _check_insertion_conditions_for_compare(
        node1: NodeType, node2: NodeType
    ) -> bool:
        """Check if node1 and node2 can be merged based on insertion info."""
        flag = True
        insertion_info1 = node1.insertion_info
        insertion_info2 = node2.insertion_info
        if insertion_info1 is None and insertion_info2 is None:
            # When None = None
            return flag
        elif insertion_info1 is not None and insertion_info2 is not None:
            if insertion_info1[0] and insertion_info2[0]:
                # 1 hit insertion that is added in the series
                return flag

            if not insertion_info1[0] and not insertion_info2[0]:
                if (
                    isinstance(insertion_info1[1], NovelInsertion)
                    and isinstance(insertion_info2[1], NovelInsertion)
                    and (
                        insertion_info1[1].query_sequence
                        == insertion_info2[1].query_sequence
                    )
                ):
                    return flag
                elif isinstance(insertion_info1[1], MicroHomology) and isinstance(
                    insertion_info2[1], MicroHomology
                ):
                    return True

        return False

    @staticmethod
    def _compare_is_merged_helper_check_condition_for_head_tail_node_mode(
        node1: NodeType,
        node2: NodeType,
        threshold: float = 0.8,
    ) -> bool:
        """Check if node1 and node2 can be merged based on overlap info.

        node1 is tail node, node2 is head node Using mean overlap ratio to
        check if they can be merged.

        :param node1:  node1
        :param node2:  node2
        :param threshold:  threshold for checking if two nodes are merged
        :return:  True if two nodes are merged, otherwise False

        .. note::

            -> [node1]
                [node2] ->
        """
        if node1.exons is None or node2.exons is None:
            raise SystemExit from ExonsNotFoundError
        node1_first_exon_start = node1.exons[0][0]
        node1_last_exon_end = node1.exons[-1][1]
        node2_first_exon_start = node2.exons[0][0]
        node2_last_exon_end = node2.exons[-1][1]
        expression1 = node1_last_exon_end >= node2_first_exon_start
        expression2 = (
            expression1
            and node1_first_exon_start <= node2_first_exon_start
            and node1_last_exon_end <= node2_last_exon_end
        )
        if expression2:
            overlap_len = node1_last_exon_end - node2_first_exon_start
            node1_mean_overlap_ratio = overlap_len / (
                node1_last_exon_end - node1_first_exon_start
            )
            node2_mean_overlap_ratio = overlap_len / (
                node1_last_exon_end - node1_first_exon_start
            )
            return (
                0.5 * (node1_mean_overlap_ratio + node2_mean_overlap_ratio) >= threshold
            )

        return False

    @staticmethod
    def _compare_is_merged_helper(node1: NodeType, node2: NodeType) -> bool:
        """Check if node1 and node2 can be merged."""
        condition = (
            node1.sv_type == node2.sv_type
            and SpliceGraph._check_insertion_conditions_for_compare(node1, node2)
        )
        if not condition:
            if (
                node1.next_breakpoint is None and node2.prev_breakpoint is None
            ):  # node1 is end node, node2 is start node
                return SpliceGraph._compare_is_merged_helper_check_condition_for_head_tail_node_mode(
                    node1, node2
                )
            elif (
                node1.next_breakpoint is None and node2.next_breakpoint is not None
            ):  # node1 is end node, node2 is middle node
                return (
                    node1.exons[0][0] == node2.exons[0][0]  # type: ignore
                    and node1.exons[-1][1] <= node2.exons[-1][1]  # type: ignore
                )

            return condition

        if (
            node1.prev_breakpoint is None and node2.prev_breakpoint is None
        ):  # both are start nodel check last exon end
            return node1.next_breakpoint == node2.next_breakpoint

        elif (
            node1.next_breakpoint is None and node2.next_breakpoint is None
        ):  # both are end nodes  # check first exon start
            return node1.prev_breakpoint == node2.prev_breakpoint

        elif (
            node1.prev_breakpoint is None and node2.prev_breakpoint is not None
        ):  # node1 is start node, node2 is middle node
            return (
                node1.exons[-1][1] == node2.exons[-1][1]  # type: ignore
                and node1.exons[0][0] >= node2.exons[0][0]  # type: ignore
            )
        elif (
            node1.prev_breakpoint is not None
            and node1.next_breakpoint is not None
            and node2.prev_breakpoint is not None
            and node2.next_breakpoint is not None
        ):  # both are middle nodes

            return (
                node1.exons[0][0] == node2.exons[0][0]  # type: ignore
                and node1.exons[-1][1] == node2.exons[-1][1]  # type: ignore
            )
        # swap node1 and node2 to check if they can be merged again
        return False

    @staticmethod
    def _compare_is_merged(node1: NodeType, node2: NodeType) -> bool:
        """Node1 is similar as node2 is precommit of the function.

         compare if node1 can merge node2
        :param node1: node1
        :param node2: node2
        :return:
        """
        flag1 = SpliceGraph._compare_is_merged_helper(node1, node2)
        if flag1:
            return flag1
        flag2 = SpliceGraph._compare_is_merged_helper(node2, node1)
        if flag2:
            return flag2
        return False

    @staticmethod
    def update_exon_coord_sr_svtype_breakpoints_name_mode(
        updated_node: NodeType, current_node: NodeType
    ) -> None:
        """Update exon coordinates of the updated node based on current node.

        :param updated_node:  node has been inserted into graph
        :param current_node: node has not been inserted into graph
        :return: None
        """
        # update exon coordinates
        updated_node.exons[0][0] = updated_node.ref_start = min(  # type: ignore
            updated_node.exons[0][0], current_node.exons[0][0]  # type: ignore
        )
        updated_node.exons[-1][1] = updated_node.ref_end = max(  # type: ignore
            updated_node.exons[-1][1], current_node.exons[-1][1]  # type: ignore
        )
        # update sr
        updated_node.update_sr()
        # update novel insertion ao
        if updated_node.insertion_info and isinstance(
            updated_node.insertion_info[1], NovelInsertion
        ):
            updated_node.insertion_info[1].increment_ao()
        # update sv_type
        updated_node.sv_type = (
            current_node.sv_type
            if current_node.sv_type is not None
            else updated_node.sv_type
        )
        # update breakpoints
        if updated_node.prev_breakpoint is None:
            updated_node.prev_breakpoint = current_node.prev_breakpoint
        if updated_node.next_breakpoint is None:
            updated_node.next_breakpoint = current_node.next_breakpoint
        # update query name
        updated_node.query_name += "," + current_node.query_name

        # update mode of the node
        if updated_node.modes is None:
            updated_node.modes = current_node.modes

    def _check_if_current_node_is_merged_in_similar_nodes_in_graph(
        self,
        current_node: NodeType,
        similar_key: str,
        merged_nodes_pool: Set[NodeType],
    ) -> None:
        """Check if current node is merged in similar nodes in graph."""
        # get similar nodes in the graph
        similar_nodes_in_graph = self.get_nodes_with_similar_key(similar_key)

        # iterate all similar nodes in the graph
        for similar_node_in_graph in similar_nodes_in_graph:
            # check if current node is merged into similar node in the graph
            if SpliceGraph._compare_is_merged(similar_node_in_graph, current_node):
                current_node.is_merged = True

                SpliceGraph.update_exon_coord_sr_svtype_breakpoints_name_mode(
                    similar_node_in_graph, current_node
                )

                merged_nodes_pool.add(current_node)

                # nodes in merged_parent_nodes are all in the graph
                current_node.merged_parent_nodes.append(similar_node_in_graph)

                # only consider nodes that has been processed: previous node in current series
                # keep in mind next node in current series is not processed yet!!!!
                # a -> b and b <- a
                similar_node_in_graph.add_predecessor(
                    current_node.previous_node_in_series
                )

    def _check_if_current_node_added_in_graph_and_update_predecessor_successor(
        self,
        current_node: NodeType,
        similar_key: str,
        merged_nodes_pool: Set[NodeType],
    ) -> None:
        """Check if current node is added in graph and update predecessor and successor."""
        if not current_node.is_merged:  # false

            current_node.is_in_graph = True  # check if node is in graph
            self.add_node_with_similar_key(current_node)

            # only consider nodes that has been processed: previous node in series
            # keep in mind next node in series is not processed yet!!!!
            current_node.add_predecessor(current_node.previous_node_in_series)

            # check merged node to see if merge node can be merged into current node
            for merge_node in [
                merge_node
                for merge_node in merged_nodes_pool
                if merge_node.similar_key == similar_key
            ]:
                if SpliceGraph._compare_is_merged(merge_node, current_node):
                    SpliceGraph.update_exon_coord_sr_svtype_breakpoints_name_mode(
                        current_node, merge_node
                    )
                    merge_node.merged_parent_nodes.append(current_node)
                    # for merge node whose previous node and next node in series
                    # has been processed
                    current_node.add_predecessor(merge_node.previous_node_in_series)
                    current_node.add_successor(merge_node.next_node_in_series)

    def construct(self) -> None:
        """Main function to construct graph."""
        # iterate all series
        merged_nodes_pool: Set[NodeType] = set()
        self.logger.trace(f"{self.series_list=}")
        for series in self.series_list:
            # iterate all nodes in series
            for index, current_node in enumerate(series):
                self.logger.trace(f"{current_node=}")
                self.logger.trace(f"{self.nodes=}")
                self.logger.trace(f"{merged_nodes_pool=}")
                # add information about  next and previous node in series to current node
                current_node.update_next_and_previous_node_in_series(index, series)
                # initialize and get unique key of current node and set node.unique_key
                # if not set when you reach node.unique_key, will return None
                _ = current_node.get_unique_key()

                # get similar key(chrom and intron) of current node
                similar_key = current_node.similar_key
                self._check_if_current_node_is_merged_in_similar_nodes_in_graph(
                    current_node, similar_key, merged_nodes_pool
                )

                self._check_if_current_node_added_in_graph_and_update_predecessor_successor(
                    current_node, similar_key, merged_nodes_pool
                )

    def _trace(
        self,
        start_node: NodeType,
        path: List[NodeType],
        group_paths: List[List[NodeType]],
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
                    self._trace(successor, path + [start_node], group_paths)
            else:
                self._trace(successors, path + [start_node], group_paths)  # type: ignore

    def trace(self) -> Any:
        """Trace through graph and find all paths."""
        result_series_list = []
        for start_node in self.get_start_nodes():
            group_paths: Any = []
            self._trace(start_node, [], group_paths)
            result_series_list.extend(group_paths)
        return result_series_list

    def _prune(self):
        """Helper function to prune graph."""

    def prune(self) -> None:
        """Prune graph."""
        self._prune()
