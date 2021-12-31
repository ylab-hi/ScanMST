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
from typing import Tuple
from typing import Union

import networkx as nx  # type: ignore
from loguru._logger import Logger
from networkx.algorithms.clique import find_cliques  # type: ignore

from ..utils import timeit
from .basicClass import Insertion
from .basicClass import MicroHomology
from .basicClass import Node
from .basicClass import Read
from .basicClass import Series

NodeType = Union[Node, Insertion]


class Ruler:
    """Calculate the similarity distance between two series.

    using longer one as the reference
    """

    def __init__(self, logger: Logger) -> None:
        """Initialize Ruler.

        :param logger: logger
        """
        self.logger = logger

    def __repr__(self):
        """Represent Ruler."""
        return f"{self.__class__.__name__}()"

    @staticmethod
    def obtain_breakpoint_pairs(series: Series) -> List:
        """Generate breakpoint pairs from a series.

        :param series: a series
        :type series: Series object
        :return: a list of breakpoint pairs
        :rtype: list
        """
        breakpoints_pairs = []
        for idx in range(len(series) - 1):
            _sv_type = series[idx].sv_type
            _bp1 = series[idx].next_breakpoint
            _bp2 = series[idx + 1].prev_breakpoint
            breakpoints_pairs.append((_sv_type, _bp1, _bp2))
        return breakpoints_pairs

    @staticmethod
    def breakpoints_distance(
        sv_type1: str, sv_type2: str, bp1: str, bp2: str
    ) -> Union[float, int]:
        """Calculate breakpoint distance sv_type1,chrA:pos1 VS sv_type2,chrB:pos2.

        :param sv_type1: sv_type of breakpoint1
        :param sv_type2: sv_type of breakpoint2
        :param bp1: breakpoint1
        :param bp2: breakpoint2
        :return: calculated breakpoint distance
        """
        if sv_type1 != sv_type2:
            return float("inf")
        else:
            chrm1, pos1 = bp1.split(":")
            chrm2, pos2 = bp2.split(":")
            if chrm1 == chrm2:
                return abs(int(pos1) - int(pos2))
            else:
                return float("inf")

    @staticmethod
    def first_node_last_node_distance(
        first_node: NodeType, last_node: NodeType
    ) -> float:
        """Calculate distance between first node of Series A and last node of Series B.

        :param first_node: the first node of Series A
        :param last_node: the last node of Series B
        :return: calculated breakpoint distance

        ..note::
                         [x]-[x]-[x]-[x]
             [x]-[x]-[x]-[x]
             * The output distance will be [0, 1]
        """
        _ft_strand = first_node.strand
        _lt_strand = last_node.strand

        if _ft_strand != _lt_strand:
            return 1.0
        else:
            distance = 1.0
            #       [xxxx]-->--
            # -->--[xxxx]
            first_node_first_exon_start = first_node.exons[0][0]  # type: ignore
            first_node_last_exon_end = first_node.exons[-1][1]  # type: ignore
            last_node_first_exon_start = last_node.exons[0][0]  # type: ignore
            last_node_last_exon_end = last_node.exons[-1][1]  # type: ignore

            if _ft_strand == _lt_strand == "+":
                if (
                    last_node_first_exon_start
                    <= first_node_first_exon_start
                    < last_node_last_exon_end
                    <= first_node_last_exon_end
                ):
                    overlapped_len = (
                        last_node_last_exon_end - first_node_first_exon_start
                    )

                    _ft_cov = overlapped_len / (
                        first_node_last_exon_end - first_node_first_exon_start
                    )
                    _lt_cov = overlapped_len / (
                        last_node_last_exon_end - last_node_first_exon_start
                    )
                    distance = 1 - (_ft_cov + _lt_cov) / 2

            #  --<--[xxxx]
            #         [xxxx]--<--
            else:
                if (
                    first_node_first_exon_start
                    <= last_node_first_exon_start
                    < first_node_last_exon_end
                    <= last_node_last_exon_end
                ):
                    overlapped_len = (
                        first_node_last_exon_end - last_node_first_exon_start
                    )

                    _ft_cov = overlapped_len / (
                        first_node_last_exon_end - first_node_first_exon_start
                    )
                    _lt_cov = overlapped_len / (
                        last_node_last_exon_end - last_node_first_exon_start
                    )
                    distance = 1 - (_ft_cov + _lt_cov) / 2
            return distance

    @staticmethod
    def breakpoint_pairs_distance(bp_pair1: List, bp_pair2: List) -> float:
        """Calculate breakpoint distance.

        :param bp_pair1: breakpoint pair list 1: [(sv_type, bp1, bp2), ...]
        :param bp_pair2: breakpoint pair list 2: [(sv_type, bp1, bp2), ...]
        :return: calculated breakpoint distance

        .. note::
            len(bp_pair1) == len(bp_pair2) should be always true
            The output distance will be [0, 1]
        """
        distance_list = []
        effect_num_pair = 0
        for i, j in zip(bp_pair1, bp_pair2):
            a_sv_type, a_bp1, a_bp2 = i
            b_sv_type, b_bp1, b_bp2 = j
            if a_sv_type != "NA" and b_sv_type != "NA":
                distance_list.append(
                    Ruler.breakpoints_distance(a_sv_type, b_sv_type, a_bp1, b_bp1)
                )
                distance_list.append(
                    Ruler.breakpoints_distance(a_sv_type, b_sv_type, a_bp2, b_bp2)
                )
                effect_num_pair += 1
        ave_distance = sum(distance_list) / (effect_num_pair * 2)
        distance = ave_distance / (max(distance_list) - min(distance_list) + 1e-6)
        return distance

    @staticmethod
    def __decide_flag(
        left_query_node: NodeType,
        right_query_node: NodeType,
        left_subject_node: NodeType,
        right_subject_node: NodeType,
    ) -> bool:  # type: ignore
        """Decide the flag.

        :param left_query_node: left query node
        :param right_query_node: right query node
        :param left_subject_node: left subject node
        :param right_subject_node: right subject node
        :return: flag
        """
        flag = False

        if not left_subject_node and not right_subject_node:
            flag = True

        if left_subject_node:
            if (
                left_query_node.strand == left_subject_node.strand == "+"
                and left_query_node.exons[0][0] >= left_subject_node.exons[0][0]  # type: ignore
            ) or (
                left_query_node.strand == left_subject_node.strand == "-"
                and left_query_node.exons[-1][1] >= left_subject_node.exons[-1][1]  # type: ignore
            ):
                flag = True
            else:
                flag = False

        if right_subject_node:
            if (
                right_query_node.strand == right_subject_node.strand == "+"
                and right_query_node.exons[-1][1] <= right_subject_node.exons[-1][1]  # type: ignore
            ) or (
                right_query_node.strand == right_subject_node.strand == "-"
                and right_query_node.exons[0][0] <= right_subject_node.exons[0][0]  # type: ignore
            ):
                flag = True
            else:
                flag = False
        return flag

    def __call__(self, series_a: Series, series_b: Series) -> float:
        """Call Ruler to calculate the distance between two series.

        :param series_a: series a
        :param series_b: series b
        :return: distance

        :Example:

        >>> ruler = Ruler()
        >>> ruler(series_a, series_b)
        0.5
        """
        if len(series_a) < len(series_b):
            series_a, series_b = series_b, series_a
        series_a_bp_pair = Ruler.obtain_breakpoint_pairs(series_a)
        series_b_bp_pair = Ruler.obtain_breakpoint_pairs(series_b)

        calculated_distance_list = []

        """
        ref:        [x]-[x]-[x]-[x]
        query: [x]-[x]-[x]
        """
        distance = Ruler.first_node_last_node_distance(series_a[0], series_b[-1])
        calculated_distance_list.append(distance)

        """
        ref:     [O]-[x]-[x]-[x]-[x]-[O]
        query1:  [x]-[x]-[x]
        query2:      [x]-[x]-[x]
        query3:          [x]-[x]-[x]
        query4:              [x]-[x]-[x]
        """
        sliding_window_size = len(series_b_bp_pair)
        for _ in range(sliding_window_size - 1):
            _chrom = "chrN:0"
            series_a_bp_pair.insert(0, ("NA", _chrom, _chrom))
            series_a_bp_pair.append(("NA", _chrom, _chrom))

        for i in range(len(series_a_bp_pair) - sliding_window_size + 1):
            _subject_bp_pair = series_a_bp_pair[i : i + sliding_window_size]
            left_query_node = series_b[0]
            right_query_node = series_b[-1]

            left_subject_node = (
                series_a[i + 1 - sliding_window_size]
                if i >= sliding_window_size
                else None
            )

            if (
                0
                <= i
                < len(series_a_bp_pair) - sliding_window_size + 1 - sliding_window_size
            ):
                right_subject_node = series_a[i + 1]
            else:
                right_subject_node = None

            flag = Ruler.__decide_flag(
                left_query_node, right_query_node, left_subject_node, right_subject_node
            )

            if flag:
                distance = Ruler.breakpoint_pairs_distance(
                    series_b_bp_pair, _subject_bp_pair
                )
                calculated_distance_list.append(distance)
        """
         ref:    [x]-[x]-[x]-[x]
         query:              [x]-[x]-[x]
        """
        distance = Ruler.first_node_last_node_distance(series_b[0], series_a[-1])
        calculated_distance_list.append(distance)
        return min(calculated_distance_list) if calculated_distance_list else 1.0


class CliqueFinder:
    """Find cliques in a graph based on series level.

     which will help to construct splice graph base on nodes level in the future.

    :param intact_series_list: list of intact series for a bam file of one sample
    :param logger: logger
    :param threshold: threshold to determine whether two series are connected

    .. note::
        :function: `networkx.algorithms.clique.find_cliques` is used to find cliques.

    :Example:

    >>> from loguru import  logger
    >>> clique_finder = CliqueFinder([], logger)
    >>> clique_finder.find_clique()
    """

    def __init__(self, intact_series_list: Any, logger: Logger, threshold: float = 0.2):
        """Initialize CliqueFinder."""
        self.ruler = Ruler(logger)
        self.intact_series_list = intact_series_list
        self.distance_dict: Dict[Any, float] = dict()
        self.graph = nx.Graph()
        self.threshold = threshold

    def _calculate_distance(self, x: Series, y: Series) -> Tuple[bool, float]:
        """Calculate distance between two series. If distance has been calculated before.

        return True and distance value. Otherwise, calculate distance and return False and
        distance value.

        :param x: series x
        :param y: series y
        :return: is_calculated, distance value
        """
        distance1 = self.distance_dict.get((x, y), None)
        if distance1 is not None:
            return True, distance1
        distance2 = self.distance_dict.get((y, x), None)
        if distance2 is not None:
            return True, distance2

        distance = self.ruler(x, y)
        self.distance_dict[(x, y)] = distance
        return False, distance

    def _add_edge_between_two_series(self, x: Series, y: Series) -> None:
        """Add edge between two series according to the distance between them.

        if the distance is less than threshold, add edge. Otherwise, do nothing.

        :param x: series x
        :param y: series y
        :return: None

        .. note::
            if `is_calculated` is True, distance value is stored in distance_dict, which
            indicates that the two series have been checked and determined if they
            should be connected in graph.
        """
        if x != y:
            is_calculated, distance = self._calculate_distance(x, y)
            if not is_calculated and distance < self.threshold:
                self.graph.add_edge(x, y)
                x.is_in_graph = True

    @timeit
    def _creat_graph_for_series(self) -> None:
        """Create graph for all series in intact_series_list.

        add edge between two series in terms of the distance value

        :return: None
        """
        for x in self.intact_series_list:
            for y in self.intact_series_list:
                self._add_edge_between_two_series(x, y)

            if not x.is_in_graph:
                self.graph.add_node(x)

    @timeit
    def find_clique(self) -> Any:
        """Find clique in graph with help of :func:`networkx.algorithms.clique.find_clique`.

        :return:  every clique in graph as a iterator (List[Series])

        :Example:

        >>> from loguru import logger
        >>> clique_finder = CliqueFinder([], logger)
        >>> cliques = clique_finder.find_clique()
        >>> for clique in cliques:
        ...     for series_list in clique:
        ...         assert isinstance(series_list, Series)
        """
        self._creat_graph_for_series()

        yield from find_cliques(self.graph)


class SpliceGraph:
    """SpliceGraph class is used to trace the path of splice graph."""

    dict_factory = dict
    list_factory = list

    def __init__(self, logger: Logger):
        """Initialize SpliceGraph."""
        self.logger = logger
        self.dict_factory = SpliceGraph.dict_factory  # type: ignore
        self.list_factory = SpliceGraph.list_factory  # type: ignore

    def __call__(self, series_list: Iterable[Series]) -> Any:
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
        self.nodes: Dict[str, List[NodeType]] = self.dict_factory()
        self.construct()
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
        for other_node in self.get_nodes_with_similar_key(node.similar_key):
            if other_node.unique_key == node.unique_key:
                return True

        return False

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
            # None == None
            return flag
        elif insertion_info1 is not None and insertion_info2 is not None:
            if insertion_info1[0] and insertion_info2[0]:
                # 1 hit insertion that is added in the series
                return flag

            if not insertion_info1[0] and not insertion_info2[0]:
                if (
                    isinstance(insertion_info1[1], Read)
                    and isinstance(insertion_info2[1], Read)
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
        node1_first_exon_start = node1.exons[0][0]  # type: ignore
        node1_last_exon_end = node1.exons[-1][1]  # type: ignore
        node2_first_exon_start = node2.exons[0][0]  # type: ignore
        node2_last_exon_end = node2.exons[-1][1]  # type: ignore
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
    def _compare_is_merged_helper(node1: NodeType, node2: NodeType) -> Any:
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
            return node1.exons[-1][1] == node2.exons[-1][1]  # type: ignore

        elif (
            node1.next_breakpoint is None and node2.next_breakpoint is None
        ):  # both are end nodes  # check first exon start
            return node1.exons[0][0] == node2.exons[0][0]  # type: ignore

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
    def update_exon_coord_sr_svtype_breakpoints(
        updated_node: NodeType, current_node: NodeType
    ) -> None:
        """Update exon coordinates of the updated node based on current node.

        :param updated_node:  node has been inserted into graph
        :param current_node: node has not been inserted into graph
        :return:
        """
        updated_node.exons[0][0] = updated_node.ref_start = min(  # type: ignore
            updated_node.exons[0][0], current_node.exons[0][0]  # type: ignore
        )
        updated_node.exons[-1][1] = updated_node.ref_end = max(  # type: ignore
            updated_node.exons[-1][1], current_node.exons[-1][1]  # type: ignore
        )
        updated_node.update_sr()
        updated_node.sv_type = (
            current_node.sv_type
            if current_node.sv_type is not None
            else updated_node.sv_type
        )

        if updated_node.prev_breakpoint is None:
            updated_node.prev_breakpoint = current_node.prev_breakpoint
        if updated_node.next_breakpoint is None:
            updated_node.next_breakpoint = current_node.next_breakpoint

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

                SpliceGraph.update_exon_coord_sr_svtype_breakpoints(
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
                    SpliceGraph.update_exon_coord_sr_svtype_breakpoints(
                        current_node, merge_node
                    )
                    merge_node.merged_parent_nodes.append(current_node)
                    # for merge node whose previous node and next node in series
                    # has been processed
                    current_node.add_predecessor(merge_node.previous_node_in_series)
                    current_node.add_successor(merge_node.next_node_in_series)

    def construct(self):
        """Main function to construct graph."""
        # iterate all series
        merged_nodes_pool = set()
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

    def _trace(self, start_node: NodeType, path: List, group_paths: List) -> None:
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
