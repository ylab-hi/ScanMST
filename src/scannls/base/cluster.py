# !/usr/bin/env python
"""Find cliques in a graph.

@Filename:    cliqueFinder.py
@license:     MIT Licence
@Time:        1/19/22 7:59 PM
"""
from itertools import combinations
from typing import Any

import networkx as nx
from loguru import logger
from networkx import connected_components

from ..utils import timeit
from .basicClass import BreakPoint
from .basicClass import Node
from .basicClass import Series
from .mergeCondition import MergeCondition


def middle_node_signature(node: Node) -> str:
    """Middle node signature using chrom, exons and strand."""
    chrom = node.chrom
    exons = node.exons
    strand = node.strand
    exons_string = map(lambda x: f"{x[0]}-{x[1]}", exons)

    return f"{chrom}:{';'.join(exons_string)};{strand}"


def _compare_is_merged_helper_check_condition_for_two_middle_nodes_list(
    node_list1: list[Node], node_list2: list[Node]
) -> bool:
    """Check if two node list of middle nodes have shared node or not.
    :param node_list1:  node_list1
    :param node_list2:  node_list2
    :return:  True if two node list have shared node, otherwise False
    """
    shared_middle_nodes = set(map(middle_node_signature, node_list1)) & set(
        map(middle_node_signature, node_list2)
    )

    return len(shared_middle_nodes) > 0


class Ruler:
    """Calculate the similarity distance between two series.

    using longer one as the reference
    """

    def __init__(self, prune_threshold: int = 10) -> None:
        """Initialize Ruler.

        :param logger: logger
        """
        self.prune_threshold = prune_threshold

    def __repr__(self):
        """Represent Ruler."""
        return f"{self.__class__.__name__}()"

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

        head_node_a = series_a[0]
        head_node_b = series_b[0]

        tail_node_a = series_a[-1]
        tail_node_b = series_b[-1]

        middle_nodes_a: list[Node] = series_a[1:-1]  # type: ignore
        middle_nodes_b: list[Node] = series_b[1:-1]  # type: ignore
        connection = False
        merge_condition: MergeCondition = MergeCondition(self.prune_threshold)

        if (
            merge_condition.head2head(head_node_a, head_node_b)
            or merge_condition.tail2tail(tail_node_a, tail_node_b)
            or merge_condition.head2tail(head_node_a, tail_node_b)
            or merge_condition.head2tail(head_node_b, tail_node_a)
            or _compare_is_merged_helper_check_condition_for_two_middle_nodes_list(
                middle_nodes_a, middle_nodes_b
            )
        ):
            connection = True

        if not connection:
            for _middle_node_b in middle_nodes_b:
                if merge_condition.head2mid(head_node_a, _middle_node_b):
                    connection = True
                    break
                if merge_condition.tail2mid(tail_node_a, _middle_node_b):
                    connection = True
                    break

            if not connection:
                for _middle_node_a in middle_nodes_a:
                    if merge_condition.head2mid(head_node_a, _middle_node_a):
                        connection = True
                        break
                    if merge_condition.tail2mid(tail_node_b, _middle_node_a):
                        connection = True
                        break

        return 0.0 if connection else 1.0


class ClusterFinder:
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

    def __init__(
        self,
        intact_series_list: Any,
        intact_series_list_len: int,
        threshold: float = 0.2,
    ) -> None:
        """Initialize CliqueFinder."""
        self.ruler = Ruler()
        self.intact_series_list_len = intact_series_list_len
        self.intact_series_list = intact_series_list
        self.distance_dict: dict[tuple[int, int], float] = {}
        self.graph = nx.Graph()
        self.threshold = threshold

    def _calculate_distance(self, x: int, y: int) -> float:
        """Calculate distance between two series. If distance has been calculated before.

        return True and distance value. Otherwise, calculate distance and return False and
        distance value.

        :param x: series x
        :param y: series y
        :return: is_calculated, distance value
        """
        return self.ruler(self.intact_series_list[x], self.intact_series_list[y])

    def _add_edge_between_two_series(self, x: int, y: int) -> None:
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
        if self._calculate_distance(x, y) < self.threshold:
            self.graph.add_edge(x, y)
            if not self.intact_series_list[x].is_in_graph:
                self.intact_series_list[x].is_in_graph = True

            if not self.intact_series_list[y].is_in_graph:
                self.intact_series_list[y].is_in_graph = True

    def _create_graph_for_series(self) -> None:
        """Create graph for all series in intact_series_list.

        add edge between two series in terms of the distance value

        :return: None
        """
        last_x = 0
        ind_x, ind_y = 0, 0

        for ind_x, ind_y in combinations(range(self.intact_series_list_len), 2):
            self._add_edge_between_two_series(ind_x, ind_y)
            if last_x != ind_x:
                if not self.intact_series_list[last_x].is_in_graph:
                    self.graph.add_node(last_x)
                last_x = ind_x

        # solve last two node
        if not self.intact_series_list[ind_x].is_in_graph:
            self.graph.add_node(last_x)

        if not self.intact_series_list[ind_y].is_in_graph:
            self.graph.add_node(ind_y)

    @timeit
    def find_cluster(self) -> Any:
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
        self._create_graph_for_series()

        for clique_index in connected_components(self.graph):
            yield (self.intact_series_list[i] for i in clique_index)

    def find_cluster_index(self):
        """Find clique in graph with help of :func:`networkx.algorithms.clique.find_clique`.

        :return:  every clique in graph as a iterator (List[int])
        """
        self._create_graph_for_series()
        yield from connected_components(self.graph)

    def creat_merge_indexs(self, cluster) -> dict[int, list[str]]:
        result = {}

        for series in cluster:
            nodes_key = []
            for node in series:
                nodes_key.append(create_merge_key_for_node(node))

            result[series.id] = nodes_key

        return result

    @staticmethod
    def check_merge(series1: Series, series2: Series, merge_keys: dict[int, list[str]]):
        """Check if series1 can merge series2."""
        series1_nodes_key = merge_keys[series1.id]
        series2_nodes_key = merge_keys[series2.id]

        if len(merge_keys[series1.id]) == len(merge_keys[series2.id]):
            # reduce duplication
            return False
        elif len(merge_keys[series1.id]) > len(merge_keys[series2.id]):
            if set(series2_nodes_key).issubset(set(series1_nodes_key)):
                start_index = series1_nodes_key.index(series2_nodes_key[0])
                if merge_same_len_node_list(series1[start_index:], series2, 1):  # type: ignore
                    merge_series(series1, series2, start_index)  # type: ignore
                    series1.merge_factor += 1
                    return True
        else:
            raise ValueError("series1 is shorter than series2")

    def merge_cluster(self):
        for cluster_index in self.find_cluster_index():
            series_list = []
            for i in cluster_index:
                current_series = self.intact_series_list[i]
                current_series.id = i
                series_list.append(current_series)

            sorted_series = sort_cluster(series_list)
            merge_keys = self.creat_merge_indexs(sorted_series)
            new_cluster = []
            ClusterFinder._merge_cluster(sorted_series, new_cluster, merge_keys)
            yield sort_cluster(
                new_cluster,
                key=lambda x: create_sort_key_by_merge_factor(x),  # type: ignore
                reverse=True,
            )

    @staticmethod
    def _merge_cluster(series_list, result, merge_keys):
        if not series_list:
            return result

        slected_series = series_list.pop()

        for current_series in series_list:
            if ClusterFinder.check_merge(slected_series, current_series, merge_keys):
                logger.warning("merge one series")
                series_list.remove(current_series)

        result.append(slected_series)

        ClusterFinder._merge_cluster(series_list, result, merge_keys)


def merge_series(series1: list[Node], series2: list[Node], start_index: int):
    assert len(series1[start_index:]) == len(series2)
    for updated_node, current_node in zip(series1[start_index:], series2):
        # update exon coordinates
        updated_node.ref_start = min(  # type: ignore
            updated_node.exons[0][0], current_node.exons[0][0]  # type: ignore
        )

        updated_node.exons[0] = updated_node.ref_start, updated_node.exons[0][1]  # type: ignore

        updated_node.ref_end = max(  # type: ignore
            updated_node.exons[-1][1], current_node.exons[-1][1]  # type: ignore
        )

        updated_node.exons[-1] = updated_node.exons[-1][0], updated_node.ref_end  # type: ignore

        updated_node.sr += current_node.sr


def merge_same_len_node_list(series1: list[Node], series2: list[Node], threashold: int):
    """seires1 is equal than series2 and series1 merge series2.

    orignial s1: [ ] - [ ] - [ ] - [ ]
    s2:                [ ] - [ ] - [ ]

    """
    assert len(series1) == len(series2)

    merge_condition = MergeCondition(threashold)

    for node1, node2 in zip(series1, series2):
        if node1.self_identity.is_head() and node2.self_identity.is_head():
            return BreakPoint.equal(
                node1.next_breakpoint, node2.next_breakpoint, threashold
            ) and merge_condition.head2head(node1, node2)

        elif node1.self_identity.is_mid() and node2.self_identity.is_head():
            return BreakPoint.equal(
                node1.next_breakpoint, node2.next_breakpoint, threashold
            ) and merge_condition.mid2head(node1, node2)

        elif node1.self_identity.is_mid() and node2.self_identity.is_mid():
            return BreakPoint.equal(
                node1.next_breakpoint, node2.next_breakpoint, threashold
            ) and merge_condition.mid2mid(node1, node2)

        elif node1.self_identity.is_mid() and node2.self_identity.is_tail():
            return merge_condition.mid2tail(node1, node2)

        elif node1.self_identity.is_tail() and node2.self_identity.is_tail():
            return merge_condition.tail2tail(node1, node2)

        else:
            raise ValueError("invalid node identity")


def create_merge_key_for_node(node: Node):
    introns = node.introns
    introns_key = "-".join([f"{i}-{j}" for i, j in introns]) if introns else "None"
    chrom = node.chrom
    sv_type = node.sv_type
    return f"{chrom}_{sv_type}_{introns_key}"


def create_merge_key_for_series(nodes_key: list[str]):
    return "_".join(nodes_key)


def create_sort_key_for_node(node: Node):
    return middle_node_signature(node)


def creat_sort_key_for_series(series: Series):
    return (
        len(series),
        *[create_sort_key_for_node(node) for node in series],
    )


def create_sort_key_by_merge_factor(series: Series):
    return (
        len(series),
        series.merge_factor,
        *[create_sort_key_for_node(node) for node in series],
    )


def sort_cluster(cluster, key=lambda x: creat_sort_key_for_series(x), reverse=False):
    return sorted(
        cluster,
        key=key,
        reverse=reverse,
    )
