# !/usr/bin/env python
"""Find cliques in a graph.

@Filename:    cliqueFinder.py
@license:     MIT Licence
@Time:        1/19/22 7:59 PM
"""
from itertools import combinations
from typing import Any

import networkx as nx
from networkx import connected_components

from ..utils import timeit
from .basicClass import Node
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

        middle_nodes_a = series_a[1:-1]
        middle_nodes_b = series_b[1:-1]
        connection = False

        if (
            _compare_is_merged_helper_check_condition_for_two_heads_nodes_mode(
                head_node_a, head_node_b, self.prune_threshold
            )
            or _compare_is_merged_helper_check_condition_for_two_tail_nodes_mode(
                tail_node_a, tail_node_b, self.prune_threshold
            )
            or _compare_is_merged_helper_check_condition_for_head_and_tail_nodes_mode(
                head_node_a, tail_node_b
            )
            or _compare_is_merged_helper_check_condition_for_head_and_tail_nodes_mode(
                head_node_b, tail_node_a
            )
            or _compare_is_merged_helper_check_condition_for_two_middle_nodes_list(
                middle_nodes_a, middle_nodes_b
            )
        ):
            connection = True

        if not connection:
            for _middle_node_b in middle_nodes_b:
                if _compare_is_merged_helper_check_condition_for_head_and_middle_nodes_mode(
                    head_node_a, _middle_node_b
                ):
                    connection = True
                    break
                if _compare_is_merged_helper_check_condition_for_tail_and_middle_nodes_mode(
                    tail_node_a, _middle_node_b
                ):
                    connection = True
                    break

            if not connection:
                for _middle_node_a in middle_nodes_a:
                    if _compare_is_merged_helper_check_condition_for_head_and_middle_nodes_mode(
                        head_node_b, _middle_node_a
                    ):
                        connection = True
                        break
                    if _compare_is_merged_helper_check_condition_for_tail_and_middle_nodes_mode(
                        tail_node_b, _middle_node_a
                    ):
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


# TODO: add duplication reduction and series mergement <04-24-23, Yangyang Li yangyang.li@northwestern.edu>
