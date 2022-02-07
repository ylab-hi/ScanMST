# !/usr/bin/env python
"""Find cliques in a graph.

@Filename:    cliqueFinder.py
@license:     MIT Licence
@Time:        1/19/22 7:59 PM
"""
from itertools import combinations
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import networkx as nx  # type: ignore
from networkx import find_cliques  # type: ignore

from ..utils import timeit
from .basicClass import Node
from .basicClass import Series
from .exception import ExonsNotFoundError
from .type import LoggerType


class Ruler:
    """Calculate the similarity distance between two series.

    using longer one as the reference
    """

    def __init__(self, logger: LoggerType) -> None:
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

        chrm1, pos1 = bp1.split(":")
        chrm2, pos2 = bp2.split(":")
        return abs(int(pos1) - int(pos2)) if chrm1 == chrm2 else float("inf")

    @staticmethod
    def first_node_last_node_distance(first_node: Node, last_node: Node) -> float:
        """Calculate distance between first node of Series A and last node of Series B.

        :param first_node: the first node of Series A
        :param last_node: the last node of Series B
        :return: calculated breakpoint distance

        ..note::
                         [x]-[x]-[x]-[x]
             [x]-[x]-[x]-[x]
             * The output distance will be [0, 1]
        """
        if first_node.exons is None or last_node.exons is None:
            raise SystemExit from ExonsNotFoundError

        _ft_strand = first_node.strand
        _lt_strand = last_node.strand

        if _ft_strand != _lt_strand:
            return 1.0

        distance = 1.0
        #       [xxxx]-->--
        # -->--[xxxx]
        first_node_first_exon_start = first_node.exons[0][0]
        first_node_last_exon_end = first_node.exons[-1][1]
        last_node_first_exon_start = last_node.exons[0][0]
        last_node_last_exon_end = last_node.exons[-1][1]

        if _ft_strand == _lt_strand == "+":
            if (
                last_node_first_exon_start
                <= first_node_first_exon_start
                < last_node_last_exon_end
                <= first_node_last_exon_end
            ):
                overlapped_len = last_node_last_exon_end - first_node_first_exon_start

                _ft_cov = overlapped_len / (
                    first_node_last_exon_end - first_node_first_exon_start
                )
                _lt_cov = overlapped_len / (
                    last_node_last_exon_end - last_node_first_exon_start
                )
                distance = 1 - (_ft_cov + _lt_cov) / 2

        #  --<--[xxxx]
        #         [xxxx]--<--
        elif (
            first_node_first_exon_start
            <= last_node_first_exon_start
            < first_node_last_exon_end
            <= last_node_last_exon_end
        ):
            overlapped_len = first_node_last_exon_end - last_node_first_exon_start

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
        normalization_factor = 100
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

        ave_distance = (
            sum(distance_list) / (effect_num_pair * 2) if effect_num_pair > 0 else 100
        )
        distance = ave_distance / normalization_factor
        return distance

    @staticmethod
    def __decide_flag(
        left_query_node: Node,
        right_query_node: Node,
        left_subject_node: Optional[Node],
        right_subject_node: Optional[Node],
    ) -> bool:
        """Decide the flag.

        :param left_query_node: left query node
        :param right_query_node: right query node
        :param left_subject_node: left subject node
        :param right_subject_node: right subject node
        :return: flag
        """
        flag = False

        if left_subject_node is None and right_subject_node is None:
            return True

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
        ref:           [x]-[x]-[x]-[x]
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
            _dummy_bp = "chrN:0"
            series_a_bp_pair.insert(0, ("NA", _dummy_bp, _dummy_bp))
            series_a_bp_pair.append(("NA", _dummy_bp, _dummy_bp))

        for i in range(len(series_a_bp_pair) - sliding_window_size + 1):
            _subject_bp_pair = series_a_bp_pair[i : i + sliding_window_size]
            left_query_node = series_b[0]
            right_query_node = series_b[-1]

            left_subject_node = (
                series_a[i + 1 - sliding_window_size]
                if i >= sliding_window_size
                else None
            )

            right_subject_node = None
            if (
                0
                <= i
                < len(series_a_bp_pair) - sliding_window_size + 1 - sliding_window_size
            ):
                right_subject_node = series_a[i + 1]

            if Ruler.__decide_flag(
                left_query_node, right_query_node, left_subject_node, right_subject_node
            ):
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

    def __init__(
        self,
        intact_series_list: Any,
        intact_series_list_len: int,
        logger: LoggerType,
        threshold: float = 0.2,
    ):
        """Initialize CliqueFinder."""
        self.ruler = Ruler(logger)
        self.intact_series_list_len = intact_series_list_len
        self.intact_series_list = intact_series_list
        self.distance_dict: Dict[Tuple[int, int], float] = {}
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
        self._create_graph_for_series()

        for clique_index in find_cliques(self.graph):
            yield (self.intact_series_list[i] for i in clique_index)
