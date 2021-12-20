# !/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Filename:    spliceGraph.py
@license:     MIT Licence
@Time:        12/15/21 10:42 AM
"""
from collections import defaultdict
from typing import Any
from typing import List
from typing import Union

from loguru._logger import Logger

from .basicClass import Node
from .basicClass import Series


class Ruler:
    """calculate the similarity distance between two series,
    using longer one as the reference
    """

    def __init__(self, logger: Logger) -> None:
        self.logger = logger

    @staticmethod
    def obtain_breakpoint_pairs(series: Series) -> List:
        """generate breakpoint pairs from a series
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
        """calculate breakpoint distance sv_type1,chrA:pos1 VS sv_type2,chrB:pos2

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
    def first_node_last_node_distance(first_node: Node, last_node: Node) -> float:
        """calculate distance between first node of Series A and last node of Series B

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
        """calculate breakpoint distance

        :param bp_pair1: breakpoint pair list 1: [(sv_type, bp1, bp2), ...]
        :param bp_pair2: breakpoint pair list 2: [(sv_type, bp1, bp2), ...]
        :return: calculated breakpoint distance

        ..note::
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
        distance = ave_distance / (max(distance_list) - min(distance_list))
        return distance

    @staticmethod
    def __decide_flag(
        left_query_node: Node,
        right_query_node: Node,
        left_subject_node: Node,
        right_subject_node: Node,
    ):

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
        for i in range(sliding_window_size - 1):
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


class Graph:
    pass


class SpliceGraph(object):
    """
    the SpliceGraph class is used to trace the path of splice graph

    :Example:

    >>> from loguru import logger
    >>> series1 = Series(blat=None, logger=logger)
    >>> series1.nodes = [ Node(prev_bp=None,next_bp='chr17:7702552',strand='+',chrom='chr17',ref_start=7701656,ref_end=7702552,
    ...                 exons=[[7701656, 7702552]], sv_type='TRA'), Node(prev_bp='chr1:15872815',next_bp=None,strand='+',chrom='chr1',ref_start=15872815,ref_end=15873800,exons=[[15872815,15873800]], sv_type=None)]
    >>> series2 = Series(blat=None, logger=logger)
    >>> series2.nodes = [ Node(prev_bp=None,next_bp='chr17:7702552',strand='+',chrom='chr17',ref_start=7701500,ref_end=7702552,exons=[[7701500, 7702552]], sv_type='TRA'), Node(prev_bp='chr1:15872815',next_bp=None,strand='+',chrom='chr1',ref_start=15872815,ref_end=15876678,exons=[[15872815,15876678]], sv_type=None)]
    >>> series3 = Series(blat=None, logger=logger)
    >>> series3.nodes = [ Node(prev_bp=None,next_bp='chr1:15873900',strand='+',chrom='chr1',ref_start=15872890,ref_end=15873900,exons=[[15872890, 15873900]], sv_type='TRA', insertion_info=(False, NovelInsertion(hit_num=1, query_sequence='ATCGATCG'))), Node(prev_bp='chr17:872815',next_bp=None,strand='+',chrom='chr17',ref_start=872815,ref_end=876678,exons=[[872815,876678]], sv_type=None)]
    >>> series4 = Series(blat=None, logger=logger)
    >>> series4.nodes = [ Node(prev_bp=None,next_bp='chr1:15873900',strand='+',chrom='chr1',ref_start=15872890,ref_end=15873900,exons=[[15872890, 15873900]], sv_type='TRA', insertion_info=(False, NovelInsertion(hit_num=1, query_sequence='ATCGATCG'))), Node(prev_bp='chr17:872815',next_bp=None,strand='+',chrom='chr17',ref_start=872815,ref_end=876678,exons=[[872815,873400], [875500,876678]], sv_type=None)]
    >>> series5 = Series(blat=None, logger=logger)
    >>> series5.nodes = [Node(prev_bp=None,next_bp='chr17:7702552',strand='+',chrom='chr17',ref_start=7701656,ref_end=7702552,exons=[[7701656, 7702552]], sv_type='TRA'), Node(prev_bp='chr1:15872815',next_bp='chr1:15873900',strand='+',chrom='chr1',ref_start=15872815,ref_end=15873900,exons=[[15872815, 15873900]], sv_type='TRA', insertion_info=(False, NovelInsertion(hit_num=1, query_sequence='ATCGATCG'))), Node(prev_bp='chr17:872815',next_bp=None,strand='+',chrom='chr17',ref_start=872815,ref_end=876678,exons=[[872815,876678]], sv_type=None)]
    >>> series6 = Series(blat=None, logger=logger)
    >>> series6.nodes = [Node(prev_bp=None,next_bp='chr17:7702552',strand='+',chrom='chr17',ref_start=7701656,ref_end=7702552,exons=[[7701656, 7702552]], sv_type='TRA'), Node(prev_bp='chr1:15872815',next_bp='chr1:15873900',strand='+',chrom='chr1',ref_start=15872815,ref_end=15873900,exons=[[15872815, 15873900]], sv_type='TRA', insertion_info=None), Node(prev_bp='chr17:872815',next_bp=None,strand='+',chrom='chr17',ref_start=872815,ref_end=876678,exons=[[872815,876678]], sv_type=None)]
    >>> series7 = Series(blat=None, logger=logger)
    >>> series7.nodes = [Node(prev_bp=None,next_bp='chr17:7702552',strand='+',chrom='chr17',ref_start=7701656,ref_end=7702552,exons=[[7701656, 7702552]], sv_type='TRA'), Node(prev_bp='chr1:15872815',next_bp='chr1:15873900',strand='+',chrom='chr1',ref_start=15872815,ref_end=15873900,exons=[[15872815, 15873900]], sv_type='TRA', insertion_info=None), Node(prev_bp='chr17:872815',next_bp=None,strand='+',chrom='chr17',ref_start=872815,ref_end=876678,exons=[[872815,873400], [875500,876678]], sv_type=None)]
    >>> splice_graph = SpliceGraph(series_list=[series1, series2, series3, series4, series5, series6, series7], logger=logger)
    >>> splice_graph.construct()
    >>> splice_graph.trace()
    >>> splice_graph.result_series_list

    """

    def __init__(self, series_list, logger):
        self.series_list = series_list
        self.logger = logger
        self.nodes = defaultdict(list)
        self.unique_nodes_map = {}
        self.result_series_list = []

    def get_start_nodes(self):
        return [
            node
            for nodes in self.nodes.values()
            for node in nodes
            if node.is_start_node
        ]

    def get_similar_nodes(self, similar_key: str) -> List[Any]:
        return self.nodes[similar_key]

    def add_similar_node(self, node):
        key = node.similar_key
        self.nodes[key].append(node)

    def __contains__(self, node):
        unique_key = node.unique_key

        for nodes in self.nodes.values():
            for other_node in nodes:
                if other_node.unique_key == unique_key:
                    return True
                else:
                    return False

    def __iter__(self):
        for nodes in self.nodes.values():
            for node in nodes:
                yield node

    @staticmethod
    def _compare_is_merged_helper(node1, node2):

        condition = node1.sv_type == node2.sv_type

        if (
            condition and node1.prev_bp is None and node2.prev_bp is None
        ):  # both are start nodel
            return node1.exons[-1][1] == node2.exons[-1][1]  # check last exon end

        elif (
            condition and node1.next_bp is None and node2.next_bp is None
        ):  # both are end nodes
            return node1.exons[0][0] == node2.exons[0][0]  # check first exon start

        elif (
            condition and node1.prev_bp is None and node2.prev_bp is not None
        ):  # node1 is start node, node2 is middle node
            return (
                node1.exons[-1][1] == node2.exons[-1][1]
                and node1.exons[0][0] >= node2.exons[0][0]
            )

        elif (
            condition and node1.next_bp is None and node2.next_bp is not None
        ):  # node1 is end node, node2 is middle node
            return (
                node1.exons[0][0] == node2.exons[0][0]
                and node1.exons[-1][1] <= node2.exons[-1][1]
            )

        elif (
            node1.next_bp is None and node2.prev_bp is None
        ):  # node1 is end node, node2 is start node
            return (
                node2.exons[0][0] >= node1.exons[0][0]
                and node2.exons[-1][1] >= node1.exons[-1][1]
            )

        else:  # both are middle nodes TODO: the condition may need to more tight

            return (
                node1.exons[0][0] == node2.exons[0][0]
                and node1.exons[-1][1] == node2.exons[-1][1]
            )

    @staticmethod
    def _compare_is_merged(node1, node2):
        """node1 is similar as node2 is precommit of the function

         compare if node1 can merge node2
        :param node1:
        :param node2:
        :return:
        """
        flag1 = SpliceGraph._compare_is_merged_helper(node1, node2)
        if flag1:
            return True
        flag2 = SpliceGraph._compare_is_merged_helper(node2, node1)
        if flag2:
            return True

    @staticmethod
    def update_exon_coord_sr(updated_node, current_node):
        updated_node.exons[0][0] = min(
            updated_node.exons[0][0], current_node.exons[0][0]
        )
        updated_node.exons[-1][1] = max(
            updated_node.exons[-1][1], current_node.exons[-1][1]
        )
        updated_node.update_sr()

    @staticmethod
    def add_predecessor(node: Node, predecessor: Node):
        if node is None:
            node.predecessor.append(predecessor)

    @staticmethod
    def add_successor(node, successor):
        if node is None:
            node.successor.append(successor)

    def construct(self):
        # iterate all series
        for series in self.series_list:
            # iterate all nodes in series
            for index, current_node in enumerate(series):
                # update next and previous node in series
                current_node.update_next_and_previsous_node_in_series(index, series)
                # get unique key of current node
                unique_key = current_node.get_unique_key()
                self.unique_nodes_map[unique_key] = current_node

                # get similar key(chrom and intron) of current node
                similar_key = current_node.similar_key
                # get similar nodes in the graph
                similar_nodes_in_graph = self.get_similar_nodes(similar_key)
                # iterate all similar nodes in the graph
                if similar_nodes_in_graph:
                    is_merged = False  # flag to check if current node is merged
                    # iterate all similar nodes in the graph
                    for similar_node_in_graph in similar_nodes_in_graph:
                        # check if current node is merged into similar node in the graph
                        if SpliceGraph._compare_is_merged(
                            similar_node_in_graph, current_node
                        ):
                            is_merged = True
                            current_node.is_merged = True

                            SpliceGraph.update_exon_coord_sr(
                                similar_node_in_graph, current_node
                            )

                            similar_node_in_graph.merge_nodes.append(current_node)

                            similar_node_in_graph.predecessors.extend(
                                current_node.predecessors
                            )

                            current_node.next_node_in_series.predecessors.append(
                                similar_node_in_graph
                            )

                        if not is_merged:  # false
                            for merge_node in similar_node_in_graph.merge_nodes:
                                if SpliceGraph._compare_is_merged(
                                    merge_node, current_node
                                ):
                                    SpliceGraph.update_exon_coord_sr(
                                        current_node, merge_node
                                    )
                                    current_node.merge_nodes.append(merge_node)
                                    current_node.successors.extend(
                                        merge_node.successors
                                    )

                                    if index <= len(series) - 1:
                                        successor_node = series[index + 1]
                                        current_node.add_successor(successor_node)
                                    if index > 0:
                                        predecessor_node = series[index - 1]
                                        current_node.add_predecessor(predecessor_node)

                    # add node to graph
                    if not is_merged:  # false
                        self.add_similar_node(current_node)

                        if index <= len(series) - 1:
                            successor_node = series[index + 1]
                            current_node.add_successor(successor_node)
                        if index > 0:
                            predecessor_node = series[index - 1]
                            current_node.add_predecessor(predecessor_node)

                else:
                    self.add_similar_node(current_node)

                    if index <= len(series) - 1:
                        successor_node = series[index + 1]
                        current_node.add_successor(successor_node)
                    if index > 0:
                        predecessor_node = series[index - 1]
                        current_node.add_predecessor(predecessor_node)

    def _trace(self, start_node, path, group_paths):

        if not start_node:
            group_paths.append(path)

        else:
            successors = start_node.successor
            if successors:
                for successor in successors:
                    self._trace(successor, path + [start_node], group_paths)
            else:
                self._trace(successors, path + [start_node], group_paths)

    def trace(self):
        for start_node in self.get_start_nodes():
            group_paths = []
            self._trace(start_node, [], group_paths)
            self.result_series_list.append(group_paths)

    def run(self):
        self.construct()
        self.trace()
