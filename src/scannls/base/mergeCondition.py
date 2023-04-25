"""Merge condition.
# !/usr/bin/env python

@Filename:    mergeCondition.py
@contact:     yangyang.li@northwestern.edu
@license:     MIT Licence
@Time:        4/18/22 7:51 PM
"""
from itertools import zip_longest

from .basicClass import Node
from .exception import ExonsNotFoundError


def _compare_is_merged_helper_check_condition_for_two_heads_nodes_mode(
    node1: Node, node2: Node, threshold: int
) -> bool:
    """Check if both head nodes can be merged.

    :param threshold:

    .. note::
        nodes with different length may be merged. []: exon -: intron
        node1: []-[]-[]
        node2:    []-[]

    """
    # no introns and breakpoints are same
    if not node1.introns and not node2.introns:
        return True

    introns_group = (
        zip_longest(node1.introns, node2.introns)
        if node1.strand == "-"
        else zip_longest(node1.introns[::-1], node2.introns[::-1])
    )

    # have introns and breakpoints are same
    for node1_intron, node2_intron in introns_group:
        if node1_intron != node2_intron:
            return node1_intron is None or node2_intron is None

    return True


def _compare_is_merged_helper_check_condition_for_head_and_tail_nodes_mode(
    node1: Node,
    node2: Node,
    threshold: float = 0.5,
) -> bool:
    """Check if node1 and node2 can be merged based on overlap info.

    node1 is tail node, node2 is head node Using mean overlap ratio to
    check if they can be merged.

    :param node1:  node1
    :param node2:  node2
    :param threshold:  threshold for checking if two nodes are merged
    :return:  True if two nodes are merged, otherwise False

    .. note::

        -> [node1]               [node1] <-
            [node2] ->      <- [node2]
    """
    if node1.exons is None or node2.exons is None:
        raise ExonsNotFoundError(f"{node1.query_name} or {node2.query_name}")

    # limit all introns
    if node1.introns != node2.introns:
        return False

    node1_first_exon_start = node1.exons[0][0]
    node1_last_exon_end = node1.exons[-1][1]
    node2_first_exon_start = node2.exons[0][0]
    node2_last_exon_end = node2.exons[-1][1]

    if node1.is_polya:
        return False

    if node1.strand == "+":
        condition = (
            node1_first_exon_start
            <= node2_first_exon_start
            < node1_last_exon_end
            <= node2_last_exon_end
        )

        if condition:
            overlap_len = node1_last_exon_end - node2_first_exon_start
            union_len = node2_last_exon_end - node1_first_exon_start
            return overlap_len / union_len >= threshold
    else:
        condition = (
            node2_first_exon_start
            <= node1_first_exon_start
            < node2_last_exon_end
            <= node1_last_exon_end
        )
        if condition:
            overlap_len = node2_last_exon_end - node1_first_exon_start
            union_len = node2_last_exon_end - node1_first_exon_start
            return overlap_len / union_len >= threshold

    return False


def _compare_is_merged_helper_check_condition_for_two_tail_nodes_mode(
    node1: Node, node2: Node
) -> bool:
    """Check if two end nodes can be merged or not.

    node1 is tail node, node2 is tail node
    check if they can be merged.

    :param threshold: threshold for checking if two nodes are merged
    :param node1:  node1
    :param node2:  node2
    :return:  True if two nodes are merged, otherwise False

    .. note::

        -> [node1]
        -> [node2]
    """
    if node1.exons is None or node2.exons is None:
        raise ExonsNotFoundError(f"{node1.query_name} or {node2.query_name}")

    # limit all introns
    if node1.introns != node2.introns:
        return False

    node1_first_exon_start = node1.exons[0][0]
    node1_last_exon_end = node1.exons[-1][1]
    node2_first_exon_start = node2.exons[0][0]
    node2_last_exon_end = node2.exons[-1][1]

    if node1.is_polya and node2.is_polya:
        return (
            node1_first_exon_start == node2_first_exon_start
            and node1_last_exon_end == node2_last_exon_end
        )

    if node1.is_polya and not node2.is_polya:
        if node1.strand == "+":
            return (
                node1_first_exon_start == node2_first_exon_start
                and node1_last_exon_end >= node2_last_exon_end
            )

        return (
            node1_last_exon_end == node2_last_exon_end
            and node1_first_exon_start <= node2_first_exon_start
        )

    if not node1.is_polya and node2.is_polya:
        if node1.strand == "+":
            return (
                node1_first_exon_start == node2_first_exon_start
                and node1_last_exon_end <= node2_last_exon_end
            )
        return (
            node1_last_exon_end == node2_last_exon_end
            and node1_first_exon_start >= node2_first_exon_start
        )

    if not node1.is_polya and not node2.is_polya:
        if node1.strand == "+":
            return node1_first_exon_start == node2_first_exon_start
        return node1_last_exon_end == node2_last_exon_end

    return False


def _compare_is_merged_helper_check_condition_for_head_and_middle_nodes_mode(
    node1: Node,
    node2: Node,
) -> bool:
    """Check if start node can be merged with a middle node or not.

    node1 is start node, node2 is middle node
    check if they can be merged.

    :param node1:  node1
    :param node2:  node2
    :return:  True if two nodes are merged, otherwise False

    .. note::

             [ node1 ] ->
        -> [  node2  ] ->
    """
    if node1.exons is None or node2.exons is None:
        raise ExonsNotFoundError(f"{node1.query_name} or {node2.query_name}")

    # limit all introns
    if node1.introns != node2.introns:
        return False

    node1_first_exon_start = node1.exons[0][0]
    node1_last_exon_end = node1.exons[-1][1]
    node2_first_exon_start = node2.exons[0][0]
    node2_last_exon_end = node2.exons[-1][1]

    if node1.strand == "+":
        return (
            node1_last_exon_end == node2_last_exon_end
            and node1_first_exon_start >= node2_first_exon_start
        )
    return (
        node1_first_exon_start == node2_first_exon_start
        and node1_last_exon_end <= node2_last_exon_end
    )


def _compare_is_merged_helper_check_condition_for_tail_and_middle_nodes_mode(
    node1: Node,
    node2: Node,
) -> bool:
    """Check if end node can be merged with a middle node or not.

    node1 is tail node, node2 is middle node
    check if they can be merged.

    :param node1:  node1
    :param node2:  node2
    :return:  True if two nodes are merged, otherwise False

    .. note::

        -> [ node1 ]
        -> [  node2  ] ->
    """
    if node1.exons is None or node2.exons is None:
        raise ExonsNotFoundError(f"{node1.query_name} or {node2.query_name}")

    # limit all introns
    if node1.introns != node2.introns:
        return False

    node1_first_exon_start = node1.exons[0][0]
    node1_last_exon_end = node1.exons[-1][1]
    node2_first_exon_start = node2.exons[0][0]
    node2_last_exon_end = node2.exons[-1][1]

    if node1.is_polya:
        return False

    if node1.strand == "+":
        return (
            node1_first_exon_start == node2_first_exon_start
            and node1_last_exon_end <= node2_last_exon_end
        )
    return (
        node1_last_exon_end == node2_last_exon_end
        and node1_first_exon_start >= node2_first_exon_start
    )
