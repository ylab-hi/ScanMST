"""Merge condition.
# !/usr/bin/env python

@Filename:    mergeCondition.py
@contact:     yangyang.li@northwestern.edu
@license:     MIT Licence
@Time:        4/18/22 7:51 PM
"""
from enum import auto
from enum import Enum
from itertools import zip_longest

from ..exception import ExonsNotFoundError
from .basicClass import Node


class MergeConditionMode(Enum):
    head2head = auto()
    head2mid = auto()
    head2tail = auto()
    mid2head = auto()
    mid2mid = auto()
    mid2tail = auto()
    tail2head = auto()
    tail2mid = auto()
    tail2tail = auto()

    @classmethod
    def get_mode(cls, node1: Node, node2: Node) -> "MergeConditionMode":
        if node1.self_identity.is_head() and node2.self_identity.is_head():
            return cls.head2head
        elif node1.self_identity.is_head() and node2.self_identity.is_tail():
            return cls.head2tail
        elif node1.self_identity.is_head() and node2.self_identity.is_mid():
            return cls.head2mid
        elif node1.self_identity.is_mid() and node2.self_identity.is_head():
            return cls.mid2head
        elif node1.self_identity.is_mid() and node2.self_identity.is_mid():
            return cls.mid2mid
        elif node1.self_identity.is_mid() and node2.self_identity.is_tail():
            return cls.mid2tail
        elif node1.self_identity.is_tail() and node2.self_identity.is_head():
            return cls.tail2head
        elif node1.self_identity.is_tail() and node2.self_identity.is_mid():
            return cls.tail2mid
        elif node1.self_identity.is_tail() and node2.self_identity.is_tail():
            return cls.tail2tail
        else:
            raise ValueError("Invalid node identity")


class MergeCondition:
    def __init__(self, threshold: int):
        self.threshold = threshold

    def head2head(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_two_heads_nodes_mode(
            node1, node2, self.threshold
        )

    def head2tail(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_head_and_tail_nodes_mode(
            node1, node2, self.threshold
        )

    def head2mid(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_head_and_middle_nodes_mode(
            node1, node2
        )

    def mid2head(self, node1: Node, node2: Node) -> bool:
        return self.head2mid(node2, node1)

    def mid2mid(self, node1: Node, node2: Node) -> bool:
        return (
            node1.exons[0][0] == node2.exons[0][0]  # type: ignore
            and node1.exons[-1][1] == node2.exons[-1][1]  # type: ignore
        )

    def mid2tail(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_tail_and_middle_nodes_mode(
            node2, node1
        )

    def tail2head(self, node1: Node, node2: Node) -> bool:
        return self.head2tail(node2, node1)

    def tail2mid(self, node1: Node, node2: Node) -> bool:
        return self.mid2tail(node2, node1)

    def tail2tail(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_two_tail_nodes_mode(
            node1, node2, self.threshold
        )

    def merged(self, node1: Node, node2: Node) -> bool:
        if node1.self_identity.is_head() and node2.self_identity.is_head():
            return self.head2head(node1, node2)
        elif node1.self_identity.is_head() and node2.self_identity.is_tail():
            return self.head2tail(node1, node2)
        elif node1.self_identity.is_head() and node2.self_identity.is_mid():
            return self.head2mid(node1, node2)
        elif node1.self_identity.is_mid() and node2.self_identity.is_head():
            return self.mid2head(node1, node2)
        elif node1.self_identity.is_mid() and node2.self_identity.is_mid():
            return self.mid2mid(node1, node2)
        elif node1.self_identity.is_mid() and node2.self_identity.is_tail():
            return self.mid2tail(node1, node2)
        elif node1.self_identity.is_tail() and node2.self_identity.is_head():
            return self.tail2head(node1, node2)
        elif node1.self_identity.is_tail() and node2.self_identity.is_mid():
            return self.tail2mid(node1, node2)
        elif node1.self_identity.is_tail() and node2.self_identity.is_tail():
            return self.tail2tail(node1, node2)
        else:
            raise ValueError("Invalid node identity")


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
    # NOTE: Do not compare breakporint here <06-08-23, Yangyang Li>
    assert node1.ref_end is not None
    assert node2.ref_end is not None

    if abs(node1.ref_end - node2.ref_end) > threshold:
        return False

    # no introns
    if not node1.introns and not node2.introns:
        return True

    introns_group = (
        zip_longest(node1.introns, node2.introns)
        if node1.strand == "-"
        else zip_longest(node1.introns[::-1], node2.introns[::-1])
    )

    # have introns
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
    node1: Node, node2: Node, threshold: float
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

    # assert node1.ref_start is not None and node2.ref_start is not None
    # if abs(node1.ref_start - node2.ref_start) > threshold:
    #     return False

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
