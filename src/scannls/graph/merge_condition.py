"""Merge condition.

@Filename:    mergeCondition.py
@contact:     yangyang.li@northwestern.edu
@Time:        4/18/22 7:51 PM
"""
from __future__ import annotations

from enum import Enum, auto
from itertools import zip_longest
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .basic_graph import Node


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
    def from_node(cls, node1: Node, node2: Node) -> MergeConditionMode:
        node1_self_identity = node1.self_identity
        node2_self_identity = node2.self_identity
        if node2_self_identity is None or node1_self_identity is None:
            msg = f"{node1} or {node2}'s self_identity is None"
            raise ValueError(msg)

        if node1_self_identity.is_head() and node2_self_identity.is_head():
            return cls.head2head
        if node1_self_identity.is_head() and node2_self_identity.is_tail():
            return cls.head2tail
        if node1_self_identity.is_head() and node2_self_identity.is_mid():
            return cls.head2mid
        if node1_self_identity.is_mid() and node2_self_identity.is_head():
            return cls.mid2head
        if node1_self_identity.is_mid() and node2_self_identity.is_mid():
            return cls.mid2mid
        if node1_self_identity.is_mid() and node2_self_identity.is_tail():
            return cls.mid2tail
        if node1_self_identity.is_tail() and node2_self_identity.is_head():
            return cls.tail2head
        if node1_self_identity.is_tail() and node2_self_identity.is_mid():
            return cls.tail2mid
        if node1_self_identity.is_tail() and node2_self_identity.is_tail():
            return cls.tail2tail

        msg = "Invalid node identity"
        raise ValueError(msg)


class MergeCondition:
    def __init__(self, threshold: int) -> None:
        self.threshold = threshold

    def head2head(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_two_heads_nodes_mode(
            node1,
            node2,
            self.threshold,
        )

    def head2tail(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_head_and_tail_nodes_mode(
            node1,
            node2,
            self.threshold,
        )

    def head2mid(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_head_and_middle_nodes_mode(
            node1,
            node2,
        )

    def mid2head(self, node1: Node, node2: Node) -> bool:
        return self.head2mid(node2, node1)

    def mid2mid(self, node1: Node, node2: Node) -> bool:
        if node1.exons is None or node2.exons is None:
            raise ValueError

        return node1.exons.first.start == node2.exons.first.start and node1.exons.last.end == node2.exons.last.end

    def mid2tail(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_tail_and_middle_nodes_mode(
            node2,
            node1,
        )

    def tail2head(self, node1: Node, node2: Node) -> bool:
        return self.head2tail(node2, node1)

    def tail2mid(self, node1: Node, node2: Node) -> bool:
        return self.mid2tail(node2, node1)

    def tail2tail(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_two_tail_nodes_mode(
            node1,
            node2,
            self.threshold,
        )

    def merged(self, node1: Node, node2: Node) -> bool:
        node1_self_identity = node1.self_identity
        node2_self_identity = node2.self_identity
        if node2_self_identity is None or node1_self_identity is None:
            msg = f"{node1} or {node2}'s self_identity is None"
            raise ValueError(msg)

        # merge will not work for nodes on the different chroms
        if node1.chrom != node2.chrom:
            return False

        # merge will not work for nodes on the different strands
        if node1.strand != node2.strand:
            return False

        if node1_self_identity.is_head() and node2_self_identity.is_head():
            return self.head2head(node1, node2)
        if node1_self_identity.is_head() and node2_self_identity.is_tail():
            return self.head2tail(node1, node2)
        if node1_self_identity.is_head() and node2_self_identity.is_mid():
            return self.head2mid(node1, node2)
        if node1_self_identity.is_mid() and node2_self_identity.is_head():
            return self.mid2head(node1, node2)
        if node1_self_identity.is_mid() and node2_self_identity.is_mid():
            return self.mid2mid(node1, node2)
        if node1_self_identity.is_mid() and node2_self_identity.is_tail():
            return self.mid2tail(node1, node2)
        if node1_self_identity.is_tail() and node2_self_identity.is_head():
            return self.tail2head(node1, node2)
        if node1_self_identity.is_tail() and node2_self_identity.is_mid():
            return self.tail2mid(node1, node2)
        if node1_self_identity.is_tail() and node2_self_identity.is_tail():
            return self.tail2tail(node1, node2)
        msg = "Invalid node identity"
        raise ValueError(msg)


def _compare_is_merged_helper_check_condition_for_two_heads_nodes_mode(
    node1: Node,
    node2: Node,
    threshold: int,
) -> bool:
    """Check if both head nodes can be merged.

    :param threshold:

    .. note::
        nodes with different length may be merged. []: exon -: intron
        node1: []-[]-[]
        node2:    []-[]

    """
    # WARN:  compare break point, and edge still compare break point <06-08-23, Yangyang Li>
    if node1.strand.is_forward():
        if abs(node1.ref_end - node2.ref_end) > threshold:
            return False
    elif abs(node1.ref_start - node2.ref_start) > threshold:
        return False

    # no introns
    if not node1.introns and not node2.introns:
        return True

    node1_introns = [] if node1.introns is None else node1.introns
    node2_introns = [] if node2.introns is None else node2.introns

    introns_group = (
        zip_longest(node1_introns, node2_introns)  # type:ignore
        if node1.strand.is_reverse()
        else zip_longest(node1_introns[::-1], node2_introns[::-1])  # type:ignore
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
    if node1.introns != node2.introns:
        return False

    if node1.is_polya:
        return False

    if node1.strand.is_forward():
        if ret := node1.exons.first.join(node2.exons.last):
            overlap, union = ret
            return len(overlap) / len(union) >= threshold

    elif ret := node2.exons.first.join(node1.exons.last):
        overlap, union = ret
        return len(overlap) / len(union) >= threshold

    return False


def _compare_is_merged_helper_check_condition_for_two_tail_nodes_mode(
    node1: Node,
    node2: Node,
    threshold: float,
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
    if node1.introns != node2.introns:
        return False

    if node1.strand.is_forward():
        if abs(node1.ref_start - node2.ref_start) > threshold:
            return False
    elif abs(node1.ref_end - node2.ref_end) > threshold:
        return False

    # For nodes with polyA, a small difference in polyA positions is allowed.
    if node1.is_polya and node2.is_polya:
        return (
            abs(node1.exons.first.start - node2.exons.first.start) <= threshold
            and abs(node1.exons.last.end - node2.exons.last.end) <= threshold
        )

    if node1.is_polya and not node2.is_polya:
        if node1.strand.is_forward():
            return node1.contains(node2, same_left=True)
        return node1.contains(node2, same_right=True)

    if not node1.is_polya and node2.is_polya:
        if node1.strand.is_forward():
            return node2.contains(node1, same_left=True)
        return node2.contains(node1, same_right=True)

    if not node1.is_polya and not node2.is_polya:
        if node1.strand.is_forward():
            return abs(node1.exons.first.start - node2.exons.first.start) <= threshold
        return abs(node1.exons.last.end - node2.exons.last.end) <= threshold

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
    # limit all introns
    if node1.introns != node2.introns:
        return False

    # WARN: first exon start = ref start, last exon end = ref end <06-30-23, Yangyang Li>
    # we save same value in different variable in which it is diffficult to change them at same time

    if node1.strand.is_forward():
        return node2.contains(node1, same_right=True)
    return node2.contains(node1, same_left=True)


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
    if node1.introns != node2.introns:
        return False

    if node1.is_polya:
        return False

    if node1.strand.is_forward():
        return node2.contains(node1, same_left=True)

    return node2.contains(node1, same_right=True)
