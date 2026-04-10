"""Merge condition."""

from __future__ import annotations

from collections import Counter
from enum import Enum, auto
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pyfaidx

    from .basic_graph import Node


SINGLE_FASTA_INSTANCE = None


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


def is_polya(node: Node, genome_fasta: pyfaidx.Fasta, ratio: float = 0.7, length: int = 20) -> bool:
    """Check whether the tail node is bona fide polyA or internal priming event."""
    if node.ref_end is None or node.ref_start is None:
        msg = f"{node} has no start or end position"
        raise SystemExit(msg)

    if node.strand.is_forward():
        seq = genome_fasta[node.chrom][node.ref_end : node.ref_end + length].seq
    else:
        seq = genome_fasta[node.chrom][node.ref_start - length : node.ref_start].reverse.complement.seq

    counter: dict[str, int] = Counter(seq)
    return counter["A"] <= ratio * len(seq)


class MergeCondition:
    def __init__(self, threshold: int, fasta: pyfaidx.Fasta) -> None:
        """Prune threshold."""
        self.threshold = threshold

        global SINGLE_FASTA_INSTANCE
        if fasta is not None:
            SINGLE_FASTA_INSTANCE = fasta
        elif SINGLE_FASTA_INSTANCE is None:
            msg = "fasta is not provided"
            raise ValueError(msg)

        self.fasta = SINGLE_FASTA_INSTANCE

    def is_polya(self, node: Node, ratio: float = 0.7, length: int = 20) -> bool:
        return is_polya(node, self.fasta, ratio, length)

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
            self.fasta,
            self.threshold,
        )

    def head2mid(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_head_and_middle_nodes_mode(
            node1,
            node2,
            self.threshold,
        )

    def mid2head(self, node1: Node, node2: Node) -> bool:
        return self.head2mid(node2, node1)

    @staticmethod
    def mid2mid(node1: Node, node2: Node) -> bool:
        if node1.chrom != node2.chrom or node1.strand != node2.strand or node1.introns != node2.introns:
            return False

        if node1.exons is None or node2.exons is None:
            raise ValueError

        return node1.exons.first.start == node2.exons.first.start and node1.exons.last.end == node2.exons.last.end

    def mid2tail(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_tail_and_middle_nodes_mode(
            node2,
            node1,
            self.fasta,
            self.threshold,
        )

    def tail2head(self, node1: Node, node2: Node) -> bool:
        return self.head2tail(node2, node1)

    def tail2mid(self, node1: Node, node2: Node) -> bool:
        return self.mid2tail(node2, node1)

    def tail2tail(self, node1: Node, node2: Node) -> bool:
        return _compare_is_merged_helper_check_condition_for_two_tail_nodes_mode(
            node1,
            node2,
            self.fasta,
            self.threshold,
        )

    def merged(self, node1: Node, node2: Node) -> bool:
        node1_self_identity = node1.self_identity
        node2_self_identity = node2.self_identity
        if node2_self_identity is None or node1_self_identity is None:
            msg = f"{node1} or {node2}'s self_identity is None"
            raise ValueError(msg)

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
        node1: [ ]-[ ]-[ ]-[ ]
        node2:      []-[ ]-[ ]

        node1: [ ]-[ ]-[ ]-[ ]
        node2:          []-[ ]

        node1: [ ]-[ ]-[ ]-[ ]
        node2:              []
    """
    # WARN:  compare break point, and edge still compare break point <06-08-23, Yangyang Li>
    if node1.chrom != node2.chrom:
        return False

    # merge will not work for nodes on the different strands
    if node1.strand != node2.strand:
        return False

    # checking if breakpoints satisfy the threshold
    if node1.strand.is_forward():
        if abs(node1.ref_end - node2.ref_end) > threshold:
            return False

    elif abs(node1.ref_start - node2.ref_start) > threshold:
        return False

    # both have no introns
    if not node1.introns and not node2.introns:
        return True

    if not node1.introns and node2.introns:
        if node1.strand.is_forward():
            return node2.exons.last.start <= node1.ref_start
        return node2.exons.first.end >= node1.ref_end
    if node1.introns and not node2.introns:
        if node1.strand.is_forward():
            return node1.exons.last.start <= node2.ref_start
        return node1.exons.first.end >= node2.ref_end
    if node1.introns and node2.introns:
        if node1.strand.is_forward():
            return __intron_lists_containment_checker(node1, node2, None, reverse_strand=False)
        return __intron_lists_containment_checker(node1, node2, None, reverse_strand=True)

    return False


def __intron_lists_containment_checker(
    node1: Node,
    node2: Node,
    ref_node: Node | None = None,
    *,
    reverse_strand: bool = False,
    control_start_or_end_when_equal_length: str = "start",
) -> bool:
    """Checks if intron_list from node1 are fully contained within intron_list of node2
       Or intron_list of node2 are fully contained within intron_list of node1 consecutively.
       When the number of introns is the same, control start position or end position based on
       reference node.

    .. note::
        nodes with different length may be merged. []: exon -: intron
       list1: [ ]-[ ]-[ ]-[ ]
       list2:      []-[ ]-[ ]

       list2: [ ]-[ ]-[ ]-[ ]
       list1: [ ]-[ ]-[]
    """
    node1_introns = node1.introns
    node2_introns = node2.introns
    node1_exons = node1.exons
    node2_exons = node2.exons
    len1 = len(node1_introns)
    len2 = len(node2_introns)

    if len1 == len2:
        have_identical_introns = node1_introns == node2_introns
        if ref_node == node1:
            if control_start_or_end_when_equal_length == "start":
                return have_identical_introns and node1.exons.first.start <= node2.exons.first.start
            return have_identical_introns and node1.exons.last.end >= node2.exons.last.end
        if ref_node == node2:
            if control_start_or_end_when_equal_length == "start":
                return have_identical_introns and node2.exons.first.start <= node1.exons.first.start
            return have_identical_introns and node2.exons.last.end >= node1.exons.last.end
        return have_identical_introns

    if len1 > len2:
        full_list = node1_introns
        sub_list = node2_introns
        full_exons = node1_exons
        sub_exons = node2_exons
    else:
        full_list = node2_introns
        sub_list = node1_introns
        full_exons = node2_exons
        sub_exons = node1_exons

    sub_length = len(sub_list)
    if sub_length < 1:
        return False

    if reverse_strand:
        # For reverse strand, check from start
        return full_list[:sub_length] == sub_list[:] and full_exons[sub_length].end >= sub_exons.last.end
    # For forward strand, check from end
    return full_list[-sub_length:] == sub_list[:] and full_exons[-(sub_length + 1)].start <= sub_exons.first.start


def find_shared_interval_indices(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> tuple[list[int], list[int]]:
    """
    Finds the indices of shared intervals between two lists of intervals using set operations.

    Args:
      a: A list of tuples, where each tuple represents an interval (start, end).
      b: A list of tuples, where each tuple represents an interval (start, end).

    Returns:
      A tuple containing two lists:
        - The indices of shared intervals in list 'a'.
        - The indices of shared intervals in list 'b'.
    """
    # Convert to sets for faster lookup
    a_set = set(a)
    b_set = set(b)

    # Find shared intervals
    shared_intervals = a_set.intersection(b_set)

    # Build interval-to-index mappings for O(1) lookup
    a_index_map = {interval: idx for idx, interval in enumerate(a)}
    b_index_map = {interval: idx for idx, interval in enumerate(b)}

    shared_indices_a = [a_index_map[interval] for interval in shared_intervals]
    shared_indices_b = [b_index_map[interval] for interval in shared_intervals]

    return shared_indices_a, shared_indices_b


def is_consecutive_from_beginning(indices: list[int]) -> bool:
    """
    Checks if the given indices are consecutive and start from the beginning of a list.

    Args:
      indices: A list of indices.

    Returns:
      True if the indices are consecutive and start from 0, False otherwise.
    """

    # Check if the indices are consecutive
    return bool(indices) and indices[0] == 0 and all(indices[i] == indices[i - 1] + 1 for i in range(1, len(indices)))


def is_consecutive_from_end(indices: list[int], list_length: int) -> bool:
    """
    Checks if the given indices are consecutive and end at the end of a list.

    Args:
      indices: A list of indices.
      list_length: The length of the list.

    Returns:
      True if the indices are consecutive and end at the last index of the list,
      False otherwise.
    """
    return bool(indices) and indices[-1] == list_length - 1 and all(indices[i] == indices[i + 1] - 1 for i in range(len(indices) - 1))


def __intron_lists_sharing_checker(
    node1: Node,
    node2: Node,
    reverse_strand: bool = False,
) -> bool:
    """Checks if introns from node1 (head node) and introns of node2 (tail node)
       are shared and if they can be merged.

    .. note::
       nodes with different length may be merged. []: exon -: intron
       node2: [ ]-[ ]-[ ]-[ ]
       node1:              []-[ ]-[ ]

       node2: [ ]-[ ]-[ ]-[]
       node1:      []-[ ]-[ ]-[ ]
    """
    node1_introns = node1.introns
    node2_introns = node2.introns
    node1_exons = node1.exons
    node2_exons = node2.exons
    len1 = len(node1_introns)
    len2 = len(node2_introns)

    shared_node1_indices, shared_node2_indices = find_shared_interval_indices(node1_introns, node2_introns)

    # no shared introns
    if len(shared_node1_indices) == 0:
        if reverse_strand:
            return node1_exons.last.end <= node2_exons.first.end and node1_exons.last.start <= node2_exons.first.start
        return node2_exons.last.end <= node1_exons.first.end and node2_exons.last.start <= node1_exons.first.start
    # shared introns on the positive strand
    if is_consecutive_from_beginning(shared_node1_indices) and is_consecutive_from_end(shared_node2_indices, len2):
        return (
            node2_exons.last.end <= node1_exons[shared_node1_indices[-1] + 1].end
            and node1_exons.first.start >= node2_exons[shared_node2_indices[0]].start
        )

    # shared introns on the negative strand
    if reverse_strand and is_consecutive_from_beginning(shared_node2_indices) and is_consecutive_from_end(shared_node1_indices, len1):
        return (
            node1_exons.last.end <= node2_exons[shared_node2_indices[-1] + 1].end
            and node2_exons.first.start >= node1_exons[shared_node1_indices[0]].start
        )

    return False


def _compare_is_merged_helper_check_condition_for_head_and_tail_nodes_mode(
    node1: Node,
    node2: Node,
    fasta,
    threshold: int = 20,
) -> bool:
    """Check if node1 and node2 can be merged based on overlap info.

    node1 is head node, node2 is tail node Using threshold to
    check if they can be merged.

    :param node1:  node1
    :param node2:  node2
    :param threshold:  threshold for checking if two nodes are merged
    :return:  True if two nodes are merged, otherwise False

    .. note::

        -> [node1]               [node1] <-
            [node2] ->      <- [node2]
    """
    if node1.chrom != node2.chrom:
        return False

    # merge will not work for nodes on the different strands
    if node1.strand != node2.strand:
        return False

    # tail node must have an internal polyA
    if is_polya(node2, fasta):
        return False

    if node1.strand.is_forward():
        if node2.ref_start <= node1.ref_start < node2.ref_end <= node1.ref_end:
            overlap = node2.ref_end - node1.ref_start
            if overlap < threshold:
                return False
            if not node1.introns and not node2.introns:
                return True
            if node1.introns and not node2.introns:
                return node1.exons.first.end >= node2.ref_end
            if not node1.introns and node2.introns:
                return node2.exons.last.start <= node1.ref_start
            if node1.introns and node2.introns:
                return __intron_lists_sharing_checker(node1, node2, reverse_strand=False)

    elif node1.ref_start <= node2.ref_start < node1.ref_end <= node2.ref_end:
        overlap = node1.ref_end - node2.ref_start
        if overlap < threshold:
            return False
        if not node1.introns and not node2.introns:
            return True
        if node1.introns and not node2.introns:
            return node1.exons.last.start <= node2.ref_start
        if not node1.introns and node2.introns:
            return node2.exons.first.end >= node1.ref_end
        if node1.introns and node2.introns:
            return __intron_lists_sharing_checker(node1, node2, reverse_strand=True)

    return False


def __obtain_the_node_with_longer_span(node1, node2) -> tuple[Node, Node]:
    """Obtain the node with longer span."""
    node1_ref_span = node1.ref_end - node1.ref_start
    node2_ref_span = node2.ref_end - node2.ref_start
    if node1_ref_span >= node2_ref_span:
        return node1, node2
    return node2, node1


def _compare_is_merged_helper_check_condition_for_two_tail_nodes_mode(
    node1: Node,
    node2: Node,
    fasta,
    threshold: int,
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

    if node1.chrom != node2.chrom:
        return False

    # merge will not work for nodes on the different strands
    if node1.strand != node2.strand:
        return False

    if node1.strand.is_forward():
        if abs(node1.ref_start - node2.ref_start) > threshold:
            return False
    elif abs(node1.ref_end - node2.ref_end) > threshold:
        return False

    # For nodes with polyA, a small difference in polyA positions is allowed.
    # introns must be the same
    if is_polya(node1, fasta) and is_polya(node2, fasta):
        return (
            node1.introns == node2.introns
            and abs(node1.exons.first.start - node2.exons.first.start) <= threshold
            and abs(node1.exons.last.end - node2.exons.last.end) <= threshold
        )

    long_node, short_node = __obtain_the_node_with_longer_span(node1, node2)

    if is_polya(long_node, fasta) and not is_polya(short_node, fasta):
        if long_node.strand.is_forward():
            if not long_node.introns and not short_node.introns:
                return long_node.ref_end >= short_node.ref_end
            if long_node.introns and not short_node.introns:
                return long_node.exons.first.end >= short_node.ref_end
            if not long_node.introns and short_node.introns:
                return False
            if long_node.introns and short_node.introns:
                return __intron_lists_containment_checker(
                    long_node,
                    short_node,
                    long_node,
                    reverse_strand=False,
                    control_start_or_end_when_equal_length="end",
                )
        elif not long_node.introns and not short_node.introns:
            return long_node.ref_start <= short_node.ref_start
        elif long_node.introns and not short_node.introns:
            return long_node.exons.last.start <= short_node.ref_start
        elif not long_node.introns and short_node.introns:
            return False
        elif long_node.introns and short_node.introns:
            return __intron_lists_containment_checker(
                long_node,
                short_node,
                long_node,
                reverse_strand=True,
                control_start_or_end_when_equal_length="start",
            )

    if not is_polya(long_node, fasta) and is_polya(short_node, fasta):
        return False

    if not is_polya(long_node, fasta) and not is_polya(short_node, fasta):
        if long_node.strand.is_forward():
            if not long_node.introns and not short_node.introns:
                return long_node.ref_end >= short_node.ref_end
            if long_node.introns and not short_node.introns:
                return long_node.exons.first.end >= short_node.ref_end
            if not long_node.introns and short_node.introns:
                return False
            if long_node.introns and short_node.introns:
                return __intron_lists_containment_checker(
                    long_node,
                    short_node,
                    long_node,
                    reverse_strand=False,
                    control_start_or_end_when_equal_length="end",
                )
        elif not long_node.introns and not short_node.introns:
            return long_node.ref_start <= short_node.ref_start
        elif long_node.introns and not short_node.introns:
            return long_node.exons.last.start <= short_node.ref_start
        elif not long_node.introns and short_node.introns:
            return False
        elif long_node.introns and short_node.introns:
            return __intron_lists_containment_checker(
                long_node,
                short_node,
                long_node,
                reverse_strand=True,
                control_start_or_end_when_equal_length="start",
            )
    return False


def _compare_is_merged_helper_check_condition_for_head_and_middle_nodes_mode(
    node1: Node,
    node2: Node,
    threshold: int,
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
    if node1.chrom != node2.chrom:
        return False

    # merge will not work for nodes on the different strands
    if node1.strand != node2.strand:
        return False

    # checking if shared breakpoints satisfy the threshold
    if node1.strand.is_forward():
        if abs(node1.ref_end - node2.ref_end) > threshold:
            return False

    elif abs(node1.ref_start - node2.ref_start) > threshold:
        return False

    # both have no introns
    if not node1.introns and not node2.introns:
        if node1.strand.is_forward():
            return node1.ref_start >= node2.ref_start
        return node1.ref_end <= node2.ref_end
    # start node has no intron, while middle node has introns
    if not node1.introns and node2.introns:
        if node1.strand.is_forward():
            return node2.exons.last.start <= node1.ref_start
        return node2.exons.first.end >= node1.ref_end

    # start node has introns, while middle node no introns
    if node1.introns and not node2.introns:
        return False

    # both have introns
    if node1.introns and node2.introns:
        if len(node1.introns) > len(node2.introns):
            return False
        if node1.strand.is_forward():
            return __intron_lists_containment_checker(
                node1,
                node2,
                node2,
                reverse_strand=False,
                control_start_or_end_when_equal_length="start",
            )
        return __intron_lists_containment_checker(
            node1,
            node2,
            node2,
            reverse_strand=True,
            control_start_or_end_when_equal_length="end",
        )

    return False


def _compare_is_merged_helper_check_condition_for_tail_and_middle_nodes_mode(
    node1: Node,
    node2: Node,
    fasta,
    threshold: int,
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
    if node1.chrom != node2.chrom:
        return False

    # merge will not work for nodes on the different strands
    if node1.strand != node2.strand:
        return False

    if is_polya(node1, fasta):
        return False

    # checking if shared breakpoints satisfy the threshold
    if node1.strand.is_forward():
        if abs(node1.ref_start - node2.ref_start) > threshold:
            return False

    elif abs(node1.ref_end - node2.ref_end) > threshold:
        return False

    # both have no introns
    if not node1.introns and not node2.introns:
        if node1.strand.is_forward():
            return node1.ref_end <= node2.ref_end
        return node1.ref_start >= node2.ref_start
    # tail node has no intron, while middle node has introns
    if not node1.introns and node2.introns:
        if node1.strand.is_forward():
            return node2.exons.first.end >= node1.ref_end
        return node2.exons.last.start <= node1.ref_start

    # tail node has introns, while middle node no introns
    if node1.introns and not node2.introns:
        return False

    # both have introns
    if node1.introns and node2.introns:
        if len(node1.introns) > len(node2.introns):
            return False
        if node1.strand.is_forward():
            return __intron_lists_containment_checker(
                node1,
                node2,
                node2,
                reverse_strand=False,
                control_start_or_end_when_equal_length="end",
            )
        return __intron_lists_containment_checker(
            node1,
            node2,
            node2,
            reverse_strand=True,
            control_start_or_end_when_equal_length="start",
        )

    return False
