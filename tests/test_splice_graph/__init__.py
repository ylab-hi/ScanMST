# !/usr/bin/env python
"""Test the splice graph module."""
from collections.abc import Sequence

from scanmst.graph import Node


def add_edge_according_order(
    nodes: Sequence[Node],
    parent_order: int,
    child_order: int,
):
    """Link two nodes by their 1-based position in `nodes`.

    This wires the predecessor/successor lists directly rather than going
    through NLGraph.add_edge, so a topology can be described without also
    inventing EdgeData for every connection.
    """
    nodes[parent_order - 1].successors.append(nodes[child_order - 1])
    nodes[child_order - 1].predecessors.append(nodes[parent_order - 1])
