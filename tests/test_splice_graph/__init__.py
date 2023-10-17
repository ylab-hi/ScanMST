# !/usr/bin/env python
"""Test the splice graph module.

@Filename:    __init__.py.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        1/20/22 6:58 PM
"""
from collections.abc import Sequence

from scannls import Node


def add_edge_according_order(
    nodes: Sequence[Node],
    parent_order: int,
    child_order: int,
):
    """Add edge according order.

    .. note::
        order is 1-based.
    """
    nodes[parent_order - 1].successors.append(nodes[child_order - 1])
    nodes[child_order - 1].predecessors.append(nodes[parent_order - 1])
