# !/usr/bin/env python
"""Test the splice graph module.

@Filename:    __init__.py.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        1/20/22 6:58 PM
"""
from typing import Any
from typing import Dict
from typing import Sequence

from scannls import NodeType


def assign_value_for_node(node: NodeType, **kwargs: Dict[str, Any]):
    """Assign value to node."""
    for key, value in kwargs.items():
        if key in node.__slots__:
            setattr(node, key, value)


def add_edge_according_order(
    nodes: Sequence[NodeType], parent_order: int, child_order: int
):
    """Add edge according order.

    .. note::
        order is 1-based.
    """
    nodes[parent_order - 1].successors.append(nodes[child_order - 1])
    nodes[child_order - 1].predecessors.append(nodes[parent_order - 1])
