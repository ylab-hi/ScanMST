"""Test suite for the ScanNLS package."""
from typing import Any
from typing import Dict

from scannls import Node


def assign_value_for_node(node: Node, **kwargs: Dict[str, Any]):
    """Assign value to node."""
    for key, value in kwargs.items():
        if key in node.__slots__:
            setattr(node, key, value)
