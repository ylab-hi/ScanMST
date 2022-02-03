"""Test suite for the ScanNLS package."""
from typing import Mapping

from scannls import Node


def assign_value_for_instance(node: Node, **kwargs: Mapping[str, object]):
    """Assign value to node."""
    for key, value in kwargs.items():
        if key in node.__slots__:
            setattr(node, key, value)
