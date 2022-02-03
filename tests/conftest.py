# !/usr/bin/env python
"""Conftest for pytest."""
from typing import List

import pytest

from . import assign_value_for_node
from scannls import Node


@pytest.fixture(scope="function")
def nodes() -> List[Node]:
    """Return a list of nodes."""
    param_dict = [
        # n1
        {
            "chrom": "chr2",
            "ref_start": 190659106,
            "ref_end": 190670461,
            "strand": "+",
            "exons": [[190659106, 190659995], [190670325, 190670461]],
            "sv_type": "TRA",
            "prev_sv_type": None,
            "modes": [1, 2],
            "query_name": "one,two",
            "sr": 2,
            "successors": [],
            "predecessors": [],
            "trace_id": -1,
            "is_traced": False,
            "prev_breakpoint": None,
            "next_breakpoint": "chr2:190670461",
        },
        # n2
        {
            "chrom": "chr17",
            "ref_start": 49502062,
            "ref_end": 49502204,
            "strand": "+",
            "exons": [[49502062, 49502204]],
            "sv_type": "TDUP",
            "prev_sv_type": "TRA",
            "modes": [1, 2],
            "query_name": "one,three,two",
            "sr": 3,
            "successors": [],
            "predecessors": [],
            "trace_id": -1,
            "is_traced": False,
            "prev_breakpoint": "chr2:190670461",
            "next_breakpoint": "chr17:49502204",
        },
    ]
    node_list = [Node() for _ in range(len(param_dict))]
    for ind, node in enumerate(node_list):
        assign_value_for_node(node, **param_dict[ind])  # type: ignore
        node.get_unique_key()
    return node_list
