# !/usr/bin/env python
"""Test the splice graph.

@Filename:    test_splice_graph.py
@license:     MIT Licence
@Time:        1/19/22 1:16 PM
"""
import pytest
from loguru import logger

from . import add_edge_according_order
from . import assign_value_for_node
from scannls import Node
from scannls import SpliceGraph


@pytest.fixture(scope="module", autouse=True)
def graph_for_prun():
    """Create a splice graph for testing pruning.

    .. note::

        Total Series: 4
        Total Nodes: 13 (n1 ... n13)
        Splice Type: Bubble, diff start, diff end

         n1    n2     n3   n4    n5
        1,2---1,2,3--1,2--1,2,3--1,2
             /  |        /   |
           /     |     /      |
        3          3           3
        n6        n7         n8

        4 --- 4  --- 4 ---- 4 --- 4
        n9    n10   n11    n12   n13

         Values:  "chrom", "ref_start", "ref_end", "strand", "exons", "sv_type",
                   "modes", "query_name", "sr", "successors", "predecessors"
    """
    param_dict = [
        # n1
        {
            "chrom": "chr2",
            "ref_start": 190659106,
            "ref_end": 190670461,
            "strand": "+",
            "exons": [[190659106, 190659995], [190670325, 190670461]],
            "sv_type": "TRA",
            "modes": [1, 2],
            "query_name": "one,two",
            "sr": 2,
            "successors": [],
            "predecessors": [],
        },
        # n2
        {
            "chrom": "chr17",
            "ref_start": 49502062,
            "ref_end": 49502204,
            "strand": "+",
            "exons": [[49502062, 49502204]],
            "sv_type": "TDUP",
            "modes": [1, 2],
            "query_name": "one,three,two",
            "sr": 3,
            "successors": [],
            "predecessors": [],
        },
        # n3
        {
            "chrom": "chr17",
            "ref_start": 49458588,
            "ref_end": 49458738,
            "strand": "+",
            "exons": [[49458588, 49458738]],
            "sv_type": "TRA",
            "modes": [1, 2],
            "query_name": "one,two",
            "sr": 2,
            "successors": [],
            "predecessors": [],
        },
        # n4
        {
            "chrom": "chr22",
            "ref_start": 49827496,
            "ref_end": 49827873,
            "strand": "+",
            "exons": [[49827496, 49827873]],
            "sv_type": "INV",
            "modes": [1, 1],
            "query_name": "one,three,two",
            "sr": 3,
            "successors": [],
            "predecessors": [],
        },
        # n5
        {
            "chrom": "chr22",
            "ref_start": 49797704,
            "ref_end": 49798116,
            "strand": "-",
            "exons": [[49797704, 49798116]],
            "sv_type": None,
            "modes": None,
            "query_name": "one,two",
            "sr": 2,
            "successors": [],
            "predecessors": [],
        },
        # n6
        {
            "chrom": "chr2",
            "ref_start": 190659176,
            "ref_end": 190670465,
            "strand": "+",
            "exons": [[190659176, 190659995], [190670325, 190670465]],
            "sv_type": "TRA",
            "modes": [1, 2],
            "query_name": "three",
            "sr": 1,
            "successors": [],
            "predecessors": [],
        },
        # n7
        {
            "chrom": "chr17",
            "ref_start": 49458583,
            "ref_end": 49458741,
            "strand": "+",
            "exons": [[49458583, 49458741]],
            "sv_type": "TRA",
            "modes": [1, 2],
            "query_name": "three",
            "sr": 1,
            "successors": [],
            "predecessors": [],
        },
        # n8
        {
            "chrom": "chr22",
            "ref_start": 49797754,
            "ref_end": 49798122,
            "strand": "-",
            "exons": [[49797754, 49798122]],
            "sv_type": None,
            "modes": None,
            "query_name": "three",
            "sr": 1,
            "successors": [],
            "predecessors": [],
        },
        # n9-n13
        {
            "chrom": "chr2",
            "ref_start": 190659106,
            "ref_end": 190670469,
            "strand": "+",
            "exons": [[190659106, 190659995], [190670325, 190670469]],
            "sv_type": "TRA",
            "modes": [1, 2],
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
        },
        {
            "chrom": "chr17",
            "ref_start": 49502052,
            "ref_end": 49502214,
            "strand": "+",
            "exons": [[49502052, 49502214]],
            "sv_type": "TDUP",
            "modes": [1, 2],
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
        },
        {
            "chrom": "chr17",
            "ref_start": 49458578,
            "ref_end": 49458726,
            "strand": "+",
            "exons": [[49458578, 49458726]],
            "sv_type": "TRA",
            "modes": [1, 2],
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
        },
        {
            "chrom": "chr22",
            "ref_start": 49827486,
            "ref_end": 49827863,
            "strand": "+",
            "exons": [[49827486, 49827863]],
            "sv_type": "INV",
            "modes": [1, 1],
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
        },
        {
            "chrom": "chr22",
            "ref_start": 49797754,
            "ref_end": 49798126,
            "strand": "-",
            "exons": [[49797754, 49798126]],
            "sv_type": None,
            "modes": None,
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
        },
    ]
    num_nodes = 13
    nodes = [Node() for _ in range(num_nodes)]
    for ind, node in enumerate(nodes):
        assign_value_for_node(node, **param_dict[ind])

    # add successors and predecessors  12 edges
    add_edge_according_order(nodes, 1, 2)  # n1 -> n2
    add_edge_according_order(nodes, 2, 3)  # n2 -> n3
    add_edge_according_order(nodes, 3, 4)  # n3 -> n4
    add_edge_according_order(nodes, 4, 5)  # n4 -> n5
    add_edge_according_order(nodes, 6, 2)  # n6 -> n2
    add_edge_according_order(nodes, 2, 7)  # n2 -> n7
    add_edge_according_order(nodes, 7, 4)  # n7 -> n4
    add_edge_according_order(nodes, 4, 8)  # n4 -> n8
    add_edge_according_order(nodes, 9, 10)  # n9 -> n10
    add_edge_according_order(nodes, 10, 11)  # n10 -> n11
    add_edge_according_order(nodes, 11, 12)  # n11 -> n12
    add_edge_according_order(nodes, 12, 13)  # n12 -> n13

    graph = SpliceGraph(logger)
    graph.nodes = graph.dict_factory()
    for node in nodes:
        graph.add_node_with_similar_key(node)

    yield graph


def test_graph_data_before_prun(graph_for_prun, num_nodes=13):
    """Test the pruning function in splice graph."""
    number = sum(1 for _ in graph_for_prun)
    assert number == num_nodes
    start_nodes = graph_for_prun.get_start_nodes()
    expected_result = {"one,two": 2, "three": 1, "four": 1}
    for start_node in start_nodes:
        assert expected_result[start_node.query_name] == start_node.sr


def test_node_num_after_prune(graph_for_prun):
    """Test number of node after the pruning function in splice graph."""
    graph_for_prun.prune()
    number = sum(1 for _ in graph_for_prun)
    assert number == 10


def test_start_node_num_after_prune(graph_for_prun):
    """Test number of start node after the pruning function in splice graph."""
    graph_for_prun.prune()
    start_nodes = graph_for_prun.get_start_nodes()
    assert len(start_nodes) == 2


def test_start_node_sr_after_prune(graph_for_prun):
    """Test start node sr after prune."""
    graph_for_prun.prune()
    start_nodes = graph_for_prun.get_start_nodes()
    expected_result = {
        "one,two": 3,
        "four": 1,
    }
    for start_node in start_nodes:
        assert expected_result[start_node.query_name] == start_node.sr


def test_bubble_after_prun(graph_for_prun):
    """Test bubble after pruning."""
    graph_for_prun.prune()
    for node in graph_for_prun:
        if (
            node.query_name == "one,two"
            and not node.is_start_node()
            and not node.is_end_node()
        ):
            assert node.sr == 3


def test_end_node_sr_after_prune(graph_for_prun):
    """Test end node sr after prune."""
    graph_for_prun.prune()
    expected_result = {
        "one,two": 3,
        "four": 1,
    }
    for node in graph_for_prun:
        if node.is_end_node():
            assert expected_result[node.query_name] == node.sr
