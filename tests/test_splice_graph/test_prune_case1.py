# !/usr/bin/env python
"""Test the splice graph.

@Filename:    test_splice_graph.py
@license:     MIT Licence
@Time:        1/19/22 1:16 PM
"""
import pytest
from loguru import logger
from scannls import NLGraph
from scannls import Node
from scannls import SpliceType

from . import add_edge_according_order
from .. import assign_value_for_instance


@pytest.fixture(scope="function", autouse=True)
def graph_for_prun():
    """Create a splice graph for testing pruning.

    .. note::

        Total Series: 4
        Total Nodes: 13 (n1 ... n13)
        Splice Type: Bubble, diff start, diff end

        After construction:
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

        After Prune:

        1,2---1,2,3--1,2--1,2,3--1,2
                            \
                              4

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
            "prev_sv_type": None,
            "modes": [1, 2],
            "query_name": "one,two",
            "sr": 2,
            "successors": [],
            "predecessors": [],
            "trace_id": -1,
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
            "prev_breakpoint": "chr2:190670461",
            "next_breakpoint": "chr17:49502204",
        },
        # n3
        {
            "chrom": "chr17",
            "ref_start": 49458588,
            "ref_end": 49458738,
            "strand": "+",
            "exons": [[49458588, 49458738]],
            "sv_type": "TRA",
            "prev_sv_type": "TDUP",
            "modes": [1, 2],
            "query_name": "one,two",
            "sr": 2,
            "successors": [],
            "predecessors": [],
            "trace_id": -1,
            "prev_breakpoint": "chr17:49458588",
            "next_breakpoint": "chr17:49458738",
        },
        # n4
        {
            "chrom": "chr22",
            "ref_start": 49827496,
            "ref_end": 49827873,
            "strand": "+",
            "exons": [[49827496, 49827873]],
            "sv_type": "INV",
            "prev_sv_type": "TRA",
            "modes": [1, 1],
            "query_name": "one,three,two",
            "sr": 3,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": "chr22:49827496",
            "next_breakpoint": "chr22:49827873",
        },
        # n5
        {
            "chrom": "chr22",
            "ref_start": 49797704,
            "ref_end": 49798116,
            "strand": "-",
            "exons": [[49797704, 49798116]],
            "sv_type": None,
            "prev_sv_type": "INV",
            "modes": None,
            "query_name": "one,two",
            "sr": 2,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": "chr22:49798116",
            "next_breakpoint": None,
        },
        # n6
        {
            "chrom": "chr2",
            "ref_start": 190659176,
            "ref_end": 190670465,
            "strand": "+",
            "exons": [[190659176, 190659995], [190670325, 190670465]],
            "sv_type": "TRA",
            "prev_sv_type": None,
            "modes": [1, 2],
            "query_name": "three",
            "sr": 1,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": None,
            "next_breakpoint": "chr2:190670465",
        },
        # n7
        {
            "chrom": "chr17",
            "ref_start": 49458583,
            "ref_end": 49458741,
            "strand": "+",
            "exons": [[49458583, 49458741]],
            "sv_type": "TRA",
            "prev_sv_type": "TDUP",
            "modes": [1, 2],
            "query_name": "three",
            "sr": 1,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": "chr17:49458583",
            "next_breakpoint": "chr17:49458741",
        },
        # n8
        {
            "chrom": "chr22",
            "ref_start": 49797754,
            "ref_end": 49798122,
            "strand": "-",
            "exons": [[49797754, 49798122]],
            "sv_type": None,
            "prev_sv_type": "INV",
            "modes": None,
            "query_name": "three",
            "sr": 1,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": "chr22:49798122",
            "next_breakpoint": None,
        },
        # n9-n13
        {
            "chrom": "chr2",
            "ref_start": 190659106,
            "ref_end": 190670469,
            "strand": "+",
            "exons": [[190659106, 190659995], [190670325, 190670469]],
            "sv_type": "TRA",
            "prev_sv_type": None,
            "modes": [1, 2],
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": None,
            "next_breakpoint": "chr2:190670469",
        },
        {
            "chrom": "chr17",
            "ref_start": 49502052,
            "ref_end": 49502214,
            "strand": "+",
            "exons": [[49502052, 49502214]],
            "sv_type": "TDUP",
            "prev_sv_type": "TRA",
            "modes": [1, 2],
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": "chr17:49502052",
            "next_breakpoint": "chr17:49502214",
        },
        {
            "chrom": "chr17",
            "ref_start": 49458578,
            "ref_end": 49458726,
            "strand": "+",
            "exons": [[49458578, 49458726]],
            "sv_type": "TRA",
            "prev_sv_type": "TDUP",
            "modes": [1, 2],
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": "chr17:49458578",
            "next_breakpoint": "chr17:49458726",
        },
        {
            "chrom": "chr22",
            "ref_start": 49827486,
            "ref_end": 49827863,
            "strand": "+",
            "exons": [[49827486, 49827863]],
            "sv_type": "INV",
            "prev_sv_type": "TRA",
            "modes": [1, 1],
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": "chr22:49827486",
            "next_breakpoint": "chr22:49827863",
        },
        {
            "chrom": "chr22",
            "ref_start": 49797754,
            "ref_end": 49798126,
            "strand": "-",
            "exons": [[49797754, 49798126]],
            "sv_type": None,
            "prev_sv_type": "INV",
            "modes": None,
            "query_name": "four",
            "sr": 1,
            "successors": [],
            "predecessors": [],
            "prev_breakpoint": "chr22:49798126",
            "next_breakpoint": None,
        },
    ]
    num_nodes = 13
    nodes = [Node() for _ in range(num_nodes)]
    for ind, node in enumerate(nodes):
        assign_value_for_instance(node, **param_dict[ind])  # type: ignore
        node.get_unique_key()

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

    graph = NLGraph(logger, None)
    graph.nodes = graph.dict_factory()
    for node in nodes:
        graph.add_node_with_similar_key(node)

    yield graph


@pytest.fixture(autouse=True)
def start_node_with_name_onetwo(graph_for_prun):
    """Start Node with name one two.

    :param graph_for_prun: graph for pruning
    :return: two nodes with name onetwo
    """
    for start_node in graph_for_prun.get_start_nodes():
        if start_node.query_name == "one,two":
            return start_node
    return None


@pytest.fixture(autouse=True)
def end_node_with_name_onetwo(graph_for_prun):
    """End Node with name one two."""
    for end_node in graph_for_prun.get_end_nodes():
        if end_node.query_name == "one,two":
            return end_node
    return None


@pytest.fixture(autouse=True)
def can_battle_nodes_onetwo_and_three(graph_for_prun):
    """Can battle nodes onetwo and three."""
    node_one_two, node_three = None, None
    for node in graph_for_prun.get_start_nodes():
        if node.query_name == "one,two":
            node_one_two = node
        if node.query_name == "three":
            node_three = node
    return node_one_two, node_three


def test_graph_data_before_prun(graph_for_prun, num_nodes=13):
    """Test the pruning function in splice graph."""
    number = sum(1 for _ in graph_for_prun)
    assert number == num_nodes
    start_nodes = graph_for_prun.get_start_nodes()
    expected_result = {"one,two": 2, "three": 1, "four": 1}
    for start_node in start_nodes:
        assert expected_result[start_node.query_name] == start_node.sr


def test_check_can_battle(graph_for_prun, can_battle_nodes_onetwo_and_three):
    """Test can battle."""
    one_node, other_node = can_battle_nodes_onetwo_and_three
    if one_node is not None or other_node is not None:
        assert graph_for_prun.check_can_battle(one_node, other_node)


def test_trace_forward(graph_for_prun, start_node_with_name_onetwo):
    """Test trace forward."""
    if start_node_with_name_onetwo is not None:
        start_node_with_name_onetwo.set_trace_id(1)
        paths = []
        graph_for_prun._trace_forward(start_node_with_name_onetwo, 2, [], paths)
        assert len(paths) == 4
        assert len(paths[0]) == 5


def test_trace_backward(graph_for_prun, end_node_with_name_onetwo):
    """Test trace backward."""
    if end_node_with_name_onetwo is not None:
        end_node_with_name_onetwo.set_trace_id(1)
        paths = []
        graph_for_prun._trace_backward(end_node_with_name_onetwo, 2, [], paths)
        assert len(paths) == 4
        assert len(paths[0]) == 5


def test_trace(graph_for_prun):
    """Test trace."""
    graph_for_prun._trace(SpliceType.forward)
    graph_for_prun._trace(SpliceType.backward)

    with pytest.raises(ValueError):
        graph_for_prun._trace(1)  # type: ignore


def test_create_same_level_node_list(graph_for_prun):
    """Test creat same level node list."""
    graph_for_prun._trace(SpliceType.forward)
    assert len(graph_for_prun.create_same_level_node_list()) == 5


def test_rule_out(graph_for_prun, can_battle_nodes_onetwo_and_three):
    """Test rule out."""
    node_one_two, node_three = can_battle_nodes_onetwo_and_three
    graph_for_prun._rule_out(node_one_two, node_three)
    assert node_one_two.sr == 3
    assert node_three not in graph_for_prun


def test_set_harmoic_mean_sr(graph_for_prun):
    """Test set harmonic mean sr."""
    for start_node in graph_for_prun.get_start_nodes():
        start_node.set_harmoic_mean_sr(start_node.sr)
        graph_for_prun._trace_forward(start_node, 2, [], [])

    for start_node in graph_for_prun.get_start_nodes():
        if start_node.query_name == "one,two":
            assert start_node.harmonic_mean_sr == 2
            assert start_node.successors[0].harmonic_mean_sr == (18 / 11)
        else:
            assert start_node.sr == 1
