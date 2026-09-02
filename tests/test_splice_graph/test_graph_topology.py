# !/usr/bin/env python
"""Test NLGraph's node/edge container behaviour.

This replaces test_prune_case1.py. That file drove `_trace(SpliceType...)`,
`_trace_backward`, `create_same_level_node_list`, `_rule_out` and the
`harmonic_mean_sr` machinery, none of which survive in NLGraph -- tracing is
now a single `trace()` entry point and support-read counts moved onto edges.

The 13-node topology below is carried over unchanged, because it is the part
worth keeping:

    n1    n2     n3   n4    n5
    1,2---1,2,3--1,2--1,2,3--1,2
         /  |        /   |
       /     |     /      |
    3          3           3
    n6        n7         n8

    4 --- 4  --- 4 ---- 4 --- 4
    n9    n10   n11    n12   n13
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pytest
from scanmst.base import Intervals
from scanmst.graph import NLGraph, Node, NodeIdentity

from . import add_edge_according_order

# chrom, ref_start, ref_end, strand, query_name, exons
NODE_SPECS: list[tuple[str, int, int, str, str, list[list[int] | tuple[int, int]]]] = [
    ("chr2", 190659106, 190670461, "+", "one,two", [[190659106, 190659995], [190670325, 190670461]]),
    ("chr17", 49502062, 49502204, "+", "one,three,two", [[49502062, 49502204]]),
    ("chr17", 49458588, 49458738, "+", "one,two", [[49458588, 49458738]]),
    ("chr22", 49827496, 49827873, "+", "one,three,two", [[49827496, 49827873]]),
    ("chr22", 49797704, 49798116, "-", "one,two", [[49797704, 49798116]]),
    ("chr2", 190659176, 190670465, "+", "three", [[190659176, 190659995], [190670325, 190670465]]),
    ("chr17", 49458583, 49458741, "+", "three", [[49458583, 49458741]]),
    ("chr22", 49797754, 49798122, "-", "three", [[49797754, 49798122]]),
    ("chr2", 190659106, 190670469, "+", "four", [[190659106, 190659995], [190670325, 190670469]]),
    ("chr17", 49502052, 49502214, "+", "four", [[49502052, 49502214]]),
    ("chr17", 49458578, 49458726, "+", "four", [[49458578, 49458726]]),
    ("chr22", 49827486, 49827863, "+", "four", [[49827486, 49827863]]),
    ("chr22", 49797754, 49798126, "-", "four", [[49797754, 49798126]]),
]

# 1-based (parent, child) pairs; 12 edges over the 13 nodes above.
EDGE_ORDERS = [
    (1, 2), (2, 3), (3, 4), (4, 5),
    (6, 2), (2, 7), (7, 4), (4, 8),
    (9, 10), (10, 11), (11, 12), (12, 13),
]


def _identity(index: int) -> NodeIdentity:
    """First node of a chain is the head, last is the tail."""
    if index == 0:
        return NodeIdentity.HEAD
    if index == len(NODE_SPECS) - 1:
        return NodeIdentity.TAIL
    return NodeIdentity.MID


@pytest.fixture()
def graph_nodes() -> list[Node]:
    """The 13 nodes, wired into the topology drawn in the module docstring."""
    nodes = [
        Node(
            query_name=query_name,
            chrom=chrom,
            strand=strand,
            ref_start=ref_start,
            ref_end=ref_end,
            identity=_identity(index),
            exons=Intervals.from_list(exons),
        )
        for index, (chrom, ref_start, ref_end, strand, query_name, exons) in enumerate(NODE_SPECS)
    ]
    for parent, child in EDGE_ORDERS:
        add_edge_according_order(nodes, parent, child)
    return nodes


@pytest.fixture()
def graph(graph_nodes, fake_logger, tmp_path) -> NLGraph:
    """An NLGraph holding the 13 nodes."""
    nlgraph = NLGraph(
        logger=fake_logger,
        rescuer=None,
        merge_threshold=10,
        support_reads=1,
        input_bam_path=Path("unused.bam"),
        output_dir=tmp_path,
        rescue_sr=False,
    )
    # NLGraph initialises its containers in __call__, not __init__, so a graph
    # built directly has to be given them.
    nlgraph.nodes = nlgraph.dict_factory()
    nlgraph.edges = defaultdict(list)
    for node in graph_nodes:
        nlgraph.add_node_with_similar_key(node)
    return nlgraph


class TestNLGraphContainer:
    """NLGraph as a node container."""

    def test_len_counts_every_node(self, graph):
        """All 13 nodes are reachable through the graph."""
        assert len(graph) == len(NODE_SPECS)

    def test_iteration_yields_every_node(self, graph, graph_nodes):
        """Iteration covers the same set of nodes that was added."""
        assert {node.unique_key for node in graph} == {node.unique_key for node in graph_nodes}

    def test_contains(self, graph, graph_nodes):
        """Membership is decided by unique key."""
        assert all(node in graph for node in graph_nodes)

    def test_nodes_are_bucketed_by_similar_key(self, graph, graph_nodes):
        """similar_key is chromosome-scoped, so nodes group per chromosome."""
        assert set(graph.nodes) == {"chr2_None", "chr17_None", "chr22_None"}
        assert len(graph.get_nodes_with_similar_key("chr2_None")) == 3
        assert len(graph.get_nodes_with_similar_key("chr17_None")) == 5
        assert len(graph.get_nodes_with_similar_key("chr22_None")) == 5

    def test_get_nodes_with_unknown_similar_key(self, graph):
        """An unseen chromosome yields no nodes rather than raising."""
        assert graph.get_nodes_with_similar_key("chrX_None") == []

    def test_get_node_with_unique_key(self, graph, graph_nodes):
        """A node can be recovered by its unique key."""
        target = graph_nodes[0]
        assert graph.get_node_with_unique_key(target.unique_key) is target
        assert graph.get_node_with_unique_key("no-such-key") is None

    def test_start_nodes(self, graph, graph_nodes):
        """n1, n6 and n9 have no predecessors in the fixture topology."""
        starts = {node.unique_key for node in graph.get_start_nodes()}
        expected = {graph_nodes[index].unique_key for index in (0, 5, 8)}
        assert starts == expected

    def test_end_nodes(self, graph, graph_nodes):
        """n5, n8 and n13 have no successors in the fixture topology."""
        ends = {node.unique_key for node in graph.get_end_nodes()}
        expected = {graph_nodes[index].unique_key for index in (4, 7, 12)}
        assert ends == expected

    def test_remove_node(self, graph, graph_nodes):
        """A removed node is no longer a member and no longer counted."""
        victim = graph_nodes[5]
        graph.remove_node(victim)
        assert victim not in graph
        assert len(graph) == len(NODE_SPECS) - 1


class TestNLGraphEdges:
    """NLGraph edge bookkeeping."""

    def test_add_and_find_edge(self, graph, graph_nodes, edge_data):
        """An added edge is retrievable between the same node pair."""
        node1, node2 = graph_nodes[0], graph_nodes[1]
        graph.add_edge(node1, node2, edge_data)
        assert graph.find_edges(node1, node2)

    def test_find_edges_without_any(self, graph, graph_nodes):
        """Looking up a pair with no edge is an error, not an empty result."""
        with pytest.raises(KeyError):
            graph.find_edges(graph_nodes[0], graph_nodes[12])
