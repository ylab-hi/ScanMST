# !/usr/bin/env python
"""Test Node class."""
import pytest
from scanmst.graph import Node, NodeIdentity


class TestNode:
    """Test Node class."""

    def test_introns(self, nodes):
        """Two exons yield one intron; a single exon yields none."""
        node1, node2 = nodes
        assert [(i.start, i.end) for i in node1.introns] == [(190659995, 190670325)]
        assert list(node2.introns) == []

    def test_similar_key(self, nodes):
        """Similar key is chromosome-scoped; introns are deliberately excluded."""
        node1, node2 = nodes
        assert node1.similar_key == "chr2_None"
        assert node2.similar_key == "chr17_None"

    def test_unique_key(self, nodes):
        """Unique key pins chrom, span, strand and query name."""
        node1, _ = nodes
        assert node1.unique_key == "chr2-190659106-190670461-+-one,two"

    def test_self_identity(self, nodes):
        """Identity is recorded per query name."""
        node1, node2 = nodes
        assert node1.self_identity is NodeIdentity.HEAD
        assert node2.self_identity is NodeIdentity.TAIL
        assert node1.identity("one,two") is NodeIdentity.HEAD

    def test_exons_length(self, nodes):
        """Exon length sums the exons, ignoring the intron."""
        node1, node2 = nodes
        assert node1.exons_length == (190659995 - 190659106) + (190670461 - 190670325)
        assert node2.exons_length == 49502204 - 49502062

    def test_is_start_node(self, nodes):
        """A node with no predecessors is a start node."""
        node1, node2 = nodes
        assert node1.is_start_node() is True
        node1.predecessors.append(node2)
        assert node1.is_start_node() is False

    def test_is_end_node(self, nodes):
        """A node with no successors is an end node."""
        node1, node2 = nodes
        assert node1.is_end_node() is True
        node1.successors.append(node2)
        assert node1.is_end_node() is False

    def test_has_predecessor(self, nodes):
        """Has predecessor tracks the predecessors list."""
        node1, node2 = nodes
        assert node1.has_predecessor() is False
        node1.predecessors.append(node2)
        assert node1.has_predecessor() is True

    def test_has_successor(self, nodes):
        """Has successor tracks the successors list."""
        node1, node2 = nodes
        assert node1.has_successor() is False
        node1.successors.append(node2)
        assert node1.has_successor() is True

    def test_add_successor_requires_node_in_graph(self, nodes):
        """A successor is only linked once it is part of the graph."""
        node1, node2 = nodes
        node1.add_successor(node2)
        assert node1.successors == []
        node2.is_in_graph = True
        node1.add_successor(node2)
        assert node1.successors == [node2]

    def test_add_successor_from_list(self, nodes):
        """Adding from a list honours the same in-graph rule."""
        node1, node2 = nodes
        node2.is_in_graph = True
        node1.add_successor_from_list([node2], None, None)
        assert node1.successors == [node2]

    def test_add_predecessor_from_list(self, nodes):
        """A predecessor not yet in the graph is not linked."""
        node1, node2 = nodes
        node2.add_predecessor_from_list([node1], None, None)
        assert node2.predecessors == []

    def test_update_next_and_previous_node_in_nlpath(self, nodes, nlpath):
        """Neighbour links are set from the node's index within the path.

        The last node also resolves its incoming edge, so a real NLPath (not a
        bare list) is required here.
        """
        node1, node2 = nodes
        node1.update_next_and_previous_node_in_nlpath(0, nlpath)
        assert node1.next_node_in_nlpath == node2

        node2.update_next_and_previous_node_in_nlpath(1, nlpath)
        assert node2.previous_node_in_nlpath == node1
        assert node2.previous_edge_in_nlapth is nlpath.edges[
            f"{node1.unique_key}-{node2.unique_key}"
        ]

        node1.clear_next_and_previous_node_in_series()
        assert node1.next_node_in_nlpath is None

    @pytest.mark.parametrize(
        ("identity", "expected"),
        [(NodeIdentity.HEAD, "HEAD"), (NodeIdentity.TAIL, "TAIL"), (NodeIdentity.MID, "MID")],
    )
    def test_node_identity_from_str(self, identity, expected):
        """NodeIdentity round-trips through its string form."""
        assert NodeIdentity.from_str(expected) is identity

    def test_equality_is_identity(self, nodes):
        """Nodes compare by object identity, not by value."""
        node1, _ = nodes
        twin = Node(
            query_name=node1.query_name,
            chrom=node1.chrom,
            strand="+",
            ref_start=node1.ref_start,
            ref_end=node1.ref_end,
            identity=NodeIdentity.HEAD,
            exons=node1.exons,
        )
        assert twin.unique_key == node1.unique_key
        assert twin != node1
