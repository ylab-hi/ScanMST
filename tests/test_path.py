# !/usr/bin/env python
"""Test the NLPath class.

NLPath replaced the old ``Path``/``Series`` class: it is constructed from a
list of nodes rather than a blat/logger pair, and support-read counts now live
on the edges instead of the nodes.
"""
import copy

import pytest
from scanmst.graph import NLPath, VariationType


class TestNLPath:
    """Test NLPath."""

    def test_len_and_iteration(self, nlpath, nodes):
        """A path is a sequence of its nodes."""
        assert len(nlpath) == len(nodes)
        assert list(nlpath) == list(nodes)
        assert nlpath[0] is nodes[0]

    def test_unique_key(self, nlpath, nodes):
        """Unique key concatenates the node keys in order."""
        assert nlpath.unique_key == "".join(node.unique_key for node in nodes)

    def test_id_is_stable(self, nlpath):
        """The generated id is deterministic and namespaced."""
        assert nlpath.id.startswith("TSP")
        assert nlpath.id == NLPath(list(nlpath.nodes)).id

    def test_sum_sr(self, nlpath, edge):
        """Support reads are summed across edges, not nodes."""
        assert nlpath.sum_sr() == edge.sr

    def test_is_all_type_del(self, nlpath, edge):
        """The fixture edge is a TDUP, so the path is not all-DEL."""
        assert nlpath.is_all_type_del() is False
        # Edge.variation_type is read-only; the value lives on the EdgeData.
        edge.edge_data.variantion_type = VariationType.DEL
        assert nlpath.is_all_type_del() is True

    def test_is_minimum_node_length_larger_than_threshold(self, nlpath):
        """Shortest node in the fixture spans 142 bases."""
        assert nlpath.is_minimum_node_length_larger_than_threshold(threshold=10) is True
        assert nlpath.is_minimum_node_length_larger_than_threshold(threshold=100_000) is False

    def test_create_node_signature(self, nodes):
        """Node signature covers chromosome, exons and strand."""
        signature = NLPath.create_node_signature(nodes[1])
        assert signature.startswith("chr17:")
        assert signature.endswith(";+")

    def test_setup_breakpoints(self, nlpath, nodes):
        """Breakpoints are derived from identity and strand."""
        nlpath.setup_breakpoints()
        # node1 is a forward-strand HEAD, so its breakpoint is its ref_end.
        assert nodes[0].breakpoints[nodes[0].ref_end] == 1
        # node2 is a forward-strand TAIL, so its breakpoint is its ref_start.
        assert nodes[1].breakpoints[nodes[1].ref_start] == 1

    def test_getitem_out_of_range(self, nlpath):
        """Indexing past the end raises, as for any sequence."""
        with pytest.raises(IndexError):
            nlpath[len(nlpath)]

    def test_reorder_event(self, event):
        """Reordering swaps the two breakpoint modes."""
        original_event = copy.deepcopy(event)
        NLPath.reorder_event(event)
        assert event.mode1 == original_event.mode2

    def test_order_events_by_trancription_direction(self, event):
        """Ordering a single event applies the same swap."""
        original_event = copy.deepcopy(event)
        result = NLPath.order_events_by_trancription_direction([event])
        assert len(result) == 1
        assert result[0].mode1 == original_event.mode2
