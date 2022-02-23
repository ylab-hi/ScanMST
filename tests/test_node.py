# !/usr/bin/env python
"""Test Node class.

@Filename:    test_node.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        2/1/22 9:50 PM
"""
import pytest

from scannls import Node


class FakePysamAlignmentFile:
    """Fake pysam.AlignmentFile."""

    @staticmethod
    def count(*args):
        """Fake count."""
        return sum(1 for _ in args)


@pytest.fixture(scope="module")
def bam():
    """Fake bam."""
    return FakePysamAlignmentFile()


class TestNode:
    """Test Node class."""

    @pytest.mark.parametrize("num", [1, 4, 5])
    def test_create_nodes(self, num):
        """Test create nodes."""
        results = Node.create_nodes(num)
        assert len(results) == num

    def test_introns(self, nodes):
        """Test intron."""
        node1, node2 = nodes
        assert node1.introns[0] == (190659995, 190670325)
        assert node2.introns == []

    def test_similar_key(self, nodes):
        """Test similar key."""
        node1, node2 = nodes
        assert node1.similar_key == "chr2-190659995-190670325"
        assert node2.similar_key == "chr17-None"

    def test_update_sr(self, nodes):
        """Test Update sr."""
        node1, node2 = nodes
        node1_original_sr = node1.sr
        node1.update_sr(node2.sr)
        assert node1.sr == node2.sr + node1_original_sr

    @pytest.mark.parametrize("sr", [1, 4, 5])
    def test_set_original_sr(self, nodes, sr):
        """Test set original sr."""
        node1, _ = nodes
        node1.set_original_sr(sr)
        assert node1.original_sr == sr

    @pytest.mark.parametrize("trace_id", [1, 4, 5])
    def test_set_trace_id(self, nodes, trace_id):
        """Test Set trace id."""
        node1, _ = nodes
        node1.is_traced = False
        node1.set_trace_id(trace_id)
        assert node1.trace_id == trace_id

    def test_reset_trace_id(self, nodes):
        """Test reset trace id."""
        node1, _ = nodes
        node1.is_traced = False
        node1.set_trace_id(1)
        assert node1.trace_id == 1
        node1.reset_trace_id()
        assert node1.trace_id == -1

    def test_is_start_node(self, nodes):
        """Test is start node."""
        node1, _ = nodes
        assert node1.is_start_node() is True
        node1.predecessors.append(1)
        assert node1.is_start_node() is False

    def test_is_end_node(self, nodes):
        """Test is end node."""
        node1, _ = nodes
        assert node1.is_end_node() is True
        node1.successors.append(1)
        assert node1.is_end_node() is False

    def test_has_predecessor(self, nodes):
        """Test has predecessor."""
        node1, _ = nodes
        assert node1.has_predecessor() is False
        node1.predecessors.append(1)
        assert node1.has_predecessor() is True

    def test_has_successor(self, nodes):
        """Test has successor."""
        node1, _ = nodes
        assert node1.has_successor() is False
        node1.successors.append(1)
        assert node1.has_successor() is True

    def test_add_successor_from_list(self, nodes):
        """Test add successor from list."""
        node1, node2 = nodes
        node2.is_in_graph = True
        node1.add_successor_from_list([node2])
        assert node1.successors == [node2]

    def test_add_predecessor_from_list(self, nodes):
        """Test add predecessor from list."""
        node1, node2 = nodes
        node2.add_predecessor_from_list([node1])
        assert node2.predecessors == []

    def test_add_successor(self, nodes):
        """Test add successor."""
        node1, node2 = nodes
        node1.add_successor(node2)
        assert node1.successors == []
        node2.is_in_graph = True
        node1.add_successor(node2)
        assert node1.successors == [node2]

    def test_update_next_and_previous_node_in_series(self, nodes):
        """Test update next and previous node in series."""
        node1, node2 = nodes
        node1.update_next_and_previous_node_in_series(0, nodes)
        assert node1.next_node_in_series == node2
        node2.update_next_and_previous_node_in_series(1, nodes)
        assert node2.previous_node_in_series == node1

        node1.clear_next_and_previous_node_in_series()
        assert node1.next_node_in_series is None
