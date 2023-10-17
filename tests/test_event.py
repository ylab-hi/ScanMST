# !/usr/bin/env python
"""Test for Event class.

@Filename:    test_event.py
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        2/2/22 6:43 PM
"""
import copy

import pytest
from scannls import ReadNotFoundError


class TestEvent:
    """Test for Event class."""

    def test_reverse(self, event):
        """Test reverse."""
        original_event = copy.deepcopy(event)
        event.reverse()
        assert event.bp1 == original_event.bp2
        assert event.mode1 == original_event.mode2
        assert event.strand1 == original_event.strand2
        assert event.read2_ref_start == original_event.read1_ref_start
        assert event.read1_exons == original_event.read2_exons

    def test_modes(self, event):
        """Test modes."""
        assert [event.mode1, event.mode2] == event.modes

    def test_chrom1(self, event):
        """Test chrom1."""
        assert event.chrom1 == "chr17"

    def test_insertion_seq1(self, event):
        """Test insertion seq."""
        assert event.insertion_seq1 == "GG"

    def test_source_s1_s2(self, event):
        """Test source s1 s2."""
        assert event.source_s1 == "left"
        assert event.source_s2 == "right"

    def test_is_type_na(self, event):
        """Test is type na."""
        assert not event.is_type_na()

    def test_has_insertion(self, event):
        """Test has insertion."""
        assert not event.has_insertion()

    def test_has_microhomology(self, event):
        """Test has microhomology."""
        assert event.has_microhomology()

    def test_is_same_strand(self, event):
        """Test is same strand."""
        assert event.is_same_strand()

    def test_get_read(self, event, reads):
        """Test get read."""
        assert event.read1(reads) == reads[0]
        assert event.read2(reads) == reads[1]
        with pytest.raises(ReadNotFoundError):
            event.read1([])

    def test_update_node_info(self, event, nodes, microhomology):
        """Test Update node info."""
        node1, _ = nodes
        event.update_node_info(flag=False, new_node=node1, insertion=microhomology)
        assert node1.sv_type == event.sv_type
        assert node1.annotation_code == event.annotation_code
        assert node1.splicing_code == event.splicing_code
        assert node1.modes == event.modes
        assert node1.genes == event.genes
        assert node1.insertion_info == (False, microhomology)

    def test_update_insertion_node_info(self, event, nodes):
        """Test update insertion node info."""
        node1, _ = nodes
        event.update_insertion_node_info(node1)
        assert node1.prev_breakpoint == "chr2:190659106"
        assert node1.annotation_code == event.annotation_code
