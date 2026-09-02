# !/usr/bin/env python
"""Test for Event class.
"""
import copy

import pytest
from scanmst import ReadNotFoundError


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
        assert event.read1(reads, 0) == reads[0]
        assert event.read2(reads, 0) == reads[1]
        with pytest.raises(ReadNotFoundError):
            event.read1([], 0)

    def test_update_specific_info_within_event(self, event):
        """Selected event attributes are copied onto the target by name.

        This replaces the old update_node_info/update_insertion_node_info
        pair. Node uses __slots__ and none of the event's info keys are Node
        slots, so the copy is exercised against a plain target object -- which
        is all the method itself promises.
        """

        class Target:
            """Minimal stand-in for the object being annotated."""

        target = Target()
        keys = ["sv_type", "annotation_code", "splicing_code", "genes"]
        result = event.update_specific_info_within_event(target, keys)

        assert result is target
        for key in keys:
            assert getattr(target, key) == getattr(event, key)

    def test_update_specific_info_within_event_rejects_unknown_key(self, event):
        """An attribute the event does not carry is an error, not a silent skip."""

        class Target:
            """Minimal stand-in for the object being annotated."""

        with pytest.raises(AttributeError):
            event.update_specific_info_within_event(Target(), ["not_an_event_field"])
