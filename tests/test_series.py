# !/usr/bin/env python
"""Test for series.

@Filename:    test_series.py
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        2/2/22 6:42 PM
"""
import copy

import pytest
from scannls import Series


@pytest.fixture()
def empty_series(fake_logger, fake_blat):
    """Empty Series fixture."""
    return Series(blat=fake_blat, logger=fake_logger)


@pytest.fixture()
def series(nodes, fake_logger, fake_blat):
    """Series fixture."""
    series = Series(blat=fake_blat, logger=fake_logger)
    for node in nodes:
        series.add_node(node)
    return series


class TestSeries:
    """Test Series."""

    def test_add_node(self, empty_series, nodes):
        """Test Add node."""
        assert empty_series.nodes == []
        for node in nodes:
            empty_series.add_node(node)
        assert len(empty_series) == len(nodes)

    def test_is_all_type_del(self, series):
        """Test is all type del."""
        assert not series.is_all_type_del()

    def test_is_all_node_sr_higher_than_threshold(self, series):
        """Test is all node sr higher than threshold."""
        assert not series.is_all_node_sr_higher_than_threshold(threshold=5)
        assert series.is_all_node_sr_higher_than_threshold(threshold=2)

    def test_get_sr_sum_for_all_node(self, series, nodes):
        """Test get sr sum for all node."""
        node1, node2 = nodes
        assert series.get_sr_sum_for_all_node() == node1.sr + node2.sr

    def test_create_series_from_node_list(self, nodes, fake_logger):
        """Test create series from node list."""
        series_instance = Series.create_series_from_node_list(
            nodes, fake_logger, set(), False
        )
        assert len(series_instance) == len(nodes)

    def test_disable_blat_logger(self, series):
        """Test disable blat logger."""
        assert series.logger is not None
        series.disable_blat_logger()
        assert series.logger is None

    def test_unique_key(self, series):
        """Test unique key."""
        assert series.unique_key is not None

    def test_reorder_event(self, event):
        """Test reorder event."""
        original_event = copy.deepcopy(event)
        Series.reorder_event(event)
        assert event.mode1 == original_event.mode2

    def test_order_events_by_trancription_direction(self, event):
        """Test order events by trancription direction."""
        original_event = copy.deepcopy(event)
        result = Series.order_events_by_trancription_direction([event])
        assert len(result) == 1
        assert result[0].mode1 == original_event.mode2

    @pytest.mark.skip(reason="Not implemented")
    def test_init(self):
        """Test init.

        .. todo:: Add test for init.
        """
