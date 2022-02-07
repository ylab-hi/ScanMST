# !/usr/bin/env python
"""Test the utils.py module."""
import subprocess

import pytest

from scannls.utils import get_softclip_length


class FakeRead:
    """Fake Pysam Read."""

    def __init__(
        self,
        query_name: str,
        reference_name: str,
        reference_start: int,
        cigarstring: str,
        query_sequence: str,
    ):
        """Initialize fake logger."""
        self.query_name = query_name
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.cigarstring = cigarstring
        self.query_sequence = query_sequence
        self.mapping_quality = 60

    def get_tag(self, *_):
        """Fake trace."""
        return 0

    def is_reverse(self, *_):
        """Fake debug."""
        return False


@pytest.mark.parametrize(
    "read, mode, expected_result",
    [
        (
            FakeRead(
                "cero",
                "chr20",
                391287,
                "20S20M10S",
                "TAAGACTAGGAATGGAGCAGTTCAGTCTAAAAAATATCACAGGTAACAGA",
            ),
            0,
            (20, "TAAGACTAGGAATGGAGCAG", 391287, 2),
        ),
        (
            FakeRead(
                "uno",
                "chr20",
                391287,
                "20S20M10S",
                "TAAGACTAGGAATGGAGCAGTTCAGTCTAAAAAATATCACAGGTAACAGA",
            ),
            1,
            (10, "AGGTAACAGA", 391307, 1),
        ),
        (
            FakeRead(
                "dos",
                "chr20",
                391287,
                "20S20M10S",
                "TAAGACTAGGAATGGAGCAGTTCAGTCTAAAAAATATCACAGGTAACAGA",
            ),
            2,
            (20, "TAAGACTAGGAATGGAGCAG", 391287, 2),
        ),
    ],
)
def test_get_softclip_length(read, mode, expected_result):
    """Test get_softclip_length func."""
    assert get_softclip_length(read, mode) == expected_result


@pytest.mark.parametrize(
    "tool, expected_result",
    [
        ("gfClient", "gfClient_output"),
        ("gfServer", "gfServer_output"),
    ],
)
def test_external_tool_checking(tool, expected_result, fake_process, fake_logger):
    """Test external_tool_checking func."""
    fake_process.register_subprocess([tool], stdout=(f"{tool}_output"))
    assert subprocess.getoutput(tool) == expected_result
