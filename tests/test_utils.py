"""Test the utils.py module."""
import shutil

import pytest
from scannls import ToolNotFoundError
from scannls.utils import external_tool_checking, get_softclip_length


class FakeRead:
    """Fake Pysam Read."""

    def __init__(
        self,
        query_name: str,
        reference_name: str,
        reference_start: int,
        cigarstring: str,
        query_sequence: str,
    ) -> None:
        """Initialize fake logger."""
        self.query_name = query_name
        self.reference_name = reference_name
        self.reference_start = reference_start
        self.cigarstring = cigarstring
        self.query_sequence = query_sequence
        self.mapping_quality = 60

    @staticmethod
    def get_tag(*_):
        """Fake trace."""
        return 0

    @staticmethod
    def is_reverse(*_):
        """Fake debug."""
        return False


@pytest.mark.parametrize(
    ("read", "mode", "expected_result"),
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


def test_external_tool_checking(monkeypatch, fake_logger):
    """Test external_tool_checking func."""
    monkeypatch.setattr(shutil, "which", lambda x: "TEST")
    assert external_tool_checking(fake_logger) is None

    monkeypatch.setattr(shutil, "which", lambda x: None)
    with pytest.raises(ToolNotFoundError):
        external_tool_checking(fake_logger)
