"""Test suite for the ScanMST package."""
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from scanmst.base import Blat
from scanmst.graph import Node

__all__ = [
    "assign_value_for_instance",
    "FakeHsp",
    "FakeLogger",
    "FakeBlat",
]


def assign_value_for_instance(node: Node, **kwargs: Mapping[str, object]):
    """Override already-initialised attributes on a constructed Node.

    Nodes are built through ``Node(...)`` and only tweaked here; the keys must
    exist in ``__slots__`` or the assignment is silently dropped, so anything
    unknown is rejected loudly instead.
    """
    for key, value in kwargs.items():
        if key not in node.__slots__:
            msg = f"Node has no slot {key!r}; it cannot be set on a Node instance"
            raise AttributeError(msg)
        setattr(node, key, value)


class FakeLogger:
    """Fake logger."""

    def __init__(self, name: str = "fake_logger") -> None:
        """Initialize fake logger."""
        self.name = name

    def trace(self, *_):
        """Fake trace."""

    def debug(self, *_):
        """Fake debug."""

    def info(self, *_):
        """Fake info."""

    def warning(self, *_):
        """Fake warning."""

    def error(self, *_):
        """Fake error."""

    def critical(self, *_):
        """Fake critical."""

    def success(self, *_):
        """Fake success."""

    def complete(self, *_):
        """Fake complete."""


class FakeBlat(Blat):
    """Fake Blat class.

    Subclasses Blat without calling its __init__: ReadsConnector gates the
    realignment path on `isinstance(self.aligner, Blat)`, and the real
    constructor needs the gfServer executable, which is downloaded on demand.
    """

    def __init__(
        self,
        name: str = "fake_blat",
        query_return: str = "query_return",
        query_insertion_return: str = "query_insertion_return",
        psl2sam_return: tuple[Any, ...] = ("psl2sam_return",),
    ) -> None:
        """Init."""
        self.name = name
        self.query_return = query_return
        self.query_insertion_return = query_insertion_return
        self.psl2sam_return = psl2sam_return

    # These deliberately diverge from Blat's signatures -- the real
    # _query_insertion and psl2sam are staticmethods -- because the double just
    # hands back canned values.
    def query(self, *_, **__) -> str:  # type: ignore[override]
        """Query."""
        return self.query_return

    def _query_insertion(self, *_, **__) -> str:  # type: ignore[override]
        """Query insertion."""
        return self.query_insertion_return

    def psl2sam(self, *_, **__) -> tuple[Any, ...]:  # type: ignore[override]
        """PSL2SAM."""
        return self.psl2sam_return


@dataclass
class FakeHsp:
    """Fake hsp."""

    query_start: int
    query_end: int
    query_seq: str
    # (start, end) of the hit on the reference; only [1] is read.
    hit_range: tuple[int, int] = (1, 2)
