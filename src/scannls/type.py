# !/usr/bin/env python
"""Type of the scannls.

@Filename:    type.py
@license:     MIT Licence
@Time:        12/30/21 2:20 PM
"""
from typing import Any
from typing import NewType
from typing import Protocol

from . import Read


class LoggerType(Protocol):
    """Logger type."""

    def trace(self, msg: str) -> None:
        """Trace."""

    def debug(self, msg: str) -> None:
        """Debug."""

    def info(self, msg: str) -> None:
        """Info."""

    def warning(self, msg: str) -> None:
        """Warning."""

    def error(self, msg: str) -> None:
        """Error."""

    def critical(self, msr: str) -> None:
        """Critical."""

    def success(self, msg: str) -> None:
        """Success."""

    def complete(self) -> Any:
        """Complete."""


ReadType = NewType("ReadType", Read)
