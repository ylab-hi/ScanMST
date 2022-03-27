# !/usr/bin/env python
"""Type of the scannls.

@Filename:    type.py
@license:     MIT Licence
@Time:        12/30/21 2:20 PM
"""
from dataclasses import dataclass
from typing import Any
from typing import List
from typing import NewType
from typing import Protocol
from typing import Tuple

from .basicRead import Read

ReadType = NewType("ReadType", Read)

EventType = Tuple[
    str,
    int,
    int,
    Tuple[str, str, int, int],
    Tuple[int, int, Any],
    Tuple[int, int, Any],
    Tuple[str, str],
    Tuple[str, str],
    List[str],
]


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


@dataclass
class Options:
    """Cli options for testing."""

    input: str
    ref: str
    gtf: str
    output: str
    two_bit: str
    support_reads: int = 1
    splice_bin: int = 5
    mapq: int = 15
    noncanonical: bool = False
    closed: bool = True
    nsleep: bool = False
    log: str = "info"
    parallel: int = 1
    port: int = 88888
    min_soft_seg_len: int = 200
    max_allowed_nm: int = 60
    ident_cutoff: float = 0.99
    tmp_dir: str = "/tmp"
    soft_len: int = 5
    mismatch: int = 3
    alignment_fraction: float = 0.8
    long_indel_length: int = 5
    substitutions_num: int = 3
    substitutions_fraction: float = 0.2
    indel_fraction: float = 0.2
