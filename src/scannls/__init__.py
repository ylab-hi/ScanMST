"""Init file for scannls package."""
__version__ = "0.0.1"
__PACKAGE_NAME__ = "scannls"

from rich.traceback import install

from . import blat, cli, cppext, graph, type, utils
from .base import (
    Blat,
    BreakPoint,
    CircRNAFilter,
    Event,
    ExonFilter,
    Insertion,
    MicroHomology,
    MyLogger,
    NovelInsertion,
    ParallelWorker,
    Read,
    ReadsConnector,
    RTSwitchingFilter,
    detect_read_read_connections_from_cigar,
    reverse_complement,
)
from .exception import (
    BreakpointNotFoundError,
    ModesNotEqualError,
    ReadNotFoundError,
    ToolNotFoundError,
)
from .writer import FastaWriter, GTFWriter, VCFWriter, Writers

__all__ = [
    "CircRNAFilter",
    "ExonFilter",
    "RTSwitchingFilter",
    "utils",
    "cli",
    "BreakPoint",
    "type",
    "graph",
    "Read",
    "Event",
    "reverse_complement",
    "Blat",
    "MyLogger",
    "ParallelWorker",
    "ReadsConnector",
    "detect_read_read_connections_from_cigar",
    "FastaWriter",
    "GTFWriter",
    "VCFWriter",
    "Writers",
    "ToolNotFoundError",
    "ReadNotFoundError",
    "BreakpointNotFoundError",
    "ModesNotEqualError",
    "NovelInsertion",
    "MicroHomology",
    "Insertion",
    "cppext",
    "blat",
]


install(show_locals=True)
