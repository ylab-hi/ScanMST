"""Init file for scannls package."""
__version__ = "0.0.1"
__PACKAGE_NAME__ = "scannls"

from rich.traceback import install

from . import cli, cppext, graph, type, utils
from .base import (
    Aligner,
    BreakPoint,
    CigarCode,
    CircRNAFilter,
    Event,
    ExonFilter,
    Insertion,
    Interval,
    Intervals,
    MicroHomology,
    MyLogger,
    NovelInsertion,
    ParallelWorker,
    Read,
    ReadsConnector,
    RTSwitchingFilter,
    Strand,
    detect_read_read_connections_from_cigar,
    reverse_complement,
)
from .exception import (
    BreakpointNotFoundError,
    ModesNotEqualError,
    ReadNotFoundError,
    ToolNotFoundError,
)
from .graph import graphvis
from .writer import FastaWriter, GTFWriter, VCFWriter, Writers

__all__ = [
    "Aligner",
    "CigarCode",
    "Strand",
    "Interval",
    "Intervals",
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
    "graphvis",
]


install(show_locals=True)
