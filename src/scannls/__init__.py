"""Init file for scannls package."""

__version__ = "0.0.1"
__PACKAGE_NAME__ = "scannls"

import sys

from rich.traceback import install

from . import aligner, blat, cli, cppext, graph, mtype, utils
from .base import (
    Blat,
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

MAX_RECURSION_LIMIT = 10000
if sys.getrecursionlimit() < MAX_RECURSION_LIMIT:
    sys.setrecursionlimit(MAX_RECURSION_LIMIT)

__all__ = [
    "Blat",
    "BreakPoint",
    "BreakpointNotFoundError",
    "CigarCode",
    "CircRNAFilter",
    "Event",
    "ExonFilter",
    "FastaWriter",
    "GTFWriter",
    "Insertion",
    "Interval",
    "Intervals",
    "MicroHomology",
    "ModesNotEqualError",
    "MyLogger",
    "NovelInsertion",
    "ParallelWorker",
    "RTSwitchingFilter",
    "Read",
    "ReadNotFoundError",
    "ReadsConnector",
    "Strand",
    "ToolNotFoundError",
    "VCFWriter",
    "Writers",
    "aligner",
    "blat",
    "cli",
    "cppext",
    "detect_read_read_connections_from_cigar",
    "graph",
    "graphvis",
    "mtype",
    "reverse_complement",
    "utils",
]


install(show_locals=True)
