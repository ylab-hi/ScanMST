# !/usr/bin/env python
"""Init file for scannls package."""
__version__ = "0.0.1"
__PACKAGE_NAME__ = "scannls"

from rich.traceback import install

from . import blat, core, cppext
from ._class.basicClass import (
    Event,
    Insertion,
    MicroHomology,
    Node,
    NovelInsertion,
    Series,
    check_end_node_is_ploya,
    reverse_complement,
)
from ._class.basicRead import Read
from ._class.blat import Blat
from ._class.cliqueFinder import CliqueFinder
from ._class.exception import (
    BreakpointNotFoundError,
    ModesNotEqualError,
    ReadNotFoundError,
    ToolNotFoundError,
)
from ._class.myLogger import MyLogger
from ._class.parallel import ParallelWorker
from ._class.readConnector import (
    ReadsConnector,
    detect_read_read_connections_from_cigar,
)
from ._class.spliceGraph import SpliceGraph, SpliceType
from ._class.srRescuer import SRRescuer
from ._class.type import LoggerType
from ._class.writer import FastaWriter, GTFWriter, VCFWriter, Writers
from .cli.arg import DefaultOptions
from .utils import (
    cigarstring2cigartuples,
    external_tool_checking,
    get_longest_insertion_sequence,
    get_softclip_length,
)

__all__ = [
    "Read",
    "Event",
    "reverse_complement",
    "Series",
    "Node",
    "Blat",
    "MyLogger",
    "ParallelWorker",
    "ReadsConnector",
    "detect_read_read_connections_from_cigar",
    "get_softclip_length",
    "get_longest_insertion_sequence",
    "cigarstring2cigartuples",
    "SpliceGraph",
    "SpliceType",
    "SRRescuer",
    "CliqueFinder",
    "FastaWriter",
    "GTFWriter",
    "VCFWriter",
    "Writers",
    "ToolNotFoundError",
    "ReadNotFoundError",
    "BreakpointNotFoundError",
    "ModesNotEqualError",
    "LoggerType",
    "external_tool_checking",
    "NovelInsertion",
    "MicroHomology",
    "Insertion",
    "core",
    "cppext",
    "check_end_node_is_ploya",
    "DefaultOptions",
    "blat",
]


install(show_locals=True)
