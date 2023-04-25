# !/usr/bin/env python
"""Init file for scannls package."""
__version__ = "0.0.1"
__PACKAGE_NAME__ = "scannls"

from rich.traceback import install

from . import blat, cppext
from .base.basicClass import (
    Event,
    Insertion,
    MicroHomology,
    Node,
    NovelInsertion,
    Series,
    check_end_node_is_ploya,
    reverse_complement,
)
from .base.basicRead import Read
from .base.blat import Blat
from .base.cluster import ClusterFinder
from .base.exception import (
    BreakpointNotFoundError,
    ModesNotEqualError,
    ReadNotFoundError,
    ToolNotFoundError,
)
from .base.myLogger import MyLogger
from .base.parallel import ParallelWorker
from .base.readConnector import (
    ReadsConnector,
    detect_read_read_connections_from_cigar,
)
from .base.spliceGraph import SpliceGraph, SpliceType
from .base.srRescuer import SRRescuer
from .base.type import LoggerType
from .base.writer import FastaWriter, GTFWriter, VCFWriter, Writers
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
    "ClusterFinder",
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
    "cppext",
    "check_end_node_is_ploya",
    "DefaultOptions",
    "blat",
]


install(show_locals=True)
