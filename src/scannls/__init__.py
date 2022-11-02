# !/usr/bin/env python
"""Init file for scannls package."""
__version__ = "0.0.1"
__PACKAGE_NAME__ = "scannls"

from rich.traceback import install


from ._class.basicRead import Read
from ._class.exception import ReadNotFoundError
from ._class.exception import ToolNotFoundError
from ._class.exception import BreakpointNotFoundError
from ._class.exception import ModesNotEqualError
from ._class.basicClass import Event
from ._class.basicClass import reverse_complement
from ._class.basicClass import Series, Node
from ._class.basicClass import NovelInsertion
from ._class.basicClass import MicroHomology
from ._class.basicClass import Insertion
from ._class.basicClass import check_end_node_is_ploya
from ._class.blat import Blat
from ._class.myLogger import MyLogger
from ._class.parallel import ParallelWorker
from ._class.readConnector import detect_read_read_connections_from_cigar
from ._class.readConnector import ReadsConnector
from .utils import get_softclip_length, external_tool_checking, get_longest_insertion_sequence, cigarstring2cigartuples
from ._class.cliqueFinder import CliqueFinder
from ._class.spliceGraph import SpliceGraph
from ._class.spliceGraph import SpliceType
from ._class.srRescuer import SRRescuer
from ._class.writer import FastaWriter
from ._class.writer import GTFWriter
from ._class.writer import VCFWriter
from ._class.writer import Writers
from ._class.type import LoggerType
from .cli.arg import DefaultOptions
from . import core
from . import cppext
from . import blat

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
    "FastaWriter",
    "GTFWriter",
    "VCFWriter",
    "Writers",
    "ToolNotFoundError",
    "ReadNotFoundError",
    "BreakpointNotFoundError",
    "ModesNotEqualError",
    "cli",
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
