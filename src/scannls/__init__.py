# !/usr/bin/env python
"""Init file for scannls package."""
__version__ = "0.0.1"

try:
    import pysam  # type: ignore
    import numpy as np
    import HTSeq  # type: ignore
except ModuleNotFoundError as e:
    raise SystemExit from e

from ._class.basicRead import Read
from ._class.exception import ReadNotFoundError
from ._class.exception import ToolNotFoundError
from ._class.exception import BreakpointNotFoundError
from ._class.exception import ModesNotEqualError
from ._class.basicClass import Event
from ._class.basicClass import reverse_complement
from ._class.basicClass import Series, NodeType, Node
from ._class.blat import Blat
from ._class.myLogger import MyLogger
from ._class.parallel import ParallelWorker
from ._class.readConnector import detect_read_read_connections_from_cigar
from .utils import get_softclip_length, external_tool_checking
from ._class.cliqueFinder import CliqueFinder
from ._class.spliceGraph import SpliceGraph
from ._class.spliceGraph import SpliceType
from ._class.srRescuer import SRRescuer
from ._class.writer import FastaWriter
from ._class.writer import GTFWriter
from ._class.writer import VCFWriter
from ._class.writer import Writers
from .type import Options, LoggerType

__all__ = [
    "Read",
    "Event",
    "reverse_complement",
    "Series",
    "Node",
    "NodeType",
    "Blat",
    "MyLogger",
    "ParallelWorker",
    "detect_read_read_connections_from_cigar",
    "get_softclip_length",
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
    "Options",
    "LoggerType",
    "external_tool_checking",
]
