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
from ._class.exception import ReadNotFoundError, ToolNotFoundError
from ._class.basicClass import Event
from ._class.basicClass import reverse_complement
from ._class.basicClass import Series
from ._class.blat import Blat
from ._class.myLogger import MyLogger
from ._class.parallel import ParallelWorker
from ._class.readConnector import detect_read_read_connections_from_cigar
from .utils import get_softclip_length

from ._class.spliceGraph import CliqueFinder
from ._class.spliceGraph import SpliceGraph
from ._class.srRescuer import SRRescuer
from ._class.writer import FastaWriter
from ._class.writer import GTFWriter

__all__ = [
    "Read",
    "Event",
    "reverse_complement",
    "Series",
    "Blat",
    "MyLogger",
    "ParallelWorker",
    "detect_read_read_connections_from_cigar",
    "get_softclip_length",
    "CliqueFinder",
    "SpliceGraph",
    "SRRescuer",
    "FastaWriter",
    "GTFWriter",
    "ToolNotFoundError",
    "ReadNotFoundError",
]
