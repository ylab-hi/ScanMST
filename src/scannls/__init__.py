"""Init file for scannls package."""
__version__ = "0.0.1"
__PACKAGE_NAME__ = "scannls"

from rich.traceback import install

from . import blat, cppext, graph, type
from .base.basic_class import (
    Event,
    Insertion,
    MicroHomology,
    NovelInsertion,
    reverse_complement,
)
from .base.basic_read import Read
from .base.blat import Blat
from .base.cluster import ClusterFinder
from .base.my_logger import MyLogger
from .base.parallel import ParallelWorker
from .base.read_connector import (
    ReadsConnector,
    detect_read_read_connections_from_cigar,
)
from .cli.arg import DefaultOptions
from .exception import (
    BreakpointNotFoundError,
    ModesNotEqualError,
    ReadNotFoundError,
    ToolNotFoundError,
)
from .utils import (
    cigarstring2cigartuples,
    external_tool_checking,
    get_longest_insertion_sequence,
    get_softclip_length,
)
from .writer import FastaWriter, GTFWriter, VCFWriter, Writers

__all__ = [
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
    "get_softclip_length",
    "get_longest_insertion_sequence",
    "cigarstring2cigartuples",
    "ClusterFinder",
    "FastaWriter",
    "GTFWriter",
    "VCFWriter",
    "Writers",
    "ToolNotFoundError",
    "ReadNotFoundError",
    "BreakpointNotFoundError",
    "ModesNotEqualError",
    "external_tool_checking",
    "NovelInsertion",
    "MicroHomology",
    "Insertion",
    "cppext",
    "DefaultOptions",
    "blat",
]


install(show_locals=True)
