"""TSG writer class."""

from __future__ import annotations

import typing
from functools import singledispatchmethod
from typing import IO, Any

from loguru import logger

from scanmst.graph import NLGraph
from scanmst.graph.graphvis import create_nxgraph

from .writer import Writer


class TSGWriter(Writer):
    """Writer for TSG files.

    .. note::
    """

    HEADER: typing.ClassVar = {
        "TSG": 1.0,
        "reference": "GRCh38",
        "PG": "scanmst",
    }

    def __init__(self, file_path: str) -> None:
        """Initialize GTFWriter object."""
        super().__init__(file_path)
        self.path_writer = False

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def open(self, mode: str = "w") -> IO:
        """Open file."""
        if self.is_opened:
            logger.warning(f"{self.__class__.__name__}: File is already opened.")
        self.io = self.file_path.open(mode)

        if hasattr(self, "write_header"):
            self.write_header()  # type: ignore
        return self.io

    def write_header(self) -> None:
        """Write header to file."""
        header = "\n".join([f"H\t{k}\t{v}" for k, v in self.HEADER.items()])
        self.write_line(header)

    def close(self) -> None:
        """Close file."""
        if self.is_opened:
            logger.trace(f"{self.__class__.__name__}: Closing file.")
            self.io.close()  # type: ignore
            self.io = None

    def write_line(self, line: str) -> None:
        """Write line to file."""
        if self.is_opened:
            self.io.write(line)  # type: ignore
        else:
            logger.warning(f"{self.__class__.__name__}: File is not opened.")

    @singledispatchmethod
    def write_data(self, data_object: Any, object_id: str) -> None:  # type: ignore
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: NLGraph, object_id: str) -> None:
        """Write Series to GTF file.

        :param data_object:
        :return:
        """
        if len(data_object.nodes) == 0:
            logger.warning(
                f"{self.__class__.__name__}: No nodes to write to file in Cluster {object_id}",
            )
        self.write_line(get_tsg_from_nlgraph(data_object, gid=object_id))


def get_tsg_from_nlgraph(nlgraph, gid=None, min_support_reads=1) -> str:
    nxgraph = create_nxgraph(nlgraph, min_support_reads=min_support_reads, possible_paths=nlgraph.possible_paths)
    result = []

    if gid:
        result.append(f"\nG\t{gid}")
    else:
        result.append(f"\nG\t{nlgraph.id}")

    for node in nxgraph.nodes(data=True):
        result.append(f"N\t{node[0]}\t{node[1]['chrom']}:{node[1]['strand']!s}:{node[1]['exons'][1:-1]!s}\t{node[1]['reads']}")

    # write edges
    for edge in nxgraph.edges(data=True):
        result.append(f"E\t{edge[2]['id']}\t{edge[0]}\t{edge[1]}\t{edge[2]['breakpoints']}")

    # write possible paths
    for path_id, path in nxgraph.graph["possible_paths"].items():
        path_str = "\t".join([f"{ele_id}+" for ele_id in path])
        result.append(f"P\t{path_id}\t{path_str}")

    # write node attributes sr
    for node in nxgraph.nodes(data=True):
        result.append(f"A\tN\t{node[0]}\tptc:i:{node[1]['ptc']}")
        result.append(f"A\tN\t{node[0]}\tptf:f:{node[1]['ptf']}")

    # write edge attributes sr
    for edge in nxgraph.edges(data=True):
        result.append(f"A\tE\t{edge[2]['id']}\tsr:i:{edge[2]['weight']}")

        # write insertion info if exists
        if "insertion_info" in edge[2]:
            insertion_info = edge[2]["insertion_info"]

            if insertion_info:
                insertion_type, insertion_seq = insertion_info.split("(")
                if insertion_type == "NovelInsertion":
                    seq = insertion_seq.strip(")").split(":")[0]
                    result.append(f"A\tE\t{edge[2]['id']}\tnovel_insertion:Z:{seq}")
                elif insertion_type == "MicroHomology":
                    seq = insertion_seq.strip(")")
                    result.append(f"A\tE\t{edge[2]['id']}\tmicrohomology:Z:{seq}")

    return "\n".join(result)
