# !/usr/bin/env python
"""FastaWriter class.

@Filename:    fastaWriter.py
@license:     MIT Licence
@Time:        1/30/22 6:18 PM
"""
from functools import singledispatchmethod
from pathlib import Path
from typing import Any
from typing import IO
from typing import Tuple

from pyfaidx import Fasta
from pyfaidx import FastaNotFoundError

from ..basicClass import MicroHomology
from ..basicClass import Node
from ..basicClass import NovelInsertion
from ..basicClass import reverse_complement
from ..basicClass import Series
from ..type import LoggerType
from .writer import Writer


class FastaWriter(Writer):
    """Writer for Fasta files."""

    def __init__(self, file_path: str, reference: str, logger: LoggerType):
        """Initialize FastaWriter object."""
        super().__init__(file_path, logger)
        self.reference = Path(reference)
        if not self.reference.exists():
            raise FastaNotFoundError(f"{self.reference} does not exist.")
        self.reference_io = Fasta(reference, sequence_always_upper=True)
        self.id = 1

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, seq_id: str, sequence: str) -> str:
        """Formatter for writing data."""
        if sequence == "":
            self.logger.warning(
                f"{self.__class__.__name__}: Sequence ID or sequence is empty."
            )
        return f">{seq_id:0>6}\n{sequence}\n"

    def open(self, mode: str = "w") -> IO:
        """Open file."""
        if self.is_opened:
            self.logger.warning(f"{self.__class__.__name__}: File is already opened.")
        self.io = self.file_path.open(mode)
        if hasattr(self, "write_header"):
            self.write_header()  # type: ignore
        return self.io

    def close(self) -> None:
        """Close file."""
        if self.is_opened:
            self.logger.trace(f"{self.__class__.__name__}: Closing file.")
            self.io.close()  # type: ignore
            self.io = None

    def write_line(self, line: str) -> None:
        """Write line to file."""
        if self.is_opened:
            self.id += 1
            self.io.write(line)  # type: ignore
        else:
            self.logger.warning(f"{self.__class__.__name__}: File is not opened.")

    @singledispatchmethod
    def write_data(self, data_object: Any, object_id: int):
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: Series, object_id: int = -1):
        """Write Series to fasta file."""
        if len(data_object.nodes) == 0:
            self.logger.warning(
                f"{self.__class__.__name__}: No nodes to write to file in Clique {object_id} Series."
            )
        sequence, node_length_str = get_nodes_sequence_from_series(
            data_object, reference_io=self.reference_io
        )

        self.write_line(self.formatter(f"{self.id} {node_length_str}", sequence))


def get_nodes_sequence_from_series(
    series: Series, reference_io: Fasta
) -> Tuple[str, str]:
    """Get sequence of nodes of series.

    :param series: Series including nodes.
    :param reference_io: ReferenceIO object.

    :return: Sequence of nodes.
    """
    sequence = ""
    node_length_str = ""
    for node in series:
        node_seq = get_exon_sequence_from_node(node, reference_io)
        sequence += node_seq
        node_length_str += f"{len(node_seq)}|"
    return sequence, node_length_str[:-1]


def get_exon_sequence_from_node(node: Node, reference_io: Fasta) -> str:
    """Get exon sequence of a node.

    remove microhomology from the sequence, and add novel insertion sequence.

    :param node: Node and InsertionType
    :param reference_io: ReferenceIO object
    :return: sequence of exon for one node
    """
    node_sequence = ""

    # positive strand sequence for novel insertion
    microhomology_sequence, novel_insertion_sequence = "", ""
    if node.insertion_info and not node.insertion_info[0]:
        insertion = node.insertion_info[1]
        if isinstance(insertion, NovelInsertion):
            novel_insertion_sequence += insertion.query_sequence
        if isinstance(insertion, MicroHomology):
            microhomology_sequence += insertion.query_sequence

    for start, end in node.exons:  # type: ignore
        node_sequence += reference_io.get_seq(node.chrom, start + 1, end).seq  # 1-based

    node_sequence = (
        node_sequence if node.strand == "+" else reverse_complement(node_sequence)
    )

    node_sequence += novel_insertion_sequence
    return node_sequence[: len(node_sequence) - len(microhomology_sequence)]
