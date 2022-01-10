"""Writer for Fasta files and GTF files for Series object.

@Filename:    writer.py
@license:     MIT Licence
@Time:        12/30/21 4:02 PM
"""
from abc import ABC
from abc import abstractmethod
from functools import singledispatchmethod
from pathlib import Path
from typing import Any
from typing import IO
from typing import List
from typing import Optional

from pyfaidx import Fasta  # type: ignore
from pyfaidx import FastaNotFoundError

from ..type import LoggerType
from .basicClass import MicroHomology
from .basicClass import NodeType
from .basicClass import NovelInsertion
from .basicClass import reverse_complement
from .basicClass import Series
from .exception import ExonsNotFoundError


# todo: add comments line
# todo: add asyncio support


class Writer(ABC):
    """Abstract class for writing object to file."""

    def __init__(self, file_path: str, logger: LoggerType):
        """Initialize Writer object."""
        self.logger = logger
        self.file_path = Path(file_path)
        if self.file_path.exists():
            self.logger.warning(f"{self.file_path} exists, will be overwritten.")
        self.io: Optional[IO] = None

    @abstractmethod
    def write_data(self, data_object: Any):
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @abstractmethod
    def write_line(self, line: str):
        """Write line to file.

        :param: line: Line to write to file.
        """

    @abstractmethod
    def open(self, mode: str = "w"):
        """Open file.

        :param: mode: Mode to open file.
        """

    @abstractmethod
    def close(self):
        """Close file."""


class FastaWriter(Writer):
    """Writer for Fasta files."""

    def __init__(self, file_path: str, reference: str, logger: LoggerType):
        """Initialize FastaWriter object."""
        super().__init__(file_path, logger)
        self.reference = Path(reference)
        if not self.reference.exists():
            raise FastaNotFoundError
        self.reference_io = Fasta(reference, sequence_always_upper=True)
        self.id = 1

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, seq_id: int, sequence: str) -> str:
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
        self.io = open(self.file_path, mode)  # add asyncio support
        return self.io

    def close(self) -> None:
        """Close file."""
        if self.is_opened:
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
    def write_data(self, data_object: Any):
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: Series):
        """Write Series to fasta file."""
        if len(data_object.nodes) == 0:
            self.logger.warning(
                f"{self.__class__.__name__}: No nodes to write to file."
            )
        self.logger.trace(f"{self.__class__.__name__}: Writing Series to file.")
        sequence = get_nodes_sequence_from_series(
            data_object, reference_io=self.reference_io
        )
        self.write_line(self.formatter(self.id, sequence))


class GTFWriter(Writer):
    """Writer for GTF files.

    .. note::
        1. seqname: chromosome
        2. source:  name of the program that generated the feature
        3. feature: feature type name, eg. gene, mRNA, exon, CDS
        4. start: start position of the feature
        5. end: end position of the feature
        6. score: a floating point value
        7. strand: defined as + (forward) - (reverse) or . (unknown)
        8. frame: one of '0', '1' or '2'. '0' indicates that the first base of the feature
                  is the first base of a codon, '1' that the second base is the first base
                  of a codon, and so on..
        9. attribute: a semicolon-separated list of tag-value pairs (separated by spaces)
    """

    num_fields: int = 9

    def __init__(self, file_path: str, logger: LoggerType) -> None:
        """Initialize GTFWriter object."""
        super().__init__(file_path, logger)
        self.id = 1

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, fields: List[str], delimiter: str = "\t") -> str:
        """Formatter for writing data."""
        if fields is None or len(fields) != GTFWriter.num_fields:
            self.logger.warning(
                f"{self.__class__.__name__}: Number of fields is not equal to 9."
            )
            raise SystemExit
        return delimiter.join(fields) + "\n"

    def open(self, mode: str = "w") -> IO:
        """Open file."""
        if self.is_opened:
            self.logger.warning(f"{self.__class__.__name__}: File is already opened.")
        self.io = open(self.file_path, mode)  # add asyncio support
        return self.io

    def close(self) -> None:
        """Close file."""
        if self.is_opened:
            self.io.close()  # type: ignore
            self.io = None

    def write_line(self, line: str) -> None:
        """Write line to file."""
        if self.is_opened:
            self.io.write(line)  # type: ignore
        else:
            self.logger.warning(f"{self.__class__.__name__}: File is not opened.")

    @singledispatchmethod
    def write_data(self, data_object: Any) -> None:
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: Series) -> None:
        """Write Series to fasta file.

        :param data_object:
        :return:
        """
        if len(data_object.nodes) == 0:
            self.logger.warning(
                f"{self.__class__.__name__}: No nodes to write to file."
            )
        self.logger.trace(f"{self.__class__.__name__}: Writing Series to file.")
        for node_gtf_feature in get_nodes_gtf_features_from_series(
            data_object, self.id
        ):
            self.write_line(self.formatter(node_gtf_feature))
        self.id += 1


def get_nodes_sequence_from_series(series: Series, reference_io) -> str:
    """Get sequence of nodes of series.

    :param series: Series including nodes.
    :param reference_io: ReferenceIO object.

    :return: Sequence of nodes.
    """
    sequence = ""
    for node in series:
        sequence += get_exon_sequence_from_node(node, reference_io)
    return sequence


def get_exon_sequence_from_node(node: NodeType, reference_io) -> str:
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


def get_nodes_gtf_features_from_series(
    series: Series, series_id: int
) -> List[List[str]]:
    """Get GTF features of nodes of series.

    :param series: Series including nodes.
    :param series_id: Series ID.

    :return: List of GTF features for node and insertions in the series.
    """
    series_gtf_features = []
    for node_id, node in enumerate(series, 1):
        series_gtf_features.extend(get_gtf_features_from_node(node, series_id, node_id))
        if node.insertion_info and isinstance(node.insertion_info[1], NovelInsertion):
            series_gtf_features.append(
                get_gtf_features_from_insertion(
                    node.insertion_info[1], series_id, node_id
                )
            )
    return series_gtf_features


def get_gtf_features_from_insertion(
    insertion: NovelInsertion, series_id: int, node_id: int
) -> List[str]:
    """Get GTF features of novel insertion."""
    return [
        ".",
        "scannls",
        "insertion",
        ".",
        ".",
        ".",
        ".",
        ".",
        f'transcript_id "{series_id:0>6}"; mega_exon_id "{node_id:0>6}"; '
        f'sequence "{insertion.query_sequence}" ',
    ]


def get_gtf_features_from_node(
    node: NodeType, series_id: int, node_id: int
) -> List[List[str]]:
    """Get exon gtf features of a node.

    :param node_id: node id
    :param series_id: series id
    :param node: Node and Insertion
    :return: list of gtf features (8 columns) except for the attribute column

    .. note::
        1. seqname: chromosome
        2. source:  name of the program that generated the feature
        3. feature: feature type name, eg. gene, mRNA, exon, CDS
        4. start: start position of the feature
        5. end: end position of the feature
        6. score: a floating point value
        7. strand: defined as + (forward) - (reverse) or . (unknown)
        8. frame: one of '0', '1' or '2'. '0' indicates that the first base of the feature
                  is the first base of a codon, '1' that the second base is the first base
                  of a codon, and so on..
        9. attribute: a semicolon-separated list of tag-value pairs (separated by spaces)
    """
    if node.exons is None:
        raise SystemExit from ExonsNotFoundError

    exons = node.exons[::-1] if node.strand == "-" else node.exons

    nodes_gtf_features = []

    for index, (start, end) in enumerate(exons, 1):
        info = [
            f'transcript_id "{series_id:0>6}"; '
            f'mega_exon_id "{node_id:0>6}"; '
            f'exon_id "{index:0>6}" '
        ]
        nodes_gtf_features.append(
            [
                f"{node.chrom}",
                "exon",
                "scannls",
                f"{start + 1}",
                f"{end}",
                ".",
                f"{node.strand}",
                ".",
            ]
            + info
        )
    return nodes_gtf_features
