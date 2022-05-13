# !/usr/bin/env python
"""GTF writer class.

@Filename:    gtfWriter.py
@Author:      YangyangLi
@license:     MIT Licence
@Time:        1/30/22 6:18 PM
"""
from functools import singledispatchmethod
from typing import Any
from typing import IO
from typing import List

from ..basicClass import Node
from ..basicClass import NovelInsertion
from ..basicClass import Series
from ..exception import ExonsNotFoundError
from ..type import LoggerType
from .writer import Writer


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
        return delimiter.join(fields) + "\n"

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
            self.io.write(line)  # type: ignore
        else:
            self.logger.warning(f"{self.__class__.__name__}: File is not opened.")

    @singledispatchmethod
    def write_data(self, data_object: Any, object_id: int) -> None:  # type: ignore
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: Series, object_id: int) -> None:
        """Write Series to GTF file.

        :param data_object:
        :return:
        """
        if len(data_object.nodes) == 0:
            self.logger.warning(
                f"{self.__class__.__name__}: No nodes to write to file in Clique {object_id} Series."
            )
        for node_gtf_feature in get_nodes_gtf_features_from_series(
            data_object, self.id
        ):
            self.write_line(self.formatter(node_gtf_feature))
        self.id += 1


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
        "+",
        ".",
        f'transcript_id "{series_id:0>6}"; mega_exon_id "{node_id:0>3}"; '
        f'sequence "{insertion.query_sequence}";',
    ]


def get_gtf_features_from_node(
    node: Node, series_id: int, node_id: int
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
        raise ExonsNotFoundError(f"{node.query_name}")

    exons = node.exons[::-1] if node.strand == "-" else node.exons

    nodes_gtf_features = []

    for index, (start, end) in enumerate(exons, 1):
        info = [
            f'transcript_id "{series_id:0>6}"; '
            f'mega_exon_id "{node_id:0>3}"; '
            f'exon_id "{index:0>3}";'
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
