"""GTF writer class.

@Filename:    gtfWriter.py
@Author:      YangyangLi
@Time:        1/30/22 6:18 PM
"""
from __future__ import annotations

import copy
from functools import singledispatchmethod
from typing import IO, TYPE_CHECKING, Any

from loguru import logger

from scannls import MicroHomology, NovelInsertion
from scannls.exception import ExonsNotFoundError
from scannls.graph import Edge, NLPath, Node

from .writer import Writer

if TYPE_CHECKING:
    from scannls.graph import Edge, NLPath, Node  # noqa: F811


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

    def __init__(self, file_path: str) -> None:
        """Initialize GTFWriter object."""
        super().__init__(file_path)
        self.id = 1

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, fields: list[str], delimiter: str = "\t") -> str:
        """Formatter for writing data."""
        if fields is None or len(fields) != GTFWriter.num_fields:
            logger.warning(
                f"{self.__class__.__name__}: Number of fields is not equal to 9.",
            )
        return delimiter.join(fields) + "\n"

    def open(self, mode: str = "w") -> IO:
        """Open file."""
        if self.is_opened:
            logger.warning(f"{self.__class__.__name__}: File is already opened.")
        self.io = self.file_path.open(mode)
        if hasattr(self, "write_header"):
            self.write_header()  # type: ignore
        return self.io

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
    def write_data(self, data_object: Any, object_id: int) -> None:  # type: ignore
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: NLPath, object_id: int) -> None:
        """Write Series to GTF file.

        :param data_object:
        :return:
        """
        if len(data_object.nodes) == 0:
            logger.warning(
                f"{self.__class__.__name__}: No nodes to write to file in Clique {object_id} Series.",
            )

        for node_gtf_feature in get_nodes_gtf_features_from_series(
            data_object,
            self.id,
        ):
            self.write_line(self.formatter(node_gtf_feature))
        self.id += 1


def get_nodes_gtf_features_from_series(
    nlpath: NLPath,
    nlpath_id: int,
) -> list[list[str]]:
    """Get GTF features of nodes of series.

    :param series: Series including nodes.
    :param series_id: Series ID.

    :return: List of GTF features for node and insertions in the series.
    """
    series_gtf_features = []

    nlpath_sr_list = []
    nlpath_originla_sr_list = []
    for node_id, node in enumerate(nlpath, 1):
        edge = nlpath.next_edge(node, node_id - 1)
        if edge is not None:
            nlpath_sr_list.append(edge.sr)
            nlpath_originla_sr_list.append(edge.original_sr)

    nlpath_sr = min(nlpath_sr_list)
    nlpath_originla_sr = min(nlpath_originla_sr_list)

    series_gtf_features.append(
        get_gtf_features_for_nlpath(nlpath_id, nlpath_sr, nlpath_originla_sr)
    )

    for node_id, node in enumerate(nlpath, 1):
        edge = nlpath.next_edge(node, node_id - 1)
        insertion_info = None if edge is None else edge.insertion_info

        series_gtf_features.extend(
            get_gtf_features_from_node(node, edge, nlpath_id, node_id),
        )
        if insertion_info and isinstance(insertion_info[1], NovelInsertion):
            series_gtf_features.append(
                get_gtf_features_from_insertion(
                    insertion_info[1],
                    nlpath_id,
                    node_id,
                ),
            )
    return series_gtf_features


def get_gtf_features_for_nlpath(
    nlpath_id: int, nlpath_sr: int, nlpath_originla_sr: int
) -> list[str]:
    """Get GTF features of transcript."""
    return [
        ".",
        "scannls",
        "transcript",
        ".",
        ".",
        ".",
        ".",
        ".",
        f'transcript_id "{nlpath_id:0>6}"; '
        f'sr "{nlpath_sr}"; '
        f'osr "{nlpath_originla_sr}";',
    ]


def get_gtf_features_from_insertion(
    insertion: NovelInsertion,
    nlpath_id: int,
    node_id: int,
) -> list[str]:
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
        f'transcript_id "{nlpath_id:0>6}"; mega_exon_id "{node_id:0>3}"; '
        f'sequence "{insertion.query_sequence}";',
    ]


def get_gtf_features_from_node(
    node: Node,
    edge: Edge | None,
    nlpath_id: int,
    node_id: int,
) -> list[list[str]]:
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
        msg = f"{node.query_name}"
        raise ExonsNotFoundError(msg)

    exons = node.exons.reversed() if node.strand.is_reverse() else node.exons

    if exons is None:
        msg = f"{node.query_name}"
        raise ExonsNotFoundError(msg)

    copy_exons = copy.deepcopy(exons)

    if edge is not None:
        # WARN: should be changed after correcting sr <08-02-23, Yangyang Li>
        # sr may be not exported in gtf file
        node_sr, node_original_sr, insertion_info = (
            edge.sr,
            edge.original_sr,
            edge.insertion_info,
        )
    else:
        # WARN: last node has no edge <08-02-23, Yangyang Li>
        node_sr, node_original_sr, insertion_info = 0, 0, None

    microhomology_sequence = ""
    if insertion_info and not insertion_info[0]:
        insertion = insertion_info[1]
        if isinstance(insertion, MicroHomology):
            microhomology_sequence += insertion.query_sequence

    # last exon end position needs a correction if there is a microhomology.
    if node.strand.is_forward() and exons.last.start < exons.last.end - len(
        microhomology_sequence,
    ):
        copy_exons.last.end -= len(microhomology_sequence)

    elif (
        node.strand.is_reverse()
        and exons.last.start + len(microhomology_sequence) < exons.last.end
    ):
        copy_exons.last.start += len(microhomology_sequence)

    nodes_gtf_features = []

    for index, (start, end) in enumerate(copy_exons, 1):
        if node_sr + node_original_sr > 0:
            info = [
                f'transcript_id "{nlpath_id:0>6}"; '
                f'mega_exon_id "{node_id:0>3}"; '
                f'exon_id "{index:0>3}"; '
                f'sr "{node_sr}"; '
                f'osr "{node_original_sr}";',
            ]
        else:
            info = [
                f'transcript_id "{nlpath_id:0>6}"; '
                f'mega_exon_id "{node_id:0>3}"; '
                f'exon_id "{index:0>3}";',
            ]
        nodes_gtf_features.append(
            [
                f"{node.chrom}",
                "scannls",
                "exon",
                f"{start + 1}",
                f"{end}",
                ".",
                f"{node.strand}",
                ".",
                *info,
            ],
        )
    return nodes_gtf_features
