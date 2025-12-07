"""GTF writer class."""

from __future__ import annotations

import copy
from functools import singledispatchmethod
from typing import IO, TYPE_CHECKING, Any

from loguru import logger

from scanmst.base import MicroHomology, NovelInsertion
from scanmst.exception import ExonsNotFoundError
from scanmst.graph import Edge, NLPath, Node

from .writer import Writer

if TYPE_CHECKING:
    from scanmst.graph import Edge, Node


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

    def __init__(self, file_path: str, rescue_sr: bool) -> None:
        """Initialize GTFWriter object."""
        super().__init__(file_path)
        self.rescue_sr = rescue_sr

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
    def write_data(self, data_object: Any, object_id: str) -> None:  # type: ignore
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: NLPath, object_id: str) -> None:
        """Write Series to GTF file.

        :param data_object:
        :return:
        """
        if len(data_object.nodes) == 0:
            logger.warning(
                f"{self.__class__.__name__}: No nodes to write to file in Clique {object_id} Series.",
            )

        for node_gtf_feature in get_nodes_gtf_features_from_nlpath(
            data_object,
            self.rescue_sr,
            object_id,
        ):
            self.write_line(self.formatter(node_gtf_feature))


def get_nodes_gtf_features_from_nlpath(
    nlpath: NLPath,
    rescue_sr: bool,
    cluster_id: str,
) -> list[list[str]]:
    """Get GTF features of nodes of series.

    :param series: Series including nodes.
    :param series_id: Series ID.

    :return: List of GTF features for node and insertions in the series.
    """
    nlpath_gtf_features = [None]

    min_nlpath_sr = float("inf")
    min_nlpath_originla_sr = float("inf")

    for idx, node in enumerate(nlpath, 1):
        edge = nlpath.next_edge(node, idx - 1)
        insertion_info = None if edge is None else edge.insertion_info

        if edge is not None:
            min_nlpath_sr = min(min_nlpath_sr, edge.sr)
            min_nlpath_originla_sr = min(min_nlpath_originla_sr, edge.original_sr)

        nlpath_gtf_features.extend(
            list(
                get_gtf_features_from_node(
                    node,
                    edge,
                    nlpath.id,
                    cluster_id,
                )
            ),
        )

        if insertion_info and isinstance(insertion_info[1], NovelInsertion):
            nlpath_gtf_features.append(
                get_gtf_features_from_insertion(
                    insertion_info[1],
                    nlpath.id,
                    node.id,
                    cluster_id=cluster_id,
                ),
            )

    nlpath_gtf_features[0] = format_gtf_features_for_nlpath(
        nlpath_id=nlpath.id,
        nlpath_sr=min_nlpath_sr,
        nlpath_originla_sr=min_nlpath_originla_sr,
        cluster_id=cluster_id,
        rescue_sr=rescue_sr,
        extend=nlpath.extension,
    )

    return nlpath_gtf_features


def add_info_to_attribute_column(col_list: list[str], add_info: str):
    """Add additional info. to the 9th column of GTF."""
    col_list[-1] += add_info
    return col_list


def format_gtf_features_for_nlpath(
    nlpath_id: str,
    nlpath_sr: float,
    nlpath_originla_sr: float,
    cluster_id: str,
    *,
    rescue_sr: bool,
    extend: bool = False,
) -> list[str]:
    """Get GTF features of transcript."""
    if rescue_sr:
        return [
            ".",
            "scanmst",
            "transcript",
            ".",
            ".",
            ".",
            ".",
            ".",
            f'sr "{nlpath_sr}"; osr "{nlpath_originla_sr}"; transcript_id "{nlpath_id}"; gene_id "{cluster_id}"; extend "{extend}";',
        ]
    return [
        ".",
        "scanmst",
        "transcript",
        ".",
        ".",
        ".",
        ".",
        ".",
        f'sr "{nlpath_sr}"; osr "{nlpath_originla_sr}"; transcript_id "{nlpath_id}"; gene_id "{cluster_id}"; extend "{extend}";',
    ]


def get_gtf_features_from_insertion(
    insertion: NovelInsertion,
    nlpath_id: str,
    node_id: str,
    cluster_id: str,
) -> list[str]:
    """Get GTF features of novel insertion."""
    return [
        ".",
        "scanmst",
        "insertion",
        ".",
        ".",
        ".",
        "+",
        ".",
        f'transcript_id "{nlpath_id}"; gene_id "{cluster_id}"; sequence "{insertion.query_sequence}";',
    ]


def get_gtf_features_from_node(
    node: Node,
    edge: Edge | None,
    nlpath_id: str,
    cluster_id: str,
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

    insertion_info = None
    if edge is not None:
        # sr may be not exported in gtf file
        insertion_info = edge.insertion_info

    microhomology_sequence = ""
    if insertion_info is not None and not insertion_info[0]:
        insertion = insertion_info[1]
        if isinstance(insertion, MicroHomology):
            microhomology_sequence += insertion.query_sequence

    nodes_gtf_features = []

    for index, (start, end) in enumerate(copy_exons, 1):
        nodes_gtf_features.append(
            [
                f"{node.chrom}",
                "scanmst",
                "exon",
                f"{start + 1}",
                f"{end}",
                ".",
                f"{node.strand}",
                ".",
                f'exon_id "{index:0>3}"; segment_id "{node.id}"; ptc "{node.ptc}"; ptf "{node.ptf}"; transcript_id "{nlpath_id}"; gene_id "{cluster_id}";',
            ],
        )

    return nodes_gtf_features
