"""FastaWriter class."""

from functools import singledispatchmethod
from pathlib import Path
from typing import IO, Any

import pyabpoa
from loguru import logger
from pyfaidx import Fasta, FastaNotFoundError

from scannls.base import MicroHomology, NovelInsertion, reverse_complement
from scannls.graph import NLPath, Node

from .writer import Writer


class FastaWriter(Writer):
    """Writer for Fasta files."""

    def __init__(
        self,
        file_path: str,
        reference: str,
        output_sequence_choice: str,
        read_name_to_seq_dict: dict,
    ) -> None:
        """Initialize FastaWriter object."""
        super().__init__(file_path)
        self.reference = Path(reference)
        if not self.reference.exists():
            msg = f"{self.reference} does not exist."
            raise FastaNotFoundError(msg)
        self.reference_io = Fasta(reference, sequence_always_upper=True)
        self.output_sequence_choice = output_sequence_choice
        if output_sequence_choice == "consensus":
            self.read_name_to_seq_dict = read_name_to_seq_dict
            self.msa_aligner = pyabpoa.msa_aligner()
        self.id = 1

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, seq_id: str, sequence: str) -> str:
        """Formatter for writing data."""
        if not sequence:
            logger.warning(
                f"{self.__class__.__name__}: Sequence ID or sequence is empty.",
            )
        return f">{seq_id:0>6}\n{sequence}\n"

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
            self.id += 1
            self.io.write(line)  # type: ignore
        else:
            logger.warning(f"{self.__class__.__name__}: File is not opened.")

    @singledispatchmethod
    def write_data(self, data_object: Any, object_id: str):  # type: ignore
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: NLPath, object_id: str):
        """Write NLPath to fasta file. object_id is cluster_id."""

        if len(data_object.nodes) == 0:
            logger.warning(
                f"{self.__class__.__name__}: No nodes to write to file in Clique {object_id} Series.",
            )
        if self.output_sequence_choice == "reference":
            sequence, node_length_str = get_nodes_sequence_from_series(
                data_object,
                reference_io=self.reference_io,
            )
        elif self.output_sequence_choice == "consensus":
            sequence, node_length_str = get_consensus_sequence_from_series(
                data_object,
                reference_io=self.reference_io,
                read_name_to_seq_dict=self.read_name_to_seq_dict,
                msa_aligner=self.msa_aligner,
            )

        self.write_line(
            self.formatter(f"{data_object.id} {node_length_str}", sequence),
        )


def get_consensus_sequence_from_series(
    nlpath: NLPath,
    reference_io: Fasta,
    read_name_to_seq_dict: dict,
    msa_aligner: pyabpoa.msa_aligner,
) -> tuple[str, str]:
    """Get sequence of nodes of series.

    :param series: Series including nodes.
    :param reference_io: ReferenceIO object.

    :return: Sequence of nodes.
    """
    node_length_str = ""
    read_names_for_nlpath = set()
    for idx, node in enumerate(nlpath):
        edge = nlpath.next_edge(node, idx)

        read_names_at_edge = set() if edge is None else set(edge.read_ids)

        if len(read_names_for_nlpath) == 0 or 0 < len(read_names_at_edge) < len(read_names_for_nlpath):
            read_names_for_nlpath = read_names_at_edge

        insertion_info = None if edge is None else edge.insertion_info
        node_seq = get_exon_sequence_from_node(node, insertion_info, reference_io)
        node_length_str += f"{len(node_seq)}|"

    read_supporting_seqs = []
    for read_name in read_names_for_nlpath:
        if read_name in read_name_to_seq_dict:
            read_supporting_seqs.append(read_name_to_seq_dict[read_name])
        else:
            logger.warning(f"{read_name=} is not found in name_to_seq_dict.")

    result = msa_aligner.msa(read_supporting_seqs, out_cons=True, out_msa=False)
    sequence = result.cons_seq[0]

    return sequence, node_length_str[:-1]


def get_nodes_sequence_from_series(
    nlpath: NLPath,
    reference_io: Fasta,
) -> tuple[str, str]:
    """Get sequence of nodes of series.

    :param series: Series including nodes.
    :param reference_io: ReferenceIO object.

    :return: Sequence of nodes.
    """
    sequence = ""
    node_length_str = ""
    for idx, node in enumerate(nlpath):
        edge = nlpath.next_edge(node, idx)
        insertion_info = None if edge is None else edge.insertion_info
        node_seq = get_exon_sequence_from_node(node, insertion_info, reference_io)
        sequence += node_seq
        node_length_str += f"{len(node_seq)}|"
    return sequence, node_length_str[:-1]


def get_exon_sequence_from_node(node: Node, insertion_info, reference_io: Fasta) -> str:
    """Get exon sequence of a node.

    remove microhomology from the sequence, and add novel insertion sequence.

    :param node: Node and InsertionType
    :param reference_io: ReferenceIO object
    :return: sequence of exon for one node
    """
    node_sequence = ""

    # positive strand sequence for novel insertion
    microhomology_sequence, novel_insertion_sequence = "", ""
    if insertion_info and not insertion_info[0]:
        insertion = insertion_info[1]
        if isinstance(insertion, NovelInsertion):
            novel_insertion_sequence += insertion.query_sequence
        if isinstance(insertion, MicroHomology):
            microhomology_sequence += insertion.query_sequence

    for start, end in node.exons:  # type: ignore
        node_sequence += reference_io.get_seq(node.chrom, start + 1, end).seq  # 1-based

    node_sequence = node_sequence if node.strand.is_forward() else reverse_complement(node_sequence)

    node_sequence += novel_insertion_sequence
    return node_sequence[: len(node_sequence)]
