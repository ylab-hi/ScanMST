"""Writer for Fasta files and GTF files for Series object.

@Filename:    writer.py
@license:     MIT Licence
@Time:        12/30/21 4:02 PM
"""
import datetime
from abc import ABC
from abc import abstractmethod
from functools import singledispatchmethod
from pathlib import Path
from typing import Any
from typing import IO
from typing import List
from typing import Optional

import pysam  # type: ignore
from pyfaidx import Fasta  # type: ignore
from pyfaidx import FastaNotFoundError

from .. import __version__
from ..type import LoggerType
from .basicClass import MicroHomology
from .basicClass import NodeType
from .basicClass import NovelInsertion
from .basicClass import reverse_complement
from .basicClass import Series
from .exception import ExonsNotFoundError
from .exception import GenesNotFoundError
from .exception import ModesNotFoundError


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
        """Write Series to GTF file.

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


class VCFWriter(Writer):
    """Writer for VCF files.

    .. note::
        1. CHROM: The name of the sequence (typically a chromosome) on which the variation
            is being called. This sequence is usually known as 'the reference sequence',
            i.e. the sequence against which the given sample varies.
        2. POS: The 1-based position of the variation on the given sequence.
        3. ID: The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown
            a ".". Multiple identifiers should be separated by semi-colons without white-space.
        4. REF:The reference base (or bases in the case of an indel) at the given position
            on the given reference sequence.
        5. ALT: The list of alternative alleles at this position.
        6. QUAL: A quality score associated with the inference of the given alleles.
        7. FILTER: A flag indicating which of a given set of filters the variation has
            failed or PASS if all the filters were passed successfully.
        8. INFO: An extensible list of key-value pairs (fields) describing the variation.
            See below for some common fields. Multiple fields are separated by semicolons
            with optional values in the format: <key>=<data>[,data].
        9. FORMAT: An (optional) extensible list of fields for describing the samples.
            See below for some common fields.
        10. SAMPLE: For each (optional) sample described in the file,
            values are given for the fields listed in FORMAT
    """

    def __init__(
        self,
        file_path: str,
        reference: str,
        bam_io: pysam.AlignmentFile,
        output_prefix: str,
        logger: LoggerType,
    ) -> None:
        """Initialize VCFWriter object."""
        super().__init__(file_path, logger)
        self.reference = Path(reference)
        if not self.reference.exists():
            raise FastaNotFoundError
        self.reference_io = Fasta(reference, sequence_always_upper=True)
        self.id = 1
        self.bam_io = bam_io
        self.bam_header = bam_io.header
        self.sample_name = output_prefix

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, fields: List[str], delimiter: str = "\t") -> str:
        """Formatter for writing data."""
        if fields is None or len(fields) != 10:
            self.logger.warning(
                f"{self.__class__.__name__}: Number of fields is not equal to 10."
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

    def write_header(self) -> None:
        """Write header to VCF file."""
        self.write_line(self.header)

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
                f"{self.__class__.__name__}: No nodes to write to VCF file."
            )
        self.logger.trace(f"{self.__class__.__name__}: Writing Series to VCF file.")
        for hop_vcf_feature in get_hops_vcf_features_from_series(
            data_object, self.id, self.bam_io
        ):
            self.write_line(self.formatter(hop_vcf_feature))
        self.id += 1

    @property
    def header(self) -> str:
        """VCF header provides metadata describing the body of the file."""
        # Metadata parsers/constants
        reserved_info = {
            "CANONICAL": "Flag",
            "NONCANONICAL": "Flag",
            "BOUNDARY": "String",
            "DP": "Integer",
            "DP1": "Integer",
            "DP2": "Integer",
            "SR": "Integer",
            "PSO": "Float",
            "AO": "Integer",
            "AF": "Float",
            "SVMETHOD": "String",
            "SVTYPE": "String",
            "SVLEN": "Integer",
            "CHR2": "String",
            "END": "Integer",
            "STRAND": "String",
            "STRAND1": "String",
            "STRAND2": "String",
            "MODE1": "String",
            "MODE2": "String",
            "GENE": "String",
            "GENE1": "String",
            "GENE2": "String",
        }
        reserved_format = {"GT": "String"}
        reserved_alt = ["INS", "DEL", "TDUP", "IDUP", "INV", "TRA"]

        description = {
            "CANONICAL": "Canonical splice site",
            "NONCANONICAL": "Noncanonical splice site",
            "BOUNDARY": "The coding exon boundary type of event, BOTH, LEFT, RIGHT, NEITHER.",
            "DP": "Total read depth at the breakpoint for insertion",
            "DP1": "Total read depth at the breakpoint1",
            "DP2": "Total read depth at the breakpoint2",
            "SR": "The number of support reads for the breakpoints",
            "AO": "Alternate allele observations, "
            "with partial observations recorded fractionally",
            "AF": "Estimated allele frequency in the range (0,1], "
            "representing the ratio of reads showing the alternative allele to all reads",
            "PSO": "Estimated Percent splice-out in the range (0,1], "
            "representing the percentage of NLS transcripts",
            "SVTYPE": "The type of event, INS, DEL, TDUP, IDUP, INV, TRA.",
            "SVLEN": "Difference in length between REF and ALT alleles",
            "CHR2": "Chromosome for END coordinate in case of a translocation",
            "END": "2nd position of the structural variant",
            "GENE": "Overlapped coding gene for insertion",
            "GENE1": "Overlapped coding gene for breakpoint1",
            "GENE2": "Overlapped coding gene for breakpoint2",
            "TRANSCRIPT_ID": "Transcript ID",
            "SVMETHOD": "Type of approach used to detect SV",
            "STRAND": "Strand for insertion",
            "STRAND1": "Strand for breakpoint1",
            "STRAND2": "Strand for breakpoint2",
            "MODE1": "Mode for softclipped reads at breakpoint1",
            "MODE2": "Mode for softclipped reads at breakpoint2",
            "GT": "Genotype",
            "INS": "Insertion",
            "DEL": "Deletion",
            "TDUP": "Tandem duplication",
            "IDUP": "Inverted duplication",
            "INV": "Inversion",
            "TRA": "Translocation",
        }

        date = datetime.datetime.today().strftime("%Y%m%d")
        source = f"ScanNLS v{__version__}"
        reference = f'<CMD={obtain_reference_from_bam_header(self.bam_header)},Description="Alignment parameters">'
        self.sample_name

        header_lines = []
        header_lines.append("##fileformat=VCFv4.3")
        header_lines.append(f"##fileDate={date}")
        header_lines.append(f"##source={source}")
        header_lines.append(f"##reference={reference}")

        for _id in reserved_info:
            if reserved_info[_id] == "Flag":
                _number = 0
            else:
                _number = 1
            header_lines.append(
                f'##INFO=<ID={_id},Number={_number},Type={reserved_info[_id]},Description="{description[_id]}">'
            )

        for _id in reserved_format:
            header_lines.append(
                f'##FORMAT=<ID={_id},Number=1,Type={reserved_format[_id]},Description="{description[_id]}">'
            )

        for _id in reserved_alt:
            header_lines.append(f'##ALT=<ID={_id},Description="{description[_id]}">')
        header_lines.append(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}"
        )

        return "\n".join(header_lines) + "\n"


def obtain_reference_from_bam_header(bam_header):
    """Obtain reference info from BAM header."""
    _aligners = {
        "CLC",
        "ContextMap2",
        "CRAC",
        "GSNAP",
        "HISAT",
        "HISAT2",
        "MapSplice2",
        "Novoalign",
        "OLego",
        "RUM",
        "SOAPsplice",
        "STAR",
        "Subread",
        "TopHat",
        "TopHap2",
        "bwa",
        "bowtie",
        "bowtie2",
        "minimap2",
    }
    avail_aligners = {x.upper() for x in _aligners}
    aln_cmd = "Unknown"
    if "PG" in bam_header:
        for j in bam_header["PG"]:
            if j["ID"].upper() in avail_aligners:
                aln_cmd = j["CL"]
                break
            else:
                aln_cmd = "Unknown"
    return aln_cmd


def get_novel_insertion_sequence_from_node(node: NodeType, reference_io) -> tuple:
    """Get novel insertion sequence of a node.

    :param node: Node and InsertionType
    :param reference_io: ReferenceIO object
    :return: sequence of exon for one node
    """
    # positive strand sequence for novel insertion
    novel_insertion_sequence = ""
    if node.insertion_info and not node.insertion_info[0]:
        insertion = node.insertion_info[1]
        if isinstance(insertion, NovelInsertion):
            novel_insertion_sequence += insertion.query_sequence

    if node.exons is None:
        raise SystemExit from ExonsNotFoundError

    _pos = node.exons[-1][1] if node.strand == "+" else node.exons[0][0]

    ref_allele = reference_io.get_seq(node.chrom, _pos + 1, _pos + 1).seq  # 1-based

    if not novel_insertion_sequence:
        return (ref_allele, None)

    alt_allele = (
        novel_insertion_sequence
        if node.strand == "+"
        else reverse_complement(novel_insertion_sequence)
    )

    return (ref_allele, alt_allele)


def get_hops_vcf_features_from_series(
    series: Series, series_id: int, reference_io
) -> List[List[str]]:
    """Obtain hop vcf features from one series."""
    series_hops_features = []
    ins_id = 1
    for event_id, idx in enumerate(range(len(series) - 1), 1):
        current_node = series[idx]
        next_node = series[idx + 1]
        _sv_type = current_node.sv_type
        _sr = current_node.sr
        annotation_code = current_node.annotation_code
        splicing_code = current_node.splicing_code

        can_field = "NONCANONICAL"
        if splicing_code == 1:
            can_field = "CANONICAL"

        if annotation_code == 0:
            anno_field = "NEITHER"
        elif annotation_code == 1:
            anno_field = "RIGHT"
        elif annotation_code == 2:
            anno_field = "LEFT"
        else:
            anno_field = "BOTH"

        if current_node.modes is None:
            raise SystemExit from ModesNotFoundError

        if current_node.genes is None:
            raise SystemExit from GenesNotFoundError

        _mode1, _mode2 = current_node.modes
        mode1 = "MS" if _mode1 == 1 else "SM"
        mode2 = "MS" if _mode2 == 1 else "SM"
        gene1, gene2 = current_node.genes
        _chrom1, _pos1 = current_node.next_breakpoint.split(":")
        _chrom2, _pos2 = next_node.prev_breakpoint.split(":")

        sv_distance = 0
        if _sv_type == "TRA":
            sv_distance = 0
        else:
            sv_distance = abs(int(_pos1) - int(_pos2))
        _dp1 = current_node.next_breakpoint_depth
        _dp2 = next_node.prev_breakpoint_depth
        if _dp1 is None or _dp2 is None:
            _pso = 0
        else:
            _pso = _sr / (_sr + (_dp1 + _dp2) / 2)
        _strand1 = current_node.strand
        _strand2 = next_node.strand

        ref_allele = "."
        alt_allele = f"<{_sv_type}>"

        position = int(_pos1) + 1  # 1 based
        end = int(_pos2) + 1  # 1 based

        fields_dict = {
            "chrom": _chrom1,
            "pos": str(position),
            "event_id": f"HOP_{event_id}",
            "ref": ref_allele,
            "alt": alt_allele,
            "quality": ".",
            "filter": ".",
            "info": (
                f"{can_field};BOUNDARY={anno_field};SVTYPE={_sv_type};CHR2={_chrom2};END={end};DP1={_dp1};"
                f"DP2={_dp2};PSO={_pso:.3g};SVLEN={sv_distance};GENE1={gene1};GENE2={gene2};"
                f"STRAND1={_strand1};STRAND2={_strand2};MODE1={mode1};MODE2={mode2};"
                f"TRANSCRIPT_ID={series_id};SVMETHOD=ScanNLS"
            ),
            "format": "GT",
            "sample": "0/1",
        }
        series_hops_features.append(get_line_from_hop(fields_dict))
        if current_node.insertion_info and isinstance(
            current_node.insertion_info[1], NovelInsertion
        ):
            _ao = current_node.insertion_info[1].ao
            ref_allele, alt_allele = get_novel_insertion_sequence_from_node(
                current_node, reference_io
            )
            _dp = _dp1
            _af = _ao / _dp
            gene = gene1
            sv_distance = len(alt_allele)
            _sv_type = "INS"
            if annotation_code in {0, 1}:
                anno_field = "NEITHER"
            else:
                anno_field = "LEFT"
            fields_dict = {
                "chrom": _chrom1,
                "pos": str(position),
                "event_id": f"INS_{ins_id}",
                "ref": ref_allele,
                "alt": alt_allele,
                "quality": ".",
                "filter": ".",
                "info": f"{can_field};BOUNDARY={anno_field};SVTYPE={_sv_type};"
                f"CHR2={_chrom2};END={end};DP={_dp};AF={_af:.3g};"
                f"SVLEN={sv_distance};GENE={gene};STRAND={_strand1};TRANSCRIPT_ID={series_id};SVMETHOD=ScanNLS",
                "format": "GT",
                "sample": "0/1",
            }
            series_hops_features.append(get_line_from_hop(fields_dict))
            ins_id += 1

    return series_hops_features


def get_line_from_hop(fields_dict: dict) -> List[str]:
    """Get the VCF line from one hop, including INS/DEL/TDUP/IDUP/INV/TRA."""
    fields_names = [
        "chrom",
        "pos",
        "event_id",
        "ref",
        "alt",
        "quality",
        "filter",
        "info",
        "format",
        "sample",
    ]
    return [fields_dict[_name] for _name in fields_names]


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
        "+",
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
