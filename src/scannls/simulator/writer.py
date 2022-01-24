# !/usr/bin/env python
"""Writer for simulator.

@license:     MIT Licence
@Time:        1/24/22 9:17 AM
"""
import datetime
from functools import singledispatchmethod
from pathlib import Path
from typing import Any
from typing import IO
from typing import List
from typing import Optional

from .oneHop import MetaExon
from scannls import __version__
from scannls import BreakpointNotFoundError
from scannls import LoggerType


class SimVCFWriter:
    """Writer for VCF files for Simulated data.

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

    num_fields = 10

    reserved_info = {
        "SVMETHOD": "String",
        "SVTYPE": "String",
        "STRAND1": "String",
        "STRAND2": "String",
        "SVLEN": "Integer",
        "CHR2": "String",
        "END": "Integer",
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

    def __init__(
        self,
        file_path: str,
        output_prefix: str,
        logger: LoggerType,
    ) -> None:
        """Initialize SimVCFWriter object."""
        self.file_path = Path(file_path)
        self.logger = logger
        if self.file_path.exists():
            self.logger.warning(f"{self.file_path} exists, will be overwritten.")
        self.id = 1
        self.sample_name = output_prefix
        self.io: Optional[IO] = None

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, fields: List[str], delimiter: str = "\t") -> str:
        """Formatter for writing data."""
        if fields is None or len(fields) != SimVCFWriter.num_fields:
            self.logger.warning(
                f"{self.__class__.__name__}: Number of fields is not equal to 10."
            )
            raise SystemExit
        return delimiter.join(fields) + "\n"

    def open(self, mode: str = "w") -> IO:
        """Open file."""
        if self.is_opened:
            self.logger.warning(f"{self.__class__.__name__}: File is already opened.")
        self.io = open(self.file_path, mode)  # noqa
        return self.io

    def close(self) -> None:
        """Close file."""
        if self.is_opened:
            self.io.close()  # type: ignore
            self.io = None  # type: ignore

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
        if len(data_object) == 0:
            self.logger.warning(
                f"{self.__class__.__name__}: No nodes to write to VCF file."
            )
        self.logger.trace(f"{self.__class__.__name__}: Writing metaexons to VCF file.")
        for hop_vcf_feature in get_vcf_features_from_metaexons(data_object, self.id):
            self.write_line(self.formatter(hop_vcf_feature))
        self.id += 1

    @property
    def header(self) -> str:
        """VCF header provides metadata describing the body of the file."""
        # Metadata parsers/constants

        date = datetime.datetime.today().strftime("%Y%m%d")
        source = f"ScanNLS v{__version__}"

        header_lines = [
            "##fileformat=VCFv4.3",
            f"##fileDate={date}",
            f"##source={source}",
        ]

        for _id in SimVCFWriter.reserved_info:
            _number = 0 if SimVCFWriter.reserved_info[_id] == "Flag" else 1
            header_lines.append(
                f"##INFO=<ID={_id},Number={_number},Type={SimVCFWriter.reserved_info[_id]},"
                f'Description="{SimVCFWriter.description[_id]}">'
            )

        for _id in SimVCFWriter.reserved_format:
            header_lines.append(
                f"##FORMAT=<ID={_id},Number=1,Type={SimVCFWriter.reserved_format[_id]},"
                f'Description="{SimVCFWriter.description[_id]}">'
            )

        for _id in SimVCFWriter.reserved_alt:
            header_lines.append(
                f'##ALT=<ID={_id},Description="{SimVCFWriter.description[_id]}">'
            )
        header_lines.append(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}"
        )

        return "\n".join(header_lines) + "\n"


def get_vcf_features_from_metaexons(
    series: List[MetaExon],
    series_id: int,
) -> List[List[str]]:
    """Obtain hop vcf features from one list of metaexons."""
    series_hops_features = []

    for event_id, current_node in enumerate(series[:-1], 1):
        next_node = series[event_id]

        _chrom1 = current_node.chrom
        _chrom2 = next_node.chrom

        _strand1, _pos1 = current_node.strand, current_node.p3_pos
        _strand2, _pos2 = next_node.strand, next_node.p5_pos

        if _pos1 is None or _pos2 is None:
            raise SystemExit from BreakpointNotFoundError

        sv_distance = abs(_pos1 - _pos2) if current_node.nls_type != "TRA" else 0

        series_hops_features.append(
            [
                _chrom1,
                f"{int(_pos1) + 1}",
                f"HOP_{event_id}",
                ".",
                f"<{current_node.nls_type}>",
                ".",
                ".",
                (
                    f"SVTYPE={current_node.nls_type};"
                    f"CHR2={_chrom2};END={int(_pos2) + 1};"
                    f"SVLEN={sv_distance};"
                    f"STRAND1={_strand1};STRAND2={_strand2};"
                    f"TRANSCRIPT_ID={series_id};SVMETHOD=ScanNLS_Simulator"
                ),
                "GT",
                "0/1",
            ]
        )
        # current_node insertion_seq: #TODO
    return series_hops_features
