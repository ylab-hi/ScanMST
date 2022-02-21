# !/usr/bin/env python
"""Writer for simulator.

@license:     MIT Licence
@Time:        1/24/22 9:17 AM
"""
import datetime
from functools import singledispatchmethod
from pathlib import Path
from typing import Any
from typing import Dict
from typing import IO
from typing import List
from typing import Optional

from .helper import real_path
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
        logger: LoggerType,
    ) -> None:
        """Initialize SimVCFWriter object."""
        file_path = real_path(file_path)
        self.file_path = Path(file_path)
        self.logger = logger
        if self.file_path.exists():
            self.logger.warning(f"{self.file_path} exists, will be overwritten.")
        self.id = 1
        self.sample_name = self.file_path.stem
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
        if isinstance(data_object, list):
            if len(data_object) == 0:
                self.logger.warning(
                    f"{self.__class__.__name__}: No nodes to write to VCF file."
                )
            self.logger.trace(
                f"{self.__class__.__name__}: Writing metaexons to VCF file."
            )
            for hop_vcf_feature in get_vcf_features_from_metaexons(
                data_object, self.id
            ):
                self.write_line(self.formatter(hop_vcf_feature))
            self.id += 1
        elif isinstance(data_object, dict):
            self.logger.trace(
                f"{self.__class__.__name__}: Writing metaexons to VCF file."
            )

            hops_feature_in_trx_list = []
            for trx_id in data_object:
                # _hop_vcf_feature is a dict, key: sv_type, chrom1|pos1, chrom2|pos2
                for _hop_vcf_feature in _get_vcf_features_from_metaexons(
                    data_object[trx_id], trx_id
                ):
                    hops_feature_in_trx_list.append(_hop_vcf_feature)

            out_vcf_dict = {}
            for hop_feature in hops_feature_in_trx_list:
                type_position_key = [*hop_feature][0]
                if type_position_key not in out_vcf_dict:
                    out_vcf_dict[type_position_key] = hop_feature[type_position_key]
                else:
                    out_vcf_dict[type_position_key][
                        "TRANSCRIPT_ID"
                    ] += f',{hop_feature[type_position_key]["TRANSCRIPT_ID"]}'

            for _idx, _out_vcf_hop in enumerate(out_vcf_dict, 1):
                hop_vcf_feature = vcf_feature_transformer(
                    out_vcf_dict[_out_vcf_hop], _idx
                )
                self.write_line(self.formatter(hop_vcf_feature))

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
            raise BreakpointNotFoundError(f"{current_node} {next_node}")

        sv_distance = abs(_pos1 - _pos2) if current_node.nls_type != "TRA" else 0
        item = [
            f"{_chrom1}",
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
        series_hops_features.append(item)
        # current_node insertion_seq: #TODO
    return series_hops_features


def _get_vcf_features_from_metaexons(
    series: List[MetaExon],
    series_id: int,
) -> List[Dict[str, Any]]:
    """Obtain hop vcf features from one list of metaexons."""
    series_hops_features = []

    for event_id, current_node in enumerate(series[:-1], 1):
        next_node = series[event_id]

        _chrom1 = current_node.chrom
        _chrom2 = next_node.chrom

        _pos1 = current_node.p3_pos
        _pos2 = next_node.p5_pos

        if _pos1 is None or _pos2 is None:
            raise BreakpointNotFoundError(f"{current_node} {next_node}")

        sv_distance = abs(_pos1 - _pos2) if current_node.nls_type != "TRA" else 0
        item = {
            f"{current_node.nls_type}_{_chrom1}|{int(_pos1) + 1}"
            f"_{_chrom2}|{int(_pos2) + 1}": {
                "CHROM": _chrom1,
                "POS": f"{int(_pos1) + 1}",
                "REF": ".",
                "ALT": f"<{current_node.nls_type}>",
                "SVTYPE": current_node.nls_type,
                "CHR2": _chrom2,
                "END": f"{int(_pos2) + 1}",
                "SVLEN": f"{sv_distance}",
                "STRAND1": f"{current_node.strand}",
                "STRAND2": f"{next_node.strand}",
                "TRANSCRIPT_ID": f"{series_id}",
                "SVMETHOD": "ScanNLS_Simulator",
            }
        }
        series_hops_features.append(item)
        # current_node insertion_seq: #TODO
    return series_hops_features


def vcf_feature_transformer(feature_dict: dict, idx: int) -> List[str]:
    """VCF feature transformer."""
    if feature_dict["SVTYPE"] == "INS":
        info_field = (
            f'SVTYPE={feature_dict["SVTYPE"]};'
            f'CHR2={feature_dict["CHR2"]};END={feature_dict["END"]};'
            f'SVLEN={feature_dict["SVLEN"]};'
            f'STRAND={feature_dict["STRAND"]};'
            f'TRANSCRIPT_ID={feature_dict["TRANSCRIPT_ID"]};SVMETHOD={feature_dict["SVMETHOD"]}'
        )
    else:
        info_field = (
            f'SVTYPE={feature_dict["SVTYPE"]};'
            f'CHR2={feature_dict["CHR2"]};END={feature_dict["END"]};'
            f'SVLEN={feature_dict["SVLEN"]};'
            f'STRAND1={feature_dict["STRAND1"]};STRAND2={feature_dict["STRAND2"]};'
            f'TRANSCRIPT_ID={feature_dict["TRANSCRIPT_ID"]};SVMETHOD={feature_dict["SVMETHOD"]}'
        )

    return [
        feature_dict["CHROM"],
        feature_dict["POS"],
        str(idx),
        feature_dict["REF"],
        feature_dict["ALT"],
        ".",
        ".",
        info_field,
        "GT",
        "0/1",
    ]


class GTFWriter:
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
        file_path = real_path(file_path)
        self.file_path = Path(file_path)
        self.logger = logger
        if self.file_path.exists():
            self.logger.warning(f"{self.file_path} exists, will be overwritten.")
        self.id = 1
        self.io: Optional[IO] = None

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
    def write_data(self, data_object: Any) -> None:
        """Write data to file.

        :param: data_object: Data to write to file.
        """
        if len(data_object) == 0:
            self.logger.warning(
                f"{self.__class__.__name__}: No nodes to write to file."
            )
        for metaexon_gtf_feature in get_nodes_gtf_features_from_metaexons(
            data_object, self.id
        ):
            self.write_line(self.formatter(metaexon_gtf_feature))
        self.id += 1


def get_nodes_gtf_features_from_metaexons(
    series: List[MetaExon], series_id: int
) -> List[List[str]]:
    """Get_nodes_gtf_features_from_metaexons."""
    series_gtf_features = []
    for metaexon_id, metaexon in enumerate(series, 1):
        series_gtf_features.extend(
            get_gtf_features_from_metaexon(metaexon, series_id, metaexon_id)
        )
    return series_gtf_features


def get_gtf_features_from_metaexon(metaexon, series_id, metaexon_id) -> List[List[str]]:
    """Get_gtf_features_from_metaexon."""
    exons = metaexon.exons[::-1] if metaexon.strand == "-" else metaexon.exons
    metaexons_gtf_features = []

    for index, exon in enumerate(exons, 1):
        info = [
            f'transcript_id "{series_id:0>6}"; '
            f'metaexon_id "{metaexon_id:0>3}"; '
            f'exon_id "{index:0>3}";'
        ]
        metaexons_gtf_features.append(
            [
                f"{exon.chrom}",
                "exon",
                "scannls-simulator",
                f"{exon.start + 1}",
                f"{exon.end}",
                ".",
                f"{metaexon.strand}",
                ".",
            ]
            + info
        )
    return metaexons_gtf_features
