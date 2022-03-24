# !/usr/bin/env python
"""VCF Writer class.

@Filename:    vcfWriter.py
@license:     MIT Licence
@Time:        1/30/22 6:19 PM
"""
import datetime
from functools import singledispatchmethod
from pathlib import Path
from typing import Any
from typing import Dict
from typing import IO
from typing import List
from typing import Tuple

from pyfaidx import Fasta  # type: ignore
from pyfaidx import FastaNotFoundError  # type: ignore

from .. import __version__
from ..basicClass import Node
from ..basicClass import NovelInsertion
from ..basicClass import reverse_complement
from ..basicClass import Series
from ..exception import AnnotationCodeNotFoundError
from ..exception import BreakpointNotFoundError
from ..exception import ExonsNotFoundError
from ..exception import GenesNotFoundError
from ..exception import ModesNotFoundError
from ..exception import SplicingCodeNotFoundError
from ..type import LoggerType
from .writer import Writer


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

    num_fields = 10

    reserved_info = {
        "CANONICAL": "Flag",
        "NONCANONICAL": "Flag",
        "BOUNDARY": "String",
        "DP": "Integer",
        "DP1": "Integer",
        "DP2": "Integer",
        "SR": "Integer",
        "PSO": "Float",
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
        "TRANSCRIPT_ID": "String",
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
        reference: str,
        bam_header: Dict[str, Any],
        logger: LoggerType,
    ) -> None:
        """Initialize VCFWriter object."""
        super().__init__(file_path, logger)
        self.reference = Path(reference)
        if not self.reference.exists():
            raise FastaNotFoundError
        self.reference_io: Fasta = Fasta(reference, sequence_always_upper=True)
        self.series_id: int = 1
        self.bam_header = bam_header
        self.sample_name: str = self.file_path.stem
        self.hops_feature_in_series_list: List[Any] = []
        self.clique_id: int = 1

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, fields: List[str], delimiter: str = "\t") -> str:
        """Formatter for writing data."""
        if fields is None or len(fields) != VCFWriter.num_fields:
            self.logger.warning(
                f"{self.__class__.__name__}: Number of fields is not equal to 10."
            )
        return delimiter.join(fields) + "\n"

    def open(self, mode: str = "w") -> IO:
        """Open file."""
        if self.is_opened:
            self.logger.warning(f"{self.__class__.__name__}: File is already opened.")
        self.io = self.file_path.open(mode)  # add asyncio support
        if hasattr(self, "write_header"):
            self.write_header()
        return self.io

    def close(self) -> None:
        """Close file."""
        if self.is_opened:
            self.write_data_helper()
            self.logger.trace(f"{self.__class__.__name__}: Closing file.")
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
    def write_data(self, data_object: Any, object_id: int) -> None:
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: Series, clique_id: int) -> None:
        """Write Series to VCF file.

        :param data_object: Series to write to file.
        """
        if clique_id != self.clique_id:
            # next clique
            # write all features in the clique
            self.write_data_helper()
            # clear all features in the clique, start a new clique
            self.hops_feature_in_series_list.clear()
            self.clique_id = clique_id

        if len(data_object.nodes) == 0:
            self.logger.warning(
                f"{self.__class__.__name__}: No nodes to write to VCF file in Clique {clique_id} Series."
            )
        # hop_vcf_feature is a dict, key: sv_type, chrom1|pos1, chrom2|pos2
        for _hop_vcf_feature in get_vcf_features_from_series(
            data_object, self.series_id, self.reference_io
        ):
            self.hops_feature_in_series_list.append(_hop_vcf_feature)
        self.series_id += 1  # series/transcript id

    def write_data_helper(self) -> None:
        """Write series data for every clique."""
        out_vcf_dict = {}
        for hop_feature in self.hops_feature_in_series_list:
            type_position_key = [*hop_feature][0]
            if type_position_key not in out_vcf_dict:
                out_vcf_dict[type_position_key] = hop_feature[type_position_key]
            else:
                if (
                    hop_feature[type_position_key]["SR"]
                    >= out_vcf_dict[type_position_key]["SR"]
                ):
                    out_vcf_dict[type_position_key]["SR"] = hop_feature[
                        type_position_key
                    ]["SR"]

                    if out_vcf_dict[type_position_key]["SVTYPE"] == "INS":
                        out_vcf_dict[type_position_key]["AF"] = hop_feature[
                            type_position_key
                        ]["AF"]
                    else:
                        out_vcf_dict[type_position_key]["PSO"] = hop_feature[
                            type_position_key
                        ]["PSO"]

                    out_vcf_dict[type_position_key][
                        "TRANSCRIPT_ID"
                    ] += f',{hop_feature[type_position_key]["TRANSCRIPT_ID"]}'

        for _idx, _out_vcf_hop in enumerate(out_vcf_dict, 1):
            hop_vcf_feature = vcf_feature_transformer(out_vcf_dict[_out_vcf_hop], _idx)
            self.write_line(self.formatter(hop_vcf_feature))

    @property
    def header(self) -> str:
        """VCF header provides metadata describing the body of the file."""
        # Metadata parsers/constants

        date = datetime.datetime.today().strftime("%Y%m%d")
        source = f"ScanNLS v{__version__}"
        reference = (
            f"<CMD={obtain_reference_from_bam_header(self.bam_header)},"
            'Description="Alignment parameters">'
        )

        header_lines = [
            "##fileformat=VCFv4.3",
            f"##fileDate={date}",
            f"##source={source}",
            f"##reference={reference}",
        ]

        for _id in VCFWriter.reserved_info:
            _number = 0 if VCFWriter.reserved_info[_id] == "Flag" else 1
            header_lines.append(
                f"##INFO=<ID={_id},Number={_number},Type={VCFWriter.reserved_info[_id]},"
                f'Description="{VCFWriter.description[_id]}">'
            )

        for _id in VCFWriter.reserved_format:
            header_lines.append(
                f"##FORMAT=<ID={_id},Number=1,Type={VCFWriter.reserved_format[_id]},"
                f'Description="{VCFWriter.description[_id]}">'
            )

        for _id in VCFWriter.reserved_alt:
            header_lines.append(
                f'##ALT=<ID={_id},Description="{VCFWriter.description[_id]}">'
            )
        header_lines.append(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}"
        )

        return "\n".join(header_lines) + "\n"


def obtain_reference_from_bam_header(bam_header: Dict[str, Any]) -> str:
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
    for item in bam_header.get("PG", {}):
        if item["ID"].upper() in avail_aligners:
            return item["CL"]
    return "Unknown"


def get_vcf_features_from_series(
    series: Series,
    series_id: int,
    reference_io: Fasta,
) -> Any:
    """Obtain hop vcf features from one series."""
    series_hops_features = []

    can_field_dict = {0: "NONCANONICAL", 1: "CANONICAL"}
    anno_field_dict = {0: "NEITHER", 1: "RIGHT", 2: "LEFT"}
    for event_id, current_node in enumerate(series.nodes[:-1], 1):
        next_node = series[event_id]

        if current_node.splicing_code is None:
            raise SplicingCodeNotFoundError(current_node.query_name)
        can_field = can_field_dict[current_node.splicing_code]  # type: ignore

        if current_node.annotation_code is None:
            raise AnnotationCodeNotFoundError(current_node.query_name)
        anno_field = anno_field_dict.get(current_node.annotation_code, "BOTH")  # type: ignore

        if current_node.genes is None:
            raise GenesNotFoundError(current_node.query_name)
        gene1, gene2 = current_node.genes

        if current_node.modes is None:
            raise ModesNotFoundError(current_node.query_name)
        _mode1, _mode2 = current_node.modes
        mode1 = "MS" if _mode1 == 1 else "SM"
        mode2 = "MS" if _mode2 == 1 else "SM"

        if current_node.next_breakpoint is None:
            raise BreakpointNotFoundError(current_node.query_name)
        if next_node.prev_breakpoint is None:
            raise BreakpointNotFoundError(next_node.query_name)
        _chrom1, _pos1 = current_node.next_breakpoint.split(":")
        _chrom2, _pos2 = next_node.prev_breakpoint.split(":")

        sv_distance = (
            abs(int(_pos1) - int(_pos2)) if current_node.sv_type != "TRA" else 0
        )
        _dp1 = (
            0
            if current_node.next_breakpoint_depth is None
            else current_node.next_breakpoint_depth
        )
        _dp2 = (
            0
            if next_node.prev_breakpoint_depth is None
            else next_node.prev_breakpoint_depth
        )
        _pso = (
            0
            if _dp1 == 0 or _dp2 == 0
            else current_node.sr / (current_node.sr + (_dp1 + _dp2) / 2)
        )

        series_hops_features.append(
            {
                f"{current_node.sv_type}_{_chrom1}|{int(_pos1) + 1}"
                f"_{_chrom2}|{int(_pos2) + 1}": {
                    "CHROM": _chrom1,
                    "POS": f"{int(_pos1) + 1}",
                    "REF": ".",
                    "ALT": f"<{current_node.sv_type}>",
                    "SVTYPE": current_node.sv_type,
                    "SR": current_node.sr,
                    "CAN": can_field,
                    "BOUNDARY": anno_field,
                    "CHR2": _chrom2,
                    "END": f"{int(_pos2) + 1}",
                    "DP1": f"{_dp1}",
                    "DP2": f"{_dp2}",
                    "PSO": f"{_pso:.3g}",
                    "SVLEN": f"{sv_distance}",
                    "GENE1": f"{gene1}",
                    "GENE2": f"{gene2}",
                    "STRAND1": f"{current_node.strand}",
                    "STRAND2": f"{next_node.strand}",
                    "MODE1": f"{mode1}",
                    "MODE2": f"{mode2}",
                    "TRANSCRIPT_ID": f"{series_id}",
                    "SVMETHOD": "ScanNLS",
                }
            }
        )
        if current_node.insertion_info and isinstance(
            current_node.insertion_info[1], NovelInsertion
        ):
            insertion = current_node.insertion_info[1]
            ref_allele, alt_allele = get_vcf_features_from_insertion(
                insertion, current_node, reference_io
            )
            _af = insertion.ao / _dp1
            sv_distance = len(alt_allele)
            _sv_type = "INS"
            anno_field = "NEITHER" if current_node.annotation_code in {0, 1} else "LEFT"
            series_hops_features.append(
                {
                    f"{_sv_type}_{_chrom1}|{int(_pos1) + 1}"
                    f"_{_chrom1}|{int(_pos1) + 1}": {
                        "CHROM": _chrom1,
                        "POS": f"{int(_pos1) + 1}",
                        "REF": f"{ref_allele}",
                        "ALT": f"{alt_allele}",
                        "SVTYPE": _sv_type,
                        "SR": insertion.ao,
                        "CAN": can_field,
                        "BOUNDARY": anno_field,
                        "CHR2": _chrom1,
                        "END": f"{int(_pos1) + 1}",
                        "DP": f"{_dp1}",
                        "AF": f"{_af:.3g}",
                        "SVLEN": f"{sv_distance}",
                        "GENE": f"{gene1}",
                        "STRAND": f"{current_node.strand}",
                        "TRANSCRIPT_ID": f"{series_id}",
                        "SVMETHOD": "ScanNLS",
                    }
                }
            )
    return series_hops_features


def vcf_feature_transformer(feature_dict: Dict[str, str], idx: int) -> List[str]:
    """VCF feature transformer."""
    if feature_dict["SVTYPE"] == "INS":
        info_field = (
            f'{feature_dict["CAN"]};BOUNDARY={feature_dict["BOUNDARY"]};'
            f'SVTYPE={feature_dict["SVTYPE"]};SR={feature_dict["SR"]};'
            f'CHR2={feature_dict["CHR2"]};END={feature_dict["END"]};DP={feature_dict["DP"]};'
            f'AF={feature_dict["AF"]};SVLEN={feature_dict["SVLEN"]};'
            f'GENE={feature_dict["GENE"]};'
            f'STRAND={feature_dict["STRAND"]};'
            f'TRANSCRIPT_ID={feature_dict["TRANSCRIPT_ID"]};SVMETHOD={feature_dict["SVMETHOD"]}'
        )
    else:
        info_field = (
            f'{feature_dict["CAN"]};BOUNDARY={feature_dict["BOUNDARY"]};'
            f'SVTYPE={feature_dict["SVTYPE"]};SR={feature_dict["SR"]};'
            f'CHR2={feature_dict["CHR2"]};END={feature_dict["END"]};DP1={feature_dict["DP1"]};'
            f'DP2={feature_dict["DP2"]};PSO={feature_dict["PSO"]};SVLEN={feature_dict["SVLEN"]};'
            f'GENE1={feature_dict["GENE1"]};GENE2={feature_dict["GENE2"]};'
            f'STRAND1={feature_dict["STRAND1"]};STRAND2={feature_dict["STRAND2"]};'
            f'MODE1={feature_dict["MODE1"]};MODE2={feature_dict["MODE2"]};'
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


def get_vcf_features_from_insertion(
    insertion: NovelInsertion, node: Node, reference_io: Fasta
) -> Tuple[str, str]:
    """Get novel insertion sequence of a node.

    :param insertion:
    :param node: Node and InsertionType
    :param reference_io: ReferenceIO object
    :return: sequence of exon for one node
    """
    # positive strand sequence for novel insertion
    novel_insertion_sequence = insertion.query_sequence

    if node.exons is None:
        raise ExonsNotFoundError(node.query_name)

    _pos = node.exons[-1][1] if node.strand == "+" else node.exons[0][0]

    ref_allele = reference_io.get_seq(node.chrom, _pos + 1, _pos + 1).seq  # 1-based

    if not novel_insertion_sequence:
        return ref_allele, ""

    alt_allele = (
        novel_insertion_sequence
        if node.strand == "+"
        else reverse_complement(novel_insertion_sequence)
    )

    return ref_allele, alt_allele
