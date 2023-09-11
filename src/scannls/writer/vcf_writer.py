"""VCF Writer class.

@Filename:    vcfWriter.py
@Time:        1/30/22 6:19 PM
"""
from __future__ import annotations

import datetime
from functools import singledispatchmethod
from pathlib import Path
from typing import IO, Any, ClassVar

from loguru import logger
from pyfaidx import Fasta, FastaNotFoundError

from scannls import MicroHomology, NovelInsertion, __version__, reverse_complement
from scannls.exception import (
    AnnotationCodeNotFoundError,
    BreakpointNotFoundError,
    GenesNotFoundError,
    ModesNotFoundError,
    SplicingCodeNotFoundError,
)
from scannls.graph import NLPath, Node  # noqa: TCH001

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

    reserved_info: ClassVar[dict[str, str]] = {
        "CANONICAL": "Flag",
        "NONCANONICAL": "Flag",
        "BOUNDARY": "String",
        "DP1": "Integer",
        "DP2": "Integer",
        "SR": "Integer",
        "OSR": "Integer",
        "PSI": "Float",
        "SVMETHOD": "String",
        "SVTYPE": "String",
        "SVLEN": "Integer",
        "CHR2": "String",
        "SVEND": "Integer",
        "END": "Integer",
        "STRAND1": "String",
        "STRAND2": "String",
        "MODE1": "String",
        "MODE2": "String",
        "GENE1": "String",
        "GENE2": "String",
        "MEGAEXON1": "String",
        "MEGAEXON2": "String",
        "HOMSEQ": "String",
        "INSSEQ": "String",
        "TRANSCRIPT_ID": "String",
        "GENE_ID": "String",
        "SR_ID": "String",
    }
    reserved_format: ClassVar[dict[str, str]] = {"GT": "String"}
    reserved_alt: ClassVar[list[str]] = [
        "DEL",
        "TDUP",
        "IDUP",
        "INV",
        "TRA",
    ]

    description: ClassVar[dict[str, str]] = {
        "CANONICAL": "Canonical splice site",
        "NONCANONICAL": "Noncanonical splice site",
        "BOUNDARY": "The coding exon boundary type of event, BOTH, LEFT, RIGHT, NEITHER.",
        "DP1": "Total read depth at the breakpoint1",
        "DP2": "Total read depth at the breakpoint2",
        "SR": "The number of support reads for the breakpoints",
        "OSR": "The number of support reads for the breakpoints before rescuer",
        "PSI": "Estimated Percent splice-in in the range (0,1], "
        "representing the percentage of NLS transcripts",
        "SVTYPE": "The type of event, DEL, TDUP, IDUP, INV, TRA.",
        "SVLEN": "Difference in length between REF and ALT alleles",
        "CHR2": "Chromosome for END coordinate in case of a translocation",
        "SVEND": "2nd position of the structural variant",  # change to SVEND in order to meet vcf standard
        "END": "A placeholder for END coordinate in case of a translocation",
        "GENE1": "Overlapped coding gene for breakpoint1",
        "GENE2": "Overlapped coding gene for breakpoint2",
        "MEGAEXON1": "ID for source mega exon",
        "MEGAEXON2": "ID for target mega exon",
        "TRANSCRIPT_ID": "Transcript ID",
        "GENE_ID": "Gene ID",
        "SR_ID": "Support read ID",
        "SVMETHOD": "Type of approach used to detect SV",
        "STRAND1": "Strand for breakpoint1",
        "STRAND2": "Strand for breakpoint2",
        "MODE1": "Mode for softclipped reads at breakpoint1",
        "MODE2": "Mode for softclipped reads at breakpoint2",
        "GT": "Genotype",
        "INSSEQ": "MicroInsertion sequence",
        "HOMSEQ": "MicroHomology sequence",
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
        bam_header: dict[str, Any],
    ) -> None:
        """Initialize VCFWriter object."""
        super().__init__(file_path)
        self.reference = Path(reference)
        if not self.reference.exists():
            raise FastaNotFoundError
        self.reference_io: Fasta = Fasta(reference, sequence_always_upper=True)
        self.nlpath_id: int = 1
        self.bam_header = bam_header
        self.sample_name: str = self.file_path.stem
        self.hops_feature_in_series_list: list[Any] = []
        self.cluster_id: str = str(1)

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, fields: list[str], delimiter: str = "\t") -> str:
        """Formatter for writing data."""
        if fields is None or len(fields) != VCFWriter.num_fields:
            logger.warning(
                f"{self.__class__.__name__}: Number of fields is not equal to 10.",
            )
        return delimiter.join(fields) + "\n"

    def open(self, mode: str = "w") -> IO:
        """Open file."""
        if self.is_opened:
            logger.warning(f"{self.__class__.__name__}: File is already opened.")
        self.io = self.file_path.open(mode)  # add asyncio support
        if hasattr(self, "write_header"):
            self.write_header()
        return self.io

    def close(self) -> None:
        """Close file."""
        if self.is_opened:
            self.write_data_helper()
            logger.trace(f"{self.__class__.__name__}: Closing file.")
            self.io.close()  # type: ignore
            self.io = None

    def write_line(self, line: str) -> None:
        """Write line to file."""
        if self.is_opened:
            self.io.write(line)  # type: ignore
        else:
            logger.warning(f"{self.__class__.__name__}: File is not opened.")

    def write_header(self) -> None:
        """Write header to VCF file."""
        self.write_line(self.header)

    @singledispatchmethod
    def write_data(self, data_object: Any, object_id: str) -> None:  # type: ignore
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: NLPath, cluster_id: str) -> None:
        """Write Series to VCF file.

        :param data_object: Series to write to file.
        """
        if cluster_id != self.cluster_id:
            # next clique
            # write all features in the clique
            self.write_data_helper()
            # clear all features in the clique, start a new clique
            self.hops_feature_in_series_list.clear()
            self.cluster_id = cluster_id

        if len(data_object.nodes) == 0:
            logger.warning(
                f"{self.__class__.__name__}: No nodes to write to VCF file in Clique {cluster_id} Series.",
            )
        # hop_vcf_feature is a dict, key: sv_type, chrom1|pos1, chrom2|pos2
        for _hop_vcf_feature in get_vcf_features_from_nlpath(
            data_object,
            self.nlpath_id,
            cluster_id,
        ):
            self.hops_feature_in_series_list.append(_hop_vcf_feature)
        self.nlpath_id += 1  # series/transcript id

    def write_data_helper(self) -> None:
        """Write series data for every clique."""
        out_vcf_dict = {}
        for hop_feature in self.hops_feature_in_series_list:
            type_position_key = [*hop_feature][0]
            if type_position_key not in out_vcf_dict:
                out_vcf_dict[type_position_key] = hop_feature[type_position_key]
            else:
                # multiple transcripts go through the same one hop
                out_vcf_dict[type_position_key][
                    "TRANSCRIPT_ID"
                ] += f',{hop_feature[type_position_key]["TRANSCRIPT_ID"]}'

        for _idx, _out_vcf_hop in enumerate(out_vcf_dict, 1):
            hop_vcf_feature = vcf_feature_transformer(out_vcf_dict[_out_vcf_hop], _idx)
            self.write_line(self.formatter(hop_vcf_feature))

    @property
    def header(self) -> str:
        """VCF header provides metadata describing the body of the file."""

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
        header_lines += self.get_contigs()

        for _id in VCFWriter.reserved_info:
            _number: str | int = 0 if VCFWriter.reserved_info[_id] == "Flag" else 1
            if _id == "TRANSCRIPT_ID":
                _number = "."
            header_lines.append(
                f"##INFO=<ID={_id},Number={_number},Type={VCFWriter.reserved_info[_id]},"
                f'Description="{VCFWriter.description[_id]}">',
            )

        for _id in VCFWriter.reserved_format:
            header_lines.append(
                f"##FORMAT=<ID={_id},Number=1,Type={VCFWriter.reserved_format[_id]},"
                f'Description="{VCFWriter.description[_id]}">',
            )

        for _id in VCFWriter.reserved_alt:
            header_lines.append(
                f'##ALT=<ID={_id},Description="{VCFWriter.description[_id]}">',
            )
        header_lines.append(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}",
        )

        return "\n".join(header_lines) + "\n"

    def get_contigs(self) -> list[str]:
        """Get contigs from BAM file header."""
        return [
            f"##contig=<ID={contig_dict['SN']},length={contig_dict['LN']}>"
            for contig_dict in self.bam_header["SQ"]
        ]


def obtain_reference_from_bam_header(bam_header: dict[str, Any]) -> str:
    """Obtain reference info from BAM header."""
    aligners = {
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
        "NGMLR",
        "minimap2",
    }
    avail_aligners = {x.upper() for x in aligners}
    for item in bam_header.get("PG", {}):
        if item["ID"].upper() in avail_aligners:
            return item["CL"]
    return "Unknown"


def get_vcf_features_from_nlpath(
    nlpath: NLPath,
    nlpath_id: int,
    cluster_id: str,
):
    """Obtain hop vcf features from one series."""
    path_hops_features = []

    can_field_dict = {0: "NONCANONICAL", 1: "CANONICAL"}
    anno_field_dict = {0: "NEITHER", 1: "RIGHT", 2: "LEFT"}
    for event_id, current_node in enumerate(nlpath.nodes[:-1], 1):
        current_edge = nlpath.next_edge(current_node, event_id - 1)
        next_node = nlpath[event_id]

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

        if current_edge is None:
            raise BreakpointNotFoundError(current_node.query_name)

        _chrom1, _pos1 = current_edge.break_point1.to_tuple()
        _chrom2, _pos2 = current_edge.break_point2.to_tuple()

        microhomology_sequence = ""
        microinsertion_sequence = ""
        if current_edge.insertion_info:
            if isinstance(current_edge.insertion_info[1], NovelInsertion):
                insertion = current_edge.insertion_info[1]
                microinsertion_sequence = obtain_sequence_from_insertion(
                    insertion,
                    current_node,
                )
            elif isinstance(current_edge.insertion_info[1], MicroHomology):
                microhomology = current_edge.insertion_info[1]
                microhomology_sequence = obtain_sequence_from_insertion(
                    microhomology,
                    current_node,
                )

        # correct the breakpoint position in order to obtain a precise "sv_distance"
        _pos1 = (
            _pos1 - len(microhomology_sequence)
            if current_node.strand == "+"
            else _pos1 + len(microhomology_sequence)
        )

        sv_distance = (
            abs(_pos1 - _pos2) if not current_edge.variation_type.is_tra() else 0
        )
        _dp1 = (
            0
            if current_edge.break_point1.depth is None
            else current_edge.break_point1.depth
        )

        _dp2 = (
            0
            if current_edge.break_point2.depth is None
            else current_edge.break_point2.depth
        )

        _pso = (
            0
            if _dp1 == 0 or _dp2 == 0
            else current_edge.sr / (current_edge.sr + (_dp1 + _dp2) / 2)
        )

        path_hops_features.append(
            {
                f"{current_edge.variation_type}_{_chrom1}|{_pos1 + 1}"
                f"_{_chrom2}|{_pos2 + 1}": {
                    "CHROM": _chrom1,
                    "POS": f"{_pos1 + 1}",
                    "REF": ".",
                    "ALT": f"<{current_edge.variation_type}>",
                    "SVTYPE": current_edge.variation_type,
                    "SR": current_edge.sr,
                    "OSR": current_edge.original_sr,
                    "CAN": can_field,
                    "BOUNDARY": anno_field,
                    "CHR2": _chrom2,
                    "SVEND": f"{_pos2 + 1}",
                    "DP1": f"{_dp1}",
                    "DP2": f"{_dp2}",
                    "PSI": f"{_pso:.3g}",
                    "SVLEN": f"{sv_distance}",
                    "GENE1": f"{gene1}",
                    "GENE2": f"{gene2}",
                    "MEGAEXON1": f"{current_node.trace_id}",
                    "MEGAEXON2": f"{next_node.trace_id}",
                    "STRAND1": f"{current_node.strand}",
                    "STRAND2": f"{next_node.strand}",
                    "MODE1": f"{mode1}",
                    "MODE2": f"{mode2}",
                    "TRANSCRIPT_ID": f"{nlpath_id}",
                    "GENE_ID": f"{cluster_id}",
                    "SR_ID": f"{','.join(current_edge.read_ids)}",
                    "SVMETHOD": "ScanNLS",
                    "HOMSEQ": microhomology_sequence if microhomology_sequence else ".",
                    "INSSEQ": microinsertion_sequence
                    if microinsertion_sequence
                    else ".",
                },
            },
        )

    return path_hops_features


def vcf_feature_transformer(feature_dict: dict[str, str], idx: int) -> list[str]:
    """VCF feature transformer."""
    info_field = (
        f'{feature_dict["CAN"]};BOUNDARY={feature_dict["BOUNDARY"]};'
        f'SVTYPE={feature_dict["SVTYPE"]};SR={feature_dict["SR"]};OSR={feature_dict["OSR"]};'
        f'CHR2={feature_dict["CHR2"]};SVEND={feature_dict["SVEND"]};DP1={feature_dict["DP1"]};'
        f'DP2={feature_dict["DP2"]};PSI={feature_dict["PSI"]};SVLEN={feature_dict["SVLEN"]};'
        f'GENE1={feature_dict["GENE1"]};GENE2={feature_dict["GENE2"]};'
        f'MEGAEXON1={feature_dict["MEGAEXON1"]};MEGAEXON2={feature_dict["MEGAEXON2"]};'
        f'STRAND1={feature_dict["STRAND1"]};STRAND2={feature_dict["STRAND2"]};'
        f'MODE1={feature_dict["MODE1"]};MODE2={feature_dict["MODE2"]};'
        f'HOMSEQ={feature_dict["HOMSEQ"]};INSSEQ={feature_dict["INSSEQ"]};'
        f'TRANSCRIPT_ID={feature_dict["TRANSCRIPT_ID"]};GENE_ID={feature_dict["GENE_ID"]};'
        f'SR_ID={feature_dict["SR_ID"]};SVMETHOD={feature_dict["SVMETHOD"]}'
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


def obtain_sequence_from_insertion(
    insertion: NovelInsertion | MicroHomology,
    node: Node,
) -> str:
    """Get sequence of novel insertion or microhomology of a node.

    :param insertion:
    :param node: Node and InsertionType
    :return: sequence of novel insertion or microhomology
    """
    # positive strand sequence for novel insertion
    novel_insertion_sequence = insertion.query_sequence

    if not novel_insertion_sequence:
        return ""

    return (
        novel_insertion_sequence
        if node.strand == "+"
        else reverse_complement(novel_insertion_sequence)
    )
