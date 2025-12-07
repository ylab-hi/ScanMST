"""VCF Writer class."""

from __future__ import annotations

import datetime
from functools import singledispatchmethod
from pathlib import Path
from typing import IO, Any, ClassVar

from loguru import logger
from pyfaidx import Fasta, FastaNotFoundError

from scanmst.base import MicroHomology, NovelInsertion, reverse_complement
from scanmst.exception import (
    BreakpointNotFoundError,
)
from scanmst.graph import NLPath, Node, Edge
from scanmst.utils import determine_mst_link_type

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
        "DP1": "Integer",
        "DP2": "Integer",
        "SR": "Integer",
        "OSR": "Integer",
        "PSI": "Float",
        "SVMETHOD": "String",
        "SVTYPE": "String",
        "LINKTYPE": "String",
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
        "SEGMENT1": "String",
        "SEGMENT2": "String",
        "HOMSEQ": "String",
        "INSSEQ": "String",
        "TRANSCRIPT_ID": "String",
        "GENE_ID": "String",
        "SR_ID": "String",
    }
    reserved_format: ClassVar[dict[str, str]] = {"GT": "String"}
    reserved_alt: ClassVar[list[str]] = [
        "ICRL",
        "ICTL",
        "ITPL",
        "ITTL",
    ]

    description: ClassVar[dict[str, str]] = {
        "CANONICAL": "Canonical splice site",
        "NONCANONICAL": "Noncanonical splice site",
        "DP1": "Total read depth at the breakpoint1",
        "DP2": "Total read depth at the breakpoint2",
        "SR": "The number of support reads for the breakpoints",
        "OSR": "The number of support reads for the breakpoints before rescuer",
        "PSI": "Estimated Percent splice-in in the range (0,1], representing the percentage of MSTs",
        "LINKTYPE": "The type of link, ICRL, ICTL, ITPL, ITTL.",
        "SVLEN": "Difference in length between REF and ALT alleles",
        "CHR2": "Chromosome for END coordinate in case of a translocation",
        "SVEND": "2nd position of the structural variant",  # change to SVEND in order to meet vcf standard
        "END": "A placeholder for END coordinate in case of a translocation",
        "SVTYPE": "The matching SV type, TDUP, INV, TRA.",
        "GENE1": "Overlapped coding gene for breakpoint1",
        "GENE2": "Overlapped coding gene for breakpoint2",
        "SEGMENT1": "ID for source mega exon",  # Given multiple transcripts, there may be multiple megaexons
        "SEGMENT2": "ID for target mega exon",
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
        "ICRL": "Intra-Chromosomal Reverse Link",
        "ICTL": "Intra-Chromosomal Trans-strand Link",
        "ITPL": "InTer-chromosomal Parallel Link",
        "ITTL": "InTer-chromosomal Trans-strand Link",
    }

    def __init__(
        self,
        file_path: str,
        rescue_sr: bool,
        reference: str,
        bam_header: dict[str, Any],
    ) -> None:
        """Initialize VCFWriter object."""
        super().__init__(file_path)
        self.rescue_sr = rescue_sr
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
            data_object, self.rescue_sr, cluster_id=cluster_id
        ):
            self.hops_feature_in_series_list.append(_hop_vcf_feature)

    def write_data_helper(self) -> None:
        """Write series data for every clique."""
        out_vcf_dict = {}
        for hop_feature in self.hops_feature_in_series_list:
            type_position_key = next(iter(hop_feature))
            if type_position_key not in out_vcf_dict:
                out_vcf_dict[type_position_key] = hop_feature[type_position_key]
            else:
                # multiple transcripts go through the same one hop
                out_vcf_dict[type_position_key][
                    "TRANSCRIPT_ID"
                ] += f",{hop_feature[type_position_key]['TRANSCRIPT_ID']}"
                out_vcf_dict[type_position_key][
                    "SEGMENT1"
                ] += f",{hop_feature[type_position_key]['SEGMENT1']}"
                out_vcf_dict[type_position_key][
                    "SEGMENT2"
                ] += f",{hop_feature[type_position_key]['SEGMENT2']}"
                out_vcf_dict[type_position_key][
                    "SR_ID"
                ] += f",{hop_feature[type_position_key]['SR_ID']}"
                # deal with 'Y' shape TSG
                if not out_vcf_dict[type_position_key]["READS"].issuperset(
                    hop_feature[type_position_key]["READS"]
                ):
                    out_vcf_dict[type_position_key]["READS"].update(
                        hop_feature[type_position_key]["READS"]
                    )
                    out_vcf_dict[type_position_key]["SR"] += hop_feature[
                        type_position_key
                    ]["SR"]
                    out_vcf_dict[type_position_key]["OSR"] += hop_feature[
                        type_position_key
                    ]["OSR"]

        for _idx, _out_vcf_hop in enumerate(out_vcf_dict, 1):
            hop_vcf_feature = vcf_feature_transformer(out_vcf_dict[_out_vcf_hop], _idx)
            self.write_line(self.formatter(hop_vcf_feature))

    @property
    def header(self) -> str:
        """VCF header provides metadata describing the body of the file."""

        date = datetime.datetime.today().strftime("%Y%m%d")
        source = "ScanMST"
        reference = f'<CMD={obtain_reference_from_bam_header(self.bam_header)},Description="Alignment parameters">'

        header_lines = [
            "##fileformat=VCFv4.3",
            f"##fileDate={date}",
            f"##source={source}",
            f"##reference={reference}",
        ]
        header_lines += self.get_contigs()

        for _id in VCFWriter.reserved_info:
            number: str | int = 0 if VCFWriter.reserved_info[_id] == "Flag" else 1
            if _id in {"TRANSCRIPT_ID", "SR_ID", "SEGMENT1", "SEGMENT2"}:
                number = "."
            header_lines.append(
                f'##INFO=<ID={_id},Number={number},Type={VCFWriter.reserved_info[_id]},Description="{VCFWriter.description[_id]}">',
            )

        for _id in VCFWriter.reserved_format:
            header_lines.append(
                f'##FORMAT=<ID={_id},Number=1,Type={VCFWriter.reserved_format[_id]},Description="{VCFWriter.description[_id]}">',
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
    rescue_sr: bool,
    cluster_id: str,
):
    """Obtain hop vcf features from one series."""
    path_hops_features = []

    can_field_dict = {0: "NONCANONICAL", 1: "CANONICAL"}
    for event_id, current_node in enumerate(nlpath.nodes[:-1], 1):
        current_edge = nlpath.next_edge(current_node, event_id - 1)
        next_node = nlpath[event_id]

        can_field = can_field_dict[current_edge.splicing_code]
        gene1, gene2 = current_edge.gene1, current_edge.gene2

        mode1_, mode2_ = current_edge.modes
        mode1 = mode1_.to_str()
        mode2 = mode2_.to_str()

        if current_edge is None:
            raise BreakpointNotFoundError(current_node.query_name)

        chrom1, pos1 = current_edge.break_point1.to_tuple()
        chrom2, pos2 = current_edge.break_point2.to_tuple()

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

        sv_distance = (
            abs(pos1 - pos2) if not current_edge.variation_type.is_tra() else 0
        )
        dp1 = (
            0
            if current_edge.break_point1.depth is None
            else current_edge.break_point1.depth
        )
        dp2 = (
            0
            if current_edge.break_point2.depth is None
            else current_edge.break_point2.depth
        )
        pso = (
            0
            if dp1 == 0 or dp2 == 0
            else current_edge.sr / (current_edge.sr + (dp1 + dp2) / 2)
        )

        if rescue_sr:
            sr = current_edge.sr
            osr = current_edge.original_sr
        else:
            sr = current_edge.sr
            osr = current_edge.sr

        mst_link_type = determine_mst_link_type(
            current_edge
        )
        path_hops_features.append(
            {
                f"{current_edge.variation_type}_{chrom1}|{pos1 + 1}_{chrom2}|{pos2 + 1}": {
                    "CHROM": chrom1,
                    "POS": f"{pos1 + 1}",
                    "REF": ".",
                    "ALT": f"<{mst_link_type}>",
                    "SVTYPE": current_edge.variation_type,
                    "LINKTYPE": mst_link_type,
                    "SR": sr,
                    "OSR": osr,
                    "CAN": can_field,
                    "CHR2": chrom2,
                    "SVEND": f"{pos2 + 1}",
                    "DP1": f"{dp1}",
                    "DP2": f"{dp2}",
                    "PSI": f"{pso:.3g}",
                    "SVLEN": f"{sv_distance}",
                    "GENE1": f"{gene1}",
                    "GENE2": f"{gene2}",
                    "SEGMENT1": current_node.id,
                    "SEGMENT2": next_node.id,
                    "STRAND1": f"{current_node.strand}",
                    "STRAND2": f"{next_node.strand}",
                    "MODE1": f"{mode1}",
                    "MODE2": f"{mode2}",
                    "TRANSCRIPT_ID": nlpath.id,
                    "GENE_ID": f"{cluster_id}",
                    "SR_ID": f"{'|'.join(current_edge.read_ids)}",
                    "READS": set(current_edge.read_ids),
                    "SVMETHOD": "ScanMST",
                    "HOMSEQ": microhomology_sequence if microhomology_sequence else ".",
                    "INSSEQ": (
                        microinsertion_sequence if microinsertion_sequence else "."
                    ),
                    "ID": current_edge.id,
                },
            },
        )

    return path_hops_features


def vcf_feature_transformer(feature_dict: dict[str, str], idx: int) -> list[str]:
    """VCF feature transformer."""
    info_field = (
        f"{feature_dict['CAN']};"
        f"LINKTYPE={feature_dict['LINKTYPE']};SR={feature_dict['SR']};OSR={feature_dict['OSR']};"
        f"CHR2={feature_dict['CHR2']};SVEND={feature_dict['SVEND']};DP1={feature_dict['DP1']};"
        f"DP2={feature_dict['DP2']};PSI={feature_dict['PSI']};SVLEN={feature_dict['SVLEN']};"
        f"GENE1={feature_dict['GENE1']};GENE2={feature_dict['GENE2']};SVTYPE={feature_dict['SVTYPE']};"
        f"SEGMENT1={feature_dict['SEGMENT1']};SEGMENT2={feature_dict['SEGMENT2']};"
        f"STRAND1={feature_dict['STRAND1']};STRAND2={feature_dict['STRAND2']};"
        f"MODE1={feature_dict['MODE1']};MODE2={feature_dict['MODE2']};"
        f"HOMSEQ={feature_dict['HOMSEQ']};INSSEQ={feature_dict['INSSEQ']};"
        f"TRANSCRIPT_ID={feature_dict['TRANSCRIPT_ID']};GENE_ID={feature_dict['GENE_ID']};"
        f"SR_ID={feature_dict['SR_ID']};SVMETHOD={feature_dict['SVMETHOD']}"
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
        if node.strand.is_forward()
        else reverse_complement(novel_insertion_sequence)
    )
