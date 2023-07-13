"""Filters based on breakpoints or circurlarRNAs."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import HTSeq

from .basic import Interval, Strand

if TYPE_CHECKING:
    from scannls.type import LoggerType

    from .basic_class import Event


@dataclass
class ExonInfo:
    """Class to store the information of one annotated exon."""

    chrom: str
    interval: Interval
    strand: Optional[Strand] = None
    trx_id: Optional[str] = None

    # fmt: off
    @property
    def start(self): return self.interval.start
    @start.setter
    def start(self, value: int): self.interval.start = value
    @property
    def end(self): return self.interval.end
    @end.setter
    def end(self, value: int): self.interval.end = value
    # fmt: on

    def __repr__(self) -> str:
        """Get a string representation of an Exon."""
        return (
            f"Exon({self.chrom}:{self.start}-{self.end}:{self.strand}, "
            f"{self.trx_id})"
        )


class ExonFilter:
    """ExonFilter is used to filter out events with both breakpoints harbored in the same exon."""

    def __init__(self, gtf_file: str, boundary_size: int, logger: LoggerType) -> None:
        """Initialize the ExonFilter class."""
        self.exons_gas = _extract_annotated_exons(
            gtf_file,
            boundary_size,
            shrink=True,
            consider_strand=False,
        )
        self.logger = logger

    def is_breakpoints_in_same_exon(self, event: Event) -> bool:
        """Annotate vcf file."""
        # only apply for TDUP or INV
        if event.sv_type in {"TRA", "DEL"}:
            return False

        gp1 = HTSeq.GenomicPosition(event.chrom1, event.junction_pos1, ".")
        gp2 = HTSeq.GenomicPosition(event.chrom2, event.junction_pos2, ".")

        exon_set1 = self.exons_gas[gp1]
        exon_set2 = self.exons_gas[gp2]
        common_exons = exon_set1.intersection(exon_set2)
        return len(common_exons) > 0


def _extract_annotated_exons(
    in_file: str,
    boundary_size: int = 10,
    minimum_exon_size: int = 30,
    *,
    shrink=False,
    consider_strand=False,
) -> HTSeq.GenomicArrayOfSets:
    """Extract annotated exons from input GTF file.

    :param in_file: gene annotation file (GTF file)
    :param boundary_size: boundary size for narrow down exon
    :param shrink: shrink or expand the annotated exons
    :param consider_strand: consider strand information
    :type in_file: str
    :type boundary_size: int
    :type shrink: bool
    :type consider_strand: bool
    :return: annotated exons
    :rtype: HTSeq.GenomicArrayOfSets
    """
    gtf_file = HTSeq.GFF_Reader(in_file)
    exons_gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    trx_to_exon = defaultdict(list)

    for feature in gtf_file:
        if feature.type == "exon":
            trx_id = feature.attr["transcript_id"]
            trx_to_exon[trx_id].append(feature.iv)

    for trx_id in trx_to_exon:
        exon_list = trx_to_exon[trx_id]
        exon_list.sort(key=lambda x: x.start)  # type: ignore
        first_exon = exon_list[0]
        strand = first_exon.strand
        for _exon in exon_list:
            chrom = _exon.chrom
            start = _exon.start
            end = _exon.end
            if end - start >= minimum_exon_size:
                exon_id = ExonInfo(chrom, Interval(start, end), strand, trx_id)
                if shrink:
                    if consider_strand:
                        iv = HTSeq.GenomicInterval(
                            chrom,
                            start + boundary_size,
                            end - boundary_size,
                            strand,
                        )
                    else:
                        iv = HTSeq.GenomicInterval(
                            chrom,
                            start + boundary_size,
                            end - boundary_size,
                            ".",
                        )
                elif consider_strand:
                    iv = HTSeq.GenomicInterval(
                        chrom,
                        start - boundary_size,
                        end + boundary_size,
                        strand,
                    )
                else:
                    iv = HTSeq.GenomicInterval(
                        chrom,
                        start - boundary_size,
                        end + boundary_size,
                        ".",
                    )

                exons_gas[iv] += exon_id
    return exons_gas


class CircRNAFilter:
    """CircRNAFilter is used to filter out suspicious circular RNAs.
    Features of suspicious circular RNAs:
    1) All the hops types are TDUP.
    2) There are full inclusions relationship between mega-exons, in terms of exons.
       For multi-hop transcripts, the middle mega-exons should be identical.
       e.g., [3][4]->[1][2][3][4]
             [4]->[1][2][3][4]->[1][2][3][4]->[1].
    """

    def __init__(self, gtf_file: str, boundary_size: int, logger: LoggerType) -> None:
        """Initialize the CircRNAFilter class."""
        self.logger = logger
        self.exons_gas = _extract_annotated_exons(
            gtf_file,
            boundary_size,
            shrink=False,
            consider_strand=True,
        )

    def is_circrna(self, series) -> bool:
        nodes = series.nodes
        # one-hop event
        if len(nodes) == 2:
            longest_node = CircRNAFilter.obtain_longest_mega_exon(nodes)
            return bool(
                nodes[0].sv_type == "TDUP"
                and (
                    len(nodes[0].introns) > 0
                    and len(nodes[1].introns) > 0
                    and (
                        set(nodes[0].introns).issuperset(set(nodes[1].introns))
                        or set(nodes[0].introns).issubset(set(nodes[1].introns))
                    )
                    or (
                        set(nodes[0].exons).issuperset(set(nodes[1].exons))
                        or set(nodes[0].exons).issubset(set(nodes[1].exons))
                    )
                    or (
                        nodes[0].ref_start == nodes[1].ref_start
                        or nodes[0].ref_end == nodes[1].ref_end
                    )
                )
                and self.is_megaexon_superpose_with_annotated_exons(longest_node),
            )

        # multi-hop event
        num_of_tdups = 0
        num_of_hops = len(nodes) - 1
        num_of_hops_satisfy_condition = 0
        for _id, current_node in enumerate(nodes[:-1], 1):
            next_node = series[_id]
            if current_node.sv_type == "TDUP":
                num_of_tdups += 1

            # first hop
            if _id == 1:
                if (
                    set(current_node.exons).issubset(set(next_node.exons))
                    or (
                        len(current_node.introns) > 0
                        and len(next_node.introns) > 0
                        and set(current_node.introns).issubset(set(next_node.introns))
                    )
                    or current_node.ref_start == next_node.ref_start
                    or current_node.ref_end == next_node.ref_end
                ) and self.is_megaexon_superpose_with_annotated_exons(next_node):
                    num_of_hops_satisfy_condition += 1
            # last hop
            elif _id == num_of_hops:
                if (
                    set(current_node.exons).issuperset(set(next_node.exons))
                    or (
                        len(current_node.introns) > 0
                        and len(next_node.introns) > 0
                        and set(current_node.introns).issuperset(set(next_node.introns))
                    )
                    or current_node.ref_start == next_node.ref_start
                    or current_node.ref_end == next_node.ref_end
                ) and self.is_megaexon_superpose_with_annotated_exons(current_node):
                    num_of_hops_satisfy_condition += 1
            # middle hops
            elif (
                set(current_node.exons) == set(next_node.exons)
                or (
                    current_node.ref_start == next_node.ref_start
                    and current_node.ref_end == next_node.ref_end
                )
            ) and self.is_megaexon_superpose_with_annotated_exons(current_node):
                num_of_hops_satisfy_condition += 1

        return num_of_hops_satisfy_condition == num_of_tdups == num_of_hops

    def is_megaexon_superpose_with_annotated_exons(
        self,
        node,
        threshold: int = 10,
    ) -> bool:
        """Check if all the exons in the megaexon can superpose with annotated exons."""
        flag = True
        chrom = node.chrom
        strand = node.strand
        for _exon in node.exons:
            _exon_start = _exon[0]
            _exon_end = _exon[1]
            exon_start = HTSeq.GenomicPosition(chrom, _exon_start, strand)
            exon_end = HTSeq.GenomicPosition(chrom, _exon_end, strand)
            exon_set1 = self.exons_gas[exon_start]
            exon_set2 = self.exons_gas[exon_end]
            common_exons = exon_set1.intersection(exon_set2)
            # no overlapping annotated exon
            if len(common_exons) == 0:
                return False

            # overlapping annotated exon does not satisfy condition
            if not CircRNAFilter.is_largest_overlapping_exon(
                common_exons,
                _exon_start,
                _exon_end,
                threshold,
            ):
                flag = False

        return flag

    @staticmethod
    def is_largest_overlapping_exon(
        overlapping_exons_set: set[ExonInfo],
        start_position: int,
        end_position: int,
        threshold: int = 10,
    ) -> bool:
        """Check if the exon in mega exon overlapped with annotated exon with largest fraction.

        [XXXXXXXXXXXXXXXXXXXXXXX]
          [XXXXXXXXXXXXXXXXXX]
            [XXXXXXXXXXXXX] <- largest_overlapping_exon
              [000000000]

        """
        largest_overlapping_exon = sorted(
            overlapping_exons_set,
            key=lambda x: (end_position - start_position) / (x.end - x.start),
            reverse=True,
        )[0]

        return bool(
            start_position - largest_overlapping_exon.start < threshold
            and largest_overlapping_exon.end - end_position < threshold,
        )

    @staticmethod
    def obtain_longest_mega_exon(nodes):
        """Obtain mega-exon with the longest exon length."""
        return sorted(
            nodes,
            key=lambda x: sum(len(j) for j in x.exons),
            reverse=True,
        )[0]


class RTSwitchingFilter:
    """RTSwitchingFilter is used to filter out events derived from potential RT switching."""

    def __init__(self, rt_switching_filter_len: int, logger: LoggerType) -> None:
        """Initialize the RTSwitchingFilter class."""
        self.filter_size = rt_switching_filter_len
        self.logger = logger

    def is_from_rt_switching(self, event: Event) -> bool:
        """The event is from RT switching."""
        return (
            event.has_microhomology()
            and event.insertion_microhomology_len > self.filter_size
        )
