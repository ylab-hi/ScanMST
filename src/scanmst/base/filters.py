"""Filters based on breakpoints or circurlarRNAs."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING

import HTSeq

from .basic import Interval, Strand

if TYPE_CHECKING:
    from .basic_class import Event


@dataclass
class ExonInfo:
    """Class to store the information of one annotated exon."""

    chrom: str
    interval: Interval
    strand: Strand | None = None
    trx_id: str | None = None

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
        return f"Exon({self.chrom}:{self.start}-{self.end}:{self.strand}, {self.trx_id})"

    def __hash__(self) -> int:
        """Hash an exon."""
        return hash(self.trx_id) ^ hash(self.chrom) ^ hash(self.start) ^ hash(self.end) ^ hash(self.strand)

    def obtain_trx_id(self):
        """Get the transcript id of an exon."""
        return self.trx_id


class ExonFilter:
    """ExonFilter is used to filter out events with both breakpoints harbored in the same exon."""

    def __init__(self, gtf_file: str, boundary_size: int) -> None:
        """Initialize the ExonFilter class."""
        self.exons_gas, _ = _extract_annotated_exons(
            gtf_file,
            boundary_size,
            shrink=True,
            consider_strand=False,
        )

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
) -> tuple[HTSeq.GenomicArrayOfSets, HTSeq.GenomicArrayOfSets]:
    """Extract annotated exons from input GTF file.

    :param in_file: gene annotation file (GTF file)
    :param boundary_size: boundary size for narrow down exon
    :param shrink: shrink or expand the annotated exons
    :param consider_strand: consider strand information
    :return: annotated exons
    """
    gtf_file = HTSeq.GFF_Reader(in_file)
    exons_gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    introns_gas = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    trx_to_exon = defaultdict(list)
    trx_to_intron = defaultdict(list)

    for feature in gtf_file:
        if feature.type == "exon" and "_" not in feature.iv.chrom:
            trx_id = feature.attr["transcript_id"]
            trx_to_exon[trx_id].append(feature.iv)

    for trx_id in trx_to_exon:
        exon_list = trx_to_exon[trx_id]
        exon_list.sort(key=lambda x: x.start)  # type: ignore
        first_exon = exon_list[0]
        strand = first_exon.strand

        tmp_list = []
        for _exon in exon_list:
            chrom = _exon.chrom
            start = _exon.start
            end = _exon.end
            tmp_list.append(start)
            tmp_list.append(end)
            tmp_list.pop(0)
            tmp_list.pop(-1)
            if len(tmp_list) >= 2:
                for intron_start, intron_end in zip(tmp_list[0::2], tmp_list[1::2]):
                    trx_to_intron[trx_id].append(
                        HTSeq.GenomicInterval(chrom, intron_start, intron_end, strand),
                    )
                    intron_id = ExonInfo(
                        chrom,
                        Interval(intron_start, intron_end),
                        strand,
                        trx_id,
                    )
                    if consider_strand:
                        iv_ = HTSeq.GenomicInterval(
                            chrom,
                            intron_start - boundary_size,
                            intron_end + boundary_size,
                            strand,
                        )
                    else:
                        iv_ = HTSeq.GenomicInterval(
                            chrom,
                            intron_start - boundary_size,
                            intron_end + boundary_size,
                            ".",
                        )
                    introns_gas[iv_] += intron_id

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
    return exons_gas, introns_gas


class CircRNAFilter:
    """CircRNAFilter is used to filter out suspicious circular RNAs.
    Features of suspicious circular RNAs:
    1) All the hops types are TDUP.
    2) There are full inclusions relationship between mega-exons, in terms of exons.
       For multi-hop transcripts, the middle mega-exons should be identical.
       e.g., [3][4]->[1][2][3][4]
             [4]->[1][2][3][4]->[1][2][3][4]->[1].
    3) there are inclusive relationship between mega-exons
       e.g., [1][2][3] -> [1]

    4) there are no inclusive relationship between mega-exons in the identical annotated transcript
             [XXXX]-[XXXX]->
             [2]       [1]
             [XXXXX]-[XXXX]->
             [3] [1]-[ 2  ]

           <-[XXXX]-[XXXX]
             [1]       [2]
           <-[XXXXX]-[XXXX]
             [1] [3]-[ 2  ]

    5) two mega-exons within annotated intron
          [XXXX]----------------------[XXXXX]->
                 [ 2 ]             [1]

         <-[XXXX]----------------------[XXXXX]
                 [ 1 ]             [2]
    6) no overlapping with annotations
         ------------------------------------------->
                                 [   1    ]
                [          2              ]
                [  3 ]
        ----------------------------------->
                                 [   1    ]
                [          2              ]
    7) false positive cases in ONT directRNA, since dorado cannot prefectly remove internal adapters
        multihop events with other junction types besides TDUP
        -------chr2>------ ---chr16------>
       [ 1  ] --- [ 2 ]
        [ 3 ] --- [ 4  ]
                              [ 5 ]

    """

    def __init__(
        self,
        gtf_file: str,
        boundary_size: int,
        breakpoint_diff_threshold: int = 10,
    ) -> None:
        """Initialize the CircRNAFilter class."""
        self.exons_gas, self.introns_gas = _extract_annotated_exons(
            gtf_file,
            boundary_size,
            shrink=False,
            consider_strand=True,
        )
        self.breakpoint_diff_threshold = breakpoint_diff_threshold

    def is_circrna(self, nlpath) -> bool:
        nodes = nlpath.nodes
        # one-hop event
        if len(nodes) == 2:
            longest_node = CircRNAFilter.obtain_longest_mega_exon(nodes)
            current_node, next_node = nodes
            current_edge = nlpath.next_edge(current_node, 0)
            # high-confidence circular RNA
            circular_condition1 = bool(
                current_edge.variation_type.is_tdup()
                and (
                    (
                        current_node.introns
                        and next_node.introns
                        and len(current_node.introns) > 0
                        and len(next_node.introns) > 0
                        and (
                            set(current_node.introns).issuperset(set(next_node.introns)) or set(current_node.introns).issubset(set(next_node.introns))
                        )
                    )
                    or (set(current_node.exons).issuperset(set(next_node.exons)) or set(current_node.exons).issubset(set(next_node.exons)))
                    or (current_node.ref_start == next_node.ref_start or current_node.ref_end == next_node.ref_end)
                )
                and self.is_megaexon_superpose_with_annotated_exons(longest_node),
            )

            # medium-confidence circular RNA
            circular_condition2 = bool(
                current_edge.variation_type.is_tdup()
                and self.is_two_megaexon_form_a_partial_loop_within_annotated_transcript(
                    current_node,
                    next_node,
                ),
            )
            # medium-confidence circular RNA
            circular_condition3 = bool(
                current_edge.variation_type.is_tdup()
                and (
                    abs(current_node.ref_start - next_node.ref_start) <= self.breakpoint_diff_threshold
                    or abs(current_node.ref_end - next_node.ref_end) <= self.breakpoint_diff_threshold
                )
                and (current_node.introns == next_node.introns)
            )

            # low-confidence circular RNA
            circular_condition4 = bool(
                current_edge.variation_type.is_tdup()
                and self.is_two_megaexon_within_annotated_intron(
                    current_node,
                    next_node,
                ),
            )

            return circular_condition1 or circular_condition2 or circular_condition3 or circular_condition4

        # multi-hop event
        num_of_tdups = 0
        num_of_hops = len(nodes) - 1
        num_of_hops_satisfy_condition = 0
        ont_condition = False
        for _id, current_node in enumerate(nodes[:-1], 1):
            current_edge = nlpath.next_edge(current_node, _id - 1)
            next_node = nlpath[_id]
            if current_edge.variation_type.is_tdup():
                num_of_tdups += 1
                if (
                    abs(current_node.ref_start - next_node.ref_start) <= self.breakpoint_diff_threshold
                    or abs(current_node.ref_end - next_node.ref_end) <= self.breakpoint_diff_threshold
                ) and (current_node.introns == next_node.introns):
                    ont_condition = True

            # first hop
            if _id == 1:
                if (
                    set(current_node.exons).issubset(set(next_node.exons))
                    or (
                        current_node.introns
                        and next_node.introns
                        and len(current_node.introns) > 0
                        and len(next_node.introns) > 0
                        and set(current_node.introns).issubset(set(next_node.introns))
                    )
                    or (
                        (not current_node.introns)
                        and (
                            abs(current_node.ref_start - next_node.ref_start) <= self.breakpoint_diff_threshold
                            or abs(current_node.ref_end - next_node.ref_end) <= self.breakpoint_diff_threshold
                        )
                    )
                ):
                    num_of_hops_satisfy_condition += 1

            # last hop
            elif _id == num_of_hops:
                if (
                    set(current_node.exons).issuperset(set(next_node.exons))
                    or (
                        current_node.introns
                        and next_node.introns
                        and len(current_node.introns) > 0
                        and len(next_node.introns) > 0
                        and set(current_node.introns).issuperset(set(next_node.introns))
                    )
                    or (
                        (not next_node.introns)
                        and (
                            abs(current_node.ref_start - next_node.ref_start) <= self.breakpoint_diff_threshold
                            or abs(current_node.ref_end - next_node.ref_end) <= self.breakpoint_diff_threshold
                        )
                    )
                ):
                    num_of_hops_satisfy_condition += 1

            # middle hops
            elif set(current_node.exons) == set(next_node.exons) or (
                abs(current_node.ref_start - next_node.ref_start) <= self.breakpoint_diff_threshold
                and abs(current_node.ref_end - next_node.ref_end) <= self.breakpoint_diff_threshold
            ):
                num_of_hops_satisfy_condition += 1

        return num_of_hops_satisfy_condition == num_of_tdups == num_of_hops or ont_condition

    def is_two_megaexon_within_annotated_intron(
        self,
        first_node,
        second_node,
    ) -> bool:
        """Check if two DUP megaexons form a loop within an annotated intron."""
        strand_first = first_node.strand
        strand_second = second_node.strand
        exons_of_first_node = set(first_node.exons)
        exons_of_second_node = set(second_node.exons)
        # rule out duplicated exons and interspersed exons
        if len(exons_of_first_node.intersection(exons_of_second_node)) > 0 or strand_first != strand_second or first_node.chrom != second_node.chrom:
            return False

        chrom = first_node.chrom

        anchor1 = None
        anchor2 = None

        if strand_first == strand_second and str(strand_first) == "+" and first_node.ref_start > second_node.ref_end:
            anchor1 = HTSeq.GenomicPosition(chrom, first_node.ref_end, "+")
            anchor2 = HTSeq.GenomicPosition(chrom, second_node.ref_start, "+")

        elif strand_first == strand_second and str(strand_first) == "-" and first_node.ref_end < second_node.ref_start:
            anchor1 = HTSeq.GenomicPosition(chrom, first_node.ref_start, "-")
            anchor2 = HTSeq.GenomicPosition(chrom, second_node.ref_end, "-")

        common_introns = set()

        if anchor1 and anchor2:
            intron_set1 = self.introns_gas[anchor1]
            intron_set2 = self.introns_gas[anchor2]
            common_introns = intron_set1.intersection(intron_set2)

        # no overlapping annotated transcript
        common_intron_condition = len(common_introns) > 0

        mono_exon_condition = None
        # mono-exon even not within an intronic region, could be indicative of circular RNA
        if len(exons_of_first_node) == len(exons_of_second_node) == 1:
            if str(strand_first) == "+":
                mono_exon_condition = first_node.ref_start > second_node.ref_end
            elif str(strand_first) == "-":
                mono_exon_condition = first_node.ref_end < second_node.ref_start

        return bool(common_intron_condition or mono_exon_condition)

    def is_two_megaexon_form_a_partial_loop_within_annotated_transcript(
        self,
        first_node,
        second_node,
    ) -> bool:
        """Check if two DUP megaexons form a loop within an annotated transcript."""
        strand_first = first_node.strand
        strand_second = second_node.strand
        exons_of_first_node = set(first_node.exons)
        exons_of_second_node = set(second_node.exons)
        # rule out duplicated exons and interspersed exons
        if len(exons_of_first_node.intersection(exons_of_second_node)) > 0 or strand_first != strand_second or first_node.chrom != second_node.chrom:
            return False

        chrom = first_node.chrom

        # outer breakpoints
        anchor1 = None
        anchor2 = None
        # inner breakpoints
        anchor3 = None
        anchor4 = None

        if strand_first == strand_second and str(strand_first) == "+":
            anchor1 = HTSeq.GenomicPosition(chrom, first_node.ref_end, "+")
            anchor2 = HTSeq.GenomicPosition(chrom, second_node.ref_start, "+")

            anchor3 = HTSeq.GenomicPosition(chrom, first_node.ref_start, "+")
            anchor4 = HTSeq.GenomicPosition(chrom, second_node.ref_end, "+")

        elif strand_first == strand_second and str(strand_first) == "-":
            anchor1 = HTSeq.GenomicPosition(chrom, first_node.ref_start, "-")
            anchor2 = HTSeq.GenomicPosition(chrom, second_node.ref_end, "-")

            anchor3 = HTSeq.GenomicPosition(chrom, first_node.ref_end, "-")
            anchor4 = HTSeq.GenomicPosition(chrom, second_node.ref_start, "-")

        if anchor1 and anchor2 and anchor3 and anchor4:
            exon_set1 = self.exons_gas[anchor1]
            exon_set2 = self.exons_gas[anchor2]

            transcript_set1 = {i.obtain_trx_id() for i in exon_set1}
            transcript_set2 = {i.obtain_trx_id() for i in exon_set2}
            common_transcripts_1 = transcript_set1.intersection(transcript_set2)

            exon_set3 = self.exons_gas[anchor3]
            exon_set4 = self.exons_gas[anchor4]

            transcript_set3 = {i.obtain_trx_id() for i in exon_set3}
            transcript_set4 = {i.obtain_trx_id() for i in exon_set4}
            common_transcripts_2 = transcript_set3.intersection(transcript_set4)

            common_transcripts = common_transcripts_1.intersection(common_transcripts_2)
            # no overlapping annotated transcript
            return len(common_transcripts) != 0

        return False

    def is_megaexon_superpose_with_annotated_exons(
        self,
        node,
        threshold: int = 10,
    ) -> bool:
        """Check if all the exons in the megaexon can superpose with annotated exons."""
        flag = True
        chrom = node.chrom
        strand = str(node.strand)
        for _exon in node.exons:
            exon_start_ = _exon.start
            exon_end_ = _exon.end
            exon_start = HTSeq.GenomicPosition(chrom, exon_start_, strand)
            exon_end = HTSeq.GenomicPosition(chrom, exon_end_, strand)
            exon_set1 = self.exons_gas[exon_start]
            exon_set2 = self.exons_gas[exon_end]
            common_exons = exon_set1.intersection(exon_set2)
            # no overlapping annotated exon
            if len(common_exons) == 0:
                return False

            # overlapping annotated exon does not satisfy condition
            if not CircRNAFilter.is_largest_overlapping_exon(
                common_exons,
                exon_start_,
                exon_end_,
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
            start_position - largest_overlapping_exon.start < threshold and largest_overlapping_exon.end - end_position < threshold,
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

    def __init__(self, rt_switching_filter_len: int) -> None:
        """Initialize the RTSwitchingFilter class."""
        self.filter_size = rt_switching_filter_len

    def is_from_rt_switching(self, event: Event) -> bool:
        """The event is from RT switching."""
        return event.has_microhomology() and event.insertion_microhomology_len > self.filter_size
