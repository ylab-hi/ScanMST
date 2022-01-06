# !/usr/bin/env python
"""SR rescuer.

@Filename:    srRescuer.py
@license:     MIT Licence
@Time:        12/30/21 15:00 PM
"""
from typing import Any
from typing import Dict
from typing import List

import skbio  # type: ignore
from loguru._logger import Logger  # type: ignore
from pysam import AlignmentFile  # type: ignore

from ..utils import get_softclip_length
from .basicClass import Node  # type: ignore
from .basicClass import Series  # type: ignore


class SRRescuer:
    """Rescue SR from softclipped non-chimeric reads."""

    def __init__(
        self,
        input_bam: AlignmentFile,
        mapq_cutoff: int,
        soft_len_cutoff: int,
        mismatch_cutoff: int,
        alignment_frac: float,
        logger: Logger,
    ) -> None:
        """Initialize Rescuer.

        :param logger: logger
        """
        self.in_bam = input_bam
        self.mapq_cutoff = mapq_cutoff
        self.soft_len_cutoff = soft_len_cutoff
        self.mismatch_cutoff = mismatch_cutoff
        self.alignment_frac = alignment_frac
        self.logger = logger

    def __repr__(self):
        """Represent Rescuer."""
        return (
            f"{self.__class__.__name__}({self.in_bam.filename}, "
            f"{self.soft_len_cutoff}, {self.mismatch_cutoff}, {self.alignment_frac})"
        )

    @staticmethod
    def mismatch_count(seq: str, seqs: list, alignment_frac: float, mode: int) -> float:
        """Local alignment."""
        mismatch = 1e6
        for each_seq in seqs:
            seq = skbio.DNA(seq)
            each_seq = skbio.DNA(each_seq)
            # if seq or each_seq is an empty string, ignore it
            if not seq or not each_seq:
                continue
            try:
                (
                    alignment,
                    score,
                    start_end_pos,
                ) = skbio.alignment.local_pairwise_align_ssw(seq, each_seq)
            # raise IndexError if SSW cannot work
            except (IndexError, ValueError):
                continue
            if (
                len(alignment[0]) / float(len(seq)) < alignment_frac
                and len(alignment[1]) / float(len(each_seq)) < alignment_frac
            ):
                continue
            if (
                (mode == 1 and start_end_pos[0][0] == 0 and start_end_pos[1][0] == 0)
                or (
                    mode == 2
                    and start_end_pos[0][1] == len(seq) - 1
                    and start_end_pos[1][1] == len(each_seq) - 1
                )
            ) and sum(alignment[0].mismatches(alignment[1])) < mismatch:
                mismatch = sum(alignment[0].mismatches(alignment[1]))
        return mismatch

    def _calculate_sr_for_reads(self, col, query_names, sr_list, sv_list, mode):
        """Calculate SR for reads.

        :param col:
        :param query_names:
        :param sr_list:
        :param sv_list:
        :param mode:
        :return:
        """
        for read in col.pileups:
            # read.alignment is an instance of pysam.AlignedSegment
            aln = read.alignment
            read_name = aln.query_name
            strand = "-" if aln.is_reverse else "+"
            if aln.mapq >= self.mapq_cutoff and read.query_position:
                # the read has soft-clipped part but not an anchor read
                if read_name not in query_names:
                    if "S" in aln.cigarstring:
                        (
                            soft_len,
                            soft_seq,
                            soft_pos,
                            soft_mode,
                        ) = get_softclip_length(aln, mode)
                        # the pileup position is equal to the soft-clipped connection point
                        # xxxxxxxxSyyyyyyyyMzzzzzS
                        #         ^      ^
                        soft_pos = soft_pos - 1 if mode == 1 else soft_pos
                        self.logger.trace(f"{col.reference_pos=}, {soft_pos=}")
                        if (
                            soft_pos == col.reference_pos
                            and soft_len >= self.soft_len_cutoff
                        ):
                            softclipped_seq = (
                                aln.query_sequence[read.query_position + 1]
                                if mode == 1
                                else aln.query_sequence[: read.query_position]
                            )
                            sr_list[strand].append(softclipped_seq)
                # the anchor read
                else:
                    (
                        _,
                        anchor_soft_seq,
                        anchor_soft_pos,
                        anchor_soft_mode,
                    ) = get_softclip_length(aln, mode)
                    anchor_soft_pos = (
                        anchor_soft_pos - 1 if mode == 1 else anchor_soft_pos
                    )
                    self.logger.trace(f"{col.reference_pos=}, {anchor_soft_pos=}")
                    if anchor_soft_pos == col.reference_pos:
                        sv_list[strand].append(anchor_soft_seq)

    def calculate_sr(self, region: str, mode: int, query_name: str) -> int:
        """Calculate SR from softclipped reads without SV tag, provided target region.

        region = 'chrm:start-end'
        ..note.
            rescued reads strand should be the same as the anchor ones
            For 'MS' mode, softclipped position should subtract by one
        """
        sr_list: Dict[str, List[str]] = {"+": [], "-": []}
        sv_list: Dict[str, List[str]] = {"+": [], "-": []}

        query_names = set(query_name.split(","))

        self.logger.trace(f"{region=}")

        for col in self.in_bam.pileup(
            region=region, truncate=True, stepper="nofilter", min_base_quality=0
        ):
            # dp = col.nsegments
            # read is an instance of pysam.PileupRead
            self._calculate_sr_for_reads(col, query_names, sr_list, sv_list, mode)

        rescued_sr = 0
        if sv_list["+"] and sr_list["+"]:
            for _soft_seq in sr_list["+"]:
                if (
                    SRRescuer.mismatch_count(
                        _soft_seq, sv_list["+"], self.alignment_frac, mode
                    )
                    <= self.mismatch_cutoff
                ):
                    rescued_sr += 1

        if sv_list["-"] and sr_list["-"]:
            for _soft_seq in sr_list["-"]:
                if (
                    SRRescuer.mismatch_count(
                        _soft_seq, sv_list["-"], self.alignment_frac, mode
                    )
                    <= self.mismatch_cutoff
                ):
                    rescued_sr += 1
        self.logger.trace(f"{sv_list=}")
        self.logger.trace(f"{sr_list=}")

        return rescued_sr

    @staticmethod
    def obtain_region_for_rescue_sr(node: Node, tgt_name: str) -> str:
        """Obtain target region (S-M boundary, M side) for rescuing SR purpose.

        ..note.
              Due to microhomology, prev_breakpoint/next_breakpoint locates inside the M side of S-M boundary
              Thus, exon start/end (S-M boundary) will be used to rescue SR.
        """
        strand = node.strand
        chrom = node.chrom
        exons = node.exons
        if strand == "+":
            if tgt_name == "next_breakpoint":
                pos = exons[-1][1]
            elif tgt_name == "prev_breakpoint":
                pos = exons[0][0]
        else:
            if tgt_name == "next_breakpoint":
                pos = exons[0][0]
            elif tgt_name == "prev_breakpoint":
                pos = exons[-1][1]
        region = f"{chrom}:{pos}-{pos + 1}"
        return region

    def update_sr(self, series: Series) -> Any:
        """Update SR for input series."""
        for idx in range(len(series) - 1):
            current_node = series[idx]
            next_node = series[idx + 1]
            mode1, mode2 = current_node.modes
            current_node.update_next_breakpoint_depth(self.in_bam, mode1)
            next_node.update_prev_breakpoint_depth(self.in_bam, mode2)
            query_name1 = current_node.query_name
            query_name2 = next_node.query_name

            _region1 = SRRescuer.obtain_region_for_rescue_sr(
                current_node, "next_breakpoint"
            )
            _region2 = SRRescuer.obtain_region_for_rescue_sr(
                next_node, "prev_breakpoint"
            )

            rescued_sr1 = self.calculate_sr(_region1, mode1, query_name1)
            rescued_sr2 = self.calculate_sr(_region2, mode2, query_name2)
            rescued_sr = rescued_sr1 + rescued_sr2
            self.logger.trace(f"{rescued_sr=}")
            if rescued_sr > 0:
                series[idx].update_sr(rescued_sr)
        return series
