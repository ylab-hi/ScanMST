# !/usr/bin/env python
"""SR rescuer.

@Filename:    srRescuer.py
@license:     MIT Licence
@Time:        12/30/21 15:00 PM
"""
import re
from typing import Dict
from typing import Iterable
from typing import List

import parasail  # type: ignore
from pysam import AlignmentFile  # type: ignore

from ..type import LoggerType
from ..utils import get_softclip_length
from .basicClass import NodeType
from .exception import ExonsNotFoundError
from .exception import ModesNotFoundError


# todo: change threshold of mismatch


class SRRescuer:
    """Rescue SR from softclipped non-chimeric reads."""

    def __init__(
        self,
        input_bam: AlignmentFile,
        mapq_cutoff: int,
        soft_len_cutoff: int,
        mismatch_cutoff: int,
        alignment_frac: float,
        logger: LoggerType,
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

    def __repr__(self) -> str:
        """Represent Rescuer."""
        return (
            f"{self.__class__.__name__}({self.in_bam.filename}, "
            f"{self.soft_len_cutoff}, {self.mismatch_cutoff}, {self.alignment_frac})"
        )

    def __call__(self, series: Iterable[NodeType]) -> None:
        """Rescue SR from softclipped non-chimeric reads.

        changed in place

        :param series: Series
        """
        for node in series:
            self.update_sr(node)

    @staticmethod
    def check_if_sr_rescued_depended_on_alignment(
        query_origin_len: int,
        align_result: parasail.bindings_v2.Result,
        alignment_frac: float,
        mismatch_threshold: int,
    ) -> bool:
        """Check if SR is rescued depended on alignment result.

        index style: [ )

        :param mismatch_threshold:
        :param query_origin_len:
        :param align_result:
        :param alignment_frac:
        :return:
        """
        flag = False
        cigar = align_result.cigar.decode.decode()
        pattern = re.compile(r"((?P<length>\d+)(?P<op>\D))")
        if not cigar or re.match(r"^\d+=", cigar) is None:
            return flag
        query_len, mismatch_count = 0, 0
        for match in re.finditer(pattern, cigar):
            length = int(match.group("length"))
            if match.group("op") == "=":
                query_len += length
            else:
                mismatch_count += length
        query_end = align_result.end_query + 1
        query_start = query_end - query_len
        target_end = align_result.end_ref + 1
        target_start = target_end - query_len

        if (
            query_len / query_origin_len >= alignment_frac
            and query_start + target_start == 0
            and mismatch_count <= mismatch_threshold
        ):
            flag = True
        return flag

    @staticmethod
    def determined_num_increment_sr(
        seq: str, seqs: list, alignment_frac: float, mode: int, mismatch_threshold: int
    ) -> int:
        """Local alignment."""
        gaps = 11
        gap_extend = 1
        increment_sr = 0
        if not seq:
            return increment_sr
        seq = seq[::-1] if mode == 2 else seq
        for each_seq in [each_seq for each_seq in seqs if each_seq]:
            each_seq = each_seq[::-1] if mode == 2 else each_seq

            align_result = parasail.sw_trace_striped_sat(
                seq, each_seq, gaps, gap_extend, parasail.dnafull
            )
            if SRRescuer.check_if_sr_rescued_depended_on_alignment(
                len(seq), align_result, alignment_frac, mismatch_threshold
            ):
                increment_sr += 1
        return increment_sr

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
            strand = "-" if aln.is_reverse else "+"
            if aln.mapq >= self.mapq_cutoff and read.query_position:
                (
                    _len,
                    _seq,
                    _pos,
                    _mode,
                ) = get_softclip_length(aln, mode)
                _pos = _pos - 1 if mode == 1 else _pos

                if aln.query_name in query_names:
                    if _pos == col.reference_pos:
                        sv_list[strand].append(_seq)
                elif (
                    "S" in aln.cigarstring
                    and _pos == col.reference_pos
                    and _len >= self.soft_len_cutoff
                ):
                    # the pileup position is equal to the soft-clipped connection point
                    # xxxxxxxxSyyyyyyyyMzzzzzS
                    #         ^      ^
                    sr_list[strand].append(_seq)

    def calculate_sr(self, region: str, mode: int, query_name: str) -> int:
        """Calculate SR from softclipped reads without SV tag, provided target region.

        region = 'chrm:start-end'
        .. note.
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
            # read is an instance of pysam.PileupRead
            self._calculate_sr_for_reads(col, query_names, sr_list, sv_list, mode)

        rescued_sr = 0

        if sv_list["+"] and sr_list["+"]:
            for _soft_seq in sr_list["+"]:
                rescued_sr += SRRescuer.determined_num_increment_sr(
                    _soft_seq,
                    sv_list["+"],
                    self.alignment_frac,
                    mode,
                    self.mismatch_cutoff,
                )

        if sv_list["-"] and sr_list["-"]:
            for _soft_seq in sr_list["-"]:
                rescued_sr += SRRescuer.determined_num_increment_sr(
                    _soft_seq,
                    sv_list["-"],
                    self.alignment_frac,
                    mode,
                    self.mismatch_cutoff,
                )

        return rescued_sr

    @staticmethod
    def obtain_region_for_rescue_sr(node: NodeType, tgt_name: str, mode: int) -> str:
        """Obtain target region (S-M boundary, M side) for rescuing SR purpose.

        ..note.
              Due to microhomology, prev_breakpoint/next_breakpoint locates inside the M side of S-M boundary
              Thus, exon start/end (S-M boundary) will be used to rescue SR.
        """
        strand = node.strand
        chrom = node.chrom
        exons = node.exons
        if exons is None:
            raise SystemExit from ExonsNotFoundError

        if strand == "+":
            pos = exons[-1][1] if tgt_name == "next_breakpoint" else exons[0][0]
        else:
            pos = exons[0][0] if tgt_name == "next_breakpoint" else exons[-1][1]
        region = f"{chrom}:{pos + 1}-{pos + 1}" if mode == 2 else f"{chrom}:{pos}-{pos}"
        return region

    def update_sr(self, current_node: NodeType) -> None:
        """Update SR for input node."""
        if current_node.is_end_node():
            return

        if current_node.modes is None:
            raise SystemExit from ModesNotFoundError

        mode1, mode2 = current_node.modes
        current_node.update_next_breakpoint_depth(self.in_bam, mode1)

        query_name_current = current_node.query_name
        _region_current = SRRescuer.obtain_region_for_rescue_sr(
            current_node, "next_breakpoint", mode1
        )
        rescued_sr = self.calculate_sr(_region_current, mode1, query_name_current)
        for next_node in current_node.successors:
            next_node.update_prev_breakpoint_depth(self.in_bam, mode2)
            query_name_next = next_node.query_name
            _region_next = SRRescuer.obtain_region_for_rescue_sr(
                next_node, "prev_breakpoint", mode2
            )
            rescued_sr_next = self.calculate_sr(_region_next, mode2, query_name_next)
            rescued_sr += rescued_sr_next
        self.logger.trace(f"{rescued_sr=}")
        if rescued_sr > 0:
            current_node.update_sr(rescued_sr)
