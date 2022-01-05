# !/usr/bin/env python
"""SR rescuer.

@Filename:    srRescuer.py
@license:     MIT Licence
@Time:        12/30/21 15:00 PM
"""
from typing import Any
from typing import Dict
from typing import List

import pysam  # type: ignore
import skbio  # type: ignore
from loguru._logger import Logger  # type: ignore [import]

from ..utils import get_softclip_length
from .basicClass import Series  # type: ignore [import]


class Rescuer:
    """Rescue SR from softclipped non-chimeric reads."""

    def __init__(
        self,
        input_bam: Any,
        logger: Logger,
        mapq_cutoff: int,
        soft_len_cutoff: int = 5,
        mismatch_cutoff: int = 3,
        alignment_frac: float = 0.8,
    ) -> None:
        """Initialize Rescuer.

        :param logger: logger
        """
        self.in_bam = pysam.AlignmentFile(input_bam, "rb")
        self.mapq_cutoff = mapq_cutoff
        self.soft_len_cutoff = soft_len_cutoff
        self.mismatch_cutoff = mismatch_cutoff
        self.alignment_frac = alignment_frac
        self.logger = logger

    def __repr__(self):
        """Represent Rescuer."""
        return f"{self.__class__.__name__}()"

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
            except IndexError:
                continue
            except ValueError:
                continue
            if (
                len(alignment[0]) / float(len(seq)) < alignment_frac
                and len(alignment[1]) / float(len(each_seq)) < alignment_frac
            ):
                continue
            if mode == 1 and start_end_pos[0][0] == 0 and start_end_pos[1][0] == 0:
                if sum(alignment[0].mismatches(alignment[1])) < mismatch:
                    mismatch = sum(alignment[0].mismatches(alignment[1]))
            elif (
                mode == 2
                and start_end_pos[0][1] == len(seq) - 1
                and start_end_pos[1][1] == len(each_seq) - 1
            ):
                if sum(alignment[0].mismatches(alignment[1])) < mismatch:
                    mismatch = sum(alignment[0].mismatches(alignment[1]))
            else:
                continue
        return mismatch

    @staticmethod
    def region_in_sv_checker(
        region: str, sv_type: str, mode: int, sv_aln_list: list
    ) -> bool:
        """Check if region in SV tag."""
        tgt_pos = region.split("-")[0]
        tgt_chrm, _tgt_pos = tgt_pos.split(":")
        tgt_pos = int(_tgt_pos) + 1  # type: ignore
        flag = False
        for sv_aln in sv_aln_list:
            _sv_type, _anno_can, _bp1, _bp2, modes, strands, genes = sv_aln.split(",")
            mode1, mode2 = map(int, modes)
            if _sv_type == sv_type:
                if _bp1 == f"{tgt_chrm}:{tgt_pos}" and mode1 == mode:
                    flag = True
                elif _bp2 == f"{tgt_chrm}:{tgt_pos}" and mode2 == mode:
                    flag = True
        return flag

    def calculate_sr(self, region: str, mode: int, query_name: str) -> int:
        """Calculate SR from softclipped reads without SV tag, provided target region.

        region = 'chrm:start-end'
        """
        sr_list: Dict[int, List[str]] = {1: [], 2: []}
        sv_list: Dict[int, List[str]] = {1: [], 2: []}
        query_names = set(query_name.split(","))
        for col in self.in_bam.pileup(
            region=region, truncate=True, stepper="nofilter", min_base_quality=0
        ):
            # dp = col.nsegments
            for read in col.pileups:
                # read is an instance of pysam.PileupRead
                aln = read.alignment
                read_name = aln.query_name
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
                            if soft_pos == col.reference_pos:
                                self.logger.trace(f"{col.reference_pos=}, {soft_pos=}")
                                if soft_len >= self.soft_len_cutoff:
                                    if mode == 1:
                                        softclipped_seq = aln.query_sequence[
                                            read.query_position + 1 :
                                        ]
                                        sr_list[mode].append(softclipped_seq)
                                    elif mode == 2:
                                        softclipped_seq = aln.query_sequence[
                                            : read.query_position
                                        ]
                                        sr_list[mode].append(softclipped_seq)
                    # the anchor read
                    else:
                        (
                            _,
                            anchor_soft_seq,
                            anchor_soft_pos,
                            anchor_soft_mode,
                        ) = get_softclip_length(aln, mode)
                        if anchor_soft_pos == col.reference_pos:
                            if mode == 1:
                                sv_list[mode].append(anchor_soft_seq)
                            elif mode == 2:
                                sv_list[mode].append(anchor_soft_seq)
            rescued_sr = 0
            if sv_list[mode]:
                for _soft_seq in sr_list[mode]:
                    if (
                        Rescuer.mismatch_count(
                            _soft_seq, sv_list[mode], self.alignment_frac, mode
                        )
                        <= self.mismatch_cutoff
                    ):
                        rescued_sr += 1
        return rescued_sr

    def update_sr(self, series: Series) -> Any:
        """Update SR for input series."""
        for idx in range(len(series) - 1):
            current_node = series[idx]
            next_node = series[idx + 1]
            mode1, mode2 = current_node.modes
            query_name1 = current_node.query_name
            query_name2 = next_node.query_name
            _bp1 = current_node.next_breakpoint
            _bp2 = next_node.prev_breakpoint
            _chrom1, _pos1 = _bp1.split(":")
            _chrom2, _pos2 = _bp2.split(":")
            _pos1 = int(_pos1)
            _pos2 = int(_pos2)
            _region1 = f"{_chrom1}:{_pos1+1}-{_pos1+1}"
            _region2 = f"{_chrom2}:{_pos2+1}-{_pos2+1}"

            rescued_sr1 = self.calculate_sr(_region1, mode1, query_name1)
            rescued_sr2 = self.calculate_sr(_region2, mode2, query_name2)
            rescued_sr = rescued_sr1 + rescued_sr2
            if rescued_sr > 0:
                series[idx].update_sr(rescued_sr)
        yield series
