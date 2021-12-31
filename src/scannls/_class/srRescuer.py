# !/usr/bin/env python
"""SR rescuer.

@Filename:    srRescuer.py
@license:     MIT Licence
@Time:        12/30/21 15:00 PM
"""
from typing import Any  # type: ignore

import pysam
from loguru._logger import Logger  # type: ignore

from ..utils import get_softclip_length
from .basicClass import Series  # type: ignore

# from .basicClass import Node  # type: ignore

"""
import copy
import types
from typing import Any  # type: ignore
from typing import Dict  # type: ignore
from typing import Iterable  # ignore
from typing import List  # ignore
from typing import Set  # ignore
from typing import Tuple  # ignore
from typing import Union  # ignore

from ..utils import timeit  # ignore
from .basicClass import Node  # ignore
from .basicClass import Read  # ignore

NodeType = Union[Node, Insertion]
"""


class Rescuer:
    """Rescue SR from softclipped non-chimeric reads."""

    def __init__(
        self, input_bam: Any, mapq_cutoff: int, soft_len_cutoff: int, logger: Logger
    ) -> None:
        """Initialize Rescuer.

        :param logger: logger
        """
        self.in_bam = pysam.AlignmentFile(input_bam, "rb")
        self.mapq_cutoff = mapq_cutoff
        self.logger = logger
        self.soft_len_cutoff = soft_len_cutoff

    def __repr__(self):
        """Represent Rescuer."""
        return f"{self.__class__.__name__}()"

    def update_sr(self, series: Series) -> Any:
        """Update SR for input series."""
        for node in series:
            chrom = node.chrom
            prev_breakpoint = node.prev_breakpoint
            next_breakpoint = node.next_breakpoint

            _region = f"{chrom}:{prev_breakpoint}-{next_breakpoint}"
            for col in self.in_bam.pileup(
                region=_region, truncate=True, stepper="nofilter", min_base_quality=0
            ):
                # dp = col.nsegments
                sr_list = {1: [], 2: []}
                for read in col.pileups:
                    # read is an instance of pysam.PileupRead
                    if read.alignment.mapq >= self.mapq_cutoff and read.query_position:
                        # the read has soft-clipped part but no SV tag
                        if not read.alignment.has_tag("SV"):
                            if "S" in read.alignment.cigarstring:
                                (
                                    soft_len,
                                    soft_seq,
                                    soft_pos,
                                    left_or_right,
                                ) = get_softclip_length(read.alignment)
                                # the pileup position is equal to the soft-clipped connection point
                                # xxxxxxxxSyyyyyyyyMzzzzzS
                                #         ^      ^
                                # if soft_pos == col.reference_pos and 'N' not in soft_seq:
                                if "N" not in soft_seq:
                                    if soft_len >= self.soft_len_cutoff:
                                        if left_or_right == 1:
                                            sr_list[left_or_right].append(
                                                read.alignment.query_sequence[
                                                    read.query_position + 1 :
                                                ]
                                            )
                                        elif left_or_right == 2:
                                            sr_list[left_or_right].append(
                                                read.alignment.query_sequence[
                                                    : read.query_position
                                                ]
                                            )
