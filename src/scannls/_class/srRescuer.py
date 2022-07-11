# !/usr/bin/env python
"""SR rescuer.

@Filename:    srRescuer.py
@license:     MIT Licence
@Time:        12/30/21 15:00 PM
"""
from typing import Any
from typing import Iterable
from typing import List
from typing import Optional
from typing import Tuple

from .basicClass import Node
from .exception import ExonsNotFoundError
from .exception import ModesNotFoundError
from .type import LoggerType
from scannls import cppext


class SRRescuer:
    """Rescue SR from soft-clipped non-chimeric reads."""

    def __init__(
        self,
        input_bam_file: str,
        mapq_cutoff: int,
        soft_len_cutoff: int,
        mismatch_cutoff: int,
        alignment_frac: float,
        average_read_depth: int,
        logger: LoggerType,
    ) -> None:
        """Initialize Rescuer.

        :param logger: logger
        """
        self.cppext_rescuer = cppext.Rescuer(
            input_bam_file,
            mapq_cutoff,
            soft_len_cutoff,
            mismatch_cutoff,
            alignment_frac,
            10,
            average_read_depth,
        )
        self.logger = logger

    def __call__(self, nodes_in_graph: Iterable[Node]) -> None:
        """Rescue SR from soft-clipped non-chimeric reads.

        changed in place

        :param nodes_in_graph: Series
        """
        query_names_in_graph = set()

        for node in nodes_in_graph:
            query_names_in_graph.update(node.query_name.split(","))

        query_names_in_graph_list = list(query_names_in_graph)

        for node in nodes_in_graph:
            node.original_sr = node.sr
            self.cppext_rescuer.reset_names_list(query_names_in_graph_list)
            self.update_sr(node, query_names_in_graph_list)

        del query_names_in_graph

    @staticmethod
    def obtain_region_for_rescue_sr(
        strand: Optional[str],
        chrom: Optional[str],
        exons: Any,
        tgt_name: str,
        mode: int,
    ) -> Tuple[str, int]:
        """Obtain target region (S-M boundary, M side) for rescuing SR purpose.

        ..note.
              Due to micro homology, prev_breakpoint/next_breakpoint locates inside the M side of S-M boundary
              Thus, exon start/end (S-M boundary) will be used to rescue SR.
        """
        if exons is None or strand is None or chrom is None:
            raise ExonsNotFoundError(f"{chrom=} {strand=} {exons=}")

        if strand == "+":
            pos = exons[-1][1] if tgt_name == "next_breakpoint" else exons[0][0]
        else:
            pos = exons[0][0] if tgt_name == "next_breakpoint" else exons[-1][1]

        if mode == 2:
            pos += 1
        else:
            pos = pos + 1 if pos == 0 else pos

        return chrom, pos

    def update_sr(self, current_node: Node, query_names_in_graph: List[str]) -> None:
        """Update SR for input node."""
        if current_node.is_end_node():
            return

        if current_node.modes is None:
            raise ModesNotFoundError(f"{current_node.query_name}")

        mode1, mode2 = current_node.modes

        chrom_n, pos_n = current_node.get_breakpoint_depth_pos(mode1, "next")
        current_node.next_breakpoint_depth = self.cppext_rescuer.count_reads(
            chrom_n, pos_n, pos_n + 1
        )

        query_name_current = current_node.query_name.split(",")
        chrom, start = SRRescuer.obtain_region_for_rescue_sr(
            current_node.strand,
            current_node.chrom,
            current_node.exons,
            "next_breakpoint",
            mode1,
        )

        self.logger.trace(
            f"{chrom=} {start=} {mode1=}  {current_node.strand=} {current_node.ref_start=} "
            f"{current_node.cigartuples_without_soft=} {query_name_current=} "
            f"{query_names_in_graph=} "
        )

        # reset query_names in graph
        if (
            current_node.ref_start is None
            or current_node.cigartuples_without_soft is None
        ):
            raise ValueError(f"{current_node.query_name} with None value")

        rescued_sr = self.cppext_rescuer.calculate_sr(
            chrom,
            start,
            start,
            mode1,
            current_node.strand,
            current_node.ref_start,
            query_name_current,
            current_node.cigartuples_without_soft,
        )

        self.logger.trace(f"current {rescued_sr=}")

        for next_node in current_node.successors:
            chrom_p, pos_p = next_node.get_breakpoint_depth_pos(mode2, "prev")
            next_node.prev_breakpoint_depth = self.cppext_rescuer.count_reads(
                chrom_p, pos_p, pos_p + 1
            )

            query_name_next = next_node.query_name.split(",")
            chrom, start = SRRescuer.obtain_region_for_rescue_sr(
                next_node.strand,
                next_node.chrom,
                next_node.exons,
                "prev_breakpoint",
                mode2,
            )
            self.logger.trace(
                f"{chrom=} {start=} {mode2=} {next_node.strand} {next_node.ref_start=} "
                f"{next_node.cigartuples_without_soft=} {query_name_next=}"
                f" {query_names_in_graph=}"
            )

            if (
                next_node.ref_start is None
                or next_node.cigartuples_without_soft is None
            ):
                raise ValueError(f"{next_node.query_name} with None value")

            rescued_sr += self.cppext_rescuer.calculate_sr(
                chrom,
                start,
                start,
                mode2,
                next_node.strand,
                next_node.ref_start,
                query_name_next,
                next_node.cigartuples_without_soft,
            )

            self.logger.trace(f"successors {rescued_sr=}")

        self.logger.trace(f"final {rescued_sr=}")
        if rescued_sr > 0:
            current_node.update_sr(rescued_sr)
