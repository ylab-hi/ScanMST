"""SR rescuer.

@Filename:    srRescuer.py
@author:      Yangyang Li
@Time:        12/30/21 15:00 PM
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Any

from loguru import logger

from scannls import cppext
from scannls.base import Mode
from scannls.exception import ExonsNotFoundError, ModesNotFoundError

if TYPE_CHECKING:
    from collections.abc import Iterable

    from scannls.graph import Node

# TODO: Rescure sr based on edge <05-01-23, Yangyang Li yangyang.li@northwestern.edu>


def is_middle_node(node: Node) -> bool:
    """Check if the node is middle node."""
    return node.predecessors and node.successors


def make_breakpoint(node: Node, mode: int) -> cppext.BreakPoint:
    """Make breakpoint."""
    return cppext.BreakPoint(
        node.ref_start,
        node.ref_end,
        mode,
        node.is_reverse(),
        is_middle_node(node),
    )


class SRRescuer:
    """Rescue SR from soft-clipped non-chimeric reads."""

    def __init__(
        self,
        input_bam_file: str,
        mapq_cutoff: int,
        soft_len_cutoff: int,
        mismatch_cutoff: int,
        alignment_frac: float,
        node_rescued_sr_maximum: int,
        average_read_depth: int | None,
    ) -> None:
        """Initialize Rescuer.

        :param logger: logger
        """
        options = (
            cppext.Options()
            .file(input_bam_file)
            .mapq(mapq_cutoff)
            .soft_len(soft_len_cutoff)
            .mismatch(mismatch_cutoff)
            .identity(alignment_frac)
            .min_seq_align_len(10)
        )

        if average_read_depth is not None:
            options = options.average_read_depth(average_read_depth)

        self.cppext_rescuer = cppext.Rescuer(options)
        self.node_rescued_sr_maximum = node_rescued_sr_maximum

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
            self.update_sr(
                node,
                query_names_in_graph_list,
                self.node_rescued_sr_maximum,
            )

        del query_names_in_graph

    @staticmethod
    def obtain_region_for_rescue_sr(
        strand: str | None,
        chrom: str | None,
        exons: Any,
        tgt_name: str,
        mode: int,
    ) -> tuple[str, int]:
        """Obtain target region (S-M boundary, M side) for rescuing SR purpose.

        ..note.
              Due to micro homology, prev_breakpoint/next_breakpoint locates inside the M side of S-M boundary
              Thus, exon start/end (S-M boundary) will be used to rescue SR.
        """
        if exons is None or strand is None or chrom is None:
            msg = f"chrom={chrom!r} strand={strand!r} exons={exons!r}"
            raise ExonsNotFoundError(msg)

        if strand == "+":
            pos = exons[-1][1] if tgt_name == "next_breakpoint" else exons[0][0]
        else:
            pos = exons[0][0] if tgt_name == "next_breakpoint" else exons[-1][1]

        if mode == 2:
            pos += 1
        else:
            pos = pos + 1 if pos == 0 else pos

        return chrom, pos

    @staticmethod
    def obtain_region_for_rescue_sr2(node: Node, mode: int, tag_name: str):
        """Obtain region from rescue."""
        if node.exons is None or node.strand is None or node.chrom is None:
            msg = f"node.chrom={node.chrom!r} node.strand={node.strand!r} node.exons={node.exons!r}"
            raise ExonsNotFoundError(msg)

        if node.strand == "+":
            pre_pos = node.exons[0][0]
            next_pos = node.exons[-1][1]
        else:
            pre_pos = node.exons[-1][1]
            next_pos = node.exons[0][0]

        if mode == Mode.Type2:
            pre_pos += 1
            next_pos += 1
        else:
            pre_pos = max(pre_pos, 1)
            next_pos = max(next_pos, 1)

        if node.prev_breakpoint is None:
            pre_pos = 0

        if node.next_breakpoint is None:
            next_pos = 0

        if tag_name == "prev_breakpoint":
            return node.chrom, pre_pos, next_pos

        return node.chrom, next_pos, pre_pos

    def update_sr(
        self,
        current_node: Node,
        query_names_in_graph: list[str],
        node_rescued_sr_maximum: int,
    ) -> None:
        """Update SR for input node."""
        if current_node.is_end_node():
            return

        if current_node.modes is None:
            msg = f"{current_node.query_name}"
            raise ModesNotFoundError(msg)

        mode1, mode2 = current_node.modes

        chrom_n, pos_n = current_node.get_breakpoint_depth_pos(mode1, "next")

        current_node.next_breakpoint_depth = self.cppext_rescuer.count_reads(
            chrom_n,
            pos_n,
            pos_n + 1,
        )

        query_name_current = current_node.query_name.split(",")

        chrom, start, check_pos = SRRescuer.obtain_region_for_rescue_sr2(
            current_node,
            mode1,
            "next_breakpoint",
        )

        logger.trace(
            f"{chrom=} {start=} {mode1=}  {current_node.strand=} {current_node.ref_start=} "
            f"{current_node.ref_end=} {is_middle_node(current_node)} "
            f"{current_node.cigartuples_without_soft=} {query_name_current=} "
            f"{query_names_in_graph=} ",
        )

        # reset query_names in graph
        if (
            current_node.ref_start is None
            or current_node.cigartuples_without_soft is None
        ):
            msg = f"{current_node.query_name} with None value"
            raise ValueError(msg)

        region = cppext.Region(chrom, start - 1, start)
        break_point = make_breakpoint(current_node, mode1)

        rescued_sr = self.cppext_rescuer.calculate_sr(
            region,
            break_point,
            query_name_current,
            current_node.cigartuples_without_soft,
        )

        logger.trace(f"current {rescued_sr=}")

        for next_node in current_node.successors:
            chrom_p, pos_p = next_node.get_breakpoint_depth_pos(mode2, "prev")
            next_node.prev_breakpoint_depth = self.cppext_rescuer.count_reads(
                chrom_p,
                pos_p,
                pos_p + 1,
            )

            query_name_next = next_node.query_name.split(",")

            chrom, start, check_pos = SRRescuer.obtain_region_for_rescue_sr2(
                next_node,
                mode2,
                "prev_breakpoint",
            )

            logger.trace(
                f"{chrom=} {start=} {mode2=} {next_node.strand=} {next_node.ref_start=} {next_node.ref_end=} "
                f"{is_middle_node(next_node)} "
                f"{next_node.cigartuples_without_soft=} {query_name_next=}"
                f" {query_names_in_graph=}",
            )

            if (
                next_node.ref_start is None
                or next_node.cigartuples_without_soft is None
            ):
                msg = f"{next_node.query_name} with None value"
                raise ValueError(msg)

            region = cppext.Region(chrom, start - 1, start)
            break_point = make_breakpoint(next_node, mode2)
            rescued_sr += self.cppext_rescuer.calculate_sr(
                region,
                break_point,
                query_name_next,
                next_node.cigartuples_without_soft,
            )

            logger.trace(f"successors {rescued_sr=}")
            if rescued_sr > node_rescued_sr_maximum:
                logger.trace(f"rescued SR has reached the higher bound, {rescued_sr=}")
                break

        logger.trace(f"final {rescued_sr=}")
        if rescued_sr > 0:
            current_node.update_sr(rescued_sr)
