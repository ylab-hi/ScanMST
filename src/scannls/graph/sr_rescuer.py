"""SR rescuer."""

from __future__ import annotations

from typing import TYPE_CHECKING

from loguru import logger

from scannls import cppext

if TYPE_CHECKING:
    from scannls.base import MappingMode
    from scannls.graph import NLGraph, Node

MIN_SEQ_ALIGN_LEN = 10


def is_middle_node(node: Node) -> bool:
    """Check if the node is middle node."""
    return node.self_identity.is_mid()


def make_breakpoint(node: Node, mode: int) -> cppext.BreakPoint:
    """Make breakpoint."""
    return cppext.BreakPoint(
        node.ref_start,
        node.ref_end,
        mode,
        node.strand.is_reverse(),
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
            .min_seq_align_len(MIN_SEQ_ALIGN_LEN)
        )

        if average_read_depth is not None:
            options = options.average_read_depth(average_read_depth)

        self.cppext_rescuer = cppext.Rescuer(options)
        self.node_rescued_sr_maximum = node_rescued_sr_maximum
        self.cache = {}

    def init(self, graph) -> None:
        """Calculate the depth on breakpoints no matter rescue SR or not.

        changed in place

        :param nodes_in_graph: Series
        """
        for node in graph:
            self.update_breakpoint_depth(
                graph,
                node,
            )

    def __call__(self, graph) -> None:
        """Rescue SR from soft-clipped non-chimeric reads.

        changed in place

        :param nodes_in_graph: Series
        """
        query_names_in_graph = set()

        for node in graph:
            query_names_in_graph.update(node.read_ids)

        query_names_in_graph_list = list(query_names_in_graph)

        for node in graph:
            self.cppext_rescuer.reset_names_list(query_names_in_graph_list)
            self.update_sr(
                graph,
                node,
                query_names_in_graph_list,
                self.node_rescued_sr_maximum,
            )

        del query_names_in_graph

    @staticmethod
    def obtain_region_for_rescue_sr(node: Node, mode: MappingMode):
        """Obtain region from rescue."""
        next_pos = node.exons.last.end if node.strand.is_forward() else node.exons.first.start

        if mode.is_sm():
            next_pos += 1
        else:
            next_pos = max(next_pos, 1)

        return node.chrom, next_pos

    def update_predecessor_sr(
        self,
        current_node: Node,
        mode1: MappingMode,
    ) -> int:
        query_name_current = current_node.read_ids
        chrom, start = SRRescuer.obtain_region_for_rescue_sr(
            current_node,
            mode1,
        )

        logger.trace(
            f"{chrom=} {start=} {mode1=}  {current_node.strand=} {current_node.ref_start=} "
            f"{current_node.ref_end=} {is_middle_node(current_node)} "
            f"{current_node.cigartuples_without_soft=} {query_name_current=} ",
        )

        region = cppext.Region(chrom, start - 1, start)
        break_point = make_breakpoint(current_node, int(mode1))

        if current_node.cigartuples_without_soft is None:
            msg = f"{current_node.query_name} with None value"
            raise ValueError(msg)

        current_node_rescued_sr = self.cppext_rescuer.calculate_sr(
            region,
            break_point,
            query_name_current,
            current_node.cigartuples_without_soft,
        )

        logger.trace(f"{current_node_rescued_sr=}")
        return current_node_rescued_sr

    def update_sr(
        self,
        graph: NLGraph,
        current_node: Node,
        query_names_in_graph: list[str],
        node_rescued_sr_maximum: int,
    ) -> None:
        """Update SR for input node."""
        current_node_rescued_sr = 0
        rescued_pre = False

        for next_node in current_node.successors:
            edges = graph.find_edges(current_node, next_node)

            if len(edges) > 1:
                logger.warning("detect multiple edges")

            for edge in edges[:1]:
                mode1, mode2 = edge.modes

                if not rescued_pre:
                    rescued_pre = True
                    current_node_rescued_sr = self.update_predecessor_sr(
                        current_node,
                        mode1,
                    )

                edge.original_sr = edge.sr
                edge.sr += current_node_rescued_sr

                query_name_next = next_node.read_ids
                chrom = edge.break_point2.chrom
                start = edge.break_point2.pos
                if mode2.is_sm():
                    start += 1

                logger.trace(
                    f"{chrom=} {start=} {mode2=} {next_node.strand=} {next_node.ref_start=} {next_node.ref_end=}"
                    f"{is_middle_node(next_node)} "
                    f"{next_node.cigartuples_without_soft=} {query_name_next=}"
                    f" {query_names_in_graph=}",
                )

                if next_node.cigartuples_without_soft is None:
                    msg = f"{next_node.query_name} with None value"
                    raise ValueError(msg)

                region = cppext.Region(chrom, start - 1, start)
                break_point = make_breakpoint(next_node, mode2)

                increased_sr = self.cppext_rescuer.calculate_sr(
                    region,
                    break_point,
                    query_name_next,
                    next_node.cigartuples_without_soft,
                )

                logger.trace(f"edge {region.to_string().strip()}, rescue sr {increased_sr}")

                edge.sr += increased_sr

    def update_breakpoint_depth(
        self,
        graph: NLGraph,
        current_node: Node,
    ) -> None:
        """Update edge depth for input node."""
        rescued_pre = False

        for next_node in current_node.successors:
            edges = graph.find_edges(current_node, next_node)

            if len(edges) > 1:
                logger.warning("detect multiple edges")

            for edge in edges[:1]:
                mode1, mode2 = edge.modes

                edge.original_sr = edge.sr
                # update depth for breakpoint1 of edge
                chrom_n, pos_n = edge.break_point1.to_tuple()

                if mode1.is_ms():
                    pos_n -= 1

                edge.break_point1.depth = self.cppext_rescuer.count_reads(
                    chrom_n,
                    pos_n,
                    pos_n + 1,
                )

                # update depth for breakpoint2 of edge
                chrom_n, pos_n = edge.break_point2.to_tuple()

                if mode2.is_ms():
                    pos_n -= 1

                edge.break_point2.depth = self.cppext_rescuer.count_reads(
                    chrom_n,
                    pos_n,
                    pos_n + 1,
                )

                if next_node.cigartuples_without_soft is None:
                    msg = f"{next_node.query_name} with None value"
                    raise ValueError(msg)
