"""Tests for `scanmst.graph.merge_condition` and `scanmst.graph.cluster`."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pyfaidx import Fasta

from scanmst.base import BreakPoint, Intervals, MappingMode, MicroHomology, NovelInsertion
from scanmst.graph import ClusterFinder, EdgeData, NLPath, Node, NodeIdentity, VariationType
from scanmst.graph.cluster import (
    Ruler,
    create_merge_key_for_nlpath,
    create_merge_key_for_node,
    create_sort_key_by_merge_factor,
    export_connected_component_to_graph,
    merge_nlpath,
    merge_same_len_node_list,
    obtain_edge_info_signature_for_nlpath,
    sort_cluster,
)
from scanmst.graph.merge_condition import (
    MergeCondition,
    MergeConditionMode,
    find_shared_interval_indices,
    is_consecutive_from_beginning,
    is_consecutive_from_end,
    is_polya,
)
from scanmst.utils import change_dir


def _node(
    chrom: str,
    strand: str,
    exons: list[list[int]],
    identity: NodeIdentity,
    qname: str = "r1",
) -> Node:
    ivs = Intervals.from_list(exons)
    return Node(
        query_name=qname,
        chrom=chrom,
        strand=strand,
        ref_start=ivs.first.start,
        ref_end=ivs.last.end,
        identity=identity,
        exons=ivs,
    )


def _path(
    nodes: list[Node],
    sv_type: VariationType = VariationType.TRA,
    insertion_info: object = None,
) -> NLPath:
    edges = [
        EdgeData(
            variantion_type=sv_type,
            break_point1=BreakPoint(nodes[i].chrom, nodes[i].ref_end if nodes[i].strand.is_forward() else nodes[i].ref_start),
            break_point2=BreakPoint(
                nodes[i + 1].chrom,
                nodes[i + 1].ref_start if nodes[i + 1].strand.is_forward() else nodes[i + 1].ref_end,
            ),
            sr=1,
            read_ids=[nodes[0].query_name],
            mode1=MappingMode.MS,
            mode2=MappingMode.SM,
            gene1="G1",
            gene2="G2",
            annotation_code=3,
            splicing_code=1,
            insertion_info=insertion_info,
        )
        for i in range(len(nodes) - 1)
    ]
    return NLPath.from_nodes_and_edges_data(nodes, edges)


@pytest.fixture
def merge_fasta_path(tmp_path: Path) -> Path:
    """Create a synthetic FASTA with both bona-fide non-A genomic regions and poly-A/T internal-priming tracts."""
    seq = list("G" * 5000)
    seq[1500:1540] = list("A" * 40)
    seq[2460:2500] = list("T" * 40)
    fa = tmp_path / "merge_ref.fa"
    fa.write_text(">chr1\n" + "".join(seq) + "\n>chr2\n" + "".join(seq) + "\n>chr3\n" + "".join(seq) + "\n")
    return fa


@pytest.fixture
def merge_fasta(merge_fasta_path: Path) -> Fasta:
    return Fasta(str(merge_fasta_path), sequence_always_upper=True)


class TestMergeCondition:
    """Test MergeConditionMode, is_polya, interval helpers, and all 9 MergeCondition pairs."""

    def test_merge_condition_mode_all_nine_pairs(self) -> None:
        h = _node("chr1", "+", [[100, 200]], NodeIdentity.HEAD)
        m = _node("chr1", "+", [[100, 200]], NodeIdentity.MID)
        t = _node("chr1", "+", [[100, 200]], NodeIdentity.TAIL)
        expected = [
            (h, h, MergeConditionMode.head2head),
            (h, m, MergeConditionMode.head2mid),
            (h, t, MergeConditionMode.head2tail),
            (m, h, MergeConditionMode.mid2head),
            (m, m, MergeConditionMode.mid2mid),
            (m, t, MergeConditionMode.mid2tail),
            (t, h, MergeConditionMode.tail2head),
            (t, m, MergeConditionMode.tail2mid),
            (t, t, MergeConditionMode.tail2tail),
        ]
        for n1, n2, exp_mode in expected:
            assert MergeConditionMode.from_node(n1, n2) is exp_mode

    def test_interval_index_helpers_and_is_polya(self, merge_fasta: Fasta) -> None:
        idx_a, idx_b = find_shared_interval_indices([(0, 10), (20, 30)], [(20, 30), (40, 50)])
        assert (idx_a, idx_b) == ([1], [0])
        assert is_consecutive_from_beginning([0, 1, 2]) is True
        assert is_consecutive_from_beginning([1, 2]) is False
        assert is_consecutive_from_beginning([]) is False
        assert is_consecutive_from_end([1, 2], 3) is True
        assert is_consecutive_from_end([0, 1], 3) is False

        tail_real_fwd = _node("chr1", "+", [[200, 400]], NodeIdentity.TAIL)
        assert is_polya(tail_real_fwd, merge_fasta) is True

        tail_ip_fwd = _node("chr1", "+", [[1300, 1500]], NodeIdentity.TAIL)
        tail_ip_rev = _node("chr1", "-", [[2500, 2700]], NodeIdentity.TAIL)
        assert is_polya(tail_ip_fwd, merge_fasta) is False
        assert is_polya(tail_ip_rev, merge_fasta) is False

    def test_head_mid_and_tail_merge_modes(self, merge_fasta: Fasta) -> None:
        mc = MergeCondition(threshold=10, fasta=merge_fasta)
        assert mc.is_polya(_node("chr1", "+", [[200, 400]], NodeIdentity.TAIL)) is True

        h1 = _node("chr1", "+", [[100, 200], [300, 400], [500, 600]], NodeIdentity.HEAD)
        h2 = _node("chr1", "+", [[320, 400], [500, 598]], NodeIdentity.HEAD)
        h3 = _node("chr1", "+", [[510, 600]], NodeIdentity.HEAD)
        h_diff_chr = _node("chr2", "+", [[510, 600]], NodeIdentity.HEAD)
        assert mc.merged(h1, h2) is True
        assert mc.merged(h1, h3) is True
        assert mc.merged(h3, h1) is True
        assert mc.merged(h1, h_diff_chr) is False

        h_rev1 = _node("chr1", "-", [[100, 200], [300, 400]], NodeIdentity.HEAD)
        h_rev2 = _node("chr1", "-", [[103, 190]], NodeIdentity.HEAD)
        assert mc.merged(h_rev1, h_rev2) is True

        mid_full = _node("chr1", "+", [[100, 200], [300, 400], [500, 600]], NodeIdentity.MID)
        assert mc.merged(h2, mid_full) is True
        assert mc.merged(mid_full, h2) is True

        mid_twin = _node("chr1", "+", [[100, 200], [300, 400], [500, 600]], NodeIdentity.MID)
        mid_shifted = _node("chr1", "+", [[105, 200], [300, 400], [500, 600]], NodeIdentity.MID)
        assert mc.merged(mid_full, mid_twin) is True
        assert mc.merged(mid_full, mid_shifted) is False

        head_ov = _node("chr1", "+", [[1450, 1700]], NodeIdentity.HEAD)
        tail_ip = _node("chr1", "+", [[1300, 1500]], NodeIdentity.TAIL)
        assert mc.merged(head_ov, tail_ip) is True
        assert mc.merged(tail_ip, head_ov) is True

        head_sp = _node("chr1", "+", [[1320, 1350], [1420, 1650]], NodeIdentity.HEAD)
        tail_sp = _node("chr1", "+", [[1200, 1350], [1420, 1500]], NodeIdentity.TAIL)
        assert mc.merged(head_sp, tail_sp) is True

        mid_ip = _node("chr1", "+", [[1300, 1700]], NodeIdentity.MID)
        assert mc.merged(mid_ip, tail_ip) is True
        assert mc.merged(tail_ip, mid_ip) is True

        t_poly1 = _node("chr1", "+", [[200, 400]], NodeIdentity.TAIL)
        t_poly2 = _node("chr1", "+", [[204, 396]], NodeIdentity.TAIL)
        assert mc.merged(t_poly1, t_poly2) is True

        t_long_poly = _node("chr1", "+", [[1300, 1800]], NodeIdentity.TAIL)
        assert mc.merged(t_long_poly, tail_ip) is True

        # Additional reverse-strand and spliced tail/mid/head branches
        tail_ip_rev = _node("chr1", "-", [[2500, 2700]], NodeIdentity.TAIL)  # 2480..2500 is T-tract -> is_polya False
        head_ov_rev = _node("chr1", "-", [[2300, 2550]], NodeIdentity.HEAD)
        mid_ip_rev = _node("chr1", "-", [[2300, 2700]], NodeIdentity.MID)
        assert mc.merged(head_ov_rev, tail_ip_rev) is True
        assert mc.merged(tail_ip_rev, mid_ip_rev) is True

        # Both internal-primed tails on '+' and '-'
        tail_ip2 = _node("chr1", "+", [[1302, 1515]], NodeIdentity.TAIL)
        assert mc.merged(tail_ip, tail_ip2) is True
        tail_ip_rev2 = _node("chr1", "-", [[2485, 2698]], NodeIdentity.TAIL)
        assert mc.merged(tail_ip_rev, tail_ip_rev2) is True

        # Long non-polyA + short polyA returns False
        t_long_nonpoly = _node("chr1", "+", [[1100, 1510]], NodeIdentity.TAIL)
        t_short_poly = _node("chr1", "+", [[1100, 1250]], NodeIdentity.TAIL)
        assert mc.merged(t_long_nonpoly, t_short_poly) is False

        # Spliced tail2mid on '+' and '-'
        tail_sp_ip = _node("chr1", "+", [[1200, 1300], [1400, 1500]], NodeIdentity.TAIL)
        mid_sp_ip = _node("chr1", "+", [[1200, 1300], [1400, 1600]], NodeIdentity.MID)
        assert mc.merged(tail_sp_ip, mid_sp_ip) is True

        tail_sp_ip_rev = _node("chr1", "-", [[2500, 2600], [2700, 2800]], NodeIdentity.TAIL)
        mid_sp_ip_rev = _node("chr1", "-", [[2400, 2600], [2700, 2800]], NodeIdentity.MID)
        assert mc.merged(tail_sp_ip_rev, mid_sp_ip_rev) is True

        # Unspliced tail into spliced mid, and spliced tail into unspliced mid (rejected)
        assert mc.merged(_node("chr1", "+", [[1205, 1290]], NodeIdentity.TAIL), mid_sp_ip) is False
        assert mc.merged(_node("chr1", "+", [[1200, 1290]], NodeIdentity.TAIL), mid_sp_ip) is False
        assert mc.merged(_node("chr1", "+", [[1200, 1500]], NodeIdentity.TAIL), mid_sp_ip) is False
        assert mc.merged(tail_sp_ip, _node("chr1", "+", [[1200, 1600]], NodeIdentity.MID)) is False

        # Spliced head2mid on '-' and unspliced head into spliced mid on '+'/'-'
        head_sp_rev = _node("chr1", "-", [[2400, 2600], [2700, 2790]], NodeIdentity.HEAD)
        assert mc.merged(head_sp_rev, mid_sp_ip_rev) is True
        assert mc.merged(_node("chr1", "+", [[1420, 1600]], NodeIdentity.HEAD), mid_sp_ip) is True
        assert mc.merged(_node("chr1", "-", [[2400, 2580]], NodeIdentity.HEAD), mid_sp_ip_rev) is True
        assert mc.merged(head_sp_rev, _node("chr1", "-", [[2400, 2800]], NodeIdentity.MID)) is False

        # Spliced head2tail on '-' (2 introns each so shared intron is at end of node1 and start of node2) and mixed spliced/unspliced
        head_sp_rev_2i = _node("chr1", "-", [[2350, 2420], [2480, 2600], [2700, 2750]], NodeIdentity.HEAD)
        tail_sp_rev_2i = _node("chr1", "-", [[2500, 2600], [2700, 2760], [2820, 2900]], NodeIdentity.TAIL)
        assert mc.merged(head_sp_rev_2i, tail_sp_rev_2i) is True
        assert mc.merged(_node("chr1", "+", [[1450, 1700]], NodeIdentity.HEAD), tail_sp_ip) is True
        assert mc.merged(_node("chr1", "+", [[1450, 1550], [1600, 1700]], NodeIdentity.HEAD), tail_ip) is True
        assert mc.merged(_node("chr1", "-", [[2300, 2550]], NodeIdentity.HEAD), tail_sp_ip_rev) is True
        assert mc.merged(_node("chr1", "-", [[2300, 2400], [2450, 2550]], NodeIdentity.HEAD), tail_ip_rev) is True

        # Spliced tail2tail (long polyA vs short non-polyA and both non-polyA on '+' and '-')
        t_long_sp_poly = _node("chr1", "+", [[1200, 1300], [1400, 1800]], NodeIdentity.TAIL)
        assert mc.merged(t_long_sp_poly, tail_sp_ip) is True
        t_long_sp_poly_unspliced_pair = _node("chr1", "+", [[1300, 1520], [1600, 1800]], NodeIdentity.TAIL)
        assert mc.merged(t_long_sp_poly_unspliced_pair, tail_ip) is True
        t_long_sp_poly_rev = _node("chr1", "-", [[2200, 2600], [2700, 2800]], NodeIdentity.TAIL)
        assert mc.merged(t_long_sp_poly_rev, tail_sp_ip_rev) is True
        t_long_sp_poly_rev_unspliced_pair = _node("chr1", "-", [[2200, 2400], [2480, 2700]], NodeIdentity.TAIL)
        assert mc.merged(t_long_sp_poly_rev_unspliced_pair, tail_ip_rev) is True
        tail_sp_ip2 = _node("chr1", "+", [[1200, 1300], [1400, 1515]], NodeIdentity.TAIL)
        assert mc.merged(tail_sp_ip, tail_sp_ip2) is True
        tail_sp_ip_rev2 = _node("chr1", "-", [[2485, 2600], [2700, 2800]], NodeIdentity.TAIL)
        assert mc.merged(tail_sp_ip_rev, tail_sp_ip_rev2) is True


class TestRulerAndClusterFinder:
    """Test Ruler distance calculation and ClusterFinder clustering/merging."""

    def test_edge_signatures_and_merge_helpers(self, merge_fasta_path: Path) -> None:
        ruler = Ruler(prune_threshold=10, fasta_file=str(merge_fasta_path))
        assert repr(ruler) == "Ruler()"

        p_blunt = _path(
            [_node("chr1", "+", [[100, 200]], NodeIdentity.HEAD, "r1"), _node("chr2", "+", [[300, 400]], NodeIdentity.TAIL, "r1")],
            insertion_info=None,
        )
        p_mh = _path(
            [_node("chr1", "+", [[105, 200]], NodeIdentity.HEAD, "r2"), _node("chr2", "+", [[300, 398]], NodeIdentity.TAIL, "r2")],
            insertion_info=(False, MicroHomology("GG")),
        )
        p_ni = _path(
            [_node("chr3", "+", [[100, 200]], NodeIdentity.HEAD, "r3"), _node("chr3", "+", [[800, 900]], NodeIdentity.TAIL, "r3")],
            insertion_info=(False, NovelInsertion(0, "ACGT")),
        )

        assert obtain_edge_info_signature_for_nlpath(p_blunt) == 2
        assert obtain_edge_info_signature_for_nlpath(p_mh) == 1
        assert obtain_edge_info_signature_for_nlpath(p_ni) == 0

        assert create_merge_key_for_node(p_blunt[0]) == "chr1_[]"
        assert create_merge_key_for_nlpath(["chr1_[]", "chr2_[]"]) == "chr1_[]_chr2_[]"
        assert len(sort_cluster([p_mh, p_blunt], key=create_sort_key_by_merge_factor)) == 2

        assert ruler(p_blunt, p_mh) == 0.0
        assert ruler(p_blunt, p_ni) == 1.0

        p_twin = _path(
            [_node("chr1", "+", [[90, 200]], NodeIdentity.HEAD, "r4"), _node("chr2", "+", [[300, 405]], NodeIdentity.TAIL, "r4")],
            insertion_info=None,
        )
        assert merge_same_len_node_list(p_blunt, p_twin, 0, ruler.merge_condition) is True
        merge_nlpath(p_blunt, p_twin, 0)
        assert p_blunt[0].ref_start == 90
        assert p_blunt.edges[next(iter(p_blunt.edges))].sr == 2

    def test_cluster_finder_precomputed_and_on_the_fly(self, merge_fasta_path: Path, tmp_path: Path) -> None:
        p1 = _path(
            [_node("chr1", "+", [[100, 200]], NodeIdentity.HEAD, "r1"), _node("chr2", "+", [[300, 400]], NodeIdentity.TAIL, "r1")],
        )
        p2 = _path(
            [_node("chr1", "+", [[104, 200]], NodeIdentity.HEAD, "r2"), _node("chr2", "+", [[300, 398]], NodeIdentity.TAIL, "r2")],
        )
        p3 = _path(
            [_node("chr3", "+", [[700, 800]], NodeIdentity.HEAD, "r3"), _node("chr3", "+", [[900, 990]], NodeIdentity.TAIL, "r3")],
        )

        cf_fly = ClusterFinder([p1, p2, p3], prune_threshold=10, fasta_file=str(merge_fasta_path))
        with pytest.raises(ValueError, match="not precomputed"):
            cf_fly.get_distance(0, 1)

        comps_fly = [set(c) for c in cf_fly.find_cluster_index(use_precomputed=False)]
        assert {0, 1} in comps_fly
        assert {2} in comps_fly

        with change_dir(str(tmp_path)):
            export_connected_component_to_graph(cf_fly._graph, [p1, p2, p3], 1, {0, 1})
            exported = json.loads((tmp_path / "cluster_1.json").read_text())
            assert len(exported["elements"]["nodes"]) == 2

        p1_b = _path(
            [_node("chr1", "+", [[100, 200]], NodeIdentity.HEAD, "r1"), _node("chr2", "+", [[300, 400]], NodeIdentity.TAIL, "r1")],
        )
        p2_b = _path(
            [_node("chr1", "+", [[104, 200]], NodeIdentity.HEAD, "r2"), _node("chr2", "+", [[300, 398]], NodeIdentity.TAIL, "r2")],
        )
        p3_b = _path(
            [_node("chr3", "+", [[700, 800]], NodeIdentity.HEAD, "r3"), _node("chr3", "+", [[900, 990]], NodeIdentity.TAIL, "r3")],
        )
        cf_pre = ClusterFinder([p1_b, p2_b, p3_b], prune_threshold=10, fasta_file=str(merge_fasta_path))
        clusters = list(cf_pre.merge_cluster(use_precomputed=True, n_jobs=1))
        assert cf_pre.get_distance(0, 1) == 0.0
        assert cf_pre.get_distance(1, 0) == 0.0
        assert len(clusters) == 2
        merged_cluster = next(c for c in clusters if c[0].nodes[0].chrom == "chr1")
        assert len(merged_cluster) == 1
        assert merged_cluster[0].merge_factor == 2
