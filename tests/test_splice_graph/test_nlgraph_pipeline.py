"""Tests for advanced `NLPath`, `NLGraph` construction/tracing/refinement, `annotate.py`, and `graphvis.py`."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

import HTSeq
import matplotlib
import pytest
from pyfaidx import Fasta

from scanmst.base import (
    BreakPoint,
    Event,
    Insertion,
    Intervals,
    MappingMode,
    MicroHomology,
    NovelInsertion,
    Read,
)
from scanmst.graph import (
    Edge,
    EdgeData,
    GraphCytoscapeExporter,
    GraphVis,
    MatplotlibVisualizeGraph,
    NLGraph,
    NLPath,
    Node,
    NodeIdentity,
    VariationType,
)
from scanmst.graph.annotate import annotate_graph, annotate_node, build_interval_tree, is_weakly_connected
from scanmst.graph.basic_graph import (
    _check_insertion_conditions_for_compare_insertion,
    merge_insertion,
    to_hash_identifier,
    to_numeric_identifier,
    update_node_with_other_node,
)
from scanmst.graph.graphvis import GraphJsonExporter, TSGraphExporter, create_nxgraph, visualize_graph_via_d3
from scanmst.graph.merge_condition import MergeCondition
from scanmst.utils import change_dir

matplotlib.use("Agg")


class _NoopRescuer:
    """Minimal rescuer double for pure graph construction/tracing tests."""

    def init(self, _graph: Any) -> None:
        return None

    def __call__(self, _graph: Any) -> None:
        return None


def _mk_node(
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


def _mk_edge_data(
    vtype: VariationType = VariationType.TRA,
    bp1: tuple[str, int] = ("chr1", 200),
    bp2: tuple[str, int] = ("chr2", 300),
    read_ids: list[str] | None = None,
    sr: int = 1,
    insertion_info: Any = None,
) -> EdgeData:
    return EdgeData(
        variantion_type=vtype,
        break_point1=BreakPoint(*bp1),
        break_point2=BreakPoint(*bp2),
        sr=sr,
        read_ids=list(read_ids or ["r1"]),
        mode1=MappingMode.MS,
        mode2=MappingMode.SM,
        gene1="G1",
        gene2="G2",
        annotation_code=3,
        splicing_code=1,
        insertion_info=insertion_info,
    )


@pytest.fixture
def graph_ref(tmp_path: Path) -> Fasta:
    fa = tmp_path / "graph_ref.fa"
    fa.write_text(">chr1\n" + ("G" * 3000) + "\n>chr2\n" + ("G" * 3000) + "\n>chr3\n" + ("G" * 3000) + "\n")
    # Seed MergeCondition's module-level SINGLE_FASTA_INSTANCE for NLGraph._compare_is_merged
    fasta = Fasta(str(fa), sequence_always_upper=True)
    MergeCondition(10, fasta)
    return fasta


class TestBasicGraphAndNLPathAdvanced:
    """Test VariationType, identifiers, insertion merging, NLPath.squeeze/is_forming_circle/new."""

    def test_variation_type_and_identifiers(self) -> None:
        for name, pred in [
            ("TRA", "is_tra"),
            ("DEL", "is_del"),
            ("TDUP", "is_tdup"),
            ("INV", "is_inv"),
            ("IDUP", "is_idup"),
        ]:
            vt = VariationType.from_str(name)
            assert str(vt) == name
            assert repr(vt) == name
            assert getattr(vt, pred)() is True
        with pytest.raises(ValueError, match="Invalid SV type"):
            VariationType.from_str("DUP")

        num_id = to_numeric_identifier("hello", 12)
        assert len(num_id) == 12
        assert num_id.isdigit()
        assert len(to_numeric_identifier("hello", None)) > 12
        with pytest.raises(TypeError):
            to_numeric_identifier(123)  # type: ignore[arg-type]
        with pytest.raises(TypeError):
            to_numeric_identifier("a", "10")  # type: ignore[arg-type]
        with pytest.raises(ValueError, match="positive"):
            to_numeric_identifier("a", 0)

        h_id = to_hash_identifier("Hello World!", 16)
        assert len(h_id) == 16
        assert h_id[0].isalpha()
        with pytest.raises(TypeError):
            to_hash_identifier(123)  # type: ignore[arg-type]
        with pytest.raises(TypeError):
            to_hash_identifier("x", 1.5)  # type: ignore[arg-type]
        with pytest.raises(ValueError, match="positive"):
            to_hash_identifier("x", -1)

    def test_node_union_merge_and_edge_insertion_rules(self) -> None:
        assert len(Node.union_of_intervals(Intervals([]))) == 0
        unioned = Node.union_of_intervals(Intervals.from_list([[100, 200], [150, 250], [300, 400]]))
        assert [(e.start, e.end) for e in unioned] == [(100, 250), (300, 400)]

        n1 = _mk_node("chr1", "+", [[100, 200]], NodeIdentity.HEAD, "r1")
        n2 = _mk_node("chr1", "+", [[120, 220]], NodeIdentity.HEAD, "r2")
        n1.set_up_breakpoints()
        n2.set_up_breakpoints()
        n1.merge(n2)
        assert (n1.ref_start, n1.ref_end) == (100, 220)
        assert set(n1.read_ids) == {"r1", "r2"}
        with pytest.raises(TypeError, match="Cannot merge"):
            n1.merge("not_a_node")  # type: ignore[arg-type]

        n1.cigartuples_without_soft = None
        n2.cigartuples_without_soft = [0, 100]
        update_node_with_other_node(n1, n2, ["cigartuples_without_soft"])
        assert n1.cigartuples_without_soft == [0, 100]

        # Insertion condition checks and merge_insertion precedence
        assert _check_insertion_conditions_for_compare_insertion(None, None) is True
        assert _check_insertion_conditions_for_compare_insertion((True, "x"), (True, "y")) is True
        assert _check_insertion_conditions_for_compare_insertion((False, MicroHomology("AA")), (False, MicroHomology("GG"))) is True
        assert _check_insertion_conditions_for_compare_insertion((False, NovelInsertion(0, "ATG")), (False, NovelInsertion(0, "ATG"))) is True
        assert _check_insertion_conditions_for_compare_insertion((False, NovelInsertion(0, "ATG")), (False, NovelInsertion(0, "CCC"))) is False

        e_ni = Edge("a", "b", _mk_edge_data(sr=1, insertion_info=(False, NovelInsertion(0, "AAAA"))))
        e_mh = Edge("a", "b", _mk_edge_data(sr=1, insertion_info=(False, MicroHomology("GG"))))
        merge_insertion(e_ni, e_mh)  # equal SR -> prefers MicroHomology over NovelInsertion
        assert isinstance(e_ni.insertion_info[1], MicroHomology)

        e_blunt = Edge("a", "b", _mk_edge_data(sr=1, insertion_info=None))
        merge_insertion(e_ni, e_blunt)  # equal SR -> prefers blunt end
        assert e_ni.insertion_info is None

    def test_nlpath_squeeze_circle_extension_and_new(self, graph_ref: Fasta) -> None:
        # 3-node path where hop 1 is DEL on '+' and hop 2 is TRA -> squeezes to 2-node TRA path
        n1 = _mk_node("chr1", "+", [[100, 200]], NodeIdentity.HEAD, "r1")
        n2 = _mk_node("chr1", "+", [[300, 400]], NodeIdentity.MID, "r1")
        n3 = _mk_node("chr2", "+", [[500, 600]], NodeIdentity.TAIL, "r1")
        ed_del = _mk_edge_data(VariationType.DEL, ("chr1", 200), ("chr1", 300), ["r1"])
        ed_tra = _mk_edge_data(VariationType.TRA, ("chr1", 400), ("chr2", 500), ["r2"])
        path = NLPath.from_nodes_and_edges_data([n1, n2, n3], [ed_del, ed_tra])

        path.anno_extension()
        assert path.extension is True  # {"r1"} & {"r2"} == set()

        path.squeeze()
        assert len(path) == 2
        assert [(e.start, e.end) for e in path[0].exons] == [(100, 200), (300, 400)]
        assert path[0].self_identity is NodeIdentity.HEAD
        assert path[1].self_identity is NodeIdentity.TAIL
        assert path.is_maximum_novel_insertion_length_valid(threshold=10) is True
        assert path.is_forming_circle(graph_ref, threshold=10) is False
        assert isinstance(hash(path), int)
        assert "NLPath(" in repr(path)

        # Circle detection via identical middle nodes
        cm1 = _mk_node("chr1", "+", [[100, 200]], NodeIdentity.HEAD, "r1")
        cm2 = _mk_node("chr2", "+", [[300, 400]], NodeIdentity.MID, "r1")
        cm3 = _mk_node("chr2", "+", [[300, 400]], NodeIdentity.MID, "r1")
        cm4 = _mk_node("chr3", "+", [[500, 600]], NodeIdentity.TAIL, "r1")
        cpath = NLPath.from_nodes_and_edges_data([cm1, cm2, cm3, cm4], [_mk_edge_data()] * 3)
        assert cpath.is_forming_circle(graph_ref, threshold=10) is True

        # NLPath.new with NovelInsertion (aligner=None) and 1-hit Insertion (aligner hit)
        r1 = Read.new("qn", "chr1", 100, "+", "50M60S", 60, 0, "A" * 110, [30] * 110)
        r2 = Read.new("qn", "chr2", 500, "+", "60S50M", 60, 0, "A" * 110, [30] * 110)
        evt_ins = Event(
            (
                "TRA",
                0,
                1,
                ("chr1:150", "chr2:500", MappingMode.MS, MappingMode.SM),
                (100, 150, [[100, 150]]),
                (500, 550, [[500, 550]]),
                ("+ACGTACGTAA", "+ACGTACGTAA"),
                ("+", "+"),
                ["G1", "G2"],
                False,
            ),
        )
        cvg = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        gene_iv = HTSeq.GenomicArrayOfSets("auto", stranded=False)
        nlpath_ni = NLPath.new([evt_ins], [r1, r2], 5, graph_ref, cvg, gene_iv, False, aligner=None)
        assert len(nlpath_ni) == 2
        assert nlpath_ni.is_maximum_novel_insertion_length_valid(threshold=5) is False
        assert nlpath_ni.is_maximum_novel_insertion_length_valid(threshold=20) is True

        class _OneHitAligner:
            def query_insertion(self, seq: str) -> tuple[bool, Insertion]:
                return True, Insertion(1, "chr3", 800, "+", f"{len(seq)}M", 60, 0, seq, [30] * len(seq))

        evt_ins2 = Event(
            (
                "TRA",
                0,
                1,
                ("chr1:150", "chr2:500", MappingMode.MS, MappingMode.SM),
                (100, 150, [[100, 150]]),
                (500, 550, [[500, 550]]),
                ("+ACGTACGTAA", "+ACGTACGTAA"),
                ("+", "+"),
                ["G1", "G2"],
                False,
            ),
        )
        nlpath_hit = NLPath.new([evt_ins2], [r1, r2], 5, graph_ref, cvg, gene_iv, False, aligner=_OneHitAligner())
        # 1-hit insertion mapped to chr3 inserts a MID node -> 3 nodes total!
        assert len(nlpath_hit) == 3
        assert nlpath_hit[1].chrom == "chr3"


class TestNLGraphAnnotateAndVis:
    """Test NLGraph construct/trace/refine/cycle handling, annotate.py, and graphvis.py."""

    def test_nlgraph_call_branching_refine_and_plot(self, graph_ref: Fasta, fake_logger, tmp_path: Path) -> None:
        # Two paths whose heads sit across a refine bin boundary (103 vs 105, refine_threshold=3)
        # sharing tail chr2:300-400
        p1_h = _mk_node("chr1", "+", [[103, 200]], NodeIdentity.HEAD, "r1")
        p1_t = _mk_node("chr2", "+", [[300, 400]], NodeIdentity.TAIL, "r1")
        p2_h = _mk_node("chr1", "+", [[105, 201]], NodeIdentity.HEAD, "r2")
        p2_t = _mk_node("chr2", "+", [[300, 400]], NodeIdentity.TAIL, "r2")

        path1 = NLPath.from_nodes_and_edges_data([p1_h, p1_t], [_mk_edge_data(read_ids=["r1"])])
        path2 = NLPath.from_nodes_and_edges_data([p2_h, p2_t], [_mk_edge_data(read_ids=["r2"])])
        path1.setup_breakpoints()
        path2.setup_breakpoints()

        nlgraph = NLGraph(
            logger=fake_logger,
            rescuer=_NoopRescuer(),
            merge_threshold=0,  # force construct() to keep p1_h and p2_h separate so refine() merges them!
            support_reads=1,
            input_bam_path=Path("sample.bam"),
            output_dir=tmp_path,
            rescue_sr=False,
            if_refine=True,
            refine_threshold=3,
        )
        traced = list(nlgraph([path1, path2], "TSG0000000001", is_plot=True))
        assert len(traced) == 1
        assert traced[0][0].ptc == 1
        assert traced[0][0].ptf == pytest.approx(1.0)
        assert (tmp_path / "graph_sample" / "sample_TSG0000000001_cy.json").exists()
        nlgraph.polish_edges()

        # Complex refine() with both common and unique predecessors & successors, plus remove_node fallbacks
        g_ref = NLGraph(
            logger=fake_logger,
            rescuer=_NoopRescuer(),
            merge_threshold=0,
            support_reads=1,
            input_bam_path=Path("ref.bam"),
            output_dir=tmp_path,
            rescue_sr=False,
            if_refine=True,
            refine_threshold=5,
        )
        g_ref.nodes = g_ref.dict_factory()
        g_ref.edges = defaultdict(list)

        h_common = _mk_node("chr1", "+", [[100, 200]], NodeIdentity.HEAD, "r1")
        h_common.read_ids.append("r2")
        h_common.identities["r2"] = NodeIdentity.HEAD
        h_uniq = _mk_node("chr1", "+", [[50, 150]], NodeIdentity.HEAD, "r3")
        m_pri = _mk_node("chr2", "+", [[300, 400]], NodeIdentity.MID, "r1")
        m_pri.read_ids.append("r4")
        m_pri.identities["r4"] = NodeIdentity.MID
        m_sec = _mk_node("chr2", "+", [[302, 401]], NodeIdentity.MID, "r2")
        m_sec.read_ids.append("r3")
        m_sec.identities["r3"] = NodeIdentity.MID
        t_common = _mk_node("chr3", "+", [[500, 600]], NodeIdentity.TAIL, "r1")
        t_common.read_ids.append("r2")
        t_common.identities["r2"] = NodeIdentity.TAIL
        t_uniq = _mk_node("chr3", "+", [[700, 800]], NodeIdentity.TAIL, "r3")

        for n in (h_common, h_uniq, m_pri, m_sec, t_common, t_uniq):
            n.is_in_graph = True
            n.set_up_breakpoints()
            g_ref.add_node_with_similar_key(n)

        h_common.add_successor(m_pri, g_ref, _mk_edge_data(bp1=("chr1", 200), bp2=("chr2", 300), read_ids=["r1"]))
        # Calling add_edge again with compatible EdgeData exercises NLGraph.add_edge merging
        g_ref.add_edge(h_common, m_pri, _mk_edge_data(bp1=("chr1", 200), bp2=("chr2", 300), read_ids=["r4"]))
        h_common.add_successor(m_sec, g_ref, _mk_edge_data(bp1=("chr1", 200), bp2=("chr2", 302), read_ids=["r2"]))
        h_uniq.add_successor(m_sec, g_ref, _mk_edge_data(bp1=("chr1", 150), bp2=("chr2", 302), read_ids=["r3"]))
        m_pri.add_successor(t_common, g_ref, _mk_edge_data(bp1=("chr2", 400), bp2=("chr3", 500), read_ids=["r1"]))
        m_sec.add_successor(t_common, g_ref, _mk_edge_data(bp1=("chr2", 401), bp2=("chr3", 500), read_ids=["r2"]))
        m_sec.add_successor(t_uniq, g_ref, _mk_edge_data(bp1=("chr2", 401), bp2=("chr3", 700), read_ids=["r3"]))

        g_ref.refine()
        assert m_sec not in g_ref
        assert h_uniq in m_pri.predecessors
        assert t_uniq in m_pri.successors

        # remove_node when node is not in graph logs warning cleanly
        ghost = _mk_node("chr9", "+", [[10, 20]], NodeIdentity.HEAD, "rg")
        g_ref.remove_node(ghost)

    def test_nlgraph_read_guided_mid_filter_and_cycle_handling(self, fake_logger, tmp_path: Path) -> None:
        # Read-guided MID check: two MID-only edges from disjoint reads ("r1" vs "r2") cannot cross
        assert (
            NLGraph.determine_edge(
                {NodeIdentity.MID: ["r1"]},
                {NodeIdentity.MID: ["r2"]},
            )
            is False
        )
        assert (
            NLGraph.determine_edge(
                {NodeIdentity.MID: ["r1", "r2"]},
                {NodeIdentity.MID: ["r2"]},
            )
            is True
        )

        # Graph with a cycle: n1 -> n2 -> n3 -> n2
        n1 = _mk_node("chr1", "+", [[100, 200]], NodeIdentity.HEAD, "r1")
        n2 = _mk_node("chr2", "+", [[300, 400]], NodeIdentity.MID, "r1")
        n3 = _mk_node("chr3", "+", [[500, 600]], NodeIdentity.TAIL, "r1")
        n1.set_up_breakpoints()
        n3.set_up_breakpoints()

        g = NLGraph(
            logger=fake_logger,
            rescuer=_NoopRescuer(),
            merge_threshold=10,
            support_reads=1,
            input_bam_path=Path("cyc.bam"),
            output_dir=tmp_path,
            rescue_sr=False,
            ignore_circle=False,
        )
        g.nodes = g.dict_factory()
        g.edges = defaultdict(list)
        for n in (n1, n2, n3):
            n.is_in_graph = True
            g.add_node_with_similar_key(n)
        n1.add_successor(n2, g, _mk_edge_data(read_ids=["r1"]))
        n2.add_successor(n3, g, _mk_edge_data(read_ids=["r1"]))
        n3.add_successor(n2, g, _mk_edge_data(read_ids=["r1"]))

        assert g.trace() == []
        assert g.has_circle is True

    def test_annotate_and_graphvis_exporters(self, fake_logger, tmp_path: Path) -> None:
        empty_g = NLGraph(fake_logger, _NoopRescuer(), 10, 1, Path("x.bam"), tmp_path, rescue_sr=False)
        empty_g.nodes = {}
        empty_g.edges = defaultdict(list)
        assert is_weakly_connected(empty_g) is True

        n1 = _mk_node("chr1", "+", [[100, 200]], NodeIdentity.HEAD, "r1")
        n2 = _mk_node("chr2", "+", [[300, 400]], NodeIdentity.TAIL, "r1")
        n3_isolated = _mk_node("chr3", "+", [[500, 600]], NodeIdentity.HEAD, "r2")
        for n in (n1, n2, n3_isolated):
            n.is_in_graph = True

        g = NLGraph(fake_logger, _NoopRescuer(), 10, 1, Path("x.bam"), tmp_path, rescue_sr=False)
        g.nodes = {}
        g.edges = defaultdict(list)
        g.add_node_with_similar_key(n1)
        g.add_node_with_similar_key(n2)
        n1.add_successor(
            n2,
            g,
            _mk_edge_data(insertion_info=(False, NovelInsertion(0, "ACGT"))),
        )
        assert is_weakly_connected(g) is True

        g.add_node_with_similar_key(n3_isolated)
        assert is_weakly_connected(g) is False
        g.remove_node(n3_isolated)

        # Gene annotation from GTF
        gtf = tmp_path / "genes.gtf"
        gtf.write_text(
            'chr1\tsrc\tgene\t50\t500\t.\t+\t.\tgene_id "GENE_ALPHA";\nchr2\tsrc\tgene\t250\t500\t.\t+\t.\tgene_id "GENE_ALPHA";\n',
        )
        trees = build_interval_tree(gtf)
        annotate_node(n1, trees)
        assert n1.gene_names == ["GENE_ALPHA"]
        top_gene = annotate_graph(g, n1, gtf)
        assert top_gene == "GENE_ALPHA"

        # GraphVis exporters and loaders
        nx_g = create_nxgraph(g, min_support_reads=1, possible_paths={"TSP1": [n1.id, g.find_edges(n1, n2)[0].id, n2.id]})
        assert len(nx_g.nodes) == 2

        with change_dir(str(tmp_path)):
            cy_exp = GraphCytoscapeExporter("demo")
            json_exp = GraphJsonExporter("demo", {"edges": "edges", "source": "source", "target": "target"})
            tsg_exp = TSGraphExporter("demo")
            vis = GraphVis(g, [cy_exp], min_support_reads=1)
            vis.register(json_exp)
            vis.visualize(visitors=[tsg_exp])

            loaded_cy = GraphVis.load_cytoscape(tmp_path / "demo_cy.json")
            assert len(loaded_cy.nodes) == 2
            assert (tmp_path / "graph_demo.json").exists()
            tsg_text = (tmp_path / "demo.tsg").read_text()
            assert "novel_insertion:Z:ACGT" in tsg_text

            # Matplotlib exporter requires 'label' attribute on edges for draw_multiedge_labels
            for _u, _v, _k, d in nx_g.edges(keys=True, data=True):
                d["label"] = d["id"]
            MatplotlibVisualizeGraph("demo_fig", support_reads=1).visit(nx_g)
            assert (tmp_path / "graph_demo_fig.png").exists()

        with pytest.raises(FileNotFoundError):
            GraphVis.load(tmp_path / "missing.json")
        with pytest.raises(FileNotFoundError):
            GraphVis.load_cytoscape(tmp_path / "missing_cy.json")
        with pytest.raises(NotImplementedError):
            visualize_graph_via_d3(nx_g, "d3")
