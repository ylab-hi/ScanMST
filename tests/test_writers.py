"""Unit tests for `scanmst.writer` (`Writers`, `FastaWriter`, `GTFWriter`, `VCFWriter`, `TSGWriter`)."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

import pytest
from pyfaidx import FastaNotFoundError

from scanmst.base import BreakPoint, Intervals, MappingMode, MicroHomology, NovelInsertion
from scanmst.exception import ExonsNotFoundError
from scanmst.graph import EdgeData, NLGraph, NLPath, Node, NodeIdentity, VariationType
from scanmst.writer import FastaWriter, GTFWriter, VCFWriter, Writers
from scanmst.writer.gtf_writer import add_info_to_attribute_column
from scanmst.writer.tsg_writer import TSGWriter
from scanmst.writer.vcf_writer import obtain_reference_from_bam_header, obtain_sequence_from_insertion


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


@pytest.fixture
def writer_ref(tmp_path: Path) -> Path:
    # chr1: 0..10 = "AAACCCGGGT", 20..30 = "TTTAAAGGGC"
    # chr2: 0..10 = "ACGTACGTAA" (reverse complement = "TTACGTACGT")
    chr1 = "AAACCCGGGT" + ("N" * 10) + "TTTAAAGGGC" + ("G" * 100)
    chr2 = "ACGTACGTAA" + ("C" * 100)
    ref = tmp_path / "writer_ref.fa"
    ref.write_text(f">chr1\n{chr1}\n>chr2\n{chr2}\n")
    return ref


@pytest.fixture
def bam_header() -> dict[str, Any]:
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 200}, {"SN": "chr2", "LN": 200}],
        "PG": [{"ID": "minimap2", "CL": "minimap2 -ax splice -Y --cs ref.fa reads.fq"}],
    }


class TestFastaAndGTFWriters:
    """Test FastaWriter (reference & consensus modes) and GTFWriter."""

    def test_fasta_writer_reference_and_consensus(self, tmp_path: Path, writer_ref: Path) -> None:
        with pytest.raises(FastaNotFoundError):
            FastaWriter(str(tmp_path / "out.fa"), str(tmp_path / "missing.fa"), "reference", {})

        n1 = _mk_node("chr1", "+", [[0, 10], [20, 30]], NodeIdentity.HEAD, "r1")
        n2 = _mk_node("chr2", "-", [[0, 10]], NodeIdentity.TAIL, "r1")
        ed = EdgeData(
            variantion_type=VariationType.TRA,
            break_point1=BreakPoint("chr1", 30),
            break_point2=BreakPoint("chr2", 10),
            sr=2,
            read_ids=["r1", "r2"],
            mode1=MappingMode.MS,
            mode2=MappingMode.MS,
            gene1="G1",
            gene2="G2",
            annotation_code=3,
            splicing_code=1,
            insertion_info=(False, NovelInsertion(0, "CCCC")),
        )
        path = NLPath.from_nodes_and_edges_data([n1, n2], [ed])

        ref_fa_out = tmp_path / "ref_mode.fa"
        fw_ref = FastaWriter(str(ref_fa_out), str(writer_ref), "reference", {})
        fw_ref.open()
        fw_ref.write_data(path, "TSG0000000001")
        fw_ref.close()

        lines = ref_fa_out.read_text().splitlines()
        # n1: "AAACCCGGGT" + "TTTAAAGGGC" + novel_insertion "CCCC" = 24bp; n2 (-): rc("ACGTACGTAA") = "TTACGTACGT" = 10bp
        assert lines[0] == f">{path.id} 24|10"
        assert lines[1] == "AAACCCGGGTTTTAAAGGGCCCCCTTACGTACGT"

        # Consensus mode via pyabpoa MSA
        cons_fa_out = tmp_path / "cons_mode.fa"
        seq_dict = {
            "r1": "AAACCCGGGTTTTAAAGGGCCCCCTTACGTACGT",
            "r2": "AAACCCGGGTTTTAAAGGGCCCCCTTACGTACGT",
        }
        fw_cons = FastaWriter(str(cons_fa_out), str(writer_ref), "consensus", seq_dict)
        fw_cons.open()
        fw_cons.write_data(path, "TSG0000000001")
        fw_cons.close()
        cons_lines = cons_fa_out.read_text().splitlines()
        assert cons_lines[0] == f">{path.id} 24|10"
        assert cons_lines[1] == "AAACCCGGGTTTTAAAGGGCCCCCTTACGTACGT"

    def test_gtf_writer_features_and_missing_exons_error(self, tmp_path: Path) -> None:
        assert add_info_to_attribute_column(["chr1", "attr;"], ' extra "1";')[-1] == 'attr; extra "1";'

        n1 = _mk_node("chr1", "-", [[0, 10], [20, 30]], NodeIdentity.HEAD, "r1")
        n2 = _mk_node("chr2", "+", [[50, 80]], NodeIdentity.TAIL, "r1")
        ed = EdgeData(
            variantion_type=VariationType.TRA,
            break_point1=BreakPoint("chr1", 0),
            break_point2=BreakPoint("chr2", 50),
            sr=3,
            read_ids=["r1"],
            mode1=MappingMode.SM,
            mode2=MappingMode.SM,
            gene1="G1",
            gene2="G2",
            annotation_code=3,
            splicing_code=1,
            insertion_info=(False, NovelInsertion(0, "ATAT")),
            original_sr=1,
        )
        path = NLPath.from_nodes_and_edges_data([n1, n2], [ed])

        gtf_out = tmp_path / "out.gtf"
        gw = GTFWriter(str(gtf_out), rescue_sr=True)
        gw.open()
        gw.write_data(path, "TSG0000000001")
        gw.close()

        rows = [line.split("\t") for line in gtf_out.read_text().splitlines()]
        # 1 transcript + 2 exons of n1 (reversed on '-') + 1 NovelInsertion + 1 exon of n2 = 5 rows
        assert len(rows) == 5
        assert all(len(r) == 9 for r in rows)
        assert rows[0][2] == "transcript"
        assert 'sr "3"; osr "1";' in rows[0][8]
        # Reverse strand exons are emitted in 5'->3' transcription order: [20,30] (-> 21..30) then [0,10] (-> 1..10)
        assert (rows[1][2], rows[1][3], rows[1][4]) == ("exon", "21", "30")
        assert (rows[2][2], rows[2][3], rows[2][4]) == ("exon", "1", "10")
        assert rows[3][2] == "insertion"
        assert 'sequence "ATAT";' in rows[3][8]

        n_bad = _mk_node("chr1", "+", [[0, 10]], NodeIdentity.HEAD, "r1")
        n_bad.exons = None  # type: ignore[assignment]
        bad_path = NLPath([n_bad])
        gw2 = GTFWriter(str(tmp_path / "bad.gtf"), rescue_sr=False)
        gw2.open()
        with pytest.raises(ExonsNotFoundError):
            gw2.write_data(bad_path, "TSG0000000001")
        gw2.close()


class TestVCFAndTSGWritersAndManager:
    """Test VCFWriter shared-hop aggregation, TSGWriter, and Writers context manager."""

    def test_vcf_header_helpers_and_multicluster_y_aggregation(
        self,
        tmp_path: Path,
        writer_ref: Path,
        bam_header: dict[str, Any],
    ) -> None:
        with pytest.raises(FastaNotFoundError):
            VCFWriter(str(tmp_path / "x.vcf"), False, str(tmp_path / "missing.fa"), bam_header)

        assert obtain_reference_from_bam_header(bam_header).startswith("minimap2")
        assert obtain_reference_from_bam_header({"PG": [{"ID": "custom_tool", "CL": "foo"}]}) == "Unknown"

        rev_node = _mk_node("chr1", "-", [[10, 50]], NodeIdentity.HEAD, "r1")
        assert obtain_sequence_from_insertion(MicroHomology("AACG"), rev_node) == "CGTT"
        assert obtain_sequence_from_insertion(NovelInsertion(0, ""), rev_node) == ""

        # Y-shape in cluster 1: two NLPaths sharing the exact same chr1:50 -> chr2:10 TDUP/TRA hop
        p1_n1 = _mk_node("chr1", "+", [[10, 50]], NodeIdentity.HEAD, "r1")
        p1_n2 = _mk_node("chr2", "+", [[10, 60]], NodeIdentity.TAIL, "r1")
        p2_n1 = _mk_node("chr1", "+", [[15, 50]], NodeIdentity.HEAD, "r2")
        p2_n2 = _mk_node("chr2", "+", [[10, 70]], NodeIdentity.TAIL, "r2")

        ed1 = EdgeData(
            VariationType.TRA,
            BreakPoint("chr1", 50, depth=10),
            BreakPoint("chr2", 10, depth=12),
            sr=2,
            read_ids=["r1"],
            mode1=MappingMode.MS,
            mode2=MappingMode.SM,
            gene1="G1",
            gene2="G2",
            annotation_code=3,
            splicing_code=1,
            insertion_info=(False, MicroHomology("GG")),
            original_sr=1,
        )
        ed2 = EdgeData(
            VariationType.TRA,
            BreakPoint("chr1", 50, depth=10),
            BreakPoint("chr2", 10, depth=12),
            sr=3,
            read_ids=["r2"],
            mode1=MappingMode.MS,
            mode2=MappingMode.SM,
            gene1="G1",
            gene2="G2",
            annotation_code=3,
            splicing_code=1,
            insertion_info=(False, MicroHomology("GG")),
            original_sr=2,
        )
        path1 = NLPath.from_nodes_and_edges_data([p1_n1, p1_n2], [ed1])
        path2 = NLPath.from_nodes_and_edges_data([p2_n1, p2_n2], [ed2])

        # Second cluster path (TDUP -> ICRL) triggers cluster-boundary flush
        c2_n1 = _mk_node("chr1", "+", [[60, 90]], NodeIdentity.HEAD, "r3")
        c2_n2 = _mk_node("chr1", "+", [[10, 40]], NodeIdentity.TAIL, "r3")
        ed_c2 = EdgeData(
            VariationType.TDUP,
            BreakPoint("chr1", 90, depth=8),
            BreakPoint("chr1", 10, depth=8),
            sr=1,
            read_ids=["r3"],
            mode1=MappingMode.MS,
            mode2=MappingMode.SM,
            gene1="G1",
            gene2="G1",
            annotation_code=3,
            splicing_code=0,
            insertion_info=(False, NovelInsertion(0, "TTAA")),
            original_sr=1,
        )
        path_c2 = NLPath.from_nodes_and_edges_data([c2_n1, c2_n2], [ed_c2])

        vcf_out = tmp_path / "sample.vcf"
        vw = VCFWriter(str(vcf_out), rescue_sr=True, reference=str(writer_ref), bam_header=bam_header)
        vw.open()
        vw.write_data(path1, "1")
        vw.write_data(path2, "1")
        vw.write_data(path_c2, "2")
        vw.close()

        text = vcf_out.read_text()
        assert "##fileformat=VCFv4.3" in text
        assert "##contig=<ID=chr1,length=200>" in text
        data_rows = [line.split("\t") for line in text.splitlines() if not line.startswith("#")]
        assert len(data_rows) == 2
        # Row 0: aggregated Y-shape TRA (ITPL) with SR = 2 + 3 = 5, OSR = 1 + 2 = 3, HOMSEQ=GG
        assert data_rows[0][4] == "<ITPL>"
        assert "SR=5;OSR=3;" in data_rows[0][7]
        assert "HOMSEQ=GG;INSSEQ=.;" in data_rows[0][7]
        assert f"{path1.id},{path2.id}" in data_rows[0][7]
        # Row 1: cluster 2 TDUP (ICRL) NONCANONICAL with INSSEQ=TTAA
        assert data_rows[1][4] == "<ICRL>"
        assert data_rows[1][7].startswith("NONCANONICAL;LINKTYPE=ICRL;")
        assert "HOMSEQ=.;INSSEQ=TTAA;" in data_rows[1][7]

    def test_writers_context_manager_and_tsg_writer(
        self,
        tmp_path: Path,
        writer_ref: Path,
        bam_header: dict[str, Any],
        fake_logger,
    ) -> None:
        n1 = _mk_node("chr1", "+", [[0, 10]], NodeIdentity.HEAD, "r1")
        n2 = _mk_node("chr2", "+", [[0, 10]], NodeIdentity.TAIL, "r1")
        n1.is_in_graph = True
        n2.is_in_graph = True
        ed = EdgeData(
            VariationType.TRA,
            BreakPoint("chr1", 10),
            BreakPoint("chr2", 0),
            sr=1,
            read_ids=["r1"],
            mode1=MappingMode.MS,
            mode2=MappingMode.SM,
            gene1="G1",
            gene2="G2",
            annotation_code=3,
            splicing_code=1,
            insertion_info=(False, MicroHomology("AC")),
        )
        path = NLPath.from_nodes_and_edges_data([n1, n2], [ed])

        nlgraph = NLGraph(fake_logger, None, 10, 1, Path("in.bam"), tmp_path, rescue_sr=False)
        nlgraph.nodes = {}
        nlgraph.edges = defaultdict(list)
        nlgraph.add_node_with_similar_key(n1)
        nlgraph.add_node_with_similar_key(n2)
        n1.add_successor(n2, nlgraph, ed)
        nlgraph.possible_paths = {path.id: [n1.id, nlgraph.find_edges(n1, n2)[0].id, n2.id]}  # type: ignore[assignment]

        fa_w = FastaWriter(str(tmp_path / "all.fa"), str(writer_ref), "reference", {})
        gtf_w = GTFWriter(str(tmp_path / "all.gtf"), rescue_sr=False)
        vcf_w = VCFWriter(str(tmp_path / "all.vcf"), rescue_sr=False, reference=str(writer_ref), bam_header=bam_header)
        tsg_w = TSGWriter(str(tmp_path / "all.tsg"))

        manager = Writers((fa_w, gtf_w, vcf_w, tsg_w))
        with manager.open():
            manager.write_path(path, "TSG0000000001")
            manager.write_graph(nlgraph, "TSG0000000001")

        tsg_content = (tmp_path / "all.tsg").read_text()
        assert "H\tTSG\t1.0" in tsg_content
        assert "G\tTSG0000000001" in tsg_content
        assert f"P\t{path.id}\t" in tsg_content
        assert "microhomology:Z:AC" in tsg_content
        assert (tmp_path / "all.fa").read_text().startswith(">")
        assert (tmp_path / "all.gtf").read_text() != ""
        assert "#CHROM\tPOS" in (tmp_path / "all.vcf").read_text()
