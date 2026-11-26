"""Unit tests for `scanmst.cli` (`arg.py`, `cli.py`, and `main.py` validation/helpers)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pysam
import pytest

from scanmst.base import Event, Intervals, MappingMode
from scanmst.cli import DefaultOptions, cli
from scanmst.cli.arg import RichArgParser, parse_args, print_banner
from scanmst.cli.cli import get_writers
from scanmst.cli.main import (
    BamScanner,
    _get_cvg_gene_iv,
    _get_genome_fasta,
    _get_read_matched_sequence,
    update_position_event_list,
)


def _write_bam(
    bam_path: Path,
    *,
    sorted_header: bool = True,
    reads: list[dict[str, Any]] | None = None,
) -> Path:
    hd: dict[str, str] = {"VN": "1.6"}
    if sorted_header:
        hd["SO"] = "coordinate"
    header = pysam.AlignmentHeader.from_dict({"HD": hd, "SQ": [{"SN": "chr1", "LN": 10000}]})
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for spec in reads or []:
            seg = pysam.AlignedSegment(header)
            seg.query_name = spec.get("qname", "r1")
            seg.reference_id = 0
            seg.reference_start = spec.get("pos", 100)
            seg.cigarstring = spec.get("cigar", "50M")
            seg.query_sequence = spec.get("seq", "A" * 50)
            seg.mapping_quality = spec.get("mapq", 60)
            seg.flag = spec.get("flag", 0)
            for tag_k, tag_v in spec.get("tags", {}).items():
                seg.set_tag(tag_k, tag_v)
            out.write(seg)
    if sorted_header:
        pysam.index(str(bam_path))
    return bam_path


class TestCliArgAndWritersSetup:
    """Test argument parser flags, RichArgParser formatting, and get_writers."""

    def test_parse_args_flags_and_banner(self) -> None:
        print_banner()
        colored = RichArgParser._color_message("Use --input and -h for help")
        assert "[bold green]--input[/]" in colored
        assert "[bold green]-h[/]" in colored

        parser = parse_args()
        ns = parser.parse_args(
            [
                "--input",
                "in.bam",
                "--ref",
                "ref.fa",
                "--gtf",
                "genes.gtf",
                "--output",
                "out/prefix",
                "--output-seq",
                "both",
                "--sr",
                "3",
                "--ncan",
                "--off-exon-filter",
                "--rescue-sr",
                "--nbound",
                "--ignore-circle",
                "--refine",
                "--refine-threshold",
                "5",
                "--species",
                "mouse",
                "--blat-nclosed",
                "--blat-nsleep",
                "--log-level",
                "info",
            ],
        )
        assert ns.output_sequence_choice == "both"
        assert ns.support_reads == 3
        assert ns.noncanonical is True
        assert ns.exon_filter is False
        assert ns.rescue_sr is True
        assert ns.bound is False
        assert ns.ignore_circle is True
        assert ns.refine is True
        assert ns.refine_threshold == 5
        assert ns.species == "mouse"
        assert ns.blat_closed is False
        assert ns.blat_sleep is False
        assert ns.log == "info"

    def test_get_writers_single_vs_both(self, tmp_path: Path) -> None:
        ref = tmp_path / "ref.fa"
        ref.write_text(">chr1\nACGTACGT\n")
        hdr = {"SQ": [{"SN": "chr1", "LN": 8}], "PG": []}

        w_single = get_writers(str(tmp_path / "s"), str(ref), False, "reference", {}, hdr)
        assert len(w_single.writers_list) == 3

        w_both = get_writers(str(tmp_path / "b"), str(ref), True, "both", {}, hdr)
        assert len(w_both.writers_list) == 4


class TestBamScannerAndMainHelpers:
    """Test BamScanner validation checks, main.py helpers, and cli error exits."""

    def test_bam_scanner_sort_softclip_and_cs_checks(self, tmp_path: Path, fake_logger) -> None:
        ref = tmp_path / "ref.fa"
        ref.write_text(">chr1\n" + ("A" * 1000) + "\n")
        gtf = tmp_path / "g.gtf"
        gtf.write_text('chr1\tsrc\texon\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')

        unsorted = _write_bam(tmp_path / "unsorted.bam", sorted_header=False)
        with pytest.raises(RuntimeError, match="not sorted"):
            BamScanner(unsorted, 20, ref, gtf, 5, fake_logger, True, 100, 200, 10, 0.05, 0.01, None)

        # Supplementary alignment (0x800 = 2048) without soft-clipping -> SystemExit (-Y required)
        no_soft = _write_bam(
            tmp_path / "no_soft.bam",
            reads=[{"cigar": "50M", "seq": "A" * 50, "flag": 2048, "tags": {"cs": ":50"}}],
        )
        scanner_no_soft = BamScanner(no_soft, 20, ref, gtf, 5, fake_logger, True, 100, 200, 10, 0.05, 0.01, None)
        with pytest.raises(SystemExit):
            scanner_no_soft.iter_bam()

        # Supplementary alignment with soft-clipping but missing cs tag -> SystemExit (--cs required)
        no_cs = _write_bam(
            tmp_path / "no_cs.bam",
            reads=[{"cigar": "10S40M", "seq": "A" * 50, "flag": 2048}],
        )
        scanner_no_cs = BamScanner(no_cs, 20, ref, gtf, 5, fake_logger, True, 100, 200, 10, 0.05, 0.01, None)
        with pytest.raises(SystemExit):
            scanner_no_cs.iter_bam()

        # Valid supplementary alignment populates representative_alignments_new_cigar
        valid_supp = _write_bam(
            tmp_path / "valid_supp.bam",
            reads=[{"qname": "q_ok", "cigar": "10S40M", "seq": "A" * 50, "flag": 2048, "tags": {"cs": ":40"}}],
        )
        scanner_ok = BamScanner(valid_supp, 20, ref, gtf, 5, fake_logger, True, 100, 200, 10, 0.05, 0.01, None)
        scanner_ok.iter_bam()
        assert scanner_ok.representative_alignments_new_cigar == {"q_ok\t10\t0": "10S40M"}

        # High-substitution supplementary alignment triggers BLAT realignment in BamScanner.iter_bam()
        noisy_supp = _write_bam(
            tmp_path / "noisy_supp.bam",
            reads=[
                {
                    "qname": "q_noisy",
                    "cigar": "10S40M",
                    "seq": "A" * 50,
                    "flag": 2048,
                    "tags": {"cs": ":10*ac*ag*at*ca*cg:25"},
                },
            ],
        )

        class _ScannerAligner:
            min_aligned_len = 20

            def query(self, in_seq: str) -> str:
                return "fake.pslx"

            def _query_insertion(self, *_a: object, **_kw: object) -> tuple[int, object, int]:
                return 1, "top_hsp", 60

            def psl2sam(self, _hsp: object, in_seq_len: int) -> tuple[str, int, str, str, int]:
                return "chr1", 200, "-", f"{in_seq_len}M", 0

            def obtain_variants_stats(self, _hsp: object, in_seq_len: int) -> tuple[float, float, float]:
                return 0.0, 0.0, 0.0

        from Bio import SearchIO

        orig_read = SearchIO.read
        SearchIO.read = lambda *_a, **_kw: "parsed_psl"  # type: ignore[assignment]
        try:
            scanner_realign = BamScanner(
                noisy_supp,
                20,
                ref,
                gtf,
                5,
                fake_logger,
                True,
                100,
                200,
                10,
                0.01,
                0.01,
                _ScannerAligner(),
            )
            scanner_realign.iter_bam()
        finally:
            SearchIO.read = orig_read
        assert scanner_realign.representative_alignments_new_record == {"q_noisy\t10\t0": "chr1,201,-,40M10S,60,0"}

    def test_main_helpers_and_update_position_event_list(self, tmp_path: Path) -> None:
        assert _get_read_matched_sequence("AAACCCTTT", 3, 3, "+") == "CCC"
        assert _get_read_matched_sequence("AAACCCTTT", 3, 3, "-") == "GGG"

        with pytest.raises(SystemExit, match="Not Found"):
            _get_genome_fasta(tmp_path / "missing_ref.fa")
        with pytest.raises(SystemExit, match="Reading GTF"):
            _get_cvg_gene_iv(tmp_path / "missing.gtf", 5)

        # Single event is returned unchanged
        ev1 = Event(
            (
                "TRA",
                0,
                1,
                ("chr1:200", "chr2:300", MappingMode.MS, MappingMode.SM),
                (100, 200, Intervals.from_list([[100, 200]])),
                (300, 380, Intervals.from_list([[300, 380]])),
                ("", ""),
                ("+", "+"),
                ["G1", "G2"],
                False,
            ),
        )
        assert update_position_event_list([ev1]) == [ev1]

        # Two-hop and three-hop co-linear middle segment updates (both forward and reversed event flags)
        ev2 = Event(
            (
                "TRA",
                0,
                1,
                ("chr2:420", "chr3:500", MappingMode.MS, MappingMode.SM),
                (310, 420, Intervals.from_list([[310, 420]])),
                (500, 600, Intervals.from_list([[500, 600]])),
                ("", ""),
                ("+", "+"),
                ["G2", "G3"],
                False,
            ),
        )
        updated = update_position_event_list([ev1, ev2])
        assert len(updated) == 2
        assert (updated[0].read2_ref_start, updated[0].read2_ref_end) == (300, 420)
        assert (updated[1].read1_ref_start, updated[1].read1_ref_end) == (300, 420)

        ev_rev1 = Event(
            (
                "TRA",
                0,
                1,
                ("chr2:300", "chr1:200", MappingMode.SM, MappingMode.MS),
                (300, 380, Intervals.from_list([[300, 380]])),
                (100, 200, Intervals.from_list([[100, 200]])),
                ("", ""),
                ("+", "+"),
                ["G2", "G1"],
                True,
            ),
        )
        ev_rev2 = Event(
            (
                "TRA",
                0,
                1,
                ("chr3:500", "chr2:420", MappingMode.SM, MappingMode.MS),
                (500, 600, Intervals.from_list([[500, 600]])),
                (310, 420, Intervals.from_list([[310, 420]])),
                ("", ""),
                ("+", "+"),
                ["G3", "G2"],
                True,
            ),
        )
        updated_rev = update_position_event_list([ev_rev1, ev_rev2])
        assert (updated_rev[0].read1_ref_start, updated_rev[0].read1_ref_end) == (300, 420)
        assert (updated_rev[1].read2_ref_start, updated_rev[1].read2_ref_end) == (300, 420)

    def test_cli_exits_on_missing_output_dir_or_zero_paths(self, tmp_path: Path) -> None:
        ref = tmp_path / "ref.fa"
        ref.write_text(">chr1\n" + ("A" * 1000) + "\n")
        gtf = tmp_path / "g.gtf"
        gtf.write_text('chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
        two_bit = tmp_path / "dummy.2bit"
        two_bit.write_bytes(b"")
        empty_bam = _write_bam(tmp_path / "empty.bam", reads=[{"cigar": "50M", "seq": "A" * 50, "tags": {"NM": 0, "cs": ":50"}}])

        # Missing parent directory
        op_bad_dir = DefaultOptions(
            input=str(empty_bam),
            output=str(tmp_path / "nonexistent_dir" / "prefix"),
            ref=str(ref),
            gtf=str(gtf),
            blat_two_bit=str(two_bit),
            aligner=("", ""),
            log="info",
        )
        with pytest.raises(SystemExit, match="does not exist"):
            cli.cli(op_bad_dir)

        # Valid directory but zero chimeric paths -> SystemExit
        op_empty = DefaultOptions(
            input=str(empty_bam),
            output=str(tmp_path / "prefix"),
            ref=str(ref),
            gtf=str(gtf),
            blat_two_bit=str(two_bit),
            aligner=("", ""),
            log="info",
        )
        with pytest.raises(SystemExit):
            cli.cli(op_empty)
