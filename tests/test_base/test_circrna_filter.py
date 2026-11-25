"""Tests for CircRNAFilter condition 4 (two TDUP mega-exons within an annotated intron).

Regression: commit 67cf6eb added a "mono-exon" rule that flagged ANY same-strand TDUP whose two
segments are unspliced, regardless of annotation. On simulated one-hop TDUPs (two different genes,
unspliced segments) it removed ~10k read paths and cut recall from 0.85 to 0.52.
"""

from types import SimpleNamespace

import pytest
import HTSeq
from scanmst.base.basic import Strand
from scanmst.base.filters import CircRNAFilter, _extract_annotated_exons

# T1 (+): exons [1000,2000) [5000,6000) [8000,9000) -> introns [2000,5000) [6000,8000)
# T2 (-): exons [20000,21000) [24000,25000)          -> intron  [21000,24000)
# T3 (+): single exon                                -> no intron
GTF = (
    'chr1\tt\texon\t1001\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tt\texon\t5001\t6000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tt\texon\t8001\t9000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tt\texon\t20001\t21000\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
    'chr1\tt\texon\t24001\t25000\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
    'chr1\tt\texon\t50001\t51000\t.\t+\t.\tgene_id "G3"; transcript_id "T3";\n'
)


def node(chrom, strand, start, end, exons=None):
    """Minimal stand-in exposing the attributes the method reads."""
    return SimpleNamespace(
        chrom=chrom,
        strand=Strand.from_str(strand),
        ref_start=start,
        ref_end=end,
        exons=exons if exons is not None else [(start, end)],
    )


@pytest.fixture(scope="module")
def circ_filter(tmp_path_factory):
    gtf = tmp_path_factory.mktemp("circ") / "mini.gtf"
    gtf.write_text(GTF)
    return CircRNAFilter(str(gtf), 10)


@pytest.mark.parametrize(
    ("first", "second"),
    [
        # + strand back-jump: 5' segment downstream of the 3' segment, both unspliced (e.g. event 166)
        (node("chr1", "+", 400_000, 403_975), node("chr1", "+", 100_000, 100_527)),
        # - strand back-jump
        (node("chr1", "-", 100_000, 103_975), node("chr1", "-", 400_000, 400_527)),
    ],
)
def test_unspliced_tdup_outside_annotated_intron_is_not_circrna(circ_filter, first, second):
    assert circ_filter.is_two_megaexon_within_annotated_intron(first, second) is False


@pytest.mark.parametrize(
    ("first", "second"),
    [
        (node("chr1", "+", 400_000, 403_975), node("chr1", "-", 100_000, 100_527)),  # strands differ
        (node("chr1", "+", 400_000, 403_975), node("chr2", "+", 100_000, 100_527)),  # chromosomes differ
    ],
)
def test_gate_rejects_mismatched_segments(circ_filter, first, second):
    assert circ_filter.is_two_megaexon_within_annotated_intron(first, second) is False


@pytest.fixture(scope="module")
def mini_gtf(tmp_path_factory):
    gtf = tmp_path_factory.mktemp("introns") / "mini.gtf"
    gtf.write_text(GTF)
    return str(gtf)


def _introns(introns_gas):
    return {(i.trx_id, i.start, i.end, str(i.strand)) for _, s in introns_gas.steps() for i in s}


def test_extract_annotated_introns_from_consecutive_exons(mini_gtf):
    _, introns_gas = _extract_annotated_exons(mini_gtf, 10, shrink=False, consider_strand=True)
    assert _introns(introns_gas) == {
        ("T1", 2000, 5000, "+"),
        ("T1", 6000, 8000, "+"),
        ("T2", 21000, 24000, "-"),
    }  # single-exon T3 contributes nothing


def test_introns_are_padded_by_boundary_size(mini_gtf):
    _, introns_gas = _extract_annotated_exons(mini_gtf, 10, shrink=False, consider_strand=True)
    assert introns_gas[HTSeq.GenomicPosition("chr1", 1990, ".")]  # 2000 - 10
    assert not introns_gas[HTSeq.GenomicPosition("chr1", 1989, ".")]
    assert introns_gas[HTSeq.GenomicPosition("chr1", 5009, ".")]  # 5000 + 10, exclusive
    assert not introns_gas[HTSeq.GenomicPosition("chr1", 5010, ".")]


@pytest.mark.parametrize(
    ("first", "second"),
    [
        # + strand back-jump, both segments inside T1's intron [2000,5000)
        (node("chr1", "+", 3500, 3800), node("chr1", "+", 2500, 2800)),
        # - strand back-jump, both segments inside T2's intron [21000,24000)
        (node("chr1", "-", 21500, 21800), node("chr1", "-", 22500, 22800)),
        # strand-agnostic by design: a - strand pair inside the + intron of T1 still counts
        (node("chr1", "-", 2500, 2800), node("chr1", "-", 3500, 3800)),
    ],
)
def test_tdup_within_annotated_intron_is_circrna(circ_filter, first, second):
    assert circ_filter.is_two_megaexon_within_annotated_intron(first, second) is True


def test_tdup_spanning_two_introns_is_not_circrna(circ_filter):
    # segments in different introns of T1 ([2000,5000) and [6000,8000)): no single shared intron
    first, second = node("chr1", "+", 6500, 6800), node("chr1", "+", 2500, 2800)
    assert circ_filter.is_two_megaexon_within_annotated_intron(first, second) is False
