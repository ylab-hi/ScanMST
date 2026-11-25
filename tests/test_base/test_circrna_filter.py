"""Tests for CircRNAFilter condition 4 (two TDUP mega-exons within an annotated intron).

Regression: commit 67cf6eb added a "mono-exon" rule that flagged ANY same-strand TDUP whose two
segments are unspliced, regardless of annotation. On simulated one-hop TDUPs (two different genes,
unspliced segments) it removed ~10k read paths and cut recall from 0.85 to 0.52.
"""

from types import SimpleNamespace

import pytest
from scanmst.base.basic import Strand
from scanmst.base.filters import CircRNAFilter

GTF = (
    'chr1\tt\texon\t1001\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tt\texon\t5001\t6000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tt\texon\t20001\t21000\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
    'chr1\tt\texon\t24001\t25000\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
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
