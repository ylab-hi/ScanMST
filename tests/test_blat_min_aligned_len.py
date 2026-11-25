"""Tests for the minimum aligned length of BLAT hits (``--min-blat-aligned-len``).

BLAT hits shorter than the cutoff are treated as unaligned, so they never become new segments.
These tests use fake PSL objects and do not need a gfServer.
"""

from types import SimpleNamespace

import pytest
from scanmst.base.blat import Blat
from scanmst.cli.arg import DefaultOptions, parse_args

from tests import FakeBlat


def hsp(aligned_len: int, score: int | None = None, hit_id: str = "chr1"):
    """A gap-free, mismatch-free HSP aligning ``aligned_len`` query bases."""
    return SimpleNamespace(
        score=aligned_len if score is None else score,
        hit_id=hit_id,
        hit_span_all=[aligned_len],
        query_span_all=[aligned_len],
        mismatch_num=0,
        hit_gap_num=0,
    )


def query(blat_result_hsps, seq_len, min_aligned_len, threshold_identity=0.9):
    return Blat._query_insertion(
        SimpleNamespace(hsps=list(blat_result_hsps)),
        "A" * seq_len,
        threshold_identity,
        top=5,
        min_aligned_len=min_aligned_len,
    )


def test_short_hit_is_unaligned_with_cutoff():
    hit, top_hsp, mapq = query([hsp(25)], seq_len=25, min_aligned_len=30)
    assert (hit, top_hsp, mapq) == (0, "", 0)


def test_short_hit_is_kept_without_cutoff():
    """min_aligned_len=0 reproduces the previous behaviour."""
    hit, _, mapq = query([hsp(25)], seq_len=25, min_aligned_len=0)
    assert hit == 1
    assert mapq == 60


def test_cutoff_is_inclusive():
    hit, _, _ = query([hsp(30)], seq_len=30, min_aligned_len=30)
    assert hit == 1


def test_only_long_enough_hits_are_kept_and_scored():
    long_hit = hsp(40, score=40)
    short_hit = hsp(25, score=45)  # better score, but too short (and below identity for a 40 bp query)
    hit, top_hsp, mapq = query([short_hit, long_hit], seq_len=40, min_aligned_len=30)
    assert hit == 1
    assert top_hsp is long_hit
    assert mapq == 60  # computed from the kept hits only


def test_short_high_identity_hits_do_not_affect_mapq():
    # two equally good 35 bp hits -> ambiguous (mapq 0); an extra 20 bp hit is ignored either way
    hit, _, mapq = query([hsp(35), hsp(35), hsp(20)], seq_len=35, min_aligned_len=30)
    assert hit == 2
    assert mapq == 0


def test_blat_default_and_test_double_have_the_attribute():
    assert Blat.min_aligned_len == 30
    assert FakeBlat().min_aligned_len == 30


@pytest.mark.parametrize(("argv", "expected"), [([], 30), (["--min-blat-aligned-len", "20"], 20)])
def test_cli_option(argv, expected):
    required = ["--input", "in.bam", "--ref", "ref.fa", "--gtf", "a.gtf", "--output", "out"]
    options = parse_args().parse_args(required + argv)
    assert options.min_blat_aligned_len == expected
    assert DefaultOptions.min_blat_aligned_len == 30
