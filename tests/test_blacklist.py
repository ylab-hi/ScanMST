"""Tests for the junction blacklist (``--blacklist``)."""

import gzip
import inspect
from types import SimpleNamespace

import pytest

from scanmst.base import BlacklistFilter
from scanmst.cli.arg import DefaultOptions, parse_args
from scanmst.cli.main import _scan_bam_helper, scanbam_run

BED = """\
track name=blacklist
browser position chr1:1-100
# comment

chr1\t100\t200\tHigh Signal Region
chr1\t500\t600
chr2 1000 2000
chr3\t50\t50
"""


@pytest.fixture
def bed(tmp_path):
    path = tmp_path / "blacklist.bed"
    path.write_text(BED)
    return path


def junction(chrom1, pos1, chrom2, pos2):
    return SimpleNamespace(chrom1=chrom1, junction_pos1=pos1, chrom2=chrom2, junction_pos2=pos2)


def test_parses_bed_skipping_headers_and_empty_intervals(bed):
    bl = BlacklistFilter(bed)
    assert bl.num_regions == 3
    assert bl.chroms == {"chr1", "chr2"}


@pytest.mark.parametrize("name", ["blacklist.bed.gz", "blacklist.bed"])
def test_parses_gzipped_bed(tmp_path, name):
    # gzip is detected from the file content, not the suffix
    path = tmp_path / name
    with gzip.open(path, "wt") as fh:
        fh.write(BED)
    bl = BlacklistFilter(path)
    assert bl.num_regions == 3
    assert bl.contains("chr2", 1500)


@pytest.mark.parametrize(
    ("chrom", "pos", "expected"),
    [
        ("chr1", 99, False),
        ("chr1", 100, True),  # BED start is inclusive
        ("chr1", 199, True),
        ("chr1", 200, False),  # BED end is exclusive
        ("chr1", 550, True),
        ("chr3", 50, False),  # zero-length interval covers nothing
        ("1", 150, False),  # chromosome names are matched verbatim
        ("chrX", 150, False),
    ],
)
def test_contains_boundaries(bed, chrom, pos, expected):
    assert BlacklistFilter(bed).contains(chrom, pos) is expected


@pytest.mark.parametrize(
    ("event", "expected"),
    [
        (junction("chr1", 150, "chr1", 5000), True),  # breakpoint 1 only
        (junction("chr1", 5000, "chr1", 150), True),  # breakpoint 2 only
        (junction("chr1", 150, "chr2", 1500), True),  # both, inter-chromosomal
        (junction("chr1", 300, "chr2", 2500), False),  # neither
        (junction("chr1", 300, "chr1", 400), False),  # span covers nothing; only breakpoints count
    ],
)
def test_is_in_blacklist(bed, event, expected):
    assert BlacklistFilter(bed).is_in_blacklist(event) is expected


@pytest.mark.parametrize("line", ["chr1\t100", "chr1\tx\t200", "chr1\t300\t200"])
def test_malformed_line_raises(tmp_path, line):
    path = tmp_path / "bad.bed"
    path.write_text(f"chr1\t1\t2\n{line}\n")
    with pytest.raises(ValueError, match=r"bad\.bed:2"):
        BlacklistFilter(path)


REQUIRED = ["--input", "a.bam", "--ref", "a.fa", "--gtf", "a.gtf", "--output", "out"]


def test_blacklist_is_optional():
    assert DefaultOptions.blacklist is None
    assert parse_args().parse_args(REQUIRED).blacklist is None


def test_blacklist_option_is_parsed():
    args = parse_args().parse_args([*REQUIRED, "--blacklist", "hg38-blacklist.v2.bed"])
    assert args.blacklist == "hg38-blacklist.v2.bed"


def test_blacklist_is_forwarded_to_scan_helper():
    # scanbam_run forwards its locals to _scan_bam_helper by keyword-only parameter name
    helper_param = inspect.signature(_scan_bam_helper).parameters["blacklist"]
    assert helper_param.kind is inspect.Parameter.KEYWORD_ONLY
    assert inspect.signature(scanbam_run).parameters["blacklist"].default is None
