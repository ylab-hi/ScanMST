# !/usr/bin/env python
"""Test for scannls/core/helper.py."""
import contextlib
import os
from dataclasses import dataclass

import HTSeq  # type: ignore
import pytest
from loguru import logger
from pyfaidx import Fasta  # type: ignore

from scannls.core.helper import cigar_validity
from scannls.core.helper import diff_chrom_diff_strand_handler
from scannls.core.helper import diff_chrom_same_strand_mode21_handler
from scannls.core.helper import extract_splice_sites
from scannls.core.helper import same_chrom_diff_strand_handler
from scannls.core.helper import same_chrom_same_strand_mode21_handler
from scannls.core.helper import splicing_confirmation

path = os.path.dirname(__file__)
os.chdir(path)


@pytest.fixture()
def fasta_setup():
    """Fasta setup."""
    yield
    with contextlib.suppress(OSError):
        os.remove("../data/dummy.fasta.fai")


@pytest.fixture()
def gtf_setup():
    """GTF setup."""
    _cvg, _gene_iv = extract_splice_sites("../data/dummy.gtf", 5)
    return _cvg, _gene_iv


@pytest.fixture()
def cigar_data():
    """Return a cigarstring."""
    return "1S2S5M3S2S"


def test_cigar_validity(cigar_data):
    """Test cigar_validity func."""
    assert cigar_validity(cigar_data) == "3S5M5S"


positions = [
    ("chr17", 247110, "AC"),
    ("chr20", 391580, "GT"),
    ("chr20", 391287, "AG"),
    ("chr20", 388301, "GT"),
]


def id_func(fixture_value):
    """A func for generating ids."""
    t = fixture_value
    return f"{t[0]}:{t[1]}"


@pytest.fixture(params=positions, ids=id_func)
def one_position(request):
    """Use id_func to generate ids."""
    return request.param


def test_extract_splice_sites(gtf_setup, one_position):
    """Test for extract_splice_sites func."""
    cvg, gene_iv = gtf_setup
    _chrm, _pos, _expect = one_position
    _result = list(cvg[HTSeq.GenomicPosition(_chrm, _pos)])[0]
    assert _result == _expect


@pytest.fixture()
def prepare_fasta_and_gtf(fasta_setup, gtf_setup):
    """Prepare Fasta and GTF."""
    _fasta = Fasta("../data/dummy.fasta", sequence_always_upper=True)
    _cvg, _gene_iv = gtf_setup
    return _fasta, _cvg, _gene_iv


@dataclass
class Breakpoints:
    """Class to store the information of one pair of breakpoints."""

    chrm1: str
    pos1: int
    chrm2: str
    pos2: int


breakpoints = [
    (Breakpoints("chr20", 391287, "chr20", 391580), (True, 3, 1)),
    (Breakpoints("chr17", 172536, "chr17", 247286), (True, 3, 1)),
]


def bp_id_func(fixture_value):
    """A func for generating ids."""
    breakpoints = fixture_value[0]
    return (
        f"{breakpoints.chrm1}:{breakpoints.pos1}-{breakpoints.chrm2}:{breakpoints.pos2}"
    )


@pytest.fixture(params=breakpoints, ids=bp_id_func)
def one_pair_pbs(request):
    """Use bp_id_func to generate ids."""
    return request.param


def test_splicing_confirmation(prepare_fasta_and_gtf, one_pair_pbs):
    """Test splicing_confirmation func."""
    _fasta, _cvg, _gene_iv = prepare_fasta_and_gtf
    bps, _expect = one_pair_pbs
    _result = splicing_confirmation(
        bps.chrm1, bps.pos1, bps.chrm2, bps.pos2, 5, _fasta, _cvg, False, True
    )

    assert _result == _expect


def test_same_chrom_same_strand_mode21_handler(prepare_fasta_and_gtf, tdup_reads):
    """Test same_chrom_same_strand_mode21_handler func."""
    genome_fasta, cvg, gene_iv = prepare_fasta_and_gtf
    read_lt, read_rt = tdup_reads
    lt_mode, rt_mode = 2, 1
    splice_bin = 5
    motif_required = True
    result = same_chrom_same_strand_mode21_handler(
        read_lt,
        read_rt,
        lt_mode,
        rt_mode,
        splice_bin,
        genome_fasta,
        cvg,
        gene_iv,
        motif_required,
        logger,
    )
    expect = (
        "TDUP",
        3,
        1,
        ("chr17:1364856", "chr17:1423648", 2, 1),
        (1364856, 1400222, [[1364856, 1365058], [1400046, 1400222]]),
        (1423349, 1423651, [[1423349, 1423651]]),
        ("-ACC", "-ACC"),
        ("-", "-"),
        ["YWHAE", "CRK"],
    )
    assert result == expect


def test_diff_chrom_same_strand_mode21_handler(
    prepare_fasta_and_gtf, trans_same_strand_reads
):
    """Test diff_chrom_same_strand_mode21_handler func (TRA)."""
    genome_fasta, cvg, gene_iv = prepare_fasta_and_gtf
    read_lt, read_rt = trans_same_strand_reads
    lt_mode, rt_mode = 1, 2
    splice_bin = 5
    motif_required = True
    result = diff_chrom_same_strand_mode21_handler(
        read_rt,
        read_lt,
        rt_mode,
        lt_mode,
        splice_bin,
        genome_fasta,
        cvg,
        gene_iv,
        motif_required,
        logger,
    )
    expect = (
        "TRA",
        3,
        1,
        ("chr20:391283", "chr17:1745408", 2, 1),
        (391283, 397559, [[391283, 391579], [396197, 397559]]),
        (
            1744991,
            1745408,
            [[1744991, 1745058], [1745174, 1745213], [1745332, 1745408]],
        ),
        ("-CAGGTG", "-CAGGTG"),
        ("+", "+"),
        ["TRIB3", "SERPINF2"],
    )
    assert result == expect


def test_diff_chrom_diff_strand_handler(prepare_fasta_and_gtf, trans_diff_strand_reads):
    """Test diff_chrom_diff_strand_handler func (TRA)."""
    genome_fasta, cvg, gene_iv = prepare_fasta_and_gtf
    read_lt, read_rt = trans_diff_strand_reads
    lt_mode, rt_mode = 1, 1
    splice_bin = 5
    motif_required = True
    result = diff_chrom_diff_strand_handler(
        read_lt,
        read_rt,
        lt_mode,
        rt_mode,
        splice_bin,
        genome_fasta,
        cvg,
        gene_iv,
        motif_required,
        logger,
    )
    expect = (
        "TRA",
        3,
        1,
        ("chr17:1745405", "chr20:439298", 1, 1),
        (
            1744991,
            1745405,
            [[1744991, 1745058], [1745174, 1745213], [1745332, 1745405]],
        ),
        (435479, 439298, [[435479, 438841], [439107, 439298]]),
        ("-CAG", "-CTG"),
        ("+", "-"),
        ["SERPINF2", "TBC1D20"],
    )
    assert result == expect


def test_same_chrom_diff_strand_handler(prepare_fasta_and_gtf, inv_reads):
    """Test same_chrom_diff_strand_handler func (INV)."""
    genome_fasta, cvg, gene_iv = prepare_fasta_and_gtf
    read_lt, read_rt = inv_reads
    lt_mode, rt_mode = 1, 1
    splice_bin = 5
    motif_required = True
    result = same_chrom_diff_strand_handler(
        read_lt,
        read_rt,
        lt_mode,
        rt_mode,
        splice_bin,
        genome_fasta,
        cvg,
        gene_iv,
        motif_required,
        logger,
    )
    expect = (
        "INV",
        3,
        1,
        ("chr17:1651554", "chr17:1730946", 1, 1),
        (
            1650628,
            1651554,
            [[1650628, 1650956], [1651107, 1651310], [1651413, 1651554]],
        ),
        (
            1727989,
            1730946,
            [[1727989, 1728626], [1730379, 1730487], [1730754, 1730946]],
        ),
        ("-CC", "-GG"),
        ("-", "+"),
        ["PRPF8", "WDR81"],
    )
    assert result == expect
