# !/usr/bin/env python
"""Test for scanmst/base/helper.py."""
import contextlib
import os
from dataclasses import dataclass

import HTSeq  # type: ignore
import pytest
from pyfaidx import Fasta  # type: ignore
from scanmst.base.helper import (
    diff_chrom_diff_strand_handler,
    diff_chrom_same_strand_mode21_handler,
    extract_splice_sites,
    gene_annotation,
    same_chrom_diff_strand_handler,
    same_chrom_same_strand_mode21_handler,
)
from scanmst.base import MappingMode
from scanmst.utils import cigar_validity

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


def assert_event(
    result,
    *,
    sv_type,
    anno,
    can,
    bp1,
    bp2,
    mode1,
    mode2,
    read1,
    read2,
    insertions,
    strands,
    genes,
):
    """Compare a handler's event tuple field by field.

    Handlers return a 10-element tuple whose members are typed: the modes are
    MappingMode, the strands are Strand and the exons are Intervals. Comparing
    against a flat literal tuple silently fails, so each field is checked
    explicitly here and exons are normalised to (start, end) pairs.
    """
    assert result is not None, "handler returned noreturn"
    (
        got_sv_type,
        got_anno,
        got_can,
        positions,
        read1_info,
        read2_info,
        insertion_info,
        got_strands,
        got_genes,
        is_read_reversed,
    ) = result

    assert got_sv_type == sv_type
    assert (got_anno, got_can) == (anno, can)
    assert positions == (bp1, bp2, mode1, mode2)

    for got_info, (exp_start, exp_end, exp_exons) in zip((read1_info, read2_info), (read1, read2)):
        got_start, got_end, got_exons = got_info
        assert (got_start, got_end) == (exp_start, exp_end)
        assert [(exon.start, exon.end) for exon in got_exons] == exp_exons

    assert insertion_info == insertions
    assert tuple(str(strand) for strand in got_strands) == strands
    assert got_genes == genes
    assert isinstance(is_read_reversed, bool)


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


# splicing_confirmation_and_correction (formerly splicing_confirmation) is no
# longer callable from breakpoints alone: it now takes the full read context
# (strand, mode, ref span and exons for both sides, 18 arguments) and returns an
# 11-tuple of corrected coordinates rather than a verdict triple. The four
# handler tests below call it through their handlers, so it stays covered.


@pytest.mark.parametrize(
    ("chrm1", "pos1", "chrm2", "pos2", "expected_result"),
    [
        ("chr20", 391287, "chr20", 410025, ("TRIB3", "RBCK1")),
        ("chr17", 172536, "chr17", 247286, ("DOC2B", "RPH3AL")),
    ],
)
def test_gene_annotation(gtf_setup, chrm1, pos1, chrm2, pos2, expected_result):
    """Test gene_annotation func."""
    cvg, gene_iv = gtf_setup
    assert gene_annotation(chrm1, pos1, chrm2, pos2, gene_iv) == expected_result


def test_same_chrom_same_strand_mode21_handler(
    prepare_fasta_and_gtf,
    tdup_reads,
    fake_logger,
):
    """Test same_chrom_same_strand_mode21_handler func."""
    genome_fasta, cvg, gene_iv = prepare_fasta_and_gtf
    read_lt, read_rt = tdup_reads
    lt_mode, rt_mode = MappingMode.SM, MappingMode.MS
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
        fake_logger,
    )
    assert_event(
        result,
        sv_type="TDUP",
        anno=3,
        can=1,
        bp1="chr17:1364856",
        bp2="chr17:1423648",
        mode1=MappingMode.SM,
        mode2=MappingMode.MS,
        read1=(1364856, 1400222, [(1364856, 1365058), (1400046, 1400222)]),
        read2=(1423349, 1423651, [(1423349, 1423651)]),
        insertions=("-ACC", "-ACC"),
        strands=("-", "-"),
        genes=["YWHAE", "CRK"],
    )


def test_diff_chrom_same_strand_mode21_handler(
    prepare_fasta_and_gtf,
    trans_same_strand_reads,
    fake_logger,
):
    """Test diff_chrom_same_strand_mode21_handler func (TRA)."""
    genome_fasta, cvg, gene_iv = prepare_fasta_and_gtf
    read_lt, read_rt = trans_same_strand_reads
    lt_mode, rt_mode = MappingMode.MS, MappingMode.SM
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
        fake_logger,
    )
    assert_event(
        result,
        sv_type="TRA",
        anno=3,
        can=1,
        bp1="chr20:391283",
        bp2="chr17:1745408",
        mode1=MappingMode.SM,
        mode2=MappingMode.MS,
        read1=(391283, 397559, [(391283, 391579), (396197, 397559)]),
        read2=(
            1744991,
            1745408,
            [(1744991, 1745058), (1745174, 1745213), (1745332, 1745408)],
        ),
        insertions=("-CAGGTG", "-CAGGTG"),
        strands=("+", "+"),
        genes=["TRIB3", "SERPINF2"],
    )


def test_diff_chrom_diff_strand_handler(
    prepare_fasta_and_gtf,
    trans_diff_strand_reads,
    fake_logger,
):
    """Test diff_chrom_diff_strand_handler func (TRA)."""
    genome_fasta, cvg, gene_iv = prepare_fasta_and_gtf
    read_lt, read_rt = trans_diff_strand_reads
    lt_mode, rt_mode = MappingMode.MS, MappingMode.MS
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
        fake_logger,
    )
    assert_event(
        result,
        sv_type="TRA",
        anno=3,
        can=1,
        bp1="chr17:1745405",
        bp2="chr20:439298",
        mode1=MappingMode.MS,
        mode2=MappingMode.MS,
        read1=(
            1744991,
            1745405,
            [(1744991, 1745058), (1745174, 1745213), (1745332, 1745405)],
        ),
        read2=(435479, 439298, [(435479, 438841), (439107, 439298)]),
        insertions=("-CAG", "-CTG"),
        strands=("+", "-"),
        genes=["SERPINF2", "TBC1D20"],
    )


def test_same_chrom_diff_strand_handler(prepare_fasta_and_gtf, inv_reads, fake_logger):
    """Test same_chrom_diff_strand_handler func (INV)."""
    genome_fasta, cvg, gene_iv = prepare_fasta_and_gtf
    read_lt, read_rt = inv_reads
    lt_mode, rt_mode = MappingMode.MS, MappingMode.MS
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
        fake_logger,
    )
    assert_event(
        result,
        sv_type="INV",
        anno=3,
        can=1,
        bp1="chr17:1651554",
        bp2="chr17:1730946",
        mode1=MappingMode.MS,
        mode2=MappingMode.MS,
        read1=(
            1650628,
            1651554,
            [(1650628, 1650956), (1651107, 1651310), (1651413, 1651554)],
        ),
        read2=(
            1727989,
            1730946,
            [(1727989, 1728626), (1730379, 1730487), (1730754, 1730946)],
        ),
        insertions=("-CC", "-GG"),
        strands=("-", "+"),
        genes=["PRPF8", "WDR81"],
    )
