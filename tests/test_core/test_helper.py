# !/usr/bin/env python
"""Test for scannls/core/helper.py."""
import contextlib
import os
from dataclasses import dataclass

import HTSeq  # type: ignore
import pytest
from pyfaidx import Fasta  # type: ignore

from scannls.core.helper import cigar_validity
from scannls.core.helper import extract_splice_sites
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
def one_postion(request):
    """Use id_func to generate ids."""
    return request.param


def test_extract_splice_sites(gtf_setup, one_postion):
    """Test for extract_splice_sites func."""
    cvg, gene_iv = gtf_setup
    _chrm, _pos, _expect = one_postion
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
