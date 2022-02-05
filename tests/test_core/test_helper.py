# !/usr/bin/env python
"""Test for scannls/core/helper.py."""
import contextlib
import os
from dataclasses import dataclass

import HTSeq  # type: ignore
import pytest
from loguru import logger
from pyfaidx import Fasta  # type: ignore

from scannls import Read
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


@pytest.fixture()
def tdup_reads():
    """Reads fixture."""
    read_init_params = [
        "query_name",
        "chrom",
        "ref_start",
        "strand",
        "cigar_str",
        "mapq",
        "nm",
        "query_seq",
    ]

    param_dict = [
        {
            "chrom": "chr17",
            "ref_start": 1364856,
            "strand": "-",
            "cigar_str": "299S202M34988N176M",
            "mapq": 60,
            "nm": 0,
            "query_seq": "TAAGACTAGGAATGGAGCAGTTCAGTCTAAAAAATATCACAGGTAACAGA"
            "ATGCTTATATAAACTAGACTGCTTTTGACATCTGTAAGAAAATTGTATAG"
            "ATGGCAGTTGGAAAAAAAAAAAAAAGATTGTTCCCATCTGTCAGCAAAAC"
            "TGTTGAACTATACTCAGCTGAAGTCCTCATCGGGATTCTGTTGATCCAGC"
            "AGACGGACATGTGTGAATGGGAAGTGACCTCGTTTGCCATTACACTCCCC"
            "TTCCCACTGACCACTCACATTAATCTTCGTAACCTTTACCAGCTCACCGA"
            "CCATTTGCCGATATTCCCGAATCATTTTTAGCTTGTCTTCTCCTCCCTTG"
            "TTTTCTTCTTTCTGTTCAATGCTGCTGATTATTCTCCAGGAGGCTCTTCT"
            "AGCTCCAATCACATTCTTATATGCAACAGATAGGAGGTTTCTTTCTTCAA"
            "CTGTCAGCTCCACATCCATCCCTGCTACTTTCTTCATTGACTCCACCATT"
            "TCGTCGTATCGCTCAGCCTGCTCGGCCAGCTTCGCCTGGTACACCAGATC"
            "CTCTCGATCATCCATAGCGGCAGCGGCTCCGGCAGGGTCTGCGCGACGGA"
            "TGGAAGCGGATAGTGTCTCCGACTCTCTCAGCCTCTCGCTCCGCGTCCGG"
            "GCAGCAAAAATGGCGGCGCCTCAATCC",
            "query_name": "one",
            "lt_soft_len": 299,
            "rt_soft_len": 0,
            "read_match_size": 378,
            "reference_match_size": 35366,
            "indel_size": 34988,
            "cigartuples_without_soft": [[0, 202], [3, 34988], [0, 176]],
            "cigartuples": [(4, 299), (0, 202), (3, 34988), (0, 176)],
            "query_length": 677,
            "adhocsms": None,
            "adhocseq": None,
            "mode": None,
            "sms": (299, 378, 0),
            "ref_end": 1400222,
        },
        {
            "chrom": "chr17",
            "ref_start": 1423349,
            "strand": "-",
            "cigar_str": "302M375S",
            "mapq": 60,
            "nm": 0,
            "query_seq": "TAAGACTAGGAATGGAGCAGTTCAGTCTAAAAAATATCACAGGTAACAGA"
            "ATGCTTATATAAACTAGACTGCTTTTGACATCTGTAAGAAAATTGTATAG"
            "ATGGCAGTTGGAAAAAAAAAAAAAAGATTGTTCCCATCTGTCAGCAAAAC"
            "TGTTGAACTATACTCAGCTGAAGTCCTCATCGGGATTCTGTTGATCCAGC"
            "AGACGGACATGTGTGAATGGGAAGTGACCTCGTTTGCCATTACACTCCCC"
            "TTCCCACTGACCACTCACATTAATCTTCGTAACCTTTACCAGCTCACCGA"
            "CCATTTGCCGATATTCCCGAATCATTTTTAGCTTGTCTTCTCCTCCCTTG"
            "TTTTCTTCTTTCTGTTCAATGCTGCTGATTATTCTCCAGGAGGCTCTTCT"
            "AGCTCCAATCACATTCTTATATGCAACAGATAGGAGGTTTCTTTCTTCAA"
            "CTGTCAGCTCCACATCCATCCCTGCTACTTTCTTCATTGACTCCACCATT"
            "TCGTCGTATCGCTCAGCCTGCTCGGCCAGCTTCGCCTGGTACACCAGATC"
            "CTCTCGATCATCCATAGCGGCAGCGGCTCCGGCAGGGTCTGCGCGACGGA"
            "TGGAAGCGGATAGTGTCTCCGACTCTCTCAGCCTCTCGCTCCGCGTCCGG"
            "GCAGCAAAAATGGCGGCGCCTCAATCC",
            "query_name": "one",
            "lt_soft_len": 0,
            "rt_soft_len": 375,
            "read_match_size": 302,
            "reference_match_size": 302,
            "indel_size": 0,
            "cigartuples_without_soft": [[0, 302]],
            "cigartuples": [(0, 302), (4, 375)],
            "query_length": 677,
            "adhocsms": None,
            "adhocseq": None,
            "mode": None,
            "sms": (0, 302, 375),
            "ref_end": 1423651,
        },
    ]

    read_instances = []
    for read_param in param_dict:
        read_instances.append(
            Read.init(**{param: read_param[param] for param in read_init_params})  # type: ignore
        )
    return read_instances


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


@pytest.fixture()
def trans_same_strand_reads():
    """Reads fixture."""
    read_init_params = [
        "query_name",
        "chrom",
        "ref_start",
        "strand",
        "cigar_str",
        "mapq",
        "nm",
        "query_seq",
    ]

    param_dict = [
        {
            "chrom": "chr17",
            "ref_start": 1744991,
            "strand": "+",
            "cigar_str": "67M116N39M119N76M1652S",
            "mapq": 60,
            "nm": 0,
            "query_seq": "GAACATGGCGCTGCTCTGGGGGCTCCTGGTGCTCAGCTGGTCCTGCCTGC"
            "AAGGCCCCTGCTCCGTGTTCTCCCCTGTGAGCGCCATGGAGCCCTTGGGC"
            "CGGCAGCTAACTAGCGGGCCGAACCAGGAGCAGGTGTCCCCACTTACCCT"
            "CCTCAAGTTGGGCAACCAGGTACAACCAGGTGTACCCCGTCCAGGAAGCC"
            "CTGGCCGTGCTGGAGCCCTATGCGCGGCTGCCCCCGCACAAGCATGTGGC"
            "TCGGCCCACTGAGGTCCTGGCTGGTACCCAGCTCCTCTACGCCTTTTTCA"
            "CTCGGACCCATGGGGACATGCACAGCCTGGTGCGAAGCCGCCACCGTATC"
            "CCTGAGCCTGAGGCTGCCGTGCTCTTCCGCCAGATGGCCACCGCCCTGGC"
            "GCACTGTCACCAGCACGGTCTGGTCCTGCGTGATCTCAAGCTGTGTCGCT"
            "TTGTCTTCGCTGACCGTGAGAGGAAGAAGCTGGTGCTGGAGAACCTGGAG"
            "GACTCCTGCGTGCTGACTGGGCCAGATGATTCCCTGTGGGACAAGCACGC"
            "GTGCCCAGCCTACGTGGGACCTGAGATACTCAGCTCACGGGCCTCATACT"
            "CGGGCAAGGCAGCCGATGTCTGGAGCCTGGGCGTGGCGCTCTTCACCATG"
            "CTGGCCGGCCACTACCCCTTCCAGGACTCGGAGCCTGTCCTGCTCTTCGG"
            "CAAGATCCGCCGCGGGGCCTACGCCTTGCCTGCAGGCCTCTCGGCCCCTG"
            "CCCGCTGTCTGGTTCGCTGCCTCCTTCGTCGGGAGCCAGCTGAACGGCTC"
            "ACAGCCACAGGCATCCTCCTGCACCCCTGGCTGCGACAGGACCCGATGCC"
            "CTTAGCCCCAACCCGATCCCATCTCTGGGAGGCTGCCCAGGTGGTCCCTG"
            "ATGGACTGGGGCTGGACGAAGCCAGGGAAGAGGAGGGAGACAGAGAAGTG"
            "GTTCTGTATGGCTAGGACCACCCTACTACACGCTCAGCTGCCAACAGTGG"
            "ATTGAGTTTGGGGGTAGCTCCAAGCCTTCTCCTGCCTCTGAACTGAGCCA"
            "AACCTTCAGTGCCTTCCAGAAGGGAGAAAGGCAGAAGCCTGTGTGGAGTG"
            "TGCTGTGTACACATCTGCTTTGTTCCACACACATGCAGTTCCTGCTTGGG"
            "TGCTTATCAGGTGCCAAGCCCTGTTCTCGGTGCTGGGAGTACAGCAGTGA"
            "GCAAAGGAGACAATATTCCCTGCTCACAGAGATGACAAACTGGCATCCTT"
            "GAGCTGACAACACTTTTCCATGACCATAGGTCACTGTCTACACTGGGTAC"
            "ACTTTGTACCAGTGTCGGCCTCCACTGATGCTGGTGCTCAGGCACCTCTG"
            "TCCAAGGACAATCCCTTTCACAAACAAACCAGCTGCCTTTGTATCTTGTA"
            "CCTTTTCAGAGAAAGGGAGGTATCCCTGTGCCAAAGGCTCCAGGCCTCTC"
            "CCCTGCAACTCAGGACCCAAGCCCAGCTCACTCTGGGAACTGTGTTCCCA"
            "GCATCTCTGTCCTCTTGATTAAGAGATTCTCCTTCCAGGCCTAAGCCTGG"
            "GATTTGGGCCAGAGATAAGAATCCAAACTATGAGGCTAGTTCTTGTCTAA"
            "CTCAAGACTGTTCTGGAATGAGGGTCCAGGCCTGTCAACCATGGGGCTTC"
            "TGACCTGAGCACCAAGGTTGAGGGACAGGATTAGGCAGGGTCTGTCCTGT"
            "GGCCACCTGGAAAGTCCCAGGTGGGACTCTTCTGGGGACACTTGGGGTCC"
            "ACAATCCCAGGTCCATACTCTAGGTTTTGGATACCATGAGTATGTATGTT"
            "TACCTGTGCCTAATAAAGGAGAATTATGAAATAA",
            "query_name": "two",
            "lt_soft_len": 0,
            "rt_soft_len": 1652,
            "read_match_size": 182,
            "reference_match_size": 417,
            "indel_size": 235,
            "cigartuples_without_soft": [[0, 67], [3, 116], [0, 39], [3, 119], [0, 76]],
            "cigartuples": [(0, 67), (3, 116), (0, 39), (3, 119), (0, 76), (4, 1652)],
            "query_length": 1834,
            "adhocsms": None,
            "adhocseq": None,
            "mode": None,
            "sms": (0, 182, 1652),
            "ref_end": 1745408,
        },
        {
            "chrom": "chr20",
            "ref_start": 391283,
            "strand": "+",
            "cigar_str": "176S296M4618N1362M",
            "mapq": 60,
            "nm": 0,
            "query_seq": "GAACATGGCGCTGCTCTGGGGGCTCCTGGTGCTCAGCTGGTCCTGCCTGC"
            "AAGGCCCCTGCTCCGTGTTCTCCCCTGTGAGCGCCATGGAGCCCTTGGGC"
            "CGGCAGCTAACTAGCGGGCCGAACCAGGAGCAGGTGTCCCCACTTACCCT"
            "CCTCAAGTTGGGCAACCAGGTACAACCAGGTGTACCCCGTCCAGGAAGCC"
            "CTGGCCGTGCTGGAGCCCTATGCGCGGCTGCCCCCGCACAAGCATGTGGC"
            "TCGGCCCACTGAGGTCCTGGCTGGTACCCAGCTCCTCTACGCCTTTTTCA"
            "CTCGGACCCATGGGGACATGCACAGCCTGGTGCGAAGCCGCCACCGTATC"
            "CCTGAGCCTGAGGCTGCCGTGCTCTTCCGCCAGATGGCCACCGCCCTGGC"
            "GCACTGTCACCAGCACGGTCTGGTCCTGCGTGATCTCAAGCTGTGTCGCT"
            "TTGTCTTCGCTGACCGTGAGAGGAAGAAGCTGGTGCTGGAGAACCTGGAG"
            "GACTCCTGCGTGCTGACTGGGCCAGATGATTCCCTGTGGGACAAGCACGC"
            "GTGCCCAGCCTACGTGGGACCTGAGATACTCAGCTCACGGGCCTCATACT"
            "CGGGCAAGGCAGCCGATGTCTGGAGCCTGGGCGTGGCGCTCTTCACCATG"
            "CTGGCCGGCCACTACCCCTTCCAGGACTCGGAGCCTGTCCTGCTCTTCGG"
            "CAAGATCCGCCGCGGGGCCTACGCCTTGCCTGCAGGCCTCTCGGCCCCTG"
            "CCCGCTGTCTGGTTCGCTGCCTCCTTCGTCGGGAGCCAGCTGAACGGCTC"
            "ACAGCCACAGGCATCCTCCTGCACCCCTGGCTGCGACAGGACCCGATGCC"
            "CTTAGCCCCAACCCGATCCCATCTCTGGGAGGCTGCCCAGGTGGTCCCTG"
            "ATGGACTGGGGCTGGACGAAGCCAGGGAAGAGGAGGGAGACAGAGAAGTG"
            "GTTCTGTATGGCTAGGACCACCCTACTACACGCTCAGCTGCCAACAGTGG"
            "ATTGAGTTTGGGGGTAGCTCCAAGCCTTCTCCTGCCTCTGAACTGAGCCA"
            "AACCTTCAGTGCCTTCCAGAAGGGAGAAAGGCAGAAGCCTGTGTGGAGTG"
            "TGCTGTGTACACATCTGCTTTGTTCCACACACATGCAGTTCCTGCTTGGG"
            "TGCTTATCAGGTGCCAAGCCCTGTTCTCGGTGCTGGGAGTACAGCAGTGA"
            "GCAAAGGAGACAATATTCCCTGCTCACAGAGATGACAAACTGGCATCCTT"
            "GAGCTGACAACACTTTTCCATGACCATAGGTCACTGTCTACACTGGGTAC"
            "ACTTTGTACCAGTGTCGGCCTCCACTGATGCTGGTGCTCAGGCACCTCTG"
            "TCCAAGGACAATCCCTTTCACAAACAAACCAGCTGCCTTTGTATCTTGTA"
            "CCTTTTCAGAGAAAGGGAGGTATCCCTGTGCCAAAGGCTCCAGGCCTCTC"
            "CCCTGCAACTCAGGACCCAAGCCCAGCTCACTCTGGGAACTGTGTTCCCA"
            "GCATCTCTGTCCTCTTGATTAAGAGATTCTCCTTCCAGGCCTAAGCCTGG"
            "GATTTGGGCCAGAGATAAGAATCCAAACTATGAGGCTAGTTCTTGTCTAA"
            "CTCAAGACTGTTCTGGAATGAGGGTCCAGGCCTGTCAACCATGGGGCTTC"
            "TGACCTGAGCACCAAGGTTGAGGGACAGGATTAGGCAGGGTCTGTCCTGT"
            "GGCCACCTGGAAAGTCCCAGGTGGGACTCTTCTGGGGACACTTGGGGTCC"
            "ACAATCCCAGGTCCATACTCTAGGTTTTGGATACCATGAGTATGTATGTT"
            "TACCTGTGCCTAATAAAGGAGAATTATGAAATAA",
            "query_name": "two",
            "lt_soft_len": 176,
            "rt_soft_len": 0,
            "read_match_size": 1658,
            "reference_match_size": 6276,
            "indel_size": 4618,
            "cigartuples_without_soft": [[0, 296], [3, 4618], [0, 1362]],
            "cigartuples": [(4, 176), (0, 296), (3, 4618), (0, 1362)],
            "query_length": 1834,
            "adhocsms": None,
            "adhocseq": None,
            "mode": None,
            "sms": (176, 1658, 0),
            "ref_end": 397559,
        },
    ]

    read_instances = []
    for read_param in param_dict:
        read_instances.append(
            Read.init(**{param: read_param[param] for param in read_init_params})  # type: ignore
        )
    return read_instances


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


@pytest.fixture()
def trans_diff_strand_reads():
    """Reads fixture."""
    read_init_params = [
        "query_name",
        "chrom",
        "ref_start",
        "strand",
        "cigar_str",
        "mapq",
        "nm",
        "query_seq",
    ]

    param_dict = [
        {
            "chrom": "chr17",
            "ref_start": 1744991,
            "strand": "+",
            "cigar_str": "67M116N39M119N73M3550S",
            "mapq": 60,
            "nm": 0,
            "query_seq": "GAACATGGCGCTGCTCTGGGGGCTCCTGGTGCTCAGCTGGTCCTGCCTGC"
            "AAGGCCCCTGCTCCGTGTTCTCCCCTGTGAGCGCCATGGAGCCCTTGGGC"
            "CGGCAGCTAACTAGCGGGCCGAACCAGGAGCAGGTGTCCCCACTTACCCT"
            "CCTCAAGTTGGGCAACCAGGTACAACCAGATTGTGTTGTATCGCGAGCAG"
            "GAAGTCCTGGACTGTGACTGTGACATGGCCTCGGTCCACCACCTGTTGTC"
            "CCAGATCCCTCAGGACTTGCCCTATGAGACACTGATCAGCAGAGCAGGAG"
            "ACCTTTTTGTTCAGTTTCCCCCATCCGAACTTGCTCGGGAGGCCGCTGCC"
            "CAACAGCAAGCTGAGAGGACGGCAGCCTCTACTTTCAAAGACTTTGAGCT"
            "GGCATCAGCCCAGCAGAGGCCTGATATGGTGCTGCGGCAGCGGTTTCGGG"
            "GACTTCTGCGGCCTGAAGATCGAACAAAAGATGTCCTGACCAAGCCAAGG"
            "ACCAACCGCTTTGTGAAATTGGCAGTGATGGGGCTGACAGTGGCACTTGG"
            "AGCGGCTGCACTGGCTGTGGTGAAAAGTGCCCTGGAATGGGCCCCTAAGT"
            "TTCAGCTGCAGCTGTTTCCCTGAAAGCCAGAGAAGACCTTCCTCTTACAT"
            "CACATTAAGGTCTCACCCTTCCATGGAAGGATTGGGAGGGGTGGAATTTC"
            "CTTTATTAAAAGGGTTTCTGTCCAGAGTTTTTTATTCCTGCCACCCTGCC"
            "CTGCCCGTTTGTGTTTGCCTTTACTTCAGAGGCCAGTAGCAGAGCCTACC"
            "TGACACAGAGGGGACCTGGGGGCCTGGCACAGAATCTTATGCCTCTCCCA"
            "GTGGGCTTTTGGATGGGCTCTCTGCCACCTCCCTGGTGTTGGAGGAAGTG"
            "GGTTTCTGAAGCTTGCAGCTGCTTTTGCAACCTGCTGGGCCTGGGCTGCC"
            "CTAGTTACCGGGAACCTCCTAGTCACCCCATAGAAGGTGAACTGGAACAT"
            "GGGATGCCCTCTTACTGCTCCTGCTTGCCTGCATCAGATGGACACAGAAG"
            "GAGGGGGCTGCTGCTGAAGGAGTCAGGCTTGGGATGGAGAGCAAGAGGCC"
            "AAAGGACTCGTCTGGGAGCATGATCGGCCCCTCTTACTAGCAGCTCCATC"
            "AGCTGATCAGAAATCAGAGACAGAAGGGAACTAGGTTACTCGTTAAATCC"
            "TACCTATTTTCAGCCTAGAAATCACATGGGCGTTTCCAGTCAAGTGGGAA"
            "GGAAGCTCATTTTTAGGGCCTAGACCTAGACTAGTTTTGCCTGCCTTTTT"
            "TTCCTCCCTGCTTTCTTCTACAAGCCACACAACATTGCTTATAAAAGGTA"
            "CTTTTATTACTATTTTTAGAATATATACACCCCATATGTAAATTACCCTC"
            "ACTTTTAAAGTAATAATAGTTCCTTCACGGAAAAAAAAATAGTTCTCATC"
            "AGAAGGATGAGCTTCTAAGATCCTGAATTTTACTTTCTACTTGATTAACA"
            "TTTTCGCTTTAGCCAGAGTTGCTATTTGCAAACAGCTTTCCACTTCCAAA"
            "GTGCCACTTTCTGACCGGTTTTTCTTCCTGTTTTCCTCAAATCATCCTAC"
            "TCCCAAATAGGCACCCACTCACTACCTTGGCGTCGTTTTTTGGGTCTGAC"
            "TTGACTCGTCAACTGCTGGTCTCTCCCACCTGGTTGGAAATGTCGTTGGA"
            "AACTTGCAAAGACTCTCCAGACCTTAGGGAACAAGAGGCATCACTCAGTC"
            "CTTCTGGGACAGCTTCCCTGTAAGCAAAGCCAGGTTTGTTTTAGATTCAC"
            "TCTGGCCTGAAAAAAAAAGTGCCTTTTGCTGCTTTAAAGAATTGGGGTAT"
            "ATGGTATGAAGCAGCCATGTACTTGTATTTTCCTGGTCTTTCCTGGGCAC"
            "TCTTCTCTCTTGGCAGATGTTTTCTTAAAGTGAACACACCAGAAGCACTC"
            "TACCCCACCATATCCAGTCCTCTTTGCATGCCTGTTGTGTGCAGGGTGGA"
            "GGGTGTCCCACGTGCCAGAGCAGCTGGGAGCCCTGTCCGCTGAGTTACAC"
            "TAAGCACTTTATTTTTATTTATTTTATTATTATTTTCTTTTTTCTGAGAT"
            "GGAGTCTCACTCTGTTGCCCAGGCTGGAGTGCAGTGGTGCAGTCTCAGCT"
            "CACTGCAACCTCCGCCTCCCAAGTTCAGGCAATTCTCGCTGCCTCAGCCT"
            "CCCGAGTAGCTGGGACTACAGGCACCTGCTACCACGCCTGACTAATTTTT"
            "GTATTTTTAGTAGAGATGGGGTTTCGCCACGTTGGCCAGGCTGGTCTCGA"
            "ACTCCTGACCTCAGGTGATCTGCCCGTCTCGTCCTCCCAAAGTGCTGGGA"
            "TTACAGGTGTGAGCCACCACGCCCAGCCACCACTACGCACTTTAAAACCA"
            "AAAAAAGCTCAAAGGGATGGATCATTGTCCATGGCTCAGGAGCCCCACTG"
            "TGGGGTAAGTGTGATGTCGATGGATTAATTTAGACAGTGGTGTATTCATC"
            "AGGGCAAAGGAAATTCTTCCAGAGGTGATGATGCGGCAGATATTCCCAGT"
            "GCTCTCCCAAAGCAAGCTGGCCCTCTTGTCATCCTCTCAGTGGGCGGCAC"
            "CAGCTGCCCCTCCTTCCACAAGTACTCAAGCCTGTTTGTAAATACTGAAG"
            "GAATTGATGGGGTTGAGGAAAGGAGGTGCATGTGACCAGGGTCCCAAGGC"
            "CACAGCTTTTCAGATCCTAGGAAGCAAGTGGCATTTGCTTGAGTTGTGGC"
            "CTCGGAAGGAGAATGTTTATCTGTTTTCTAACTTTGCTGACACCAGGATT"
            "CTCCCTGTCATTGAGAAGAAAGCATTATCTAATTACCTTCAGGTGGTTTA"
            "CTTATTCTGTAAAGAATATGTGTAAATATTTTGTACAGAGCCCTGTATCA"
            "AATAAACAGCCATATGTGGTTACTAATCACCTCTTCTGTCATTCCGTCCT"
            "TGGCCACCGCTCAGTGGGAATGGTCTCTGATCTGGATGCTCCCACCTTCC"
            "ATGTCAGGCCCAGAACTGTGCCATGGTCTGTGGACTCCTGGTCAGCCTTG"
            "ACTGGCTAGGAGACCTTGGGCAGTACCTACAGTCTTGCTGTTTCTGTTTC"
            "ATCTGCAAGAATTATGACCCACACACTCCAGCTGCAGCCCAGGGCACTGT"
            "GATATTTTATACGTGTGTAGATGTTTTTGTCCACAGTTCCTGGTTCATCA"
            "CTCCCATAACCCTTTGTTATAATGTTGGGACACTGCAGGCCTCAGAAAAC"
            "GGAATCTCTGTCTGTGACCTTCTCCTGCCCCATTTCACTTGCTCAACACC"
            "AGACTTTAATCTGACTGTAGCTCATAAGACCCTCATTCCAGAGAGGGTGC"
            "TGCCCCATACCCGGAAGGAGGAACGCTGCACAGAGAGGCCAAGAAGCATC"
            "TGGACAGACAGGCCTTGCTGGGTTTAGACCTTATGCTTTTTGTCCAGTTT"
            "CATCTCAACACAGCTGCCATGCTTCAGCCATGCCTATCCAATGACGTCTC"
            "CATAAAAGGCCCAGGAACACGGGAGCTTCTGAAGAGCTGAACATGTGGAG"
            "GGAGGGGAACGAGAACTTGTCCATGTGCCAAGAGGGTGGCGCACCCCCAC"
            "TCCATGGGGACAGAAGCTCCAGCATTTGCCCAGGACCCGTCCAGACCTCA"
            "CCCTGTGTGTATCTTCATCTGGCTGTTTACTTATTTGTATCCTTTTCTAA"
            "TAATGTTTGTAATAAACTGGTAAACATAA",
            "query_name": "three",
            "lt_soft_len": 0,
            "rt_soft_len": 3550,
            "read_match_size": 179,
            "reference_match_size": 414,
            "indel_size": 235,
            "cigartuples_without_soft": [[0, 67], [3, 116], [0, 39], [3, 119], [0, 73]],
            "cigartuples": [(0, 67), (3, 116), (0, 39), (3, 119), (0, 73), (4, 3550)],
            "query_length": 3729,
            "adhocsms": None,
            "adhocseq": None,
            "mode": None,
            "sms": (0, 179, 3550),
            "ref_end": 1745405,
        },
        {
            "chrom": "chr20",
            "ref_start": 435479,
            "strand": "-",
            "cigar_str": "3362M266N191M176S",
            "mapq": 60,
            "nm": 0,
            "query_seq": "TTATGTTTACCAGTTTATTACAAACATTATTAGAAAAGGATACAAATAAG"
            "TAAACAGCCAGATGAAGATACACACAGGGTGAGGTCTGGACGGGTCCTGG"
            "GCAAATGCTGGAGCTTCTGTCCCCATGGAGTGGGGGTGCGCCACCCTCTT"
            "GGCACATGGACAAGTTCTCGTTCCCCTCCCTCCACATGTTCAGCTCTTCA"
            "GAAGCTCCCGTGTTCCTGGGCCTTTTATGGAGACGTCATTGGATAGGCAT"
            "GGCTGAAGCATGGCAGCTGTGTTGAGATGAAACTGGACAAAAAGCATAAG"
            "GTCTAAACCCAGCAAGGCCTGTCTGTCCAGATGCTTCTTGGCCTCTCTGT"
            "GCAGCGTTCCTCCTTCCGGGTATGGGGCAGCACCCTCTCTGGAATGAGGG"
            "TCTTATGAGCTACAGTCAGATTAAAGTCTGGTGTTGAGCAAGTGAAATGG"
            "GGCAGGAGAAGGTCACAGACAGAGATTCCGTTTTCTGAGGCCTGCAGTGT"
            "CCCAACATTATAACAAAGGGTTATGGGAGTGATGAACCAGGAACTGTGGA"
            "CAAAAACATCTACACACGTATAAAATATCACAGTGCCCTGGGCTGCAGCT"
            "GGAGTGTGTGGGTCATAATTCTTGCAGATGAAACAGAAACAGCAAGACTG"
            "TAGGTACTGCCCAAGGTCTCCTAGCCAGTCAAGGCTGACCAGGAGTCCAC"
            "AGACCATGGCACAGTTCTGGGCCTGACATGGAAGGTGGGAGCATCCAGAT"
            "CAGAGACCATTCCCACTGAGCGGTGGCCAAGGACGGAATGACAGAAGAGG"
            "TGATTAGTAACCACATATGGCTGTTTATTTGATACAGGGCTCTGTACAAA"
            "ATATTTACACATATTCTTTACAGAATAAGTAAACCACCTGAAGGTAATTA"
            "GATAATGCTTTCTTCTCAATGACAGGGAGAATCCTGGTGTCAGCAAAGTT"
            "AGAAAACAGATAAACATTCTCCTTCCGAGGCCACAACTCAAGCAAATGCC"
            "ACTTGCTTCCTAGGATCTGAAAAGCTGTGGCCTTGGGACCCTGGTCACAT"
            "GCACCTCCTTTCCTCAACCCCATCAATTCCTTCAGTATTTACAAACAGGC"
            "TTGAGTACTTGTGGAAGGAGGGGCAGCTGGTGCCGCCCACTGAGAGGATG"
            "ACAAGAGGGCCAGCTTGCTTTGGGAGAGCACTGGGAATATCTGCCGCATC"
            "ATCACCTCTGGAAGAATTTCCTTTGCCCTGATGAATACACCACTGTCTAA"
            "ATTAATCCATCGACATCACACTTACCCCACAGTGGGGCTCCTGAGCCATG"
            "GACAATGATCCATCCCTTTGAGCTTTTTTTGGTTTTAAAGTGCGTAGTGG"
            "TGGCTGGGCGTGGTGGCTCACACCTGTAATCCCAGCACTTTGGGAGGACG"
            "AGACGGGCAGATCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACG"
            "TGGCGAAACCCCATCTCTACTAAAAATACAAAAATTAGTCAGGCGTGGTA"
            "GCAGGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGCGAGAATTG"
            "CCTGAACTTGGGAGGCGGAGGTTGCAGTGAGCTGAGACTGCACCACTGCA"
            "CTCCAGCCTGGGCAACAGAGTGAGACTCCATCTCAGAAAAAAGAAAATAA"
            "TAATAAAATAAATAAAAATAAAGTGCTTAGTGTAACTCAGCGGACAGGGC"
            "TCCCAGCTGCTCTGGCACGTGGGACACCCTCCACCCTGCACACAACAGGC"
            "ATGCAAAGAGGACTGGATATGGTGGGGTAGAGTGCTTCTGGTGTGTTCAC"
            "TTTAAGAAAACATCTGCCAAGAGAGAAGAGTGCCCAGGAAAGACCAGGAA"
            "AATACAAGTACATGGCTGCTTCATACCATATACCCCAATTCTTTAAAGCA"
            "GCAAAAGGCACTTTTTTTTTCAGGCCAGAGTGAATCTAAAACAAACCTGG"
            "CTTTGCTTACAGGGAAGCTGTCCCAGAAGGACTGAGTGATGCCTCTTGTT"
            "CCCTAAGGTCTGGAGAGTCTTTGCAAGTTTCCAACGACATTTCCAACCAG"
            "GTGGGAGAGACCAGCAGTTGACGAGTCAAGTCAGACCCAAAAAACGACGC"
            "CAAGGTAGTGAGTGGGTGCCTATTTGGGAGTAGGATGATTTGAGGAAAAC"
            "AGGAAGAAAAACCGGTCAGAAAGTGGCACTTTGGAAGTGGAAAGCTGTTT"
            "GCAAATAGCAACTCTGGCTAAAGCGAAAATGTTAATCAAGTAGAAAGTAA"
            "AATTCAGGATCTTAGAAGCTCATCCTTCTGATGAGAACTATTTTTTTTTC"
            "CGTGAAGGAACTATTATTACTTTAAAAGTGAGGGTAATTTACATATGGGG"
            "TGTATATATTCTAAAAATAGTAATAAAAGTACCTTTTATAAGCAATGTTG"
            "TGTGGCTTGTAGAAGAAAGCAGGGAGGAAAAAAAGGCAGGCAAAACTAGT"
            "CTAGGTCTAGGCCCTAAAAATGAGCTTCCTTCCCACTTGACTGGAAACGC"
            "CCATGTGATTTCTAGGCTGAAAATAGGTAGGATTTAACGAGTAACCTAGT"
            "TCCCTTCTGTCTCTGATTTCTGATCAGCTGATGGAGCTGCTAGTAAGAGG"
            "GGCCGATCATGCTCCCAGACGAGTCCTTTGGCCTCTTGCTCTCCATCCCA"
            "AGCCTGACTCCTTCAGCAGCAGCCCCCTCCTTCTGTGTCCATCTGATGCA"
            "GGCAAGCAGGAGCAGTAAGAGGGCATCCCATGTTCCAGTTCACCTTCTAT"
            "GGGGTGACTAGGAGGTTCCCGGTAACTAGGGCAGCCCAGGCCCAGCAGGT"
            "TGCAAAAGCAGCTGCAAGCTTCAGAAACCCACTTCCTCCAACACCAGGGA"
            "GGTGGCAGAGAGCCCATCCAAAAGCCCACTGGGAGAGGCATAAGATTCTG"
            "TGCCAGGCCCCCAGGTCCCCTCTGTGTCAGGTAGGCTCTGCTACTGGCCT"
            "CTGAAGTAAAGGCAAACACAAACGGGCAGGGCAGGGTGGCAGGAATAAAA"
            "AACTCTGGACAGAAACCCTTTTAATAAAGGAAATTCCACCCCTCCCAATC"
            "CTTCCATGGAAGGGTGAGACCTTAATGTGATGTAAGAGGAAGGTCTTCTC"
            "TGGCTTTCAGGGAAACAGCTGCAGCTGAAACTTAGGGGCCCATTCCAGGG"
            "CACTTTTCACCACAGCCAGTGCAGCCGCTCCAAGTGCCACTGTCAGCCCC"
            "ATCACTGCCAATTTCACAAAGCGGTTGGTCCTTGGCTTGGTCAGGACATC"
            "TTTTGTTCGATCTTCAGGCCGCAGAAGTCCCCGAAACCGCTGCCGCAGCA"
            "CCATATCAGGCCTCTGCTGGGCTGATGCCAGCTCAAAGTCTTTGAAAGTA"
            "GAGGCTGCCGTCCTCTCAGCTTGCTGTTGGGCAGCGGCCTCCCGAGCAAG"
            "TTCGGATGGGGGAAACTGAACAAAAAGGTCTCCTGCTCTGCTGATCAGTG"
            "TCTCATAGGGCAAGTCCTGAGGGATCTGGGACAACAGGTGGTGGACCGAG"
            "GCCATGTCACAGTCACAGTCCAGGACTTCCTGCTCGCGATACAACACAAT"
            "CTGGTTGTACCTGGTTGCCCAACTTGAGGAGGGTAAGTGGGGACACCTGC"
            "TCCTGGTTCGGCCCGCTAGTTAGCTGCCGGCCCAAGGGCTCCATGGCGCT"
            "CACAGGGGAGAACACGGAGCAGGGGCCTTGCAGGCAGGACCAGCTGAGCA"
            "CCAGGAGCCCCCAGAGCAGCGCCATGTTC",
            "query_name": "three",
            "lt_soft_len": 0,
            "rt_soft_len": 176,
            "read_match_size": 3553,
            "reference_match_size": 3819,
            "indel_size": 266,
            "cigartuples_without_soft": [[0, 3362], [3, 266], [0, 191]],
            "cigartuples": [(0, 3362), (3, 266), (0, 191), (4, 176)],
            "query_length": 3729,
            "adhocsms": None,
            "adhocseq": None,
            "mode": None,
            "sms": (0, 3553, 176),
            "ref_end": 439298,
        },
    ]

    read_instances = []
    for read_param in param_dict:
        read_instances.append(
            Read.init(**{param: read_param[param] for param in read_init_params})  # type: ignore
        )
    return read_instances


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


@pytest.fixture()
def inv_reads():
    """Reads fixture."""
    read_init_params = [
        "query_name",
        "chrom",
        "ref_start",
        "strand",
        "cigar_str",
        "mapq",
        "nm",
        "query_seq",
    ]

    param_dict = [
        {
            "chrom": "chr17",
            "ref_start": 1650628,
            "strand": "-",
            "cigar_str": "328M151N203M103N141M935S",
            "mapq": 60,
            "nm": 0,
            "query_seq": "TCCTATTTATACAAAATTTATTATTATATTTTATTCAGGATGACAAGCCA"
            "TCAGGAGGTCAACAACACAAGCACAGACAGAGGGAAAGAGGCCAAACTGC"
            "TGAATGTCAGCGGCCTGTCTGGAGGGGCTGAGGCTTCGGCCTCGGGAGGC"
            "TGAAGCAGGAGGCAGGGAAACGGTCAGGCATACAGGTCCTCCCGATCCGC"
            "AGAGTAAACCTCCCCCTCCTGCAGGAGAGCAAAGTTGAGGAAGTGAGAGG"
            "GCCTGTGCACCTCGTGGTAGAACTCTTTGGGGTTCGCCAGCTGTAGCTCA"
            "TATTTCATGTTGGGGTCATGCCGAACACCCATGAAGTTGTAGTTCCACGA"
            "GGACTGGGCAGGGACCATGAAGAAGCCAAGGAAACGGTCCGACAGCAGCA"
            "TCTGCACCCTCTCATAGTGTGAAGGCAGGTAGCCCTTGGGGTTGTTGCCC"
            "TTGTCTGTGTTCTGGCGGCCCCATTCGTAGCCACTGGGGGTCAGCTTGTA"
            "GGCCGTCAGTGTACAGGAGCCTGGCGTGAAGCTGCATGTGATGATAATGG"
            "TCTTCTCGCCATCCCAAGATGGGTTGTCAGCCATGATCTTGGCATGGGTG"
            "GTGACATCCTGGGGTGATAACTGCGGGGACTCATTGGGCTGAGTGTGGAT"
            "CCAACCTAAGGGTTCCATCTCCAGGTAGCTGATGTAGGGCAGGTACTGGT"
            "AGGTGAGGACAGGCTCCCCATACAGGCGGGCGATGTGGAGGAGGCAGCTG"
            "AGCACAGGCCCTGACACGATGTCGCCCAGGACCGGCCTCTTCTGGTAGAT"
            "GTTGCCGGCGCTCAGCGGTGGGCTCTCGCCACTGCTCACTGTGAACTGCT"
            "GCCGAGTGGGTCCAACATAACAAGACGTCAGCAGGCGGAGCAGGTTCCGG"
            "GCCACGTGGCGAGAGGCCACTGTGGGGCCGAGCTTGGCAGACAGCCAGCG"
            "GACCATCTTGCAGGCTGTATCAAGGAGGATCTTCTGTTCTTTGCCTTCTG"
            "ACTGCTCAGGCAGTGCCTCCTCCTCTTCTTCTCCATCTCCCCCACTGCCG"
            "CCGGCCACAACCGTCTCCATGGACAGCACCGTGTCAGACAGAGTGAGCTC"
            "AGATGCCCCGGTGACCTCCTCCTGCTCCCCCTCCTCCTCCTCTAGCACCA"
            "CGCAGCTGTCCTCCTCCTCCTCTTCCTCCTCGGAGCCCTCGCTTTGCTTC"
            "AAGTCCTGGCTGCTGTCACCTGATCGAAGGCTGCTCTTGTCCACGGGGGC"
            "ACCCCCCTCGTCTGGAGCCCGCTCCTCACCCAGGGAGGTCTCGCTGGTGC"
            "TGCTCTTGTCACTCAGCCGGCCCAGGCTCACAGCCTCAGCCTCCTGGGGC"
            "TGGGGAGACTCAGTCACATAGAGCCCGGCTTGGAAGTCCTCTGTCTCAGG"
            "GAGGTCAGCCTGGTCGTGGAAGCTGACGCCAGACGTGTAGTCTGGGAGCC"
            "CCAGGCCCGAGGAGGCAGGAGGCTCCCCATCCATGGGAATCTCCTCCCCA"
            "AAAGCACAGGAGCCAGGCCCGGCCCCGGGCAGCCCGCTCTCCTCCTCCTC"
            "AGCAGCCCCTGCCAGGTCCTTGCTCTCCTCCTGGGAGGCCTCTGCGCCCG"
            "CCAGCAC",
            "query_name": "four",
            "lt_soft_len": 0,
            "rt_soft_len": 935,
            "read_match_size": 672,
            "reference_match_size": 926,
            "indel_size": 254,
            "cigartuples_without_soft": [
                [0, 328],
                [3, 151],
                [0, 203],
                [3, 103],
                [0, 141],
            ],
            "cigartuples": [(0, 328), (3, 151), (0, 203), (3, 103), (0, 141), (4, 935)],
            "query_length": 1607,
            "adhocsms": None,
            "adhocseq": None,
            "mode": None,
            "sms": (0, 672, 935),
            "ref_end": 1651554,
        },
        {
            "chrom": "chr17",
            "ref_start": 1727989,
            "strand": "+",
            "cigar_str": "637M1753N108M267N192M670S",
            "mapq": 60,
            "nm": 0,
            "query_seq": "GTGCTGGCGGGCGCAGAGGCCTCCCAGGAGGAGAGCAAGGACCTGGCAGG"
            "GGCTGCTGAGGAGGAGGAGAGCGGGCTGCCCGGGGCCGGGCCTGGCTCCT"
            "GTGCTTTTGGGGAGGAGATTCCCATGGATGGGGAGCCTCCTGCCTCCTCG"
            "GGCCTGGGGCTCCCAGACTACACGTCTGGCGTCAGCTTCCACGACCAGGC"
            "TGACCTCCCTGAGACAGAGGACTTCCAAGCCGGGCTCTATGTGACTGAGT"
            "CTCCCCAGCCCCAGGAGGCTGAGGCTGTGAGCCTGGGCCGGCTGAGTGAC"
            "AAGAGCAGCACCAGCGAGACCTCCCTGGGTGAGGAGCGGGCTCCAGACGA"
            "GGGGGGTGCCCCCGTGGACAAGAGCAGCCTTCGATCAGGTGACAGCAGCC"
            "AGGACTTGAAGCAAAGCGAGGGCTCCGAGGAGGAAGAGGAGGAGGAGGAC"
            "AGCTGCGTGGTGCTAGAGGAGGAGGAGGGGGAGCAGGAGGAGGTCACCGG"
            "GGCATCTGAGCTCACTCTGTCTGACACGGTGCTGTCCATGGAGACGGTTG"
            "TGGCCGGCGGCAGTGGGGGAGATGGAGAAGAAGAGGAGGAGGCACTGCCT"
            "GAGCAGTCAGAAGGCAAAGAACAGAAGATCCTCCTTGATACAGCCTGCAA"
            "GATGGTCCGCTGGCTGTCTGCCAAGCTCGGCCCCACAGTGGCCTCTCGCC"
            "ACGTGGCCCGGAACCTGCTCCGCCTGCTGACGTCTTGTTATGTTGGACCC"
            "ACTCGGCAGCAGTTCACAGTGAGCAGTGGCGAGAGCCCACCGCTGAGCGC"
            "CGGCAACATCTACCAGAAGAGGCCGGTCCTGGGCGACATCGTGTCAGGGC"
            "CTGTGCTCAGCTGCCTCCTCCACATCGCCCGCCTGTATGGGGAGCCTGTC"
            "CTCACCTACCAGTACCTGCCCTACATCAGCTACCTGGAGATGGAACCCTT"
            "AGGTTGGATCCACACTCAGCCCAATGAGTCCCCGCAGTTATCACCCCAGG"
            "ATGTCACCACCCATGCCAAGATCATGGCTGACAACCCATCTTGGGATGGC"
            "GAGAAGACCATTATCATCACATGCAGCTTCACGCCAGGCTCCTGTACACT"
            "GACGGCCTACAAGCTGACCCCCAGTGGCTACGAATGGGGCCGCCAGAACA"
            "CAGACAAGGGCAACAACCCCAAGGGCTACCTGCCTTCACACTATGAGAGG"
            "GTGCAGATGCTGCTGTCGGACCGTTTCCTTGGCTTCTTCATGGTCCCTGC"
            "CCAGTCCTCGTGGAACTACAACTTCATGGGTGTTCGGCATGACCCCAACA"
            "TGAAATATGAGCTACAGCTGGCGAACCCCAAAGAGTTCTACCACGAGGTG"
            "CACAGGCCCTCTCACTTCCTCAACTTTGCTCTCCTGCAGGAGGGGGAGGT"
            "TTACTCTGCGGATCGGGAGGACCTGTATGCCTGACCGTTTCCCTGCCTCC"
            "TGCTTCAGCCTCCCGAGGCCGAAGCCTCAGCCCCTCCAGACAGGCCGCTG"
            "ACATTCAGCAGTTTGGCCTCTTTCCCTCTGTCTGTGCTTGTGTTGTTGAC"
            "CTCCTGATGGCTTGTCATCCTGAATAAAATATAATAATAAATTTTGTATA"
            "AATAGGA",
            "query_name": "four",
            "lt_soft_len": 0,
            "rt_soft_len": 670,
            "read_match_size": 937,
            "reference_match_size": 2957,
            "indel_size": 2020,
            "cigartuples_without_soft": [
                [0, 637],
                [3, 1753],
                [0, 108],
                [3, 267],
                [0, 192],
            ],
            "cigartuples": [
                (0, 637),
                (3, 1753),
                (0, 108),
                (3, 267),
                (0, 192),
                (4, 670),
            ],
            "query_length": 1607,
            "adhocsms": None,
            "adhocseq": None,
            "mode": None,
            "sms": (0, 937, 670),
            "ref_end": 1730946,
        },
    ]

    read_instances = []
    for read_param in param_dict:
        read_instances.append(
            Read.init(**{param: read_param[param] for param in read_init_params})  # type: ignore
        )
    return read_instances


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
    expect = None
    assert result == expect
