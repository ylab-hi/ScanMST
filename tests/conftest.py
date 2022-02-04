# !/usr/bin/env python
"""Conftest for pytest."""
from typing import List

import pytest
from tests import assign_value_for_instance
from tests import FakeBlat
from tests import FakeLogger

from scannls import Event
from scannls import Insertion
from scannls import MicroHomology
from scannls import Node
from scannls import NovelInsertion
from scannls import Read


@pytest.fixture(scope="session")
def fake_logger():
    """Fake logger."""
    return FakeLogger()


@pytest.fixture(scope="session")
def fake_blat():
    """Fake blat."""
    return FakeBlat()


@pytest.fixture()
def novel_insertion():
    """Novel insertion."""
    return NovelInsertion(hit_num=0, query_sequence="ATCA")


@pytest.fixture()
def microhomology():
    """Microhomology."""
    return MicroHomology(query_sequence="ATCA")


@pytest.fixture()
def insertion():
    """Insertion."""
    return Insertion(
        hit_num=1,
        chrom="1",
        ref_start=1,
        strand="+",
        cigarstring="1S2M1S",
        mapq=60,
        nm=0,
        query_sequence="AGCT",
    )


@pytest.fixture(scope="function")
def nodes() -> List[Node]:
    """Return a list of nodes."""
    param_dict = [
        # n1
        {
            "chrom": "chr2",
            "ref_start": 190659106,
            "ref_end": 190670461,
            "strand": "+",
            "exons": [[190659106, 190659995], [190670325, 190670461]],
            "sv_type": "TRA",
            "prev_sv_type": None,
            "modes": [1, 2],
            "query_name": "one,two",
            "sr": 2,
            "successors": [],
            "predecessors": [],
            "trace_id": -1,
            "is_traced": False,
            "prev_breakpoint": None,
            "next_breakpoint": "chr2:190670461",
        },
        # n2
        {
            "chrom": "chr17",
            "ref_start": 49502062,
            "ref_end": 49502204,
            "strand": "+",
            "exons": [[49502062, 49502204]],
            "sv_type": "TDUP",
            "prev_sv_type": "TRA",
            "modes": [1, 2],
            "query_name": "one,three,two",
            "sr": 3,
            "successors": [],
            "predecessors": [],
            "trace_id": -1,
            "is_traced": False,
            "prev_breakpoint": "chr2:190670461",
            "next_breakpoint": "chr17:49502204",
        },
    ]
    node_list = [Node() for _ in range(len(param_dict))]
    for ind, node in enumerate(node_list):
        assign_value_for_instance(node, **param_dict[ind])  # type: ignore
        node.get_unique_key()
    return node_list


@pytest.fixture()
def read_param_dict():
    """Return a dict of read parameters."""
    return [
        {
            "chrom": "chr17",
            "ref_start": 7701655,
            "strand": "+",
            "cigar_str": "711S734M",
            "mapq": 60,
            "nm": 0,
            "query_seq": "CGCCCGGTCCCCGGCTCCCCCAGTCCCCCACTTAGGCGGGCTCACAGATCCCGGGGTG"
            "CTGGCGCGTGGGCCGGGGGCGCGTAGGGCGCCTGCAGACGGCCCCTGGAAGGGCTCTG"
            "GTGGGGCTGAGCGCTCTGCCGCGGGGGCGCGGGCACAGCAGGAAGCAGGTCCGCGTGG"
            "GCGCTGGGGGCATCAGCTACCGGGGTGGTCCGGGCTGAAGAGCCAGGCAGCCAAGGCA"
            "GCCACCCCGGGGGGTGGGCGACTTTGGGGGAGTTGGTGCCCCGCCCCCCAGGCCTTGG"
            "CGGGGTCATGGGGCCCCCCCATTCTGGGCCGGGGGGCGTGCGAGTCGGGGCCCTGCTG"
            "CTGCTGGGGGTTTTGGGGCTGGTGTCTGGGCTCAGCCTGGAGCCTGTCTACTGGAACT"
            "CGGCGAATAAGAGGTTCCAGGCAGAGGGTGGTTATGTGCTGTACCCTCAGATCGGGGA"
            "CCGGCTAGACCTGCTCTGCCCCCGGGCCCGGCCTCCTGGCCCTCACTCCTCTCCTAAT"
            "TATGAGTTCTACAAGCTGTACCTGGTAGGGGGTGCTCAGGGCCGGCGCTGTGAGGCAC"
            "CCCCTGCCCCAAACCTCCTTCTCACTTGTGATCGCCCAGACCTGGATCTCCGCTTCAC"
            "CATCAAGTTCCAGGAGTATAGCCCTAATCTCTGGGGCCACGAGTTCCGCTCGCACCAC"
            "GATTACTACATCATTGGATGAGCTGACGGCAGCCCATTCGCTCTGCTTCTCCCCGGAT"
            "GGCTCCCAGCTCTTCTGTGGCTTCAACCGGACTGTGCGTGTTTTTTCCACGGCCCGGC"
            "CTGGCCGAGACTGCGAGGTCCGAGCCACATTTGGTAAGCATCTGTGCCTCCAAGGGAG"
            "GAGGAGAGGGAAGGGCACTGCCACCTGCACAGGGGCCTTTTGTGAGCCGGGGGCCACC"
            "TGTGGGGGTTCACGCCGTCCTCTGTACGGCCCCGGGAGCAGGTGCAGCCCAGTCGGCA"
            "GAGGAGCAAACAGGCTCAGAGCAGGTAGGAAACCTTCCCAAGGCCAACCAGCTGGTCAA"
            "AGGACTGCTTCCTTCCTGAACTCATACCCTGTCAGCTGTGGAGCTTTTGGTCTCTGAAA"
            "TCTTTCTAGAAAATTGTTGATAAAGCTGATTCCGTTTTCCTGTAGGCCTTCAACTTGCA"
            "TCTCTCCAAGGAAGAACTGGGATTTGAGAGGGATGAAGTGGGGCTTGGGCATTTAGGTC"
            "CTTTGGGAGGATAGATGTGGGGAGCATCAGAGGTCTTTGTCCTGCTTGTGACAGACAGC"
            "ATGGGGGGGATGTTGAGTCCAAGCATGTTGGTGCTGGGACGGGAGACAGACCTCTGCTT"
            "AGCCTGGTTAGTGCCAGGAGCCATTGCCCCCTCCCCCACTTTGTTCCTTCCCTCTCTAG"
            "CAAAAAAGCAGGGCCAGAGCGGCATCATCTCCTGCATAGCCTTCAG",
            "query_name": "B",
            "lt_soft_len": 711,
            "rt_soft_len": 0,
            "read_match_size": 734,
            "reference_match_size": 734,
            "indel_size": 0,
            "cigartuples_without_soft": [[0, 734]],
            "cigartuples": [(4, 711), (0, 734)],
            "query_length": 1445,
            "mode": 2,
            "sms": (711, 734, 0),
            "ref_end": 7702389,
        },
        {
            "chrom": "chr17",
            "ref_start": 7705301,
            "strand": "+",
            "cigar_str": "419M2237N294M732S",
            "mapq": 60,
            "nm": 0,
            "query_seq": "CGCCCGGTCCCCGGCTCCCCCAGTCCCCCACTTAGGCGGGCTCACAGATCCCGGGGTG"
            "CTGGCGCGTGGGCCGGGGGCGCGTAGGGCGCCTGCAGACGGCCCCTGGAAGGGCTCTG"
            "GTGGGGCTGAGCGCTCTGCCGCGGGGGCGCGGGCACAGCAGGAAGCAGGTCCGCGTGG"
            "GCGCTGGGGGCATCAGCTACCGGGGTGGTCCGGGCTGAAGAGCCAGGCAGCCAAGGCA"
            "GCCACCCCGGGGGGTGGGCGACTTTGGGGGAGTTGGTGCCCCGCCCCCCAGGCCTTGG"
            "CGGGGTCATGGGGCCCCCCCATTCTGGGCCGGGGGGCGTGCGAGTCGGGGCCCTGCTG"
            "CTGCTGGGGGTTTTGGGGCTGGTGTCTGGGCTCAGCCTGGAGCCTGTCTACTGGAACT"
            "CGGCGAATAAGAGGTTCCAGGCAGAGGGTGGTTATGTGCTGTACCCTCAGATCGGGGA"
            "CCGGCTAGACCTGCTCTGCCCCCGGGCCCGGCCTCCTGGCCCTCACTCCTCTCCTAAT"
            "TATGAGTTCTACAAGCTGTACCTGGTAGGGGGTGCTCAGGGCCGGCGCTGTGAGGCAC"
            "CCCCTGCCCCAAACCTCCTTCTCACTTGTGATCGCCCAGACCTGGATCTCCGCTTCAC"
            "CATCAAGTTCCAGGAGTATAGCCCTAATCTCTGGGGCCACGAGTTCCGCTCGCACCAC"
            "GATTACTACATCATTGGATGAGCTGACGGCAGCCCATTCGCTCTGCTTCTCCCCGGAT"
            "GGCTCCCAGCTCTTCTGTGGCTTCAACCGGACTGTGCGTGTTTTTTCCACGGCCCGGC"
            "CTGGCCGAGACTGCGAGGTCCGAGCCACATTTGGTAAGCATCTGTGCCTCCAAGGGAG"
            "GAGGAGAGGGAAGGGCACTGCCACCTGCACAGGGGCCTTTTGTGAGCCGGGGGCCACC"
            "TGTGGGGGTTCACGCCGTCCTCTGTACGGCCCCGGGAGCAGGTGCAGCCCAGTCGGCA"
            "GAGGAGCAAACAGGCTCAGAGCAGGTAGGAAACCTTCCCAAGGCCAACCAGCTGGTCA"
            "AAGGACTGCTTCCTTCCTGAACTCATACCCTGTCAGCTGTGGAGCTTTTGGTCTCTGA"
            "AATCTTTCTAGAAAATTGTTGATAAAGCTGATTCCGTTTTCCTGTAGGCCTTCAACTT"
            "GCATCTCTCCAAGGAAGAACTGGGATTTGAGAGGGATGAAGTGGGGCTTGGGCATTTA"
            "GGTCCTTTGGGAGGATAGATGTGGGGAGCATCAGAGGTCTTTGTCCTGCTTGTGACAG"
            "ACAGCATGGGGGGGATGTTGAGTCCAAGCATGTTGGTGCTGGGACGGGAGACAGACCT"
            "CTGCTTAGCCTGGTTAGTGCCAGGAGCCATTGCCCCCTCCCCCACTTTGTTCCTTCCC"
            "TCTCTAGCAAAAAAGCAGGGCCAGAGCGGCATCATCTCCTGCATAGCCTTCAG",
            "query_name": "B",
            "lt_soft_len": 0,
            "rt_soft_len": 732,
            "read_match_size": 713,
            "reference_match_size": 2950,
            "indel_size": 2237,
            "cigartuples_without_soft": [[0, 419], [3, 2237], [0, 294]],
            "cigartuples": [(0, 419), (3, 2237), (0, 294), (4, 732)],
            "query_length": 1445,
            "adhocsms": None,
            "adhocseq": None,
            "mode": 1,
            "sms": (0, 713, 732),
            "ref_end": 7708251,
        },
    ]


@pytest.fixture()
def reads(read_param_dict):
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

    read_instances = []
    for read_param in read_param_dict:
        temp = Read.init(**{param: read_param[param] for param in read_init_params})  # type: ignore
        read_instances.append(temp)
    return read_instances


@pytest.fixture()
def event():
    """Event fixture."""
    event_tuple = (
        "TDUP",  # svtype
        0,  # annotation_code
        1,  # splicing code
        ("chr17:7701655", "chr17:7708249", 2, 1),  # bp1, bp2, mode1, mode2
        (7701655, 7702389, [[7701655, 7702389]]),  # r1 start, r1 end, r1 exons
        (7705301, 7708251, [[7705301, 7705720], [7707957, 7708251]]),  # r2
        ("-GG", "-GG"),  # insertions info
        ("+", "+"),  # strand1, strand2
        ["INTERGENIC", "INTERGENIC"],
    )  # gene1, gene2

    return Event(event_tuple)  # type: ignore
