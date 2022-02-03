# !/usr/bin/env python
"""Test Insertion.

@Filename:    test_basicClass_novelinsertion.py
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        2/1/22 9:45 PM
"""
import pytest

from scannls import Insertion
from scannls import MicroHomology
from scannls import NovelInsertion


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


class TestInsertion:
    """Test Insertion."""

    def test_update_cigarstring_sms(
        self, insertion, sms=(5, 10, 2), source_s="left", source_strand="+"
    ):
        """Test update cirgarstring sms."""
        insertion.update_cigarstring_sms(
            sms=sms, source_s=source_s, source_strand=source_strand
        )
        assert insertion.sms == (1, 4, 12)
        assert insertion.cigarstring == "2S2M13S"

    def test_reverse_completement_query(self, insertion):
        """Test reverse completement."""
        insertion.reverse_completement_query()
        assert insertion.query_sequence == "AGCT"

    def test_reverse_strand(self, insertion):
        """Test reverse strand."""
        assert insertion.strand == "+"
        insertion.reverse_strand()
        assert insertion.strand == "-"


class TestNovelInsertion:
    """Test Novel Insertion."""

    def test_reverse_completement_query(self, novel_insertion):
        """Test reverse completement query."""
        novel_insertion.reverse_completement_query()
        assert novel_insertion.query_sequence == "TGAT"

    @pytest.mark.parametrize("num", [1, 4, 5])
    def test_increment_ao(self, novel_insertion, num):
        """Test increment ao."""
        assert novel_insertion.ao == 1
        novel_insertion.increment_ao(num)
        assert novel_insertion.ao == num + 1


class TestMicroHomology:
    """Test MicroHomology."""

    def test_reverse_completement_query(self, microhomology):
        """Test reverse completement query."""
        microhomology.reverse_completement_query()
        assert microhomology.query_sequence == "TGAT"
