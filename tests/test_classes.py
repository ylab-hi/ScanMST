#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@version: 0.0.1
@license: MIT Licence
@file: test_classes.py
@time: 16/11/2021 16:03
"""
import os
from pathlib import Path

import loguru
import pytest
from loguru import logger

from scannls.classes import Blat
from scannls.classes import Node
from scannls.classes import NovelInsertion
from scannls.classes import Series
from scannls.classes import SpliceGraph


class TestBlat:
    @pytest.fixture(scope="class")
    def blat(self):
        return Blat(ref_2bit=".", logger=logger, port=88888, output_dir=".")

    @pytest.fixture(scope="class")
    def process(self):
        names = ["gfServer", "test"]

        class _Process:
            def __init__(self, name):
                self._name = name

            def name(self):
                return self._name

            def cmdline(self):
                return True

        return [_Process(name) for name in names]

    def test_ref_dir(self, blat):
        ref_dir = Path(blat.ref_dir)

        assert ref_dir.is_absolute()

    def test_log_file(self, blat):
        log_file = Path(blat.log_file)
        assert log_file.is_absolute()

    def test_is_ready(self, blat, mocker):
        spy = mocker.spy(loguru.logger, "debug")
        assert blat.is_ready() is False

        with open(blat.log_file, "a") as log:
            log.write("Server ready")

        assert blat.is_ready() is True
        os.remove(blat.log_file)

        assert spy.call_count == 2

    def test_is_running(self, blat, process, mocker):
        process_1, process_2 = process

        mocker.patch("psutil.process_iter", return_value=[process_1])
        assert blat.is_running() is True

        mocker.patch("psutil.process_iter", return_value=[process_2])
        assert blat.is_running() is False


class TestSpliceGraph:
    @pytest.fixture(scope="class")
    def series_list(self):
        series1 = Series(blat=None, logger=None)
        series1.nodes = [
            Node(
                prev_bp=None,
                next_bp="chr17:7702552",
                strand="+",
                chrom="chr17",
                ref_start=7701656,
                ref_end=7702552,
                exons=[[7701656, 7702552]],
                sv_type="TRA",
            ),
            Node(
                prev_bp="chr1:15872815",
                next_bp=None,
                strand="+",
                chrom="chr1",
                ref_start=15872815,
                ref_end=15873800,
                exons=[[15872815, 15873800]],
                sv_type=None,
            ),
        ]
        series2 = Series(blat=None, logger=None)
        series2.nodes = [
            Node(
                prev_bp=None,
                next_bp="chr17:7702552",
                strand="+",
                chrom="chr17",
                ref_start=7701500,
                ref_end=7702552,
                exons=[[7701500, 7702552]],
                sv_type="TRA",
            ),
            Node(
                prev_bp="chr1:15872815",
                next_bp=None,
                strand="+",
                chrom="chr1",
                ref_start=15872815,
                ref_end=15876678,
                exons=[[15872815, 15876678]],
                sv_type=None,
            ),
        ]
        series3 = Series(blat=None, logger=None)
        series3.nodes = [
            Node(
                prev_bp=None,
                next_bp="chr1:15873900",
                strand="+",
                chrom="chr1",
                ref_start=15872890,
                ref_end=15873900,
                exons=[[15872890, 15873900]],
                sv_type="TRA",
                insertion_info=(
                    False,
                    NovelInsertion(hit_num=1, query_sequence="ATCGATCG"),
                ),
            ),
            Node(
                prev_bp="chr17:872815",
                next_bp=None,
                strand="+",
                chrom="chr17",
                ref_start=872815,
                ref_end=876678,
                exons=[[872815, 876678]],
                sv_type=None,
            ),
        ]
        series4 = Series(blat=None, logger=None)
        series4.nodes = [
            Node(
                prev_bp=None,
                next_bp="chr1:15873900",
                strand="+",
                chrom="chr1",
                ref_start=15872890,
                ref_end=15873900,
                exons=[[15872890, 15873900]],
                sv_type="TRA",
                insertion_info=(
                    False,
                    NovelInsertion(hit_num=1, query_sequence="ATCGATCG"),
                ),
            ),
            Node(
                prev_bp="chr17:872815",
                next_bp=None,
                strand="+",
                chrom="chr17",
                ref_start=872815,
                ref_end=876678,
                exons=[[872815, 873400], [875500, 876678]],
                sv_type=None,
            ),
        ]
        series5 = Series(blat=None, logger=None)
        series5.nodes = [
            Node(
                prev_bp=None,
                next_bp="chr17:7702552",
                strand="+",
                chrom="chr17",
                ref_start=7701656,
                ref_end=7702552,
                exons=[[7701656, 7702552]],
                sv_type="TRA",
            ),
            Node(
                prev_bp="chr1:15872815",
                next_bp="chr1:15873900",
                strand="+",
                chrom="chr1",
                ref_start=15872815,
                ref_end=15873900,
                exons=[[15872815, 15873900]],
                sv_type="TRA",
                insertion_info=(
                    False,
                    NovelInsertion(hit_num=1, query_sequence="ATCGATCG"),
                ),
            ),
            Node(
                prev_bp="chr17:872815",
                next_bp=None,
                strand="+",
                chrom="chr17",
                ref_start=872815,
                ref_end=876678,
                exons=[[872815, 876678]],
                sv_type=None,
            ),
        ]
        series6 = Series(blat=None, logger=None)
        series6.nodes = [
            Node(
                prev_bp=None,
                next_bp="chr17:7702552",
                strand="+",
                chrom="chr17",
                ref_start=7701656,
                ref_end=7702552,
                exons=[[7701656, 7702552]],
                sv_type="TRA",
            ),
            Node(
                prev_bp="chr1:15872815",
                next_bp="chr1:15873900",
                strand="+",
                chrom="chr1",
                ref_start=15872815,
                ref_end=15873900,
                exons=[[15872815, 15873900]],
                sv_type="TRA",
                insertion_info=None,
            ),
            Node(
                prev_bp="chr17:872815",
                next_bp=None,
                strand="+",
                chrom="chr17",
                ref_start=872815,
                ref_end=876678,
                exons=[[872815, 876678]],
                sv_type=None,
            ),
        ]
        series7 = Series(blat=None, logger=None)
        series7.nodes = [
            Node(
                prev_bp=None,
                next_bp="chr17:7702552",
                strand="+",
                chrom="chr17",
                ref_start=7701656,
                ref_end=7702552,
                exons=[[7701656, 7702552]],
                sv_type="TRA",
            ),
            Node(
                prev_bp="chr1:15872815",
                next_bp="chr1:15873900",
                strand="+",
                chrom="chr1",
                ref_start=15872815,
                ref_end=15873900,
                exons=[[15872815, 15873900]],
                sv_type="TRA",
                insertion_info=None,
            ),
            Node(
                prev_bp="chr17:872815",
                next_bp=None,
                strand="+",
                chrom="chr17",
                ref_start=872815,
                ref_end=876678,
                exons=[[872815, 873400], [875500, 876678]],
                sv_type=None,
            ),
        ]

        return [series1, series2, series3, series4, series5, series6, series7]

    def test_run(self, series_list):
        from loguru import logger

        splice_graph = SpliceGraph(series_list, logger=logger)
        splice_graph.run()
