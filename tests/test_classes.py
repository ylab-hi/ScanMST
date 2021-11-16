#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@version: 0.0.1
@license: MIT Licence
@file: test_classes.py
@time: 16/11/2021 16:03
"""
import os
import tempfile
from pathlib import Path

import loguru
import pytest
from loguru import logger

from scannls.classes import Blat


class TestBlat:
    @pytest.fixture(scope="class")
    def blat(self):
        return Blat(ref_2bit=".", logger=logger, port=88888, output_dir=".")

    def test_ref_dir(self, blat):
        ref_dir = Path(blat.ref_dir)

        assert ref_dir.is_absolute()

    def test_log_file(self, blat):
        log_file = Path(blat.log_file)
        assert log_file.is_absolute()

    def test_is_ready(self, blat, mocker):
        # mocker.patch("loguru.logger.debug")
        spy = mocker.spy(loguru.logger, "debug")
        assert blat.is_ready() is False

        with open(blat.log_file, "a") as log:
            log.write("Server ready")

        assert blat.is_ready() is True
        os.remove(blat.log_file)

        assert spy.call_count == 2

    def test_is_running(self, blat, mocker):

        assert False

    def test_start_server(self):
        assert False

    def test_stop_server(self):
        assert False

    def test__query(self):
        assert False

    def test_query(self):
        assert False
