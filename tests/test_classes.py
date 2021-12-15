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

from scannls._class.blat import Blat  # type: ignore


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
