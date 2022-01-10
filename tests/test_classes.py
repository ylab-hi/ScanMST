#!/usr/bin/env python3
"""Test classes.

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

from scannls import Blat  # type: ignore


class TestBlat:
    """Test Blat class."""

    @pytest.fixture(scope="class")
    def blat(self) -> Blat:
        """Create Blat instance."""
        return Blat(ref_2bit=".", logger=logger, port=88888, output_dir=".")

    @pytest.fixture(scope="class")
    def process(self):
        """Create fake process."""
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
        """Test ref_dir."""
        ref_dir = Path(blat.ref_dir)

        assert ref_dir.is_absolute()

    def test_log_file(self, blat):
        """Test log file."""
        log_file = Path(blat.log_file_path)
        assert log_file.is_absolute()

    def test_is_ready(self, blat, mocker):
        """Test is_ready."""
        spy = mocker.spy(loguru.logger, "debug")
        assert blat.is_ready() is False

        with open(blat.log_file_path, "a") as log:
            log.write("Server ready")

        assert blat.is_ready() is True
        os.remove(blat.log_file_path)

        assert spy.call_count == 2

    def test_is_running(self, blat, process, mocker):
        """Test is running."""
        process_1, process_2 = process

        mocker.patch("psutil.process_iter", return_value=[process_1])
        assert blat.is_running() is True

        mocker.patch("psutil.process_iter", return_value=[process_2])
        assert blat.is_running() is False
