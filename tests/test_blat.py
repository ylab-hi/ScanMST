#!/usr/bin/env python3
"""Test classes.

@version: 0.0.1
@license: MIT Licence
@file: test_classes.py
@time: 16/11/2021 16:03
"""
from pathlib import Path

import pytest
from loguru import logger

from scannls import Blat


@pytest.fixture(scope="module")
def blat() -> Blat:
    """Create Blat instance."""
    return Blat(ref_2bit=".", logger=logger, port=88888, output_dir=".")


@pytest.fixture()
def process():
    """Create fake process."""
    names = ["gfServer", "test"]

    class _Process:
        def __init__(self, name: str, status: str):
            self._status = status
            self._name = name

        def name(self) -> str:
            return self._name

        def cmdline(self) -> bool:
            return True

        def status(self) -> str:
            return self._status

    return [_Process(name, "running") for name in names]


@pytest.mark.usefixtures("blat")
class TestBlat:
    """Test Blat class."""

    def test_ref_dir(self, blat):
        """Test ref_dir."""
        ref_dir = Path(blat.ref_dir)

        assert ref_dir.is_absolute()

    def test_log_file(self, blat):
        """Test log file."""
        log_file = Path(blat.log_file_path)
        assert log_file.is_absolute()

    def test_is_ready_when_self_open_server(self, blat):
        """Test is_ready."""
        assert blat.is_ready() is False

        with open(blat.log_file_path, "a") as log:
            log.write("Server ready")
        blat.is_start_server = True

        assert blat.is_ready() is True
        Path(blat.log_file_path).unlink()

    def test_is_running(self, blat, process, mocker):
        """Test is running."""
        process_1, process_2 = process

        mocker.patch("psutil.process_iter", return_value=[process_1])
        assert blat.is_running() is True

        mocker.patch("psutil.process_iter", return_value=[process_2])
        assert blat.is_running() is False
