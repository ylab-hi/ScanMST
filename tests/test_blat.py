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
        def __init__(self, name):
            self._name = name

        def name(self) -> str:
            return self._name

        def cmdline(self) -> bool:
            return True

    return [_Process(name) for name in names]


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

    def test_is_ready_when_self_open_server(self, blat, mocker):
        """Test is_ready."""
        blat.stop_server()
        assert blat.is_ready() is False

        with open(blat.log_file_path, "a") as log:
            log.write("Server ready")
        blat.is_start_server = True

        assert blat.is_ready() is True
        blat.stop_server()

    def test_is_ready_when_others_open_server(self, blat):
        """Test is_ready when others open server."""
        blat.is_start_server = False
        blat.set_env(True)
        assert blat.is_ready() is True
        blat.stop_server()

    def test_is_running(self, blat, process, mocker):
        """Test is running."""
        process_1, process_2 = process

        mocker.patch("psutil.process_iter", return_value=[process_1])
        assert blat.is_running() is True

        mocker.patch("psutil.process_iter", return_value=[process_2])
        assert blat.is_running() is False
