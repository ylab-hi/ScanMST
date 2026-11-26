# !/usr/bin/env python
"""Test for blat executables.

@Filename:    test_load_blat.py
"""

import subprocess

import pytest

from scanmst import blat

pytestmark = pytest.mark.skipif(
    not blat.load_blat().exists(),
    reason="BLAT executables are downloaded from UCSC on demand and are not installed",
)


def test_load_blat():
    """Test for blat path."""
    blat_path = blat.load_blat()
    assert blat_path.exists()


def test_load_gfserver():
    """Test for gfserver."""
    gfserver_path = blat.load_gfserver()
    assert gfserver_path.exists()
    res = subprocess.run([gfserver_path], capture_output=True, text=True, check=False)
    assert res.returncode in {0, 255}
    assert "gfserver" in (res.stdout + res.stderr).lower()


def test_load_gfclient():
    """Test for gfclient."""
    gfclient_path = blat.load_gfclient()
    assert gfclient_path.exists()
    res = subprocess.run([gfclient_path], capture_output=True, text=True, check=False)
    assert res.returncode in {0, 255}
    assert "gfclient" in (res.stdout + res.stderr).lower()


def test_load_fa2bit():
    """Test for fa2bit."""
    fa2bit_path = blat.load_fa2bit()
    assert fa2bit_path.exists()
    res = subprocess.run([fa2bit_path], capture_output=True, text=True, check=False)
    assert res.returncode in {0, 255}
    assert "fatotwobit" in (res.stdout + res.stderr).lower()
