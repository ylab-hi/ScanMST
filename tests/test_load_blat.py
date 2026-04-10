# !/usr/bin/env python
"""Test for blat executables.

@Filename:    test_load_blat.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        5/23/22 11:34 AM
"""
import subprocess

from scanmst import blat


def test_load_blat():
    """Test for blat path."""
    blat_path = blat.load_blat()
    assert blat_path.exists()


def test_load_gfserver():
    """Test for gfserver."""
    gfserver_path = blat.load_gfserver()
    assert gfserver_path.exists()
    subprocess.check_call(gfserver_path)


def test_load_gfclient():
    """Test for gfclient."""
    gfclient_path = blat.load_gfclient()
    assert gfclient_path.exists()
    subprocess.check_call(gfclient_path)


def test_load_fa2bit():
    """Test for fa2bit."""
    fa2bit_path = blat.load_fa2bit()
    assert fa2bit_path.exists()
    subprocess.check_call(fa2bit_path)
