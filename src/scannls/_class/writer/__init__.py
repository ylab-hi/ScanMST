# !/usr/bin/env python
"""Module for writing data to a file.

@Filename:    __init__.py.py
@Author:      YangyangLi
@license:     MIT Licence
@Time:        1/30/22 6:15 PM
"""
from .fastaWriter import FastaWriter
from .gtfWriter import GTFWriter
from .vcfWriter import VCFWriter
from .writer import Writers

__all__ = ["FastaWriter", "GTFWriter", "VCFWriter", "Writers"]
