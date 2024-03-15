"""Module for writing data to a file.

@Filename:    __init__.py.py
@Author:      YangyangLi
@Time:        1/30/22 6:15 PM
"""

from .fasta_writer import FastaWriter
from .gtf_writer import GTFWriter
from .vcf_writer import VCFWriter
from .writer import Writers

__all__ = ["FastaWriter", "GTFWriter", "VCFWriter", "Writers"]
