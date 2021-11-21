import copy
import re
import sys
from collections import defaultdict

from align import aligner
from Bio.Seq import Seq

from . import __version__


try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)

__funcs__ = {"reverse_complement"}


def reverse_complement(in_str):
    """
    obtain reverse complement sequence
    """
    my_dna = Seq(in_str)
    return str(my_dna.reverse_complement())
