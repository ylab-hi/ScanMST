#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ===========================================================
import argparse
import os
import re
import subprocess
import sys
import textwrap
import time
from collections import defaultdict

from Bio.Seq import Seq
from loguru import logger
from pyfaidx import Fasta
from pyfaidx import FastaNotFoundError

from .. import __version__
from ..classes import Blat
from ..classes import LengthAction
from ..classes import Read
from ..classes import ReadsConnecter
from ..classes import Series
from ..common import get_softclip_length

try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)

__funcs__ = {}
