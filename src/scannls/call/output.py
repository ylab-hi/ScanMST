#!/usr/bin/env python3


try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)

__funcs__ = {}
