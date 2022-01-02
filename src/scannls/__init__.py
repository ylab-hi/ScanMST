# !/usr/bin/env python
"""Init file for scannls package."""
__version__ = "0.0.1"

try:
    import pysam  # type: ignore
    import numpy as np
    import HTSeq  # type: ignore
except ModuleNotFoundError as e:
    raise SystemExit from e
