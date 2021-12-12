__version__ = "0.0.1"

try:
    import pysam  # type: ignore
    import numpy as np
    import HTSeq  # type: ignore
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)
#
# from rich.traceback import install
#
# install()
