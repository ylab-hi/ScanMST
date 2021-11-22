__version__ = "0.0.1"

try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)


from rich.traceback import install

install()
