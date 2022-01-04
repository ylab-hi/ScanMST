# !/usr/bin/env python
"""Useful functions for scannls."""
import subprocess
import time
from functools import wraps
from typing import Any
from typing import Callable
from typing import Tuple

import pysam.libcalignedsegment  # type: ignore
from loguru import logger
from loguru._logger import Logger

from ._class.basicClass import Read  # type: ignore
from ._class.exception import ToolNotFoundError  # type: ignore

__funcs__ = {"reverse_complement", "external_tool_checking", "get_softclip_length"}


def external_tool_checking(logger: Logger) -> None:
    """Checking dependencies are installed."""
    software = ["samtools", "gfClient", "gfServer", "gapmis"]
    for tool in software:
        output = subprocess.getoutput(tool)
        if "command not found" in output:
            raise ToolNotFoundError(tool)
        else:
            logger.success("Checking for '" + tool + "': found ")  # type: ignore


def get_softclip_length(
    read: pysam.libcalignedsegment.AlignedSegment, mode: int
) -> Tuple[int, str, int, int]:
    """Extract softclipped sequence information from input read.

    :param mode: read mode
    :param read: reads from pysam
    :return: length of soft-clipped part, sequence of soft-clipped part,
     the connection point of soft-clipped part (left/right),
     mode of soft-clipped part: 0:other; 2:left[SM]; 1:right[MS]
    """
    _cigar = read.cigarstring
    _mapq = read.mapping_quality
    _nm = read.get_tag("NM")
    _seq = read.query_sequence
    _strand = "-" if read.is_reverse else "+"
    _chrm = read.reference_name
    _pos = read.reference_start
    read_obj = Read.init(
        read.query_name, _chrm, _pos, _strand, _cigar, _mapq, _nm, _seq
    )

    if not mode:
        if read_obj.lt_soft_len > read_obj.rt_soft_len:
            return (
                read_obj.lt_soft_len,
                read_obj.query_sequence[: read_obj.lt_soft_len],
                read_obj.ref_start,
                2,
            )
        elif read_obj.lt_soft_len < read_obj.rt_soft_len:
            return (
                read_obj.rt_soft_len,
                read_obj.query_sequence[read_obj.query_length - read_obj.rt_soft_len :],
                read_obj.ref_end,
                1,
            )
        else:
            return 0, "", -1, 0
    else:
        if mode == 1:
            return (
                read_obj.rt_soft_len,
                read_obj.query_sequence[read_obj.query_length - read_obj.rt_soft_len :],
                read_obj.ref_end,
                1,
            )
        elif mode == 2:
            return (
                read_obj.lt_soft_len,
                read_obj.query_sequence[: read_obj.lt_soft_len],
                read_obj.ref_start,
                2,
            )
        else:
            return 0, "", -1, 0


def write_series_to_file(file_name: str, series: Any) -> None:
    """Write series to file."""
    with open(file_name, "w") as f:
        for item in series:
            f.write(str(item) + "\n")


def timeit(func: Callable[..., Any]) -> Callable[..., Any]:
    """Time the function execution.

    :param func: the function to be timed
    """

    @wraps(func)
    def wrapped(*args, **kwargs):  # type: ignore
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        logger.debug("Function '{}' executed in {:f} s", func.__name__, end - start)
        return result

    return wrapped
