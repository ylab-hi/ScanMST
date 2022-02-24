# !/usr/bin/env python
"""Useful functions for scannls."""
import os
import secrets
import shutil
import time
from functools import wraps
from typing import Any
from typing import Callable
from typing import Tuple

import pysam  # type: ignore
from loguru import logger

from . import ToolNotFoundError
from ._class.type import LoggerType
from scannls import cppext

__all__ = ["external_tool_checking", "get_softclip_length", "timeit", "sleep"]


def external_tool_checking(log_handler: LoggerType) -> None:
    """Checking dependencies are installed."""
    software = ["gfClient", "gfServer"]
    for tool in software:
        output = shutil.which(tool)
        if not output:
            raise ToolNotFoundError(tool)
        log_handler.success(f"Checking for {tool} found ")


def sleep(input_file: str, max_time: int = 30) -> None:
    """Sleep random time."""
    file_size = os.stat(input_file).st_size
    secrets.SystemRandom().seed(file_size)
    time.sleep(secrets.randbelow(max_time))


def get_softclip_length(
    read: pysam.libcalignedsegment.AlignedSegment,
    mode: int,
) -> Tuple[int, str, int, int]:
    """Extract softclipped sequence information from input read.

    :param mode: read mode
    :param read: reads from pysam
    :return: length of soft-clipped part, sequence of soft-clipped part,
     the connection point of soft-clipped part (left/right),
     mode of soft-clipped part: 0:other; 2:left[SM]; 1:right[MS]
    """
    parse_result = cppext.parseCigar(read.cigarstring)
    ref_end = read.reference_start + parse_result.ref_match

    if mode == 0:
        if parse_result.lt_soft_len > parse_result.rt_soft_len:
            return (
                parse_result.lt_soft_len,
                read.query_sequence[: parse_result.lt_soft_len],
                read.reference_start,
                2,
            )
        elif parse_result.lt_soft_len < parse_result.rt_soft_len:
            return (
                parse_result.rt_soft_len,
                read.query_sequence[
                    parse_result.query_len - parse_result.rt_soft_len :
                ],
                ref_end,
                1,
            )
        else:
            return 0, "", -1, 0

    if mode == 1:
        return (
            parse_result.rt_soft_len,
            read.query_sequence[parse_result.query_len - parse_result.rt_soft_len :],
            ref_end,
            1,
        )
    elif mode == 2:
        return (
            parse_result.lt_soft_len,
            read.query_sequence[: parse_result.lt_soft_len],
            read.reference_start,
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
