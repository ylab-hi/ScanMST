"""Useful functions for scannls."""
import os
import secrets
import shutil
import subprocess
import time
from collections.abc import Callable
from collections.abc import Generator
from contextlib import contextmanager
from functools import wraps
from pathlib import Path
from typing import Any

import pysam
import rscannls
from loguru import logger
from rscannls import cigarstring2cigartuple  # type: ignore

from .base.exception import ToolNotFoundError
from .base.type import LoggerType
from .blat import load_fa2bit

__all__ = [
    "external_tool_checking",
    "get_softclip_length",
    "get_longest_insertion_sequence",
    "cigarstring2cigartuple",
    "timeit",
    "sleep",
    "find_2bit_file",
]


def external_tool_checking(software: list[str], log_handler: LoggerType) -> None:
    """Checking dependencies are installed."""
    for tool in software:
        output = shutil.which(tool)
        if not output:
            raise ToolNotFoundError(tool)
        log_handler.success(f"Checking for {tool} found ")


def find_2bit_file(
    fasta_path: str, log_handler: LoggerType, parameter: list[str] | None
) -> str:
    """Create 2bit file from fasta file.

     fa2bit usage:
      faToTwoBit in.fa [in2.fa in3.fa ...] out.2bit
     options:

    -long          use 64-bit offsets for index.   Allow for twoBit to contain more than 4Gb of sequence.
                   NOT COMPATIBLE WITH OLDER CODE.
    -noMask        Ignore lower-case masking in fa file.
    -stripVersion  Strip off version number after '.' for GenBank accessions.
    -ignoreDups    Convert first sequence only if there are duplicate sequence
                   names.  Use 'twoBitDup' to find duplicate sequences.
    """
    if parameter is None:
        parameter = []
    bit_file = Path(fasta_path).with_suffix(".2bit")
    if not bit_file.exists():
        log_handler.info(f"{bit_file.as_posix()} Not Found Creating...")
        subprocess.check_call(
            [load_fa2bit(), " ".join(parameter), fasta_path, bit_file.as_posix()]
        )
    return bit_file.as_posix()


def sleep(input_file: str, max_time: int = 30) -> None:
    """Sleep random time."""
    file_size = os.stat(input_file).st_size
    secrets.SystemRandom().seed(file_size)
    time.sleep(secrets.randbelow(max_time))


def get_softclip_length(
    read: pysam.libcalignedsegment.AlignedSegment,
    mode: int,
) -> tuple[int, str, int, int]:
    """Extract softclipped sequence information from input read.

    :param mode: read mode
    :param read: reads from pysam

    :return: length of soft-clipped part, sequence of soft-clipped part,
     the connection point of soft-clipped part (left/right),
     mode of soft-clipped part: 0:other; 2:left[SM]; 1:right[MS]
    """
    if read.query_sequence is None or read.cigarstring is None:
        raise ValueError(f"{read.query_name}'s query sequence or cigar is None")

    return rscannls.get_softclip_length(
        read.query_sequence, read.cigarstring, read.reference_start, mode
    )


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


@contextmanager
def change_dir(path: str) -> Generator:
    """Change current working directory.

    :param path: the path to be changed
    """
    old_dir = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(old_dir)


def change_dir_decorator(path: str):
    """Change current working directory.

    :param path: the path to be changed
    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(func)
        def wrapped(*args, **kwargs):  # type: ignore
            with change_dir(path):
                return func(*args, **kwargs)

        return wrapped

    return decorator


def get_longest_insertion_sequence(
    read: pysam.libcalignedsegment.AlignedSegment,
    insertion_length_cutoff: int = 50,
) -> tuple[int, str, int]:
    """Extract longest insertion sequences information from input read.

    :param read: reads from pysam
    :param insertion_length_cutoff: minimum insertion length considering for duplicated exon
    :return: the reference start position of insertion, the read start position of insertion, length of the insertion
    """
    if read.query_sequence is None or read.cigarstring is None:
        raise ValueError(f"{read.query_name}'s query sequence or cigar is None")

    return rscannls.get_longest_insertion_sequence(
        read.query_sequence,
        read.cigarstring,
        read.reference_start,
        insertion_length_cutoff,
    )
