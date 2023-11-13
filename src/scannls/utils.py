"""Useful functions for scannls."""
from __future__ import annotations

import os
import re
import secrets
import shutil
import subprocess
import time
from contextlib import contextmanager
from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING, Any

from loguru import logger

from scannls.base import MappingMode

from . import cppext
from .exception import ToolNotFoundError

if TYPE_CHECKING:
    from collections.abc import Callable, Generator

    import pysam

    from .type import LoggerType

__all__ = [
    "cigar_validity",
    "external_tool_checking",
    "get_softclip_length",
    "get_longest_insertion_sequence",
    "cigarstring2cigartuples",
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


def find_2bit_file(fasta_path: str, parameter: list[str] | None = None) -> str:
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
        logger.info(f"{bit_file.as_posix()} Not Found Creating...")
        subprocess.check_call(
            [load_fa2bit(), " ".join(parameter), fasta_path, bit_file.as_posix()],
        )
    return bit_file.as_posix()


def sleep(input_file: str, max_time: int = 30) -> None:
    """Sleep random time."""
    file_size = Path.stat(Path(input_file)).st_size
    secrets.SystemRandom().seed(file_size)
    time.sleep(secrets.randbelow(max_time))


def get_softclip_length(
    read: pysam.libcalignedsegment.AlignedSegment,
    mode: MappingMode,
) -> tuple[int, str, int, MappingMode] | None:
    """Extract softclipped sequence information from input read.

    :param mode: read mode
    :param read: reads from pysam
    :return: length of soft-clipped part, sequence of soft-clipped part,
     the connection point of soft-clipped part (left/right),
     mode of soft-clipped part: 0:other; 2:left[SM]; 1:right[MS]
    """
    if read.query_sequence is None or read.cigarstring is None:
        msg = f"{read.query_name}'s query sequence or cigar is None"
        raise ValueError(msg)

    parse_result = cppext.parseCigar(read.cigarstring)
    ref_end = read.reference_start + parse_result.ref_match

    if mode == MappingMode.Type0:
        if parse_result.lt_soft_len > parse_result.rt_soft_len:
            return (
                parse_result.lt_soft_len,
                read.query_sequence[: parse_result.lt_soft_len],
                read.reference_start,
                MappingMode.SM,
            )

        if parse_result.lt_soft_len < parse_result.rt_soft_len:
            return (
                parse_result.rt_soft_len,
                read.query_sequence[parse_result.query_len - parse_result.rt_soft_len :],
                ref_end,
                MappingMode.MS,
            )
        return None

    if mode == MappingMode.MS:
        return (
            parse_result.rt_soft_len,
            read.query_sequence[parse_result.query_len - parse_result.rt_soft_len :],
            ref_end,
            mode,
        )

    if mode == MappingMode.SM:
        return (
            parse_result.lt_soft_len,
            read.query_sequence[: parse_result.lt_soft_len],
            read.reference_start,
            mode,
        )

    return None


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
        msg = f"{read.query_name}'s query sequence or cigar is None"
        raise ValueError(msg)

    parse_result = cppext.parseCigar(read.cigarstring)
    ref_start = read.reference_start
    lt_soft_len = parse_result.lt_soft_len

    cigartuples_without_soft = parse_result.cigartuples_without_soft

    current_pos = ref_start
    current_len = lt_soft_len
    insertion_list = []
    for idx in range(0, len(cigartuples_without_soft), 2):
        op_code = cigartuples_without_soft[idx]
        _len = cigartuples_without_soft[idx + 1]

        if op_code == 0:  # M
            current_pos = current_pos + _len
            current_len = current_len + _len
        elif op_code in {2, 3}:  # D or N
            current_pos = current_pos + _len
        elif op_code == 1:  # I
            insertion_list.append((current_pos, current_len, _len))
            current_len = current_len + _len

    if len(insertion_list) == 0:
        return 0, "", 0

    # sorted by insertion length
    # if multiple insertions with the same size, choose the one with smallest reference position
    ins_ref_pos, ins_read_pos, ins_length = sorted(
        insertion_list,
        key=lambda x: x[2],
        reverse=True,
    )[0]
    ins_seq = read.query_sequence[ins_read_pos : (ins_read_pos + ins_length)]
    # update `ins_ref_pos` if insertion has adjacent N (100I500N)
    pattern = re.compile(re.escape(f"{ins_length}I") + r"(\d+)N")

    if ins_length >= insertion_length_cutoff:
        mat = pattern.search(read.cigarstring)
        if mat:
            ins_ref_pos += int(mat.group(1))
        return ins_ref_pos, ins_seq, ins_length
    return 0, "", 0


def cigarstring2cigartuples(cigarstring: str) -> list[tuple[int, int]]:
    """Convert cigarstring to cigartuples.

    :param cigarstring: cigarstring from reads
    :return: cigartuples is a list of (operation, length) tuples, such as [(0, 30), (1, 20), (4, 5)]
    """
    cigar_dict = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5}
    cigartuples = []
    _cigartuples = re.findall(r"(\d+)(\w)", cigarstring)
    for _, (length, operation) in enumerate(_cigartuples):
        op_code = cigar_dict[operation]
        _len = int(length)
        cigartuples.append((op_code, _len))
    return cigartuples


def cigar_validity(cigar_str: str) -> str:
    """Merge the first two OR last two 'same' operations in the CIGAR string generated by BLAT.

    :param cigar_str: BLAT generated cigarstring from 'softclipped_seq2SA_tag'
    :return: valid cigarstring

    :Example:

    >>> cigarstring =  '1S2S5M3S2S'
    >>> cigar_validity(cigarstring)
    '3S5M5S'
    """
    pattern = re.compile(r"((?P<length>\d+)(?P<op>\D))")
    items_list = pattern.findall(cigar_str)
    stack = [items_list[0]]
    for item in items_list[1:]:  # [('1S', '1', 'S'),..]
        last_item = stack[-1]
        if last_item[2] == item[2]:
            length = int(last_item[1]) + int(item[1])
            stack[-1] = (f"{length}{last_item[2]}", f"{length}", last_item[2])
        else:
            stack.append(item)

    return "".join(i[0] for i in stack)
