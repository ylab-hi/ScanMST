import subprocess
import time

from Bio.Seq import Seq
from loguru import logger

__funcs__ = {"reverse_complement", "external_tool_checking", "get_softclip_length"}

from scannls.exception import ToolNotFoundError


def reverse_complement(in_str):
    """
    obtain reverse complement sequence
    """
    my_dna = Seq(in_str)
    return str(my_dna.reverse_complement())


def external_tool_checking(logger) -> None:
    """checking dependencies are installed"""
    software = ["samtools", "gfClient", "gfServer"]
    for tool in software:
        output = subprocess.getoutput(tool)
        if "command not found" in output:
            raise ToolNotFoundError(tool)
        else:
            logger.success("Checking for '" + tool + "': found ")


def get_softclip_length(read) -> tuple:
    """Extract softclipped sequence information from input read
    :param read: reads from pysam
    :type read: pysam.libcalignedsegment.AlignedSegment
    :return: length of soft-clipped part, sequence of soft-clipped part, the connection point of soft-clipped part (left/right), mode of soft-clipped part: 0:other; 2:left[SM]; 1:right[MS]
    :rtype: tuple
    """
    if read.cigartuples[0][0] == 4:
        # there are soft-clipped segments in left and right both
        if read.cigartuples[-1][0] == 4:
            # length of left soft-clipped segment is bigger
            if read.cigartuples[0][1] > read.cigartuples[-1][1]:
                return (
                    read.cigartuples[0][1],
                    read.query_sequence[: read.cigartuples[0][1]],
                    read.reference_start,
                    2,
                )
            # length of right soft-clipped segment is bigger
            else:
                return (
                    read.cigartuples[-1][1],
                    read.query_sequence[read.query_length - read.cigartuples[-1][1] :],
                    read.reference_end - 1,
                    1,
                )
        # there are soft-clipped segments in left only
        else:
            return (
                read.cigartuples[0][1],
                read.query_sequence[: read.cigartuples[0][1]],
                read.reference_start,
                2,
            )
    # there are soft-clipped segments in right only
    elif read.cigartuples[-1][0] == 4:
        return (
            read.cigartuples[-1][1],
            read.query_sequence[read.query_length - read.cigartuples[-1][1] :],
            read.reference_end - 1,
            1,
        )
    else:
        return 0, "", -1, 0


def write_series_to_file(file_name, series):
    """
    write series to file
    """
    with open(file_name, "w") as f:
        for item in series:
            f.write(str(item) + "\n")


def timeit(func):
    def wrapped(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        logger.debug("Function '{}' executed in {:f} s", func.__name__, end - start)
        return result

    return wrapped
