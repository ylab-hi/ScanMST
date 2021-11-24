import subprocess

from Bio.Seq import Seq

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


def get_softclip_length(read):
    """0 => M
    1 => I
    2 => D
    3 => N
    4 => S
    5 => H
    Name changes:
    reference_start  == pos
    reference_end    == aend
    query_length     == rlen
    reference_length == alen
    query_sequence   == seq
    return: length of soft-clipped part, sequence of soft-clipped part, the connection point of soft-clipped part (left/right), left/right soft-clipped part: 0:other; 2:left[SM]; 1:right[MS]
    """
    if read.cigartuples[0][0] == 4:
        # there are soft-clipped segments in left and right both
        if read.cigartuples[-1][0] == 4:
            # length of left soft-clipped segment is bigger
            if read.cigartuples[0][1] > read.cigartuples[-1][1]:
                return (
                    read.cigartuples[0][1],
                    read.query_sequence[: read.cigartuples[0][1]],
                    read.ref_start,
                    2,
                )
            # length of right soft-clipped segment is bigger
            else:
                return (
                    read.cigartuples[-1][1],
                    read.query_sequence[read.query_length - read.cigartuples[-1][1] :],
                    read.ref_end - 1,
                    1,
                )
        # there are soft-clipped segments in left only
        else:
            return (
                read.cigartuples[0][1],
                read.query_sequence[: read.cigartuples[0][1]],
                read.ref_start,
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
