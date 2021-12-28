# !/usr/bin/env python
"""Aligner class based on Gapmis.

@Filename:    aligner.py
@license:     MIT Licence
@Time:        12/27/21 6:01 PM
"""
import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import List
from typing import Tuple


@dataclass
class AlignerResult:
    """Class to store the result of aligner."""

    seq1: str
    seq2: str
    start1: int
    start2: int
    end1: int
    end2: int
    score: float
    n_gaps: float
    n_mismatches: float


class Aligner:
    """Aligner class for aligning two sequences by semi-global algorithm.

    :param seqa: the sequence to be aligned
    :param seqb:  another sequence to be aligned

    .. note::
        This class is wrapper of `gapmis` based on C implementation so that
        we can simply call it from Python.
    """

    def __init__(self, seqa: str, seqb: str):
        """Initialize Aligner class."""
        self.seqa = seqa
        self.seqb = seqb
        self.cmd = "gapmis -a {seqa} -b {seqb} -o {out}".format

    def __repr__(self):
        """Represent Aligner class."""
        return f"Aligner(seqa={self.seqa}, seqb={self.seqb})"

    def run(self) -> AlignerResult:
        """Run the aligner with temp file and return the result.

        .. note::
           1. Using gapmis to align two sequences
           2. Using `module::subprocess` to run the shell command of gapmis
           3. All files are stored in temp directory and will be deleted after aligning
           4. Return the result of aligning including the aligned sequence, start, end,
               score, number of gaps, number of mismatches for two sequences
        """
        with tempfile.TemporaryDirectory() as tmpdirname:
            tempfile_seq1 = os.path.join(tmpdirname, "seq1.fa")
            tempfile_seq2 = os.path.join(tmpdirname, "seq2.fa")
            with open(tempfile_seq1, "w") as f1, open(tempfile_seq2, "w") as f2:
                f1.write(f">seq1\n{self.seqa}\n")
                f2.write(f">seq2\n{self.seqb}\n")
            tempfile_name = os.path.join(tmpdirname, "tempfile.txt")
            subprocess.check_call(
                self.cmd(seqa=tempfile_seq1, seqb=tempfile_seq2, out=tempfile_name),
                shell=True,
            )
            align_result = self.parse_gapmis_result(tempfile_name)
        return align_result

    def parse_gapmis_result(self, result_file: str) -> AlignerResult:
        """Parse the result of gapmis.

        .. note::
            1. Return the result of aligning including the aligned sequence, start, end,
                score, number of gaps, number of mismatches for two sequences
            2. start and end coordinates are 0-based and end is not included, same as Python
        """
        seqa_coords: List[Tuple] = []
        seqb_coords: List[Tuple] = []
        with open(result_file) as f:
            for line in [line.strip() for line in f if not line.startswith("#")]:
                if line.startswith("seq1"):
                    # (1, 50)
                    seqa_coords.append(
                        (int(line.split()[1]) - 1, int(line.split()[-1]))
                    )
                elif line.startswith("seq2"):
                    seqb_coords.append(
                        (int(line.split()[1]) - 1, int(line.split()[-1]))
                    )
                elif line.startswith("Alignment"):
                    score = float(line.split()[-1])
                elif line.startswith("Number"):
                    mismatches = float(line.split()[-1])
                elif line.startswith("Length"):
                    gaps = float(line.split()[-1])
        return AlignerResult(
            seq1=self.seqa[seqa_coords[0][0] : seqa_coords[-1][-1]],
            seq2=self.seqb[seqb_coords[0][0] : seqb_coords[-1][-1]],
            start1=seqa_coords[0][0],
            start2=seqb_coords[0][0],
            end1=seqa_coords[-1][-1],
            end2=seqb_coords[-1][-1],
            score=score,
            n_gaps=gaps,
            n_mismatches=mismatches,
        )
