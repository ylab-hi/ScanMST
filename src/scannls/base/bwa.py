"""Module for BLAT.

@Filename:    bwa.py
@Author:      YangyangLi
@license:     MIT Licence
@Time:        12/15/23 2:00 PM
"""
from __future__ import annotations

import secrets
from pathlib import Path

import psutil
from Bio.Sequencing.Applications import BwaIndexCommandline, BwaMemCommandline
from loguru import logger


class Aligner:
    MIN_MEMORY = 8

    def __init__(self, reference=Path) -> None:
        self.reference = Path(reference)
        if not self.reference.exists():
            msg = f"{self.reference} not found."
            raise FileNotFoundError(msg)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.reference=})"

    __str__ = __repr__

    def index_exist(self) -> bool:
        """Check if index exists."""
        return all(self.reference.with_suffix(ext).exists() for ext in ["amb", "ann", "bwt", "pac", "sa"])

    def build_index(self) -> None:
        index_cmd = BwaIndexCommandline(reference=self.reference)
        logger.trace(f"build index: {index_cmd}")
        index_cmd()

    @staticmethod
    def enough_memory() -> bool:
        """Check if there is enough memory to build index."""
        return psutil.virtual_memory().available >> 30 > Aligner.MIN_MEMORY

    def query(self, query: str, output: Path | None) -> None:
        if not self.index_exist():
            if not self.enough_memory():
                msg = "Not enough memory to build index."
                raise Exception(msg)

            # Build index
            self.build_index()

        output = self.mem(query, output)

    def mem(self, query: str, output: Path | None) -> Path:
        """Align query to reference."""
        ran_id = secrets.randbits(42)
        in_fastq = self.reference.parent / f"{ran_id}.fastq"
        with in_fastq.open("w") as fastq_file:
            fastq_file.write(f"@{ran_id}\n")
            fastq_file.write(f"{query}\n")
            fastq_file.write("+\n")
            fastq_file.write("I" * len(query) + "\n")

        if output is None:
            output = self.reference.parent / f"{ran_id}.sam"

        mem_cmd = BwaMemCommandline(reference=self.reference, read_file1=in_fastq)
        logger.trace(f"bwa mem: {mem_cmd}")
        mem_cmd(stdout=output.as_posix())
        in_fastq.unlink()

        return output

    @staticmethod
    def record_identity(record):
        """Calculate alignment identity for every record in sam file."""
        return (record.get_tag("NM") - record.get_tag("AS")) / record.query_alignment_length

    @staticmethod
    def alignment_identity(sam_file: Path | str):
        """Calculate alignment identity for every record in sam file."""
        if isinstance(sam_file, str):
            sam_file = Path(sam_file)

        import pysam

        with sam_file.open() as sam:
            for record in pysam.AlignmentFile(sam):
                print(record.query_alignment_length, record.query_length, record.reference_length)
                print(record.query_name, record.get_tag("NM"), record.get_tag("AS"), record.get_tag("XS"))
