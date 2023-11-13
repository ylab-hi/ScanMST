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
from Bio import SearchIO
from Bio.Sequencing.Applications import BwaIndexCommandline, BwaMemCommandline
from loguru import logger


class Aligner:
    MIN_MEMORY = 8

    def __init__(self, reference=Path) -> None:
        if not reference.exists():
            msg = f"{self.reference} not found."
            raise FileNotFoundError(msg)

        self.reference = reference

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.reference=})"

    __str__ = __repr__

    def index_exist(self) -> bool:
        """Check if index exists."""
        return all(self.reference.with_suffix(ext).exists() for ext in ["amb", "ann", "bwt", "pac", "sa"])

    def build_index(self) -> None:
        index_cmd = BwaIndexCommandline(reference=self.reference)
        index_cmd()
        logger.trace(f"build index: {index_cmd}")

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
        in_fasta = self.reference.parent / f"{ran_id}.fasta"
        with in_fasta.open("w") as fasta_file:
            fasta_file.write(f">{ran_id}\n")
            fasta_file.write(f"{query}\n")

        if output is None:
            output = self.reference.parent / f"{ran_id}.sam"

        mem_cmd = BwaMemCommandline(reference=self.reference, query=in_fasta)
        mem_cmd(stdout=output.as_posix())
        logger.trace(f"bwa mem: {mem_cmd}")

        return output

    @staticmethod
    def alignment_identity(sam_file: Path) -> float:
        """Calculate alignment identity for every record in sam file."""
        with sam_file.open() as sam:
            sam_records = SearchIO.parse(sam, "sam")
            for sam_record in sam_records:
                for hit in sam_record.hits:
                    for hsp in hit.hsps:
                        return hsp.ident_pct
            return None
