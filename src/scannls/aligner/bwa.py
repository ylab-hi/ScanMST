"""Module for Bwa.

@Filename:    bwa.py
@Author:      YangyangLi
@Time:        12/15/23 2:00 PM
"""

from __future__ import annotations

import secrets
from pathlib import Path
from subprocess import SubprocessError
from typing import TYPE_CHECKING

from Bio.Sequencing.Applications import BwaIndexCommandline, BwaMemCommandline
from loguru import logger

from scannls import Insertion, NovelInsertion

from .aligner import Aligner

if TYPE_CHECKING:
    from collections.abc import Iterator

    import pysam


class Bwa(Aligner):
    MIN_MEMORY = 8

    def __init__(self, reference=Path, min_mapq: int = 20, threshold_identity: float = 0.99) -> None:
        self.reference = Path(reference)
        self.min_mapq = min_mapq
        self.threshold_identity = threshold_identity
        if not self.reference.exists():
            msg = f"{self.reference} not found."
            raise FileNotFoundError(msg)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.reference=})"

    __str__ = __repr__

    def index_exist(self) -> bool:
        """Check if index exists."""
        return all(self.reference.with_suffix(".fa" + ext).exists() for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"])

    def build_index(self) -> None:
        index_cmd = BwaIndexCommandline(infile=self.reference)
        logger.trace(f"build index: {index_cmd}")
        index_cmd()

    def query(self, query: str, output: Path | None = None) -> Iterator[pysam.AlignedSegment]:
        if not self.index_exist():
            if not self.enough_memory(Bwa.MIN_MEMORY):
                msg = "Not enough memory to build index."
                raise Exception(msg)

            try:
                # Build index
                self.build_index()
            except SubprocessError as e:
                msg = f"Failed to build index: {e}"
                raise Exception(msg) from e

        output = self.mem(query, output)
        yield from self.alignment_records(output)
        output.unlink()

    @staticmethod
    def filters(
        records: Iterator[pysam.AlignedSegment],
        threshold_identity,
        min_mapq,
        query_sequence_length,
    ) -> Iterator[pysam.AlignedSegment]:
        return (
            record
            for record in records
            if Bwa.record_identity(record, query_sequence_length) > threshold_identity and record.mapping_quality > min_mapq
        )

    def query_insertion(
        self,
        query: str,
        output: Path | None = None,
    ):
        records = self.query(query, output)
        keep_records = list(Bwa.filters(records, self.threshold_identity, self.min_mapq, len(query)))

        logger.trace(f"alginer: keep_records: {len(keep_records)}")

        if not keep_records:
            return False, NovelInsertion(hit_num=0, query_sequence=query)

        if len(keep_records) == 1:
            top_record = keep_records[0]
            return True, Insertion(
                hit_num=1,
                chrom=top_record.reference_name,
                ref_start=top_record.reference_start,
                strand="+" if top_record.is_reverse else "-",
                cigarstring=top_record.cigarstring,
                mapq=top_record.mapping_quality,
                nm=top_record.get_tag("NM") if top_record.has_tag("NM") else 0,
                query_sequence=query,
                query_qualities=top_record.query_qualities,
            )

        return False, NovelInsertion(
            hit_num=len(keep_records),
            query_sequence=query,
        )

    def mem(self, query: str, output: Path | None) -> Path:
        """Align query to reference."""
        ran_id = secrets.token_hex(8)
        in_fastq = self.reference.parent / f"{ran_id}.fq"
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
