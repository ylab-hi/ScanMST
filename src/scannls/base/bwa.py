"""Module for BLAT.

@Filename:    bwa.py
@Author:      YangyangLi
@license:     MIT Licence
@Time:        12/15/23 2:00 PM
"""
from __future__ import annotations

import secrets
from pathlib import Path
from subprocess import SubprocessError
from typing import TYPE_CHECKING

import psutil
import pysam
from Bio.Sequencing.Applications import BwaIndexCommandline, BwaMemCommandline
from loguru import logger

from .basic_class import Insertion, NovelInsertion

if TYPE_CHECKING:
    from collections.abc import Iterator


class Aligner:
    MIN_MEMORY = 8

    def __init__(
        self, reference=Path, min_mapq: int = 20, threshold_identity: float = 0.99
    ) -> None:
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
        return all(
            self.reference.with_suffix(".fa" + ext).exists()
            for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]
        )

    def build_index(self) -> None:
        index_cmd = BwaIndexCommandline(infile=self.reference)
        logger.trace(f"build index: {index_cmd}")
        index_cmd()

    @staticmethod
    def enough_memory() -> bool:
        """Check if there is enough memory to build index."""
        return psutil.virtual_memory().available >> 30 > Aligner.MIN_MEMORY

    def query(
        self, query: str, output: Path | None = None
    ) -> Iterator[pysam.AlignedSegment]:
        if not self.index_exist():
            if not self.enough_memory():
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
            if Aligner.record_identity(record, query_sequence_length)
            > threshold_identity
            and record.mapping_quality > min_mapq
        )

    def query_insertion(
        self,
        query: str,
        output: Path | None = None,
    ):
        records = self.query(query, output)
        keep_records = list(
            Aligner.filters(records, self.threshold_identity, self.min_mapq, len(query))
        )

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

    @staticmethod
    def obtain_variants_stats(cigartuples, md_string) -> tuple[int, int]:
        """Calculate alignment matched length and number of substitution in it."""
        deletion_len = 0
        substitution_len = 0
        match_len = 0
        for _operation, _len in cigartuples:
            if _operation == 2:
                deletion_len += _len
            elif _operation == 0:
                match_len += _len

        sum_of_subs_dels = 0
        for _letter in md_string:
            if ord(_letter) >= 65 and ord(_letter) <= 90:
                sum_of_subs_dels += 1

        substitution_len = sum_of_subs_dels - deletion_len
        return match_len, substitution_len

    @staticmethod
    def record_identity(record, query_sequence_length):
        """Calculate alignment identity for every record in sam file."""
        md_str = record.get_tag("MD") if record.has_tag("MD") else ""
        match_length, substitution_len = Aligner.obtain_variants_stats(
            record.cigartuples, md_str
        )
        identity = (match_length - substitution_len) / query_sequence_length
        logger.trace(
            f"record_identity: {identity} record mapq: {record.mapping_quality} record length: {record.query_length}"
        )
        return identity

    @staticmethod
    def alignment_records(sam_file: Path | str):
        """Calculate alignment identity for every record in sam file."""
        if isinstance(sam_file, str):
            sam_file = Path(sam_file)

        with sam_file.open() as sam:
            yield from pysam.AlignmentFile(sam)
