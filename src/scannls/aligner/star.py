from __future__ import annotations

import secrets
import shlex
import subprocess
from pathlib import Path
from subprocess import SubprocessError
from typing import TYPE_CHECKING

from loguru import logger

from .aligner import Aligner

if TYPE_CHECKING:
    from collections.abc import Iterator

    import pysam


class Star(Aligner):
    MIN_MEMORY = 20

    INDEX_TEMPLATE = "STAR --runThreadN {thread}  --runMode genomeGenerate --genomeDir {genomeDir} --genomeFastaFiles {genomeFastaFiles} --sjdbGTFfile {GTFfile}  ".format
    MAPPING_TEMPLATE = "STAR --runThreadN {thread} --genomeDir {genomeDir} --readFilesIn {readFilesIn} --outFileNamePrefix {outFileNamePrefix} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI AS NM MD".format

    def __init__(self, reference=Path, gtf=Path, threads: int = 2, min_mapq: int = 20, threshold_identity: float = 0.99) -> None:
        super().__init__()
        self.reference = reference
        self.gtf = gtf
        self.min_mapq = min_mapq
        self.threshold_identity = threshold_identity
        self.threads = threads

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.reference=})"

    def index_exist(self):
        pass

    def build_index(self):
        pass

    def _query(self, query: str) -> Path:
        ran_id = secrets.token_hex(8)
        in_fastq = self.reference.parent / f"{ran_id}.fq"
        with in_fastq.open("w") as fastq_file:
            fastq_file.write(f"@{ran_id}\n")
            fastq_file.write(f"{query}\n")
            fastq_file.write("+\n")
            fastq_file.write("I" * len(query) + "\n")

        output = self.reference.parent / f"{ran_id}.bam"

        cmd = self.MAPPING_TEMPLATE(thread=self.threads, genomeDir=self.reference, readFilesIn=in_fastq, outFileNamePrefix=output)

        logger.trace(f"alinger cmd: {cmd}")
        result = subprocess.check_output(shlex.split(cmd))
        with output.open("wb") as output_file:
            output_file.write(result)
        in_fastq.unlink()
        return output

    def query(self, query: str) -> Iterator[pysam.AlignedSegment]:
        if not self.index_exist():
            if not self.enough_memory(self.MIN_MEMORY):
                msg = "Not enough memory to build index."
                raise Exception(msg)

            try:
                # Build index
                self.build_index()
            except SubprocessError as e:
                msg = f"Failed to build index: {e}"
                raise Exception(msg) from e

        output = self._query(query)
        yield from self.alignment_records(output)
        output.unlink()

    __str__ = __repr__
