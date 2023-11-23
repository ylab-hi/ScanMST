from __future__ import annotations

import secrets
from pathlib import Path

from loguru import logger

from .aligner import Aligner


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

    def _query(self, query: str, output: Path):
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

    def query(self, query: str, output: Path | None = None) -> Iterator[pysam.AlignedSegment]:
        pass

    __str__ = __repr__
