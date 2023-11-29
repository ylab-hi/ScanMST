from __future__ import annotations

import secrets
import shlex
import subprocess
from typing import TYPE_CHECKING

from loguru import logger

from scannls.base import Insertion, NovelInsertion

from .aligner import Aligner

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    import pysam


class STAR(Aligner):
    MIN_MEMORY = 20

    INDEX_TEMPLATE = (
        "STAR --runThreadN {thread}  --runMode genomeGenerate --genomeDir {genomeDir} "
        "--genomeFastaFiles {genomeFastaFiles} --sjdbGTFfile {GTFfile}  "
    ).format

    MAPPING_TEMPLATE = (
        "STAR --runThreadN {thread} --genomeDir {genomeDir} --readFilesIn {readFilesIn}"
        "--outFileNamePrefix {outFileNamePrefix} --outSAMtype BAM SortedByCoordinate "
        " --outSAMunmapped Within --outSAMattributes NH HI AS NM MD"
    ).format

    def __init__(
        self, reference: Path, index: Path, threads: int = 1, min_mapq: int = 20, threshold_identity: float = 0.99
    ) -> None:
        super().__init__()
        self.reference = reference
        self.index = index
        self.min_mapq = min_mapq
        self.threshold_identity = threshold_identity
        self.threads = threads

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.reference=})"

    __str__ = __repr__

    def index_exist(self):
        return self.index.exists()

    def build_index(self):
        raise NotImplementedError

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
            url = "https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf"
            msg = f"Index does not exist. Please reference {url} to build index first"
            raise Exception(msg)

        output = self._query(query)
        yield from self.alignment_records(output)
        output.unlink()

    def query_insertion(
        self,
        query: str,
    ):
        records = self.query(query)
        keep_records = list(self.filters(records, self.threshold_identity, self.min_mapq, len(query)))

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
