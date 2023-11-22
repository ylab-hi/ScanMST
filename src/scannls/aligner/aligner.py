from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

import logger
import psutil
import pysam


class Aligner(ABC):
    @staticmethod
    def enough_memory(min_memory=8) -> bool:
        """Check if there is enough memory to build index."""
        return psutil.virtual_memory().available >> 30 > min_memory

    @abstractmethod
    def index_exist(self) -> bool:
        """Check if index exists."""
        raise NotImplementedError

    @abstractmethod
    def build_index(self) -> None:
        """Build index."""
        raise NotImplementedError

    @abstractmethod
    def query(self, query: str, output: Path | None = None):
        """Query."""
        raise NotImplementedError

    @abstractmethod
    def filters(self, *args):
        """Filter."""
        raise NotImplementedError

    @abstractmethod
    def query_insertion(self, insertion: str, output: Path | None = None):
        """Query insertion."""
        raise NotImplementedError

    @staticmethod
    def record_identity(record, query_sequence_length):
        """Calculate alignment identity for every record in sam file."""
        md_str = record.get_tag("MD") if record.has_tag("MD") else ""
        match_length, substitution_len = Aligner.obtain_variants_stats(record.cigartuples, md_str)
        identity = (match_length - substitution_len) / query_sequence_length
        logger.trace(f"record_identity: {identity} record mapq: {record.mapping_quality} record length: {record.query_length}")
        return identity

    @staticmethod
    def alignment_records(sam_file: Path | str):
        """Calculate alignment identity for every record in sam file."""
        if isinstance(sam_file, str):
            sam_file = Path(sam_file)

        with sam_file.open() as sam:
            yield from pysam.AlignmentFile(sam)

    @staticmethod
    def obtain_variants_stats(cigartuples, md_string) -> tuple[int, int]:
        """Calculate alignment matched length and number of substitution in it."""
        deletion_len = 0
        substitution_len = 0
        match_len = 0
        # unmapped reads does not have CIGAR
        if cigartuples:
            for operation, _len in cigartuples:
                if operation == 2:
                    deletion_len += _len
                elif operation == 0:
                    match_len += _len

        sum_of_subs_dels = 0
        if md_string:
            for letter in md_string:
                if ord(letter) >= 65 and ord(letter) <= 90:
                    sum_of_subs_dels += 1

        substitution_len = sum_of_subs_dels - deletion_len
        return match_len, substitution_len
