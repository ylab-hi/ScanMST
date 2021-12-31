# !/usr/bin/env python
"""Writer for Fasta files and GTF files for Series object.

@Filename:    writer.py
@license:     MIT Licence
@Time:        12/30/21 4:02 PM
"""
from abc import ABC
from abc import abstractmethod
from functools import singledispatch
from pathlib import Path
from typing import Any

from basicClass import Series
from loguru._logger import Logger


class Writer(ABC):
    def __init__(self, file_path: Path, logger: Logger):
        self.logger = logger
        self.file_path = Path(file_path)
        if self.file_path.exists():
            self.logger.warning(f"{self.file_path} exists, will be overwritten.")

    @abstractmethod
    def write_data(self, data_object: Any):
        """Write data to file.

        Args:
            data_object: Data to write to file.
        """

    @abstractmethod
    def write_line(self, line: str):
        """Write line to file.

        Args:
            line: Line to write to file.
        """


class FastaWriter(Writer):
    def __init__(self, file_path: Path, logger: Logger):
        super().__init__(file_path, logger)

    def write_data(self, data_object: Any):
        pass

    def write_line(self, line: str):
        pass


class GTFWriter(Writer):
    def __init__(self, file_path: Path, logger: Logger):
        super().__init__(file_path, logger)

    def write_data(self, data_object: Any):
        pass

    def write_line(self, line: str):
        pass
