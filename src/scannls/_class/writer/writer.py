"""Writer for Fasta files and GTF files for Series object.

@Filename:    writer.py
@license:     MIT Licence
@Time:        12/30/21 4:02 PM
"""
from abc import ABC
from abc import abstractmethod
from contextlib import contextmanager
from pathlib import Path
from typing import Any
from typing import IO
from typing import List
from typing import Optional
from typing import Tuple

from ..basicClass import Series
from scannls import LoggerType


# todo: add asyncio support


class Writer(ABC):
    """Abstract class for writing object to file."""

    def __init__(self, file_path: str, logger: LoggerType):
        """Initialize Writer object."""
        self.logger = logger
        self.file_path = Path(file_path)
        if self.file_path.exists():
            self.logger.warning(f"{self.file_path} exists, will be overwritten.")
        self.io: Optional[IO] = None

    @abstractmethod
    def write_data(self, data_object: Any):
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @abstractmethod
    def write_line(self, line: str):
        """Write line to file.

        :param: line: Line to write to file.
        """

    @abstractmethod
    def open(self, mode: str = "w"):
        """Open file.

        :param: mode: Mode to open file.
        """

    @abstractmethod
    def close(self):
        """Close file."""


class Writers:
    """Writers."""

    def __init__(self, writers: Tuple["Writer", ...]):
        """Init writers."""
        self.writers_list = writers

    def write_series(self, series: Series) -> None:
        """Write series.

        .. note::
             This method need all writers to be opened.
        """
        for writer in self.writers_list:
            writer.write_data(series)

    def open_writers(self, mode: str = "w") -> List[IO]:
        """Open writers."""
        writers_list = []
        for writer in self.writers_list:
            writers_list.append(writer.open(mode))
        return writers_list

    def close_writers(self) -> None:
        """Close writers."""
        for writer in self.writers_list:
            writer.close()

    @contextmanager
    def open(self, mode: str = "w") -> Any:
        """Open writers contextmanager."""
        try:
            yield self.open_writers(mode)
        finally:
            self.close_writers()
