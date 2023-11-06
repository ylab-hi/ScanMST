"""Writer for Fasta files and GTF files for Series object.

@Filename:    writer.py
@Author:      YangyangLi
@Time:        12/30/21 4:02 PM
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from contextlib import contextmanager
from pathlib import Path
from typing import IO, TYPE_CHECKING, Any

from loguru import logger

if TYPE_CHECKING:
    from scannls.graph import NLPath


class Writer(ABC):
    """Abstract class for writing object to file."""

    def __init__(self, file_path: str) -> None:
        """Initialize Writer object."""
        self.file_path = Path(file_path)
        if self.file_path.exists():
            logger.warning(f"{self.file_path} exists, will be overwritten.")
        self.io: IO | None = None

    @abstractmethod
    def write_data(self, data_object: Any, object_id: str):
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

    def __init__(self, writers: tuple[Writer, ...]) -> None:
        """Init writers."""
        self.writers_list = writers

    def write_series(self, series: NLPath, clique_id: str) -> None:
        """Write series.

        .. note::
             This method need all writers to be opened.
        """
        for writer in self.writers_list:
            writer.write_data(series, clique_id)

    def open_writers(self, mode: str = "w") -> list[IO]:
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
