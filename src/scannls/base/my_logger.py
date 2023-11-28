# !/usr/bin/env python
"""Wrapper for loguru.logger.

@Filename:    mylogger.py
@Author:      YangyangLi
@Time:        12/15/21 2:08 PM
"""
from loguru._logger import Logger

from scannls.mtype import LoggerType


class MyLogger(LoggerType):
    """Wrapper for logger in order to use in multiprocessing.

    to show contig name in logging information before message.
    However, there is no way to get concrete line number the code is running.
    Hence, it is difficult to debug in parallel mode
    """

    def __init__(self, contig: str, logger: Logger) -> None:
        """Initialize logger with contig name."""
        self.logger = logger
        self.contig = contig

    def debug(self, msg: str) -> None:
        """Wrapper for debug method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.debug(msg)

    def info(self, msg: str) -> None:
        """Wrapper for info method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.info(msg)

    def warning(self, msg: str) -> None:
        """Wrapper for warning method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.warning(msg)

    def error(self, msg: str) -> None:
        """Wrapper for error method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.error(msg)

    def critical(self, msg: str) -> None:
        """Wrapper for critical method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.critical(msg)

    def trace(self, msg: str) -> None:
        """Wrapper for trace method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.trace(msg)

    def success(self, msg: str) -> None:
        """Wrapper for success method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.success(msg)

    def complete(self) -> None:
        """Wrapper for complete method in logger."""
        self.logger.complete()
