# !/usr/bin/env python
# -*- coding:utf-8 -*-
"""Wrapper for loguru.logger.

@Filename:    mylogger.py
@license:     MIT Licence
@Time:        12/15/21 2:08 PM
"""


class MyLogger(object):
    """Wrapper for logger in order to use in multiprocessing.

    to show contig name in logging information before message.
    However, there is no way to get concrete line number the code is running.
    Hence, it is difficult to debug in parallel mode
    """

    def __init__(self, contig, logger):
        """Initialize logger with contig name."""
        self.logger = logger
        self.contig = contig

    def debug(self, msg: str):
        """Wrapper for debug method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.debug(msg)

    def info(self, msg: str):
        """Wrapper for info method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.info(msg)

    def warning(self, msg: str):
        """Wrapper for warning method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.warning(msg)

    def error(self, msg: str):
        """Wrapper for error method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.error(msg)

    def critical(self, msg: str):
        """Wrapper for critical method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.critical(msg)

    def trace(self, msg: str):
        """Wrapper for trace method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.trace(msg)

    def success(self, msg: str):
        """Wrapper for success method in logger."""
        msg = f"{self.contig}: {msg}"
        self.logger.success(msg)

    def complete(self):
        """Wrapper for complete method in logger."""
        self.logger.complete()
