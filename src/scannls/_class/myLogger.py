# !/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Filename:    mylogger.py
@license:     MIT Licence
@Time:        12/15/21 2:08 PM
"""


class MyLogger(object):
    """
    wrapper for logger in order to use in multiprocessing
    to show contig name in logging information before message.
    However, there is no way to get concrete line number the code is running.
    Hence, it is difficult to debug in parallel mode
    """

    def __init__(self, contig, logger):
        self.logger = logger
        self.contig = contig

    def debug(self, msg):
        msg = f"{self.contig}: {msg}"
        self.logger.debug(msg)

    def info(self, msg):
        msg = f"{self.contig}: {msg}"
        self.logger.info(msg)

    def warning(self, msg):
        msg = f"{self.contig}: {msg}"
        self.logger.warning(msg)

    def error(self, msg):
        msg = f"{self.contig}: {msg}"
        self.logger.error(msg)

    def critical(self, msg):
        msg = f"{self.contig}: {msg}"
        self.logger.critical(msg)

    def trace(self, msg):
        msg = f"{self.contig}: {msg}"
        self.logger.trace(msg)

    def success(self, msg):
        msg = f"{self.contig}: {msg}"
        self.logger.success(msg)

    def complete(self):
        self.logger.complete()
