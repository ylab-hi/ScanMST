#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module for the exception class.

@version: 0.0.1
@license: MIT Licence
@file: exception.py
@time: 16/11/2021 11:28
"""


class ScannlsExceptionError(Exception):
    """Base class for exceptions in this module."""

    pass


class ToolNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when external tool not found."""

    def __init__(self, tool: str) -> None:
        """Initialize the exception."""
        super(ToolNotFoundError, self).__init__(
            f"external tool: {tool} not found, please install that!"
        )
        self.tool = tool


class ReadNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when read not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super(ReadNotFoundError, self).__init__(
            "Current read cannot found in read_chains"
        )


class ReadNotConnectedError(ScannlsExceptionError):
    """Exception raised for errors when read not connected."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super(ReadNotConnectedError, self).__init__(
            "Start read cannot connect all reads in candidate_nodes"
        )


class SeqNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when sequence not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super(SeqNotFoundError, self).__init__(
            "Gapmis: Sequence not found for semi-global alignment"
        )
