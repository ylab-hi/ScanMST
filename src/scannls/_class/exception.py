#!/usr/bin/env python3
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
        super().__init__(f"external tool: {tool} not found, please install that!")


class ReadNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when read not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super().__init__("Current read cannot found in read_chains")


class ReadNotConnectedError(ScannlsExceptionError):
    """Exception raised for errors when read not connected."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super().__init__("Start read cannot connect all reads in candidate_nodes")


class SeqNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when sequence not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super().__init__("Gapmis: Sequence not found for semi-global alignment")


class ExonsNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when exons not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super().__init__("Cannot find exons in current node")


class ModesNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when modes not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super().__init__("Cannot find modes in current node")


class GenesNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when genes not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super().__init__("Cannot find genes in current node")


class BreakpointNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when breakpoint not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super().__init__("Cannot find breakpoint in current node")


class NumberOfHopIsNotValidError(ScannlsExceptionError):
    """Exception raised for errors when breakpoint not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super().__init__("The number of hop is not valid")
