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

    def __init__(self, msg: str) -> None:
        """Initialize the exception."""
        super().__init__(f"Cannot find exons in current node {msg}")


class ModesNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when modes not found."""

    def __init__(self, msg: str) -> None:
        """Initialize the exception."""
        super().__init__(f"Cannot find modes in current node {msg}")


class ModesNotEqualError(ScannlsExceptionError):
    """Exception raised for errors when modes not equal."""

    def __init__(self, msg: str) -> None:
        """Initialize the exception."""
        super().__init__(
            f"modes not equal in the same chrom and the different strands from reads {msg}"
        )


class GenesNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when genes not found."""

    def __init__(self, msg: str) -> None:
        """Initialize the exception."""
        super().__init__(f"Cannot find genes in current node {msg}")


class SplicingCodeNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when splicing code not found."""

    def __init__(self, msg: str) -> None:
        """Initialize the exception."""
        super().__init__(f"Cannot find splicing code in current node {msg}")


class AnnotationCodeNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when annotation code not found."""

    def __init__(self, msg: str) -> None:
        """Initialize the exception."""
        super().__init__(f"Cannot find annotation code in current node {msg}")


class BreakpointNotFoundError(ScannlsExceptionError):
    """Exception raised for errors when breakpoint not found."""

    def __init__(self, msg: str) -> None:
        """Initialize the exception."""
        super().__init__(f"Cannot find breakpoint in current node {msg}")
