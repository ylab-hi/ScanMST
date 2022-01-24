# !/usr/bin/env python
"""Exception for the simulator.

@Filename:    exception.py
@license:     MIT Licence
@Time:        1/24/22 9:26 AM
"""


class NumberOfHopIsNotValidError(Exception):
    """Exception raised for errors when breakpoint not found."""

    def __init__(self) -> None:
        """Initialize the exception."""
        super().__init__("The number of hop is not valid")
