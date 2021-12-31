# !/usr/bin/env python
"""SR rescuer.

@Filename:    srRescuer.py
@license:     MIT Licence
@Time:        12/30/21 15:00 PM
"""
from loguru._logger import Logger

"""
import copy
import types
from typing import Any  # ignore
from typing import Dict  # ignore
from typing import Iterable  # ignore
from typing import List  # ignore
from typing import Set  # ignore
from typing import Tuple  # ignore
from typing import Union  # ignore

from ..utils import timeit  # ignore
from .basicClass import Node  # ignore
from .basicClass import Read  # ignore
from .basicClass import Series  # ignore

NodeType = Union[Node, Insertion]
"""


class Rescuer:
    """Rescue SR from softclipped non-chimeric reads."""

    def __init__(self, logger: Logger) -> None:
        """Initialize Rescuer.

        :param logger: logger
        """
        self.logger = logger

    def __repr__(self):
        """Represent Rescuer."""
        return f"{self.__class__.__name__}()"
