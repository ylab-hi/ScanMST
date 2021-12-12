#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@version: 0.0.1
@license: MIT Licence
@file: exception.py
@time: 16/11/2021 11:28
"""


class ScannlsException(Exception):
    pass


class ToolNotFoundError(ScannlsException):
    def __init__(self, tool: str) -> None:
        super(ToolNotFoundError, self).__init__(
            f"external tool: {tool} not found, please install that!"
        )
        self.tool = tool


class ReadNotFoundError(ScannlsException):
    def __init__(self) -> None:
        super(ReadNotFoundError, self).__init__(
            f"Current read cannot found in read_chains"
        )


class ReadNotConnectedError(ScannlsException):
    def __init__(self) -> None:
        super(ReadNotConnectedError, self).__init__(
            "Start read cannot connect all reads in candidate_nodes"
        )


class SeqNotFoundError(ScannlsException):
    def __init__(self) -> None:
        super(SeqNotFoundError, self).__init__(
            "Gapmis: Sequence not found for semi-global alignment"
        )
