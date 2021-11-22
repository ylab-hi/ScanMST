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
    def __init__(self, tool):
        super(ToolNotFoundError, self).__init__(
            f"external tool: {tool} not found, please install that!"
        )
        self.tool = tool


class ReadNotFoundError(ScannlsException):
    def __init__(self):
        super(ReadNotFoundError, self).__init__(
            f"Current read cannot found in read_chains"
        )
