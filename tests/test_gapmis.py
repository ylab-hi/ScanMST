# !/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@Filename:    test_gapmis.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        12/10/21 6:32 PM
"""
from scannls.draft.nls_inference import Aligner
from scannls.draft.nls_inference import AlignerResult


def test_run():
    seq1 = "ATCGACGTGCAG"
    seq2 = "ATCGACGATCGA"

    aligner = Aligner(seq1, seq2)
    alignment_result = aligner.run()

    assert isinstance(alignment_result, AlignerResult)
