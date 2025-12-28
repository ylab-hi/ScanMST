#!/usr/bin/env python3
"""Test classes.
"""
import os

import pytest
from scanmst import cppext as cpp


@pytest.mark.skip(reason="TODO")
@pytest.mark.parametrize(
    (
        "mapq_threshold",
        "soft_len_threshold",
        "mismatch_threshold",
        "align_len_threshold",
        "chrom",
        "start",
        "mode",
        "current_names",
        "all_names",
    ),
    [(15, 5, 3, 0.8, "chr17", 7708250, 1, ["one_hope_back"], ["one_hope_back"])],
)
def test_calculate_sr(
    mapq_threshold,
    soft_len_threshold,
    mismatch_threshold,
    align_len_threshold,
    chrom,
    start,
    mode,
    current_names,
    all_names,
) -> None:
    """Test calculate_sr."""
    path = os.path.dirname(__file__)
    os.chdir(path)
    #
    cppext_rescuer = cpp.Rescuer(
        "../data/test_calculate_sr.bam",
        mapq_threshold,
        soft_len_threshold,
        mismatch_threshold,
        align_len_threshold,
    )

    sr = cppext_rescuer.calculate_sr(
        chrom,
        start,
        start,
        mode,
        current_names,
        all_names,
    )
    assert sr == 4
