# !/usr/bin/env python
"""End-to-end tests for the CLI.
"""
import os
from pathlib import Path

import pytest
from scanmst.cli import DefaultOptions, cli


@pytest.mark.xfail(
    reason=(
        "cli/main.py iter_bam() exits on the bundled data. Its -Y check lives "
        "inside the read loop and compares a running count of supplementary "
        "alignments lacking soft clips against the running total, so a single "
        "bad first supplementary alignment gives 1 >= 1 and raises SystemExit. "
        "Either the check belongs after the loop or these BAMs need "
        "regenerating with -Y; both are calls for a maintainer."
    ),
    raises=SystemExit,
    strict=True,
)
@pytest.mark.parametrize(
    ("data_name", "parallel"),
    [("INV_TDUP", 1), ("TDUP_TRA", 1), ("TDUP_TRA", 2), ("INV_TDUP", 2)],
)
def test_cli(tmpdir, data_name, parallel):
    """Test the CLI."""
    path = os.path.dirname(__file__)
    os.chdir(path)

    data_dir = Path("../data/")
    out_dir = tmpdir.mkdir("out")
    op = DefaultOptions(
        input=f"{data_dir}/{data_name}.bam",
        output=f"{out_dir}/{data_name}",
        ref=f"{data_dir}/dummy.fasta",
        gtf=f"{data_dir}/dummy.gtf",
        blat_two_bit=f"{data_dir}/dummy.2bit",
        noncanonical=True,
        log="WARNING",
        thread=parallel,
    )
    # Run the pipeline and compare the outputs against the checked-in expected
    # files. These comparisons used to sit inside a pytest.raises block whose
    # call raised first, and the block ended with a hand-written
    # `raise ToolNotFoundError("This is a test")` to satisfy it -- so the suite
    # could not tell a passing run from a failing one. See the xfail above for
    # why the run currently stops early.
    cli.cli(op)

    with open(f"{out_dir}/{data_name}.fasta") as of, open(
        f"{out_dir}/{data_name}.gtf",
    ) as og, open(f"{data_dir}/{data_name}.fasta") as ef, open(
        f"{data_dir}/{data_name}.gtf",
    ) as eg:
        out_fasta = of.readlines()
        out_gtf = og.readlines()
        expect_fasta = ef.readlines()
        expect_gtf = eg.readlines()

    assert expect_fasta[1].strip() == out_fasta[1].strip()
    assert len(expect_gtf) == len(out_gtf)
    # chrom
    assert expect_gtf[0].split()[0].strip() == out_gtf[0].split()[0].strip()
    # start of the first exon
    assert expect_gtf[0].split()[3].strip() == out_gtf[0].split()[3].strip()
    # end of the first exon
    assert expect_gtf[0].split()[4].strip() == out_gtf[0].split()[4].strip()
    # start of the last exon
    assert expect_gtf[-1].split()[3].strip() == out_gtf[-1].split()[3].strip()
