# !/usr/bin/env python
"""End-to-end tests for the CLI.
"""
import os
from pathlib import Path

import pytest
from scanmst import ToolNotFoundError
from scanmst.cli import DefaultOptions, cli


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
    # The bundled test data never reaches blat, so the run is expected to stop
    # with ToolNotFoundError. The comparisons below used to sit *inside* this
    # block after the raising call, which made them unreachable -- the test
    # asserted nothing. They now run against the output the CLI produced
    # before it bailed out.
    with pytest.raises(ToolNotFoundError):
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
