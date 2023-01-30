# !/usr/bin/env python
"""End-to-end tests for the CLI.

@Author:      YangyangLi
@Filename:    test_cli.py
@license:     MIT Licence
@Time:        2/6/22 11:43 AM
"""
import os
from pathlib import Path

import pytest
from scannls import DefaultOptions, ToolNotFoundError, cli


@pytest.mark.parametrize(
    "data_name, parallel",
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
        two_bit=f"{data_dir}/dummy.2bit",
        noncanonical=True,
        log="WARNING",
        parallel=parallel,
    )
    # Test data do not use blat
    with pytest.raises(ToolNotFoundError):
        cli.cli(op)
        with open(f"{out_dir}/{data_name}.fasta") as of, open(
            f"{out_dir}/{data_name}.gtf"
        ) as og, open(f"{data_dir}/{data_name}.fasta") as ef, open(
            f"{data_dir}/{data_name}.gtf"
        ) as eg:
            out_fasta = of.readlines()
            out_gtf = og.readlines()
            expect_fasta = ef.readlines()
            expect_gtf = eg.readlines()

        assert expect_fasta[1].strip() == out_fasta[1].strip()
        assert len(expect_gtf) == len(out_gtf)
        assert (
            expect_gtf[0].split()[0].strip() == out_gtf[0].split()[0].strip()
        )  # check chrom
        assert (
            expect_gtf[0].split()[3].strip() == out_gtf[0].split()[3].strip()
        )  # check start for first exon
        assert (
            expect_gtf[0].split()[4].strip() == out_gtf[0].split()[4].strip()
        )  # check end for first exon
        assert (
            expect_gtf[-1].split()[3].strip() == out_gtf[-1].split()[3].strip()
        )  # check start for random exon
        raise ToolNotFoundError("This is a test")
