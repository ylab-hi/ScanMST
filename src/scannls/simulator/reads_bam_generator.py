#!/usr/bin/env python
"""A user-friendly script to generate simulated reads using PBSIM2 (https://github.com/ylab-hi/pbsim2)."""
import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO  # type: ignore
from loguru import logger


class ProfileNotFoundError(Exception):
    """Exception raised for errors when external tool not found."""

    def __init__(self, file: str) -> None:
        """Initialize the exception."""
        super().__init__(
            f"profile file: {file} not found, please rerun sampling-based PBSIM2!"
        )


def status_message(msg) -> None:
    """Print status message."""
    print(msg)
    sys.stdout.flush()


def remove(infile):
    """Remove files."""
    if os.path.isfile(infile):
        os.remove(infile)


def run_cmd(cmd, logger=logger):
    """Run cmd with message."""
    logger.info(cmd)
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as err:
        logger.warning(f"Error happend!: {err}\n{err.output}")
        raise SystemExit from None
    else:
        return True


def fastq_simulation(in_fa_file, out_prefix, model, logger, depth=10):
    """Generate simulated fastq using pbsim_rna."""

    def combine_fastq(out_prefix):
        with open(f"{out_prefix}.fastq", "w") as out:
            for _filename in glob.iglob(f"{out_prefix}_*.fastq"):
                with open(_filename) as f:
                    out.write(f.read())
                _prefix = os.path.splitext(os.path.basename(_filename))[0]
                remove(_filename)
                remove(f"{_prefix}.ref")
                remove(f"{_prefix}.maf")

    profile_checker(model, logger)

    cmd = f"pbsim_rna --depth {depth} --prefix {out_prefix} --sample-profile-id {model} {in_fa_file}"
    run_cmd(cmd, logger)
    combine_fastq(out_prefix)
    return f"{out_prefix}.fastq"


def profile_checker(model, logger) -> None:
    """Check profile files of PBSIM2."""
    profile_files = [f"sample_profile_{model}.stats", f"sample_profile_{model}.fastq"]
    for _file in profile_files:
        file = Path(_file)
        if not (file.exists() and file.stat().st_size > 0):
            raise ProfileNotFoundError(_file)
        else:
            logger.success(f"Checking for {_file} found ")


def combine_fastq(in_mt_fq, in_wt_fq, out_fq) -> str:
    """Combine MT and WT fastq into one."""
    with open(out_fq, "w") as out:
        for r in SeqIO.parse(in_mt_fq, "fastq"):
            SeqIO.write(r, out, "fastq")

        if in_wt_fq:
            for r in SeqIO.parse(in_wt_fq, "fastq"):
                r.id = f"{r.id}-WT"
                SeqIO.write(r, out, "fastq")
    return out_fq


def alignment_runner(in_fq, ref_fa, bigbed, data_type, thread_num, out_prefix, logger):
    """Warpper for minimap2."""
    if data_type == "pacbio":
        cmd1 = (
            f"minimap2 -t {thread_num} -Y -ax splice:hq -uf -R "
            f'"@RG\\tID:{out_prefix}\\tSM:hs\\tLB:ga\\tPL:PacBio" '
            f"--MD --secondary=no --junc-bed {bigbed} {ref_fa} {in_fq} > {out_prefix}.tmp.sam"
        )
    elif data_type == "nanopore":
        cmd1 = (
            f"minimap2 -t {thread_num} -Y -ax splice -uf -k14 -R "
            f'"@RG\\tID:{out_prefix}\\tSM:hs\\tLB:ga\\tPL:ONT" '
            f"--MD --junc-bed {bigbed} {ref_fa} {in_fq} > {out_prefix}.tmp.sam"
        )
    cmd2 = (
        f"samtools sort -@ {thread_num} -O SAM {out_prefix}.tmp.sam -o {out_prefix}.sam"
    )

    cmd3 = "samtools view -b {0}.sam -o {0}.bam".format(out_prefix)
    cmd4 = f"samtools index {out_prefix}.bam"

    step1 = run_cmd(cmd1, logger)
    if step1:
        step2 = run_cmd(cmd2, logger)
        if step2:
            remove(f"{out_prefix}.tmp.sam")
            run_cmd(f"{cmd3} && {cmd4}")
            return f"{out_prefix}.sam", f"{out_prefix}.bam"

    return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "-m",
        "--mt_input",
        action="store",
        dest="mt_input",
        help="MT fasta file",
        required=True,
    )
    parser.add_argument(
        "-w",
        "--wt_input",
        action="store",
        dest="wt_input",
        help="WT fasta file",
        required=True,
    )
    parser.add_argument(
        "--model", action="store", dest="model", help="error model name"
    )
    parser.add_argument(
        "-t",
        "--thread",
        action="store",
        dest="thread",
        type=int,
        help="The number of threads (default: %(default)s)",
        default=1,
    )
    parser.add_argument(
        "--mt_dp",
        action="store",
        dest="mt_depth",
        type=int,
        help="MT depth (default: %(default)s)",
        default=10,
    )
    parser.add_argument(
        "--wt_dp",
        action="store",
        dest="wt_depth",
        type=int,
        help="WT depth (default: %(default)s)",
        default=90,
    )
    parser.add_argument(
        "-r",
        "--ref",
        action="store",
        dest="ref",
        help="reference FASTA (default: %(default)s)",
        default="/panfs/home/yang4414/tywang/database/genome/hg38_canon.fa",
    )
    parser.add_argument(
        "-a",
        "--annot",
        action="store",
        dest="annot",
        help="gene annotation bigbed (default: %(default)s)",
        default="/panfs/home/yang4414/tywang/database/gencode/gencode.v38.bigbed",
    )
    parser.add_argument(
        "-l",
        "--library",
        action="store",
        dest="library",
        help="output type (default: %(default)s)",
        choices=["pacbio", "nanopore"],
        default="pacbio",
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        help="output prefix (default: %(default)s)",
        default="output",
    )
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")
    args = parser.parse_args()

    mt_fa = args.mt_input
    wt_fa = args.wt_input
    ref_fa = args.ref
    mt_depth = args.mt_depth
    wt_depth = args.wt_depth
    out_prefix = args.output
    thread_num = args.thread
    model = args.model

    mt_fq = fastq_simulation(
        in_fa_file=mt_fa,
        out_prefix=f"{out_prefix}.MT",
        model=model,
        logger=logger,
        depth=mt_depth,
    )

    wt_fq = None
    if wt_depth > 0:
        wt_fq = fastq_simulation(
            in_fa_file=wt_fa,
            out_prefix=f"{out_prefix}.WT",
            model=model,
            logger=logger,
            depth=wt_depth,
        )

    out_fq = f"{out_prefix}.fastq"

    combined_fq = combine_fastq(mt_fq, wt_fq, out_fq)

    alignment_runner(
        in_fq=combined_fq,
        ref_fa=ref_fa,
        bigbed=args.annot,
        data_type=args.library,
        thread_num=thread_num,
        out_prefix=out_prefix,
        logger=logger,
    )
