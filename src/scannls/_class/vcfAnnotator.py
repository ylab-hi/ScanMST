# !/usr/bin/env python
"""VcfAnnotation class using genome SV.

@Filename:    vcfAnnotator.py
@Author:      YangyangLi
@contact:     li002252@umn.edu
@license:     MIT Licence
@Time:        9/6/22 10:37 AM
"""
import shlex
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict
from typing import List

from pysam import VariantFile

from .type import LoggerType


class VcfAnnotator:
    """VcfAnnotator is used to annotate vcf file given genome structural variation info."""

    def __init__(self, vcf_file: str, sv_file: str, logger: LoggerType):
        """Initialize the VcfAnnotator class."""
        self._vcf = vcf_file
        self._out = vcf_file + ".annotated.vcf"

        self._sv = sv_file
        self.sv_map: Dict[str, List[str]] = defaultdict(list)

        self.logger = logger
        self._sv2nl_output = vcf_file + ".sv2nl"
        self.sv2nl = f"sv2nl {self._sv} {vcf_file} -o {self._sv2nl_output}"

    def annotate(self):
        """Annotate vcf file."""
        self._construct_sv_map()
        self._annotate()
        self.logger.info("Annotate vcf file finished")

    def _annotate(self):
        """Annotate vcf file and Add SSV info field."""
        vcf_handle = VariantFile(self._sv2nl_output)
        out_handle = VariantFile(self._out, "w", header=vcf_handle.header)

        for rec in vcf_handle.fetch():
            key = f"{rec.chrom}-{rec.pos}-{rec.info['SVEND']}-{rec.info['SVTYPE']}"
            if sv_info := self.sv_map.get(key):
                rec.info["SSV"] = sv_info
            out_handle.write(rec)
        vcf_handle.close()
        out_handle.close()

    def _run_sv2nl(self):
        """Run sv2nl to annotate vcf file."""
        try:
            subprocess.check_call(shlex.split(self.sv2nl))
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"Annotate sv failed due to sv2nl failed with {e}.")
        else:
            self.logger.info("Execute sv2nl successfully.")

    def _add_map(self, sv_type_file: str):
        """Add sv type file to sv map."""
        with open(sv_type_file) as sv_file:
            next(sv_file)  # skip header
            for line in sv_file:
                line = line.strip().split("\t")
                key = "-".join(line[:4])
                self.sv_map[key].append("-".join(line[4:]))

    def _construct_sv_map(self):
        """Construct sv map."""
        self._run_sv2nl()
        self._add_map(self._sv2nl_output + ".dup")
        self._add_map(self._sv2nl_output + ".inv")
        self._add_map(self._sv2nl_output + ".tra")

    def _clean(self):
        """Clean up the intermediate files."""
        Path(self._vcf).unlink()
