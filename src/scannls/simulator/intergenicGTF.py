#!/usr/bin/env python
"""INTERGENIC GTF FILE."""
import os
import shutil
import subprocess

from scannls import LoggerType

__license__ = "MIT"
__revision__ = " $Id: actor.py 1586 2009-01-30 15:56:25Z cokelaer $ "
__docformat__ = "reStructuredText"


class Intergenic:
    """Obtain gene upstream/downstream intergenic regions.

    The first line is brief explanation, which may be completed with
    a longer one. For instance to discuss about its methods. The only
    method here is :func:`function1`'s. The main idea is to document
    the class and methods's arguments with

    - **parameters**, **types**, **return** and **return types**::

          :param arg1: description
          :param arg2: description
          :type arg1: type description
          :type arg1: type description
          :return: return description
          :rtype: the return type description

    - and to provide sections such as **Example** using the double commas syntax::

          :Example:


      which appears as follow:

      :Example:

      followed by a blank line

    - Finally special sections such as **See Also**, **Warnings**, **Notes**
      use the sphinx syntax (*paragraph directives*)::

          .. seealso:: blabla
          .. warnings also:: blabla
          .. note:: blabla
          .. todo:: blabla

    .. note::
        There are many other Info fields but they may be redundant:
            * param, parameter, arg, argument, key, keyword: Description of a
              parameter.
            * type: Type of a parameter.
            * raises, raise, except, exception: That (and when) a specific
              exception is raised.
            * var, ivar, cvar: Description of a variable.
            * returns, return: Description of the return value.
            * rtype: Return type.

    .. note::
        There are many other directives such as versionadded, versionchanged,
        rubric, centered, ... See the sphinx documentation for more details.

    Here below is the results of the :func:`function1` docstring.

    """

    def __init__(
        self,
        input_gtf: str,
        output_gtf: str,
        logger: LoggerType,
        protein_coding_only: bool = True,
    ) -> None:
        """Input annotation GTF, output gene-intergenic GTF."""
        self.in_gtf = input_gtf
        self.out_gtf = output_gtf
        self.logger = logger
        self.protein_coding_only = protein_coding_only

    def _run_cmd(self, cmd: str) -> None:
        """Function is used to run the command in the system.

        :param cmd: the command to be run
        """
        subprocess.check_call(cmd, shell=True)

    @staticmethod
    def _remove(infile):
        if os.path.isfile(infile):
            os.remove(infile)
        elif os.path.isdir(infile):
            shutil.rmtree(infile)

    @staticmethod
    def prepare_genome_file(ref="hg38"):
        """Prepare genome size file for complementBed."""
        hg38_dict = {
            "chr1": 248956422,
            "chr10": 133797422,
            "chr11": 135086622,
            "chr12": 133275309,
            "chr13": 114364328,
            "chr14": 107043718,
            "chr15": 101991189,
            "chr16": 90338345,
            "chr17": 83257441,
            "chr18": 80373285,
            "chr19": 58617616,
            "chr2": 242193529,
            "chr20": 64444167,
            "chr21": 46709983,
            "chr22": 50818468,
            "chr3": 198295559,
            "chr4": 190214555,
            "chr5": 181538259,
            "chr6": 170805979,
            "chr7": 159345973,
            "chr8": 145138636,
            "chr9": 138394717,
            "chrM": 16569,
            "chrX": 156040895,
            "chrY": 57227415,
        }
        hg19_dict = {
            "chr1": 249250621,
            "chr10": 135534747,
            "chr11": 135006516,
            "chr12": 133851895,
            "chr13": 115169878,
            "chr14": 107349540,
            "chr15": 102531392,
            "chr16": 90354753,
            "chr17": 81195210,
            "chr18": 78077248,
            "chr19": 59128983,
            "chr2": 243199373,
            "chr20": 63025520,
            "chr21": 48129895,
            "chr22": 51304566,
            "chr3": 198022430,
            "chr4": 191154276,
            "chr5": 180915260,
            "chr6": 171115067,
            "chr7": 159138663,
            "chr8": 146364022,
            "chr9": 141213431,
            "chrM": 16571,
            "chrX": 155270560,
            "chrY": 59373566,
        }
        file_name = f"{ref}.genome"
        if ref == "hg38":
            with open(file_name, "w") as f:
                for i in hg38_dict:
                    f.write(f"{i}\t{hg38_dict[i]}\n")
        else:
            with open(file_name, "w") as f:
                for i in hg19_dict:
                    f.write(f"{i}\t{hg19_dict[i]}\n")
        return file_name

    def _closest_parser(self, upstream_bed, downstream_bed) -> str:
        """Parser of closestBed for gene upstream/downstream intergenic."""
        upstream_dict = {}
        with open(upstream_bed) as f:
            for line in f:
                _up_l = line.rstrip().split("\t")
                key = "\t".join(_up_l[:9])
                intergenic = f"{_up_l[-3]}-{_up_l[-2]}"
                distance = _up_l[-1]
                upstream_dict[key] = {"intergenic": intergenic, "distance": distance}

        with open(downstream_bed) as f, open(self.out_gtf, "w") as output:
            for line in f:
                _down_l = line.rstrip().split("\t")
                key = "\t".join(_down_l[:9])
                intergenic = f"{_down_l[-3]}-{_down_l[-2]}"
                distance = _down_l[-1]
                additional_field = (
                    f""" upstream_intergenic "{upstream_dict[key]['intergenic']}"; """
                    f"""upstream_distance "{upstream_dict[key]['distance']}"; """
                    f"""downstream_intergenic "{intergenic}"; downstream_distance "{distance}";"""
                )
                _down_l[8] = _down_l[8] + additional_field
                new_line = "\t".join(_down_l[:9])
                output.write(f"{new_line}\n")
        return self.out_gtf

    def run(self) -> str:
        """Run the gene to intergenic.

        :param arg1: the first value
        :param arg2: the first value
        :param arg3: the first value
        :type arg1: int, float,...
        :type arg2: int, float,...
        :type arg3: int, float,...
        :returns: arg1/arg2 +arg3
        :rtype: int, float

        :Example:

        >>> import template
        >>> a = template.MainClass1()
        >>> a.function1(1,1,1)
        2

        .. note:: can be useful to emphasize
            important feature
        .. seealso:: :class:`MainClass2`
        .. warning:: arg2 must be non-zero.
        .. todo:: check that arg2 is non zero.
        """
        name = os.path.splitext(os.path.basename(self.in_gtf))[0]

        genome_file = Intergenic.prepare_genome_file()
        _gene_flag = ""
        if self.protein_coding_only:
            _gene_flag = ' && /type "protein_coding/ '
        # refseq: /gene_biotype "protein_coding/
        # gencode: /gene_type "protein_coding/
        # step1: intergenic bed
        cmd = (
            f"""cat {self.in_gtf} | awk 'BEGIN{{OFS="\\t";}} $3=="gene"{_gene_flag}{{print $1,$4-1,$5}}' """
            f"""| sortBed | mergeBed -i - | complementBed -i stdin -g {genome_file} > {name}.intergenic.bed"""
        )
        self.logger.info(f"{cmd}")
        self._run_cmd(cmd)
        Intergenic._remove(f"{genome_file}")

        # step2: gene region gtf
        cmd = (
            f"""cat {self.in_gtf} | awk 'BEGIN{{OFS="\\t";}} $3=="gene"{_gene_flag}{{print $0}}' """
            f"""| sortBed > {name}.genic.gtf"""
        )
        self.logger.info(f"{cmd}")
        self._run_cmd(cmd)
        # step3: obtain closet intergenic of genes
        # https://bedtools.readthedocs.io/en/latest/content/tools/closest.html
        cmd_upstream = (
            f"closestBed -D a -id -a {name}.genic.gtf "
            f"-b {name}.intergenic.bed > {name}.upstream.bed"
        )
        cmd_downstream = (
            f"closestBed -D a -iu -a {name}.genic.gtf "
            f"-b {name}.intergenic.bed > {name}.downstream.bed"
        )
        self.logger.info(f"{cmd_upstream}")
        self.logger.info(f"{cmd_downstream}")
        self._run_cmd(cmd_upstream)
        self._run_cmd(cmd_downstream)

        Intergenic._remove(f"{name}.intergenic.bed")
        Intergenic._remove(f"{name}.genic.gtf")

        _upstream_file = f"{name}.upstream.bed"
        _downstream_file = f"{name}.downstream.bed"

        output_file = self._closest_parser(_upstream_file, _downstream_file)

        Intergenic._remove(_upstream_file)
        Intergenic._remove(_downstream_file)

        return output_file


if __name__ == "__main__":
    from loguru import logger

    intergenic = Intergenic(
        input_gtf="gencode.v37.annotation.gtf", output_gtf="test.gtf", logger=logger
    )
    intergenic.run()
