import copy
import re
import sys
from collections import defaultdict

from align import aligner

from .. import __version__


try:
    import pysam
    import numpy as np
    import HTSeq
except ModuleNotFoundError as e:
    raise SystemExit(e.msg)

__funcs__ = {"vcf_header", "sv_checker"}


def vcf_header(output_prefix, bam_header):
    """current_output VCF header using information from BAM header
    :param output_prefix: file prefix for VCF file, usually uses sample name
    :param bam_header: header of BAM file
    :type output_prefix: str
    :type bam_header: str
    :return: header of VCF file
    :rtype: str
    .. note::
        GROUP field is important,
        multiple NLS events with the same 'group' information will form one transcript
    """

    _aligners = {
        "CLC",
        "ContextMap2",
        "CRAC",
        "GSNAP",
        "HISAT",
        "HISAT2",
        "MapSplice2",
        "Novoalign",
        "OLego",
        "RUM",
        "SOAPsplice",
        "STAR",
        "Subread",
        "TopHat",
        "TopHap2",
        "bwa",
        "bowtie",
        "bowtie2",
        "minimap2",
    }
    avail_aligners = {x.upper() for x in _aligners}
    if "PG" in bam_header:
        for j in bam_header["PG"]:
            if j["ID"].upper() in avail_aligners:
                aln_cmd = j["CL"]
                break
            else:
                aln_cmd = "Unknown"
    else:
        aln_cmd = "Unknown"

    header = ["##fileformat=VCFv4.1"]
    header.append("##source=ScanNLS " + __version__)
    header.append(f'##reference=<CMD={aln_cmd},Description="Alignment parameters">')
    header.append('##ALT=<ID=TDUP,Description="Tandem duplication">')
    header.append('##ALT=<ID=INV,Description="Inversion">')
    header.append('##ALT=<ID=TRA,Description="Translocation">')
    header.append(
        '##INFO=<ID=CANONICAL,Number=0,Type=Flag,Description="Canonical splice site">'
    )
    header.append(
        '##INFO=<ID=NONCANONICAL,Number=0,Type=Flag,Description="Noncanonical splice site">'
    )
    header.append(
        '##INFO=<ID=BOUNDARY,Number=1,Type=String,Description="The coding exon boundary type of event, BOTH, LEFT, RIGHT, NEITHER.">'
    )
    header.append(
        '##INFO=<ID=DP1,Number=1,Type=Integer,Description="Total read depth at the breakpoint1">'
    )
    header.append(
        '##INFO=<ID=DP2,Number=1,Type=Integer,Description="Total read depth at the breakpoint2">'
    )
    header.append(
        '##INFO=<ID=SR,Number=1,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">'
    )
    header.append(
        '##INFO=<ID=PSO,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1], representing the ratio of reads showing the alternative allele to all reads">'
    )
    header.append(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The type of event, INS, DEL, TDUP, INV, TRA.">'
    )
    header.append(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'
    )
    header.append(
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">'
    )
    header.append(
        '##INFO=<ID=GROUP,Number=1,Type=String,Description="events in the group from the same transcript">'
    )
    header.append(
        '##INFO=<ID=GENE1,Number=1,Type=String,Description="Overlapped coding gene for breakpoint1">'
    )
    header.append(
        '##INFO=<ID=GENE2,Number=1,Type=String,Description="Overlapped coding gene for breakpoint2">'
    )
    header.append(
        '##INFO=<ID=STRAND1,Number=1,Type=String,Description="Strand for breakpoint1">'
    )
    header.append(
        '##INFO=<ID=STRAND2,Number=1,Type=String,Description="Strand for breakpoint2">'
    )
    header.append(
        '##INFO=<ID=MODE1,Number=1,Type=String,Description="mode for softclipped reads at breakpoint1">'
    )
    header.append(
        '##INFO=<ID=MODE2,Number=1,Type=String,Description="mode for softclipped reads at breakpoint2">'
    )
    header.append(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="nd position of the structural variant">'
    )
    header.append(
        '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">'
    )
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + output_prefix
    )
    return "\n".join(header)


# TODO
def sv_checker(read, ref_site, mapq_cutoff, sv_len_cutoff):
    """"""
    chimeric_aln = read.alignment.get_tag("SA")[:-1].split(";")
    if len(chimeric_aln) == 1:
        chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = chimeric_aln[0].split(",")
    else:
        chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = multiple_sa_tag_selector(
            chimeric_aln
        ).split(",")

    if int(mapq_sa) < mapq_cutoff:
        return False
    (
        sv_type,
        sv_anno_and_can,
        sv_pos,
        sv_length,
        read_mode,
        sa_mode,
        strands,
        genes,
    ) = read.alignment.get_tag("SV").split(",")

    if sv_type == "TRA":
        if int(sv_pos) - 2 + (int(read_mode) - 1) == ref_site:
            return True
        else:
            return False
    # TDUP and INV
    else:
        if int(sv_length) < sv_len_cutoff:
            return False
        if int(sv_pos) - 2 + (int(read_mode) - 1) == ref_site:
            return True
        elif int(sv_pos) - 2 + int(sv_length) + (int(read_mode) - 1) == ref_site:
            return True
        else:
            return False
