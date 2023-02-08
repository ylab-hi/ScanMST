# !/usr/bin/env python
"""Helper functions."""
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import HTSeq  # type: ignore
import pyfaidx  # type: ignore
import pysam  # type: ignore
import yaml  # type: ignore

from scannls import cppext

from .. import __PACKAGE_NAME__
from ..base.exception import ModesNotEqualError

__all__ = [
    "extract_splice_sites",
    "gene_annotation",
    "splicing_confirmation",
    "cigar_validity",
    "blat2chimeric_alignment",
    "insertion2chimeric_alignment",
    "same_chrom_same_strand_mode21_handler",
    "same_chrom_same_strand_handler",
    "same_chrom_diff_strand_handler",
    "diff_chrom_same_strand_mode21_handler",
    "diff_chrom_same_strand_handler",
    "diff_chrom_diff_strand_handler",
    "softclipped_length_and_event_size_checker",
    "obtain_bp_region_seq",
    "obtain_variants_stats",
    "strand_mode_checker",
]


def reverse_complement(seq: str) -> str:
    """Obtain reverse complement sequence."""
    rctrans = str.maketrans("ACGT", "TGCA")
    return str.translate(seq, rctrans)[::-1]


def extract_splice_sites(in_file: str, bin_size: int) -> Any:
    """Extract splice sites and gene regions from input GTF file.

    :param in_file: gene annotation file (GTF file)
    :param bin_size: bin size to search splice site
    :type in_file: str
    :type bin_size: int
    :return: annotated splice sites (HTSeq.GenomicArrayOfSets) and
        annotated gene regions (HTSeq.GenomicArrayOfSets)
    :rtype: tuple
    """
    # todo using real splice sites from reference genome
    gtf_file = HTSeq.GFF_Reader(in_file)
    cvg = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gene_iv = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    trx_to_exon = defaultdict(list)

    for feature in gtf_file:
        gene_name = feature.attr.get("gene_name") or feature.attr.get("gene")

        if feature.type == "exon":
            trx_id = feature.attr["transcript_id"]
            trx_to_exon[trx_id].append(feature.iv)
        if feature.type == "gene":
            gene_iv[
                HTSeq.GenomicInterval(
                    feature.iv.chrom,
                    feature.iv.start - bin_size,
                    feature.iv.end + bin_size,
                    ".",
                )
            ] += str(gene_name)

    for trx_id in trx_to_exon:
        exon_list = trx_to_exon[trx_id]
        exon_list.sort(key=lambda x: x.start)  # type: ignore
        exon_num = len(exon_list)
        first_exon = exon_list[0]
        last_exon = exon_list[-1]
        strand = first_exon.strand
        if exon_num == 1:
            cvg[
                HTSeq.GenomicInterval(
                    first_exon.chrom,
                    first_exon.start - bin_size,
                    first_exon.start + bin_size,
                    ".",
                )
            ] += "XX"
            cvg[
                HTSeq.GenomicInterval(
                    first_exon.chrom,
                    first_exon.end - bin_size,
                    first_exon.end + bin_size,
                    ".",
                )
            ] += "XX"
        elif exon_num == 2:
            if strand == "+":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin_size,
                        first_exon.start + bin_size,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin_size,
                        first_exon.end + bin_size,
                        ".",
                    )
                ] += "GT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin_size,
                        last_exon.start + bin_size,
                        ".",
                    )
                ] += "AG"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.end - bin_size,
                        last_exon.end + bin_size,
                        ".",
                    )
                ] += "XX"
            elif strand == "-":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin_size,
                        first_exon.start + bin_size,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin_size,
                        first_exon.end + bin_size,
                        ".",
                    )
                ] += "CT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin_size,
                        last_exon.start + bin_size,
                        ".",
                    )
                ] += "AC"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.end - bin_size,
                        last_exon.end + bin_size,
                        ".",
                    )
                ] += "XX"
        # exon_num > 2
        else:
            if strand == "+":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin_size,
                        first_exon.start + bin_size,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin_size,
                        first_exon.end + bin_size,
                        ".",
                    )
                ] += "GT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin_size,
                        last_exon.start + bin_size,
                        ".",
                    )
                ] += "AG"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.end - bin_size,
                        last_exon.end + bin_size,
                        ".",
                    )
                ] += "XX"
            elif strand == "-":
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.start - bin_size,
                        first_exon.start + bin_size,
                        ".",
                    )
                ] += "XX"
                cvg[
                    HTSeq.GenomicInterval(
                        first_exon.chrom,
                        first_exon.end - bin_size,
                        first_exon.end + bin_size,
                        ".",
                    )
                ] += "CT"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.start - bin_size,
                        last_exon.start + bin_size,
                        ".",
                    )
                ] += "AC"
                cvg[
                    HTSeq.GenomicInterval(
                        last_exon.chrom,
                        last_exon.end - bin_size,
                        last_exon.end + bin_size,
                        ".",
                    )
                ] += "XX"
            for _exon in exon_list[1:-1]:
                iv1 = HTSeq.GenomicInterval(
                    _exon.chrom, _exon.start - bin_size, _exon.start + bin_size, "."
                )
                iv2 = HTSeq.GenomicInterval(
                    _exon.chrom, _exon.end - bin_size, _exon.end + bin_size, "."
                )
                if strand == "+":
                    cvg[iv1] += "AG"
                    cvg[iv2] += "GT"
                elif strand == "-":
                    cvg[iv1] += "AC"
                    cvg[iv2] += "CT"
    return cvg, gene_iv


def gene_annotation(
    chrm1: str, pos1: int, chrm2: str, pos2: int, gene_iv: HTSeq.GenomicArrayOfSets
) -> tuple[str, str]:
    """Obtain gene annotations for breakpoints.

    :param chrm1: chromosome for breakpoint1
    :param chrm2: chromosome for breakpoint2
    :param pos1: position for breakpoint1
    :param pos2: position for breakpoint2
    :param gene_iv: gene annotations in HTSeq.GenomicArrayOfSets
    :type chrm1: str
    :type chrm2: str
    :type pos1: int
    :type pos2: int
    :type gene_iv: HTSeq.GenomicArrayOfSets
    :return: overlapped genes for breakpoints
    :rtype: tuple
    """
    gene1, gene2 = None, None
    try:
        gene1 = "&".join(list(gene_iv[HTSeq.GenomicPosition(chrm1, pos1)]))
    except IndexError:
        gene1 = ""
    except TypeError:
        pass
    try:
        gene2 = "&".join(list(gene_iv[HTSeq.GenomicPosition(chrm2, pos2)]))
    except IndexError:
        gene2 = ""
    except TypeError:
        pass

    if not gene1:
        gene1 = "INTERGENIC"
    if not gene2:
        gene2 = "INTERGENIC"
    return gene1, gene2


def splicing_confirmation(
    chrm1: str,
    pos1: int,
    strand1: str,
    mode1: int,
    chrm2: str,
    pos2: int,
    strand2: str,
    mode2: int,
    splice_bin: int,
    genome_fasta: pyfaidx.Fasta,
    cvg: HTSeq.GenomicArrayOfSets,
    motif_required: bool = True,
) -> tuple[bool, int, int]:
    """Judge whether the breakpoints are NLS events or not.

    if motif_required is ON: it will only report NLS events with 'canonical
    splice sites'; otherwise: it will report NLS events whatever the splice
    sites they used

    :param chrm1: chromosome for breakpoint1
    :param chrm2: chromosome for breakpoint2
    :param pos1: position for breakpoint1
    :param pos2: position for breakpoint2
    :param mode1: SM/MS mode for breakpoint1
    :param mode2: SM/MS mode for breakpoint2
    :param strand1: strand for breakpoint1
    :param strand2: strand for breakpoint2
    :param splice_bin: bin size for splice sites searching
    :param genome_fasta: reference genome (pyfaidx.Fasta object)
    :param cvg: splice site annotations (HTSeq.GenomicArrayOfSets)
    :param motif_required: canonical splice sites required;
           if True: considering canonical splice sites only;
           else: considering canonical and noncanonical splice sites both
    :return: report/not report, overlapping boundary in bits, canonical
        splice site/noncanonical splice site

    .. note::
        Possible current_output scenarios
        * True,  3(11), 1 => reported, both breakpoints overlap with annotated
            coding exons boundary, using canonical splice motif
        * True,  2(10), 0 => reported, one breakpoint overlap with annotated
            coding exons boundary, using noncanonical splice motif
        * True,  1(01), 0 => reported, one breakpoint overlap with annotated
            coding exons boundary, using noncanonical splice motif
        * False, 0(00), 0 => not reported, none breakpoint overlap with annotated
            coding exons boundary, using noncanonical splice motif

    """

    def matched_candidate_sites_checker(
        donor_seq: str, acceptor_seq, splice_motif_dict: dict[str, str]
    ) -> bool:
        """Find canonical splice sites in the input sequence.

        :param donor_seq: adjacent sequence at donor breakpoint
        :param acceptor_seq: adjacent sequence at acceptor breakpoint
        :param acceptor_seq: sequence at the other pair end
        :param splice_motif_dict: paired splice sites
        :return: canonical splice sites are paired or not
        :rtype: bool
        """
        paired = False
        for _donor, _acceptor in splice_motif_dict.items():
            if _donor in donor_seq and _acceptor in acceptor_seq:
                paired = True
                break
        return paired

    def donor_accepter_breakpoint_determintor(
        chrm1: str,
        pos1: int,
        strand1: str,
        mode1: int,
        chrm2: str,
        pos2: int,
        strand2: str,
        mode2: int,
    ) -> tuple[tuple[str, int, str], tuple[str, int, str]]:
        """Determine the donor breakpoint and the accepter breakpoint.

        :param chrm1: chromosome for breakpoint1
        :param chrm2: chromosome for breakpoint2
        :param pos1: position for breakpoint1
        :param pos2: position for breakpoint2
        :param mode1: SM/MS mode for breakpoint1
        :param mode2: SM/MS mode for breakpoint2
        :param strand1: strand for breakpoint1
        :param strand2: strand for breakpoint2
        """
        _breakpoint1 = (chrm1, pos1, strand1)
        _breakpoint2 = (chrm2, pos2, strand2)
        donor_accepter_dict = {
            "++21": (_breakpoint2, _breakpoint1),
            "++12": (_breakpoint1, _breakpoint2),
            "--21": (_breakpoint1, _breakpoint2),
            "--12": (_breakpoint2, _breakpoint1),
            "+-22": (_breakpoint2, _breakpoint1),
            "+-11": (_breakpoint1, _breakpoint2),
            "-+22": (_breakpoint1, _breakpoint2),
            "-+11": (_breakpoint2, _breakpoint1),
        }

        return donor_accepter_dict.get(f"{strand1}{strand2}{mode1}{mode2}", None)

    # key: strand of donor site, strand of accepter site
    # values: possible matched donor site and accepter site (>99% splice site using GT-AG)
    canonical_splice_dict = {
        "++": {"GT": "AG"},
        "+-": {"GT": "CT"},
        "-+": {"AC": "AG"},
        "--": {"AC": "CT"},
    }

    donor_bp, acceptor_bp = donor_accepter_breakpoint_determintor(
        chrm1, pos1, strand1, mode1, chrm2, pos2, strand2, mode2
    )

    chrm_do, pos_do, strand_do = donor_bp
    chrm_ac, pos_ac, strand_ac = acceptor_bp

    splice_motif_dict = canonical_splice_dict.get(f"{strand_do}{strand_ac}")

    if splice_motif_dict is None:
        raise ValueError("Invalid strand combination")

    possible_donors = {_v: _k for _k, _v in splice_motif_dict.items()}

    # motif_do/motif_ac will be available if chrm_do:pos_do/chrm_ac:pos_ac overlapped with annotated exon boundary
    try:
        motif_do = list(cvg[HTSeq.GenomicPosition(chrm_do, pos_do)])[0]
    except IndexError:
        motif_do = ""
    try:
        motif_ac = list(cvg[HTSeq.GenomicPosition(chrm_ac, pos_ac)])[0]
    except IndexError:
        motif_ac = ""

    # Non-annotated coding exon boundary
    if motif_do not in splice_motif_dict and motif_ac not in splice_motif_dict.values():
        donor_seq = genome_fasta[chrm_do][pos_do - splice_bin : pos_do + splice_bin].seq
        acceptor_seq = genome_fasta[chrm_ac][
            pos_ac - splice_bin : pos_ac + splice_bin
        ].seq
        if matched_candidate_sites_checker(donor_seq, acceptor_seq, splice_motif_dict):
            return True, 0, 1
        else:
            return (False, 0, 0) if motif_required else (True, 0, 0)

    # pos1 in annotated coding exon boundary, pos2 not.
    elif motif_do in splice_motif_dict and motif_ac not in splice_motif_dict.values():
        acceptor_seq = genome_fasta[chrm_ac][
            pos_ac - splice_bin : pos_ac + splice_bin
        ].seq
        if splice_motif_dict[motif_do] in acceptor_seq:
            return True, 2, 1
        else:
            return (False, 2, 0) if motif_required else (True, 2, 0)

    # pos2 in annotated coding exon boundary, pos1 not.
    elif motif_do not in splice_motif_dict and motif_ac in splice_motif_dict.values():
        donor_seq = genome_fasta[chrm_do][pos_do - splice_bin : pos_do + splice_bin].seq

        if possible_donors[motif_ac] in donor_seq:
            return True, 1, 1
        else:
            return (False, 1, 0) if motif_required else (True, 1, 0)

    # pos1 and pos2 both in annotated coding exon boundary
    elif motif_do in splice_motif_dict and motif_ac in splice_motif_dict.values():
        if splice_motif_dict[motif_do] == motif_ac:
            return True, 3, 1
        else:
            return (False, 3, 0) if motif_required else (True, 3, 0)


def cigar_validity(cigar_str: str) -> str:
    """Merge the first two OR last two 'same' operations in the CIGAR string generated by BLAT.

    :param cigar_str: BLAT generated cigarstring from 'softclipped_seq2SA_tag'
    :type cigar_str: str
    :return: valid cigarstring
    :rtype: str

    :Example:

    >>> cigarstring =  '1S2S5M3S2S'
    >>> cigar_validity(cigarstring)
    '3S5M5S'
    """
    pattern = re.compile(r"((?P<length>\d+)(?P<op>\D))")
    items_list = pattern.findall(cigar_str)
    stack = [items_list[0]]
    for item in items_list[1:]:  # [('1S', '1', 'S'),..]
        last_item = stack[-1]
        if last_item[2] == item[2]:
            length = int(last_item[1]) + int(item[1])
            stack[-1] = (f"{length}{last_item[2]}", f"{length}", last_item[2])
        else:
            stack.append(item)

    return "".join(i[0] for i in stack)


def blat2chimeric_alignment(
    in_seq: str,
    read_length: int,
    read_strand: str,
    read_mode: int,
    blat: Any,
    mapq_cutoff: int,
    max_allowed_nm: int,
    blat_ident_pct_cutoff: float = 0.95,
) -> str:
    """Create chimeric alignments from the alignments.

    the alignments which has a long softclipped segment but without SA tag.

    :param blat:
    :param in_seq: softclipped segment of the aligned read
    :param read_length: the length of the aligned read
    :param read_strand: the strand of the aligned read (-/+)
    :param read_mode: mode of the aligned read (1/2)
    :param mapq_cutoff: MAPQ cutoff
    :param max_allowed_nm: mismatches cutoff used for discarding supplementary alignments
    :param blat_ident_pct_cutoff: BLAT HSP identity cutoff
    :return: putative supplementary alignment of the alignment which is ready for put in the SA tag
    """
    chimeric_aln_str = ""
    in_seq_len = len(in_seq)

    top_hsp, __mapq = blat.fetch_mapq(in_seq, blat_ident_pct_cutoff)
    if top_hsp is None:
        return ""
    if (
        top_hsp.ident_pct / 100 >= blat_ident_pct_cutoff
        and top_hsp.query_span / in_seq_len >= blat_ident_pct_cutoff
    ):
        chrom_sa, pos_sa, strand_sa, cigar_sa_partial, nm_sa = blat.psl2sam(
            top_hsp, in_seq_len
        )
        if read_strand == strand_sa:
            # same strand: different reads mode
            # MS(1) ~ SM(2) or SM(2) ~ MS(1)
            # SM
            cigar_sa = (
                f"{read_length - in_seq_len}S{cigar_sa_partial}"
                if read_mode == 1
                else f"{cigar_sa_partial}{read_length - in_seq_len}S"
            )  # MS
        else:
            # opposite strand: same reads mode
            # MS(1) ~ MS(1) or SM(2) ~ SM(2)
            # MS
            cigar_sa = (
                f"{cigar_sa_partial}{read_length - in_seq_len}S"
                if read_mode == 1
                else f"{read_length - in_seq_len}S{cigar_sa_partial}"
            )  # SM

        valid_cigar_sa = cigar_validity(cigar_sa)

        if __mapq >= mapq_cutoff and int(nm_sa) < max_allowed_nm:
            return "{},{},{},{},{},{};".format(
                chrom_sa, pos_sa, strand_sa, valid_cigar_sa, __mapq, nm_sa
            )

    return chimeric_aln_str


def obtain_insertion_surrouding_cigarstrings(
    cigar_str: str,
    insertion_length: int,
) -> tuple[str, str]:
    """Obtain the upstream and downstream CIGAR strings of the targeted insertion."""
    insertion_str = f"{insertion_length}I"
    try:
        ins_idx = cigar_str.index(insertion_str)
    except ValueError:
        surrounding_cigar = "", ""
    else:
        surrounding_cigar = (
            cigar_str[:ins_idx],
            cigar_str[ins_idx + len(insertion_str) :],
        )
    return surrounding_cigar


def obtain_read_segment_length_from_cigar_string(cigar_str: str) -> int:
    """Obtain the read segment length providing CIGAR string."""
    parse_result = cppext.parseCigar(cigar_str)
    cigartuples = parse_result.cigartuples
    read_seg_len = 0

    for idx in range(0, len(cigartuples), 2):
        op_code = cigartuples[idx]
        _len = cigartuples[idx + 1]

        if op_code in {0, 1, 4}:  # M, I or S
            read_seg_len = read_seg_len + _len
    return read_seg_len


def insertion2chimeric_alignment(
    read: pysam.libcalignedsegment.AlignedSegment,
    insertion_ref_pos: int,
    insertion_seq: str,
    read_length: int,
    read_strand: str,
    max_allowed_nm: int,
    blat: Any,
    blat_ident_pct_cutoff: float = 0.95,
    top: int = 3,
    align_len_threshold: int = 50,
) -> tuple[str, str]:
    """Create chimeric alignments from the alignment with long insertion.

    the alignment which has a long insertion segment but without SA tag.

    :param read: pysam.libcalignedsegment.AlignedSegment
    :param insertion_ref_pos: the insertion reference start position from the original read
    :param insertion_seq: insertion segment sequence
    :param read_length: the length of the aligned read
    :param read_strand: the strand of the aligned read (-/+)
    :param max_allowed_nm: the maximum allowed NM
    :param blat: instantiated blat object
    :param blat_ident_pct_cutoff: BLAT HSP identity cutoff
    :param top: the top number of the alignments
    :param align_len_threshold: the threshold of the insertion sequence length
    :return: putative supplementary alignment of the alignment which is ready for put in the SA tag
    """
    read_strand = "-" if read.is_reverse else "+"
    if read_strand == "-":
        insertion_seq = reverse_complement(insertion_seq)

    insertion_seq_len = len(insertion_seq)

    original_cigar_str = read.cigarstring
    original_ref_start = read.reference_start
    nm_read = read.get_tag("NM")

    left_cigar_str, right_cigar_str = obtain_insertion_surrouding_cigarstrings(
        original_cigar_str, insertion_seq_len
    )
    left_cigar_read_seg_len = obtain_read_segment_length_from_cigar_string(
        left_cigar_str
    )
    right_cigar_read_seg_len = obtain_read_segment_length_from_cigar_string(
        right_cigar_str
    )

    flag, insertion_info = blat.query_insertion(
        insert_seq=insertion_seq,
        threshold_identity=blat_ident_pct_cutoff,
        top=top,
        align_len_threshold=align_len_threshold,
    )

    chimeric_aln_str = ""
    primary_aln_cigarstring = ""
    if flag:
        # BLAT unique HSP
        chrom_blat = insertion_info.chrom
        strand_blat = insertion_info.strand
        pos_blat = insertion_info.ref_start
        ref_end_blat = insertion_info.ref_end
        cigar_blat = insertion_info.cigarstring
        mapq_blat = insertion_info.mapq
        nm_blat = insertion_info.nm
        # Insertion sequence BLAT HSP and original insertion have the same chrom and reference start position
        if (
            chrom_blat == read.reference_name
            and strand_blat == read_strand
            and (
                abs(pos_blat - insertion_ref_pos) <= 10
                and ref_end_blat <= read.reference_end
            )
            or (
                abs(ref_end_blat - insertion_ref_pos) <= 10
                and pos_blat >= read.reference_start
            )
        ):
            # SM
            cigar_ra = f"{read_length - right_cigar_read_seg_len}S{right_cigar_str}"
            # MS
            cigar_sa = f"{left_cigar_str}{cigar_blat}{read_length - left_cigar_read_seg_len - insertion_seq_len}S"

            valid_cigar_ra = cigar_validity(cigar_ra)
            valid_cigar_sa = cigar_validity(cigar_sa)

            nm_sa = nm_read - insertion_seq_len + nm_blat
            if nm_sa < max_allowed_nm:
                chimeric_aln_str = (
                    f"{chrom_blat},{original_ref_start + 1},"
                    f"{read_strand},{valid_cigar_sa},{mapq_blat},{nm_sa};"
                )
                primary_aln_cigarstring = valid_cigar_ra

    return primary_aln_cigarstring, chimeric_aln_str


def strand_mode_checker(strand1: str, strand2: str, mode1: int, mode2: int) -> bool:
    """Check if the two strands are compatible with the two modes."""
    return (strand1 == strand2 and mode1 != mode2) or (
        strand1 != strand2 and mode1 == mode2
    )


def softclipped_length_and_event_size_checker(
    read, mode, event_size, bp_region_seq_len
) -> bool:
    """When read length > predicted tandem duplication size.

    check whether the softclipped length is less than the inferred event size

    :param read: a chimeric read
    :param mode: mode for the chimeric read
    :param event_size: event size inferred from 'query_offset - target_offset'

    :type read : Read
    :type mode: int
    :type event_size: int
    :return: whether event_size > softclipped_length (If it is True, it will be a TDUP event)
    :rtype: bool
    """
    return (
        read.lt_soft_len < event_size + bp_region_seq_len
        if mode == 2
        else read.rt_soft_len < event_size + bp_region_seq_len
    )


def obtain_bp_region_seq(read, mode, bp_region_seq_len, genome_fasta) -> str:
    """Obtain breakpoint region sequence from read.

    :param bp_region_seq_len: the length of the breakpoint region sequence
    :param read:  the chimeirc read
    :param mode: mode for the chimeirc read
    :param genome_fasta: reference genome (pyfaidx.Fasta object)
    :type mode: int
    :return: (putative insertion/microhomology sequence from the read; + means insertion,
        - means microhomology, mode)

    ..note:
        * inserted sequence:
          S-----SM---M    M---MS-----S
          SSSSSXXMMMMM    MMMMMXXSSSSS
        * microhomology:
          S-----SM---M    M---MS-----S
          SSSSSSSXXMMM    MMMXXSSSSSSS
    """
    read_seq = read.query_sequence
    chrom = read.chrom
    bp_region_seq = ""
    # inserted sequence
    if bp_region_seq_len > 0:
        if mode == 2:  # SM
            bp_region_seq = read_seq[: read.lt_soft_len][-bp_region_seq_len:]
        elif mode == 1:  # MS
            bp_region_seq = read_seq[-read.rt_soft_len :][:bp_region_seq_len]
        bp_region_seq = "+" + bp_region_seq
    # microhomology
    elif bp_region_seq_len < 0:
        if mode == 2:  # SM
            bp_region_seq = genome_fasta[chrom][
                read.ref_start : read.ref_start - bp_region_seq_len
            ].seq

        elif mode == 1:  # MS
            bp_region_seq = genome_fasta[chrom][
                read.ref_end + bp_region_seq_len : read.ref_end
            ].seq

        bp_region_seq = "-" + bp_region_seq

    return bp_region_seq


noreturn = "NA", 0, 0, (), (), (), (), (), []  # type: ignore


def same_chrom_same_strand_mode21_handler(
    read_lt,
    read_rt,
    lt_mode,
    rt_mode,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    logger,
    microinsertion_cutoff=20,
    is_reverse=False,
):
    """Same chrom same strand mode 21 handler."""
    lt_chrm = read_lt.chrom
    lt_exons, lt_introns = read_lt.get_exons_and_introns()
    rt_exons, rt_introns = read_rt.get_exons_and_introns()

    if lt_mode == 2 and rt_mode == 1:
        target_start = read_rt.ref_start
        target_end = read_lt.ref_end
        target_offset = target_end - target_start

        bp_region_seq_len = (
            read_lt.query_length
            - read_lt.rt_soft_len
            - read_rt.lt_soft_len
            - read_lt.read_match_size
            - read_rt.read_match_size
        )

        logger.trace(f"{bp_region_seq_len=}")

        if bp_region_seq_len > microinsertion_cutoff:
            logger.trace(f"{bp_region_seq_len=} > {microinsertion_cutoff=}")
            return noreturn

        if bp_region_seq_len > 0:
            query_offset = read_lt.reference_match_size + read_rt.reference_match_size
        else:
            query_offset = (
                read_lt.reference_match_size
                + read_rt.reference_match_size
                + bp_region_seq_len
            )

        lt_bp_seq = obtain_bp_region_seq(
            read_lt, lt_mode, bp_region_seq_len, genome_fasta
        )
        rt_bp_seq = obtain_bp_region_seq(
            read_rt, rt_mode, bp_region_seq_len, genome_fasta
        )

        evt_size = query_offset - target_offset

        logger.trace(f"{evt_size=}, {query_offset=}")
        if evt_size <= 0:  # deletion
            del_start = read_rt.ref_end
            del_end = del_start + abs(evt_size)
            _, _anno, _can = splicing_confirmation(
                lt_chrm,
                del_start,
                read_rt.strand,
                rt_mode,
                lt_chrm,
                del_end,
                read_lt.strand,
                lt_mode,
                splice_bin,
                genome_fasta,
                cvg,
                motif_required,
            )
            _genes = gene_annotation(lt_chrm, del_start, lt_chrm, del_end, gene_iv)
            # 1 => 2
            if is_reverse:
                if _anno == 1:
                    _anno = 2
                elif _anno == 2:
                    _anno = 1
                _genes = _genes[::-1]
            return (
                "DEL",
                _anno,
                _can,
                (
                    f"{lt_chrm}:{del_start}",
                    f"{lt_chrm}:{del_end}",
                    1,
                    2,
                ),
                (read_rt.ref_start, read_rt.ref_end, rt_exons),
                (read_lt.ref_start, read_lt.ref_end, lt_exons),
                (rt_bp_seq, lt_bp_seq),
                (read_rt.strand, read_lt.strand),
                [*_genes],
            )

        # reads length < tandem duplication size
        elif evt_size >= query_offset:  # large tandem duplication
            chrm_start = lt_chrm
            junc_start = read_lt.ref_start
            chrm_end = lt_chrm
            junc_end = junc_start + evt_size
            _nls, _anno, _can = splicing_confirmation(
                chrm_start,
                junc_start,
                read_lt.strand,
                lt_mode,
                chrm_end,
                junc_end,
                read_rt.strand,
                rt_mode,
                splice_bin,
                genome_fasta,
                cvg,
                motif_required,
            )
            _genes = gene_annotation(
                chrm_start, junc_start, chrm_end, junc_end, gene_iv
            )
            if is_reverse:
                if _anno == 1:
                    _anno = 2
                elif _anno == 2:
                    _anno = 1
                _genes = _genes[::-1]
            if _nls:
                if not is_reverse:
                    return (
                        "TDUP",
                        _anno,
                        _can,
                        (
                            f"{lt_chrm}:{junc_start}",
                            f"{lt_chrm}:{junc_end}",
                            2,
                            1,
                        ),
                        (read_lt.ref_start, read_lt.ref_end, lt_exons),
                        (read_rt.ref_start, read_rt.ref_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (read_lt.strand, read_rt.strand),
                        [*_genes],
                    )
                return (
                    "TDUP",
                    _anno,
                    _can,
                    (
                        f"{lt_chrm}:{junc_end}",
                        f"{lt_chrm}:{junc_start}",
                        1,
                        2,
                    ),
                    (read_rt.ref_start, read_rt.ref_end, rt_exons),
                    (read_lt.ref_start, read_lt.ref_end, lt_exons),
                    (rt_bp_seq, lt_bp_seq),
                    (read_rt.strand, read_lt.strand),
                    [*_genes],
                )
            else:
                return noreturn
        # read length > tandem duplication size
        else:
            # softclipped length < tandem duplication size (check chimeric read [SM])
            if softclipped_length_and_event_size_checker(
                read_lt, lt_mode, evt_size, bp_region_seq_len
            ):
                logger.trace("softclipped length < event size: TDUP")
                is_dup = True
            # softclipped length >= tandem duplication size
            # TDUP; Novel Insertion feature: evt_size=0 and bp_region_seq_len>0
            else:
                is_dup = True
                logger.trace("softclipped length >= event size: TDUP")
            if is_dup:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start
                chrm_end = lt_chrm
                junc_end = junc_start + evt_size
                _nls, _anno, _can = splicing_confirmation(
                    chrm_start,
                    junc_start,
                    read_lt.strand,
                    lt_mode,
                    chrm_end,
                    junc_end,
                    read_rt.strand,
                    rt_mode,
                    splice_bin,
                    genome_fasta,
                    cvg,
                    motif_required,
                )
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                if is_reverse:
                    if _anno == 1:
                        _anno = 2
                    elif _anno == 2:
                        _anno = 1
                    _genes = _genes[::-1]
                # 2 => 1
                if _nls:
                    if not is_reverse:
                        return (
                            "TDUP",
                            _anno,
                            _can,
                            (
                                f"{lt_chrm}:{junc_start}",
                                f"{lt_chrm}:{junc_end}",
                                2,
                                1,
                            ),
                            (read_lt.ref_start, read_lt.ref_end, lt_exons),
                            (read_rt.ref_start, read_rt.ref_end, rt_exons),
                            (lt_bp_seq, rt_bp_seq),
                            (read_lt.strand, read_rt.strand),
                            [*_genes],
                        )
                    return (
                        "TDUP",
                        _anno,
                        _can,
                        (
                            f"{lt_chrm}:{junc_end}",
                            f"{lt_chrm}:{junc_start}",
                            1,
                            2,
                        ),
                        (read_rt.ref_start, read_rt.ref_end, rt_exons),
                        (read_lt.ref_start, read_lt.ref_end, lt_exons),
                        (rt_bp_seq, lt_bp_seq),
                        (read_rt.strand, read_lt.strand),
                        [*_genes],
                    )
                else:
                    return noreturn


def same_chrom_same_strand_handler(
    read_lt,
    read_rt,
    lt_mode,
    rt_mode,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    logger,
    microinsertion_cutoff=20,
):
    """Handler for same chrom and same strand."""
    logger.trace("same_chrom_same_strand_handler takes over the task.")
    if lt_mode == 2 and rt_mode == 1:
        return same_chrom_same_strand_mode21_handler(
            read_lt,
            read_rt,
            lt_mode,
            rt_mode,
            splice_bin,
            genome_fasta,
            cvg,
            gene_iv,
            motif_required,
            logger,
            microinsertion_cutoff,
        )
    elif lt_mode == 1 and rt_mode == 2:
        return same_chrom_same_strand_mode21_handler(
            read_rt,
            read_lt,
            rt_mode,
            lt_mode,
            splice_bin,
            genome_fasta,
            cvg,
            gene_iv,
            motif_required,
            logger,
            microinsertion_cutoff,
            is_reverse=True,
        )


def same_chrom_diff_strand_handler(
    read_lt,
    read_rt,
    lt_mode,
    rt_mode,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    logger,
    microinsertion_cutoff=20,
):
    """Handler for same chrom and different strand."""
    logger.trace("same_chrom_diff_strand_handler takes over the task.")
    # lt_mode must be equal to rt_mode
    if lt_mode != rt_mode:
        logger.warning(
            f"ModesNotEqualError: {read_lt.query_name}, read_lt:{read_lt} read_rt:{read_rt}"
        )
        return noreturn

    lt_chrm = read_lt.chrom
    lt_exons, lt_introns = read_lt.get_exons_and_introns()
    rt_exons, rt_introns = read_rt.get_exons_and_introns()

    same_mode = lt_mode
    if same_mode == 1:
        ra_bp = read_lt.ref_start + read_lt.reference_match_size
        sa_bp = read_rt.ref_start + read_rt.reference_match_size
        bp_region_seq_len = (
            read_lt.query_length
            - read_lt.lt_soft_len
            - read_rt.lt_soft_len
            - read_lt.read_match_size
            - read_rt.read_match_size
        )
    elif same_mode == 2:
        ra_bp = read_lt.ref_start
        sa_bp = read_rt.ref_start
        bp_region_seq_len = (
            read_lt.query_length
            - read_lt.rt_soft_len
            - read_rt.rt_soft_len
            - read_lt.read_match_size
            - read_rt.read_match_size
        )
    logger.trace(f"{bp_region_seq_len=}")

    if bp_region_seq_len > microinsertion_cutoff:
        logger.trace(f"{bp_region_seq_len=} > {microinsertion_cutoff=}")
        return noreturn

    if ra_bp == sa_bp:  # inverted duplication (IDUP)
        # allow one read with noncanonical splice site for IDUP
        if not read_lt.splice_site_checker(
            genome_fasta
        ) and not read_rt.splice_site_checker(genome_fasta):
            logger.debug(
                f"Splice site checking[IDUP]: {read_lt.query_name=}, "
                f"{read_lt.cigarstring=}, {read_rt.cigarstring=}"
            )
            return noreturn
        chrm_start = lt_chrm
        junc_start = ra_bp
        chrm_end = lt_chrm
        junc_end = ra_bp
        _nls, _anno, _can = splicing_confirmation(
            chrm_start,
            junc_start,
            read_lt.strand,
            same_mode,
            chrm_end,
            junc_end,
            read_rt.strand,
            same_mode,
            splice_bin,
            genome_fasta,
            cvg,
            motif_required,
        )
        strands = (read_lt.strand, read_rt.strand)
        lt_start_end_exons = (read_lt.ref_start, read_lt.ref_end, lt_exons)
        rt_start_end_exons = (read_rt.ref_start, read_rt.ref_end, rt_exons)
        lt_bp_seq = obtain_bp_region_seq(
            read_lt, lt_mode, bp_region_seq_len, genome_fasta
        )
        rt_bp_seq = obtain_bp_region_seq(
            read_rt, rt_mode, bp_region_seq_len, genome_fasta
        )

        if _nls:
            _genes = gene_annotation(
                chrm_start, junc_start, chrm_end, junc_end, gene_iv
            )

            return (
                "IDUP",
                _anno,
                _can,
                (
                    f"{lt_chrm}:{junc_start}",
                    f"{lt_chrm}:{junc_end}",
                    same_mode,
                    same_mode,
                ),
                lt_start_end_exons,
                rt_start_end_exons,
                (lt_bp_seq, rt_bp_seq),
                (*strands,),
                [*_genes],
            )
        else:
            return noreturn

    else:  # conventional INV
        # If using noncanonical splice site, return NA
        if not read_lt.splice_site_checker(
            genome_fasta
        ) or not read_rt.splice_site_checker(genome_fasta):
            logger.debug(
                f"Splice site checking[INV]: {read_lt.query_name=}, {read_lt.cigarstring=}, {read_rt.cigarstring=}"
            )
            return noreturn

        chrm_start = lt_chrm
        junc_start = min(ra_bp, sa_bp)
        chrm_end = lt_chrm
        junc_end = junc_start + abs(ra_bp - sa_bp)

        if junc_start == ra_bp:
            strands = (read_lt.strand, read_rt.strand)
            lt_start_end_exons = (read_lt.ref_start, read_lt.ref_end, lt_exons)
            rt_start_end_exons = (read_rt.ref_start, read_rt.ref_end, rt_exons)
            lt_bp_seq = obtain_bp_region_seq(
                read_lt, lt_mode, bp_region_seq_len, genome_fasta
            )
            rt_bp_seq = obtain_bp_region_seq(
                read_rt, rt_mode, bp_region_seq_len, genome_fasta
            )
        elif junc_start == sa_bp:
            strands = (read_rt.strand, read_lt.strand)
            lt_start_end_exons = (read_rt.ref_start, read_rt.ref_end, rt_exons)
            rt_start_end_exons = (read_lt.ref_start, read_lt.ref_end, lt_exons)
            lt_bp_seq = obtain_bp_region_seq(
                read_rt, rt_mode, bp_region_seq_len, genome_fasta
            )
            rt_bp_seq = obtain_bp_region_seq(
                read_lt, lt_mode, bp_region_seq_len, genome_fasta
            )

        _nls, _anno, _can = splicing_confirmation(
            chrm_start,
            junc_start,
            strands[0],
            same_mode,
            chrm_end,
            junc_end,
            strands[1],
            same_mode,
            splice_bin,
            genome_fasta,
            cvg,
            motif_required,
        )
        _genes = gene_annotation(chrm_start, junc_start, chrm_end, junc_end, gene_iv)
        if _nls:
            return (
                "INV",
                _anno,
                _can,
                (
                    f"{lt_chrm}:{junc_start}",
                    f"{lt_chrm}:{junc_end}",
                    same_mode,
                    same_mode,
                ),
                lt_start_end_exons,
                rt_start_end_exons,
                (lt_bp_seq, rt_bp_seq),
                (*strands,),
                [*_genes],
            )
        else:
            return noreturn


def diff_chrom_same_strand_mode21_handler(
    read_lt,
    read_rt,
    lt_mode,
    rt_mode,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    logger,
    microinsertion_cutoff=20,
    is_reverse=False,
):
    """Different chrom same stand mode 21 handler."""
    lt_exons, lt_introns = read_lt.get_exons_and_introns()
    rt_exons, rt_introns = read_rt.get_exons_and_introns()

    chrm_start = read_lt.chrom
    junc_start = read_lt.ref_start
    chrm_end = read_rt.chrom
    junc_end = read_rt.ref_start + read_rt.reference_match_size

    bp_region_seq_len = (
        read_lt.query_length
        - read_lt.rt_soft_len
        - read_rt.lt_soft_len
        - read_lt.read_match_size
        - read_rt.read_match_size
    )
    logger.trace(f"{bp_region_seq_len=}")

    if bp_region_seq_len > microinsertion_cutoff:
        logger.trace(f"{bp_region_seq_len=} > {microinsertion_cutoff=}")
        return noreturn

    lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len, genome_fasta)
    rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len, genome_fasta)
    _nls, _anno, _can = splicing_confirmation(
        chrm_start,
        junc_start,
        read_lt.strand,
        lt_mode,
        chrm_end,
        junc_end,
        read_rt.strand,
        rt_mode,
        splice_bin,
        genome_fasta,
        cvg,
        motif_required,
    )
    _genes = gene_annotation(chrm_start, junc_start, chrm_end, junc_end, gene_iv)
    if is_reverse:
        if _anno == 1:
            _anno = 2
        elif _anno == 2:
            _anno = 1
        _genes = _genes[::-1]

    if _nls:
        if not is_reverse:
            return (
                "TRA",
                _anno,
                _can,
                (f"{chrm_start}:{junc_start}", f"{chrm_end}:{junc_end}", 2, 1),
                (read_lt.ref_start, read_lt.ref_end, lt_exons),
                (read_rt.ref_start, read_rt.ref_end, rt_exons),
                (lt_bp_seq, rt_bp_seq),
                (read_lt.strand, read_rt.strand),
                [*_genes],
            )
        return (
            "TRA",
            _anno,
            _can,
            (f"{chrm_end}:{junc_end}", f"{chrm_start}:{junc_start}", 1, 2),
            (read_rt.ref_start, read_rt.ref_end, rt_exons),
            (read_lt.ref_start, read_lt.ref_end, lt_exons),
            (rt_bp_seq, lt_bp_seq),
            (read_rt.strand, read_lt.strand),
            [*_genes],
        )

    return noreturn


def diff_chrom_same_strand_handler(
    read_lt,
    read_rt,
    lt_mode,
    rt_mode,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    logger,
    microinsertion_cutoff=20,
):
    """Diff chrom same strand handler."""
    logger.trace("diff_chrom_same_strand_handler takes over the task.")
    if lt_mode == 2 and rt_mode == 1:
        return diff_chrom_same_strand_mode21_handler(
            read_lt,
            read_rt,
            lt_mode,
            rt_mode,
            splice_bin,
            genome_fasta,
            cvg,
            gene_iv,
            motif_required,
            logger,
            microinsertion_cutoff,
        )
    elif lt_mode == 1 and rt_mode == 2:
        return diff_chrom_same_strand_mode21_handler(
            read_rt,
            read_lt,
            rt_mode,
            lt_mode,
            splice_bin,
            genome_fasta,
            cvg,
            gene_iv,
            motif_required,
            logger,
            microinsertion_cutoff,
            is_reverse=True,
        )


def diff_chrom_diff_strand_handler(
    read_lt,
    read_rt,
    lt_mode,
    rt_mode,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    logger,
    microinsertion_cutoff=20,
):
    """Diff chrom different strand handler."""
    logger.trace("diff_chrom_diff_strand_handler takes over the task.")
    # lt_mode must be equal to rt_mode
    if lt_mode != rt_mode:
        raise ModesNotEqualError(
            f"read_lt:{read_lt.query_name} read_rt:{read_rt.query_name}"
        )

    lt_exons, lt_introns = read_lt.get_exons_and_introns()
    rt_exons, rt_introns = read_rt.get_exons_and_introns()
    same_mode = lt_mode

    if same_mode == 1:
        chrm_start = read_lt.chrom
        junc_start = read_lt.ref_start + read_lt.reference_match_size
        chrm_end = read_rt.chrom
        junc_end = read_rt.ref_start + read_rt.reference_match_size
        bp_region_seq_len = (
            read_lt.query_length
            - read_lt.lt_soft_len
            - read_rt.lt_soft_len
            - read_lt.read_match_size
            - read_rt.read_match_size
        )
    elif same_mode == 2:
        chrm_start = read_lt.chrom
        junc_start = read_lt.ref_start
        chrm_end = read_rt.chrom
        junc_end = read_rt.ref_start
        bp_region_seq_len = (
            read_lt.query_length
            - read_lt.rt_soft_len
            - read_rt.rt_soft_len
            - read_lt.read_match_size
            - read_rt.read_match_size
        )

    logger.trace(f"{bp_region_seq_len=}")

    if bp_region_seq_len > microinsertion_cutoff:
        logger.trace(f"{bp_region_seq_len=} > {microinsertion_cutoff=}")
        return noreturn

    lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len, genome_fasta)
    rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len, genome_fasta)
    _nls, _anno, _can = splicing_confirmation(
        chrm_start,
        junc_start,
        read_lt.strand,
        same_mode,
        chrm_end,
        junc_end,
        read_rt.strand,
        same_mode,
        splice_bin,
        genome_fasta,
        cvg,
        motif_required,
    )
    _genes = gene_annotation(chrm_start, junc_start, chrm_end, junc_end, gene_iv)
    if _nls:
        return (
            "TRA",
            _anno,
            _can,
            (
                f"{chrm_start}:{junc_start}",
                f"{chrm_end}:{junc_end}",
                same_mode,
                same_mode,
            ),
            (read_lt.ref_start, read_lt.ref_end, lt_exons),
            (read_rt.ref_start, read_rt.ref_end, rt_exons),
            (lt_bp_seq, rt_bp_seq),
            (read_lt.strand, read_rt.strand),
            [*_genes],
        )
    return noreturn


def obtain_variants_stats(
    cigar_str: str, md_tag: str, indel_len_cutoff: int = 4
) -> tuple[int, float, float]:
    """Obtain variants stats from read matched part.

    :param cigar_str: CIGAR string
    :param md_tag: MD tag
    :param indel_len_cutoff: INDEL length threshold
    :return: number of substitutions, fraction of long insertions and fraction of long deletions.

    .. note::
        'A': 65
        'Z': 90
        '^': 94
        '0': 48
        https://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag
    """
    parsed_cigar_result = cppext.parseCigar(cigar_str)
    cigartuples_without_soft: list[int] = parsed_cigar_result.cigartuples_without_soft

    del_num = 0
    ins_num = 0
    del_outlier_num = 0
    ins_outlier_num = 0
    dels_len_total = 0

    for idx in range(0, len(cigartuples_without_soft), 2):
        op_code = cigartuples_without_soft[idx]
        _len = cigartuples_without_soft[idx + 1]
        if op_code == 2:  # DEL
            del_num += 1
            dels_len_total += _len
            if _len >= indel_len_cutoff:
                del_outlier_num += 1
        elif op_code == 1:  # INS
            ins_num += 1
            if _len >= indel_len_cutoff:
                ins_outlier_num += 1

    sum_of_subs_dels = 0
    for _letter in md_tag:
        if ord(_letter) >= 65 and ord(_letter) <= 90:
            sum_of_subs_dels += 1

    num_of_subs = sum_of_subs_dels - dels_len_total
    ins_fraction = 0 if ins_num == 0 else ins_outlier_num / ins_num
    del_fraction = 0 if del_num == 0 else del_outlier_num / del_num

    return num_of_subs, ins_fraction, del_fraction


def get_transcriptome_length(species: str) -> int:
    """Obtain transcriptome size from BAM header.

    :param species:  species name
    :return: reference transcriptome size.
    """
    blat_path = Path(sys.modules[__PACKAGE_NAME__].__file__).parent / "blat"

    path = blat_path / "transcriptome_length.yaml"
    transcript_length_dict = yaml.safe_load(path.open())
    return transcript_length_dict["species"][species]
