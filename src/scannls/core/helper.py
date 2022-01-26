# !/usr/bin/env python
"""Helper functions."""
import re
from collections import defaultdict

from scannls import ModesNotEqualError

__funcs__ = {
    "extract_splice_sites",
    "gene_annotation",
    "splicing_confirmation",
    "update_breakpoints",
    "cigar_validity",
    "blat2chimeric_alignment",
    "same_chrom_same_strand_mode21_handler",
    "same_chrom_same_strand_handler",
    "same_chrom_diff_strand_handler",
    "diff_chrom_same_strand_mode21_handler",
    "diff_chrom_same_strand_handler",
    "diff_chrom_diff_strand_handler",
    "softclipped_length_and_event_size_checker",
    "obtain_bp_region_seq",
}

from typing import Tuple, List, Dict, Any

import HTSeq  # type: ignore
import pyfaidx  # type: ignore


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
        biotype = feature.attr["gene_type"]
        gene_name = feature.attr["gene_name"]
        if biotype == "protein_coding":
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
) -> Tuple[str, str]:
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
        print(chrm1, pos1)
    try:
        gene2 = "&".join(list(gene_iv[HTSeq.GenomicPosition(chrm2, pos2)]))
    except IndexError:
        gene2 = ""
    except TypeError:
        print(chrm2, pos2)

    if not gene1:
        gene1 = "INTERGENIC"
    if not gene2:
        gene2 = "INTERGENIC"
    return gene1, gene2


def splicing_confirmation(
    chrm1: str,
    pos1: int,
    chrm2: str,
    pos2: int,
    splice_bin: int,
    genome_fasta: pyfaidx.Fasta,
    cvg: HTSeq.GenomicArrayOfSets,
    strand_changed: bool,
    motif_required: bool = True,
) -> Tuple[bool, int, int]:
    """Judge whether the breakpoints are NLS events or not.

    if motif_required is ON: it will only report NLS events with 'canonical
    splice sites'; otherwise: it will report NLS events whatever the splice
    sites they used

    :param chrm1: chromosome for breakpoint1
    :param chrm2: chromosome for breakpoint2
    :param pos1: position for breakpoint1
    :param pos2: position for breakpoint2
    :param splice_bin: bin size for splice sites searching
    :param genome_fasta: reference genome (pyfaidx.Fasta object)
    :param cvg: splice site annotations (HTSeq.GenomicArrayOfSets)
    :param strand_changed: whether breakpoint1 and breakpoint2 use the same
        strand or not
    :param motif_required: canonical splice sites required;
           if True: considering canonical splice sites only;
           else: considering canonical and noncanonical splice sites both
    :return: report/not report, overlapping boundary in bits, canonical
        splice site/noncanonical splice site
    :rtype: tuple

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

    def canonical_site_finder(in_seq: str) -> List[str]:
        """Find canonical splice sites in the input sequence.

        :param in_seq: input sequence (usually sequence nearby the breakpoints)
        :type in_seq: str
        :return: a list of canonical splice sites in the input sequence,
            it can be a empty list
        :rtype: list
        """
        candidate_sites = ["GT", "AG", "CT", "AC"]
        matches = (i in in_seq for i in candidate_sites)
        hit_sites = [j for i, j in zip(matches, candidate_sites) if i]
        return hit_sites

    def splice_paired_checker(
        hit_sites: List[str], pair_seq: str, splice_motif_dict: Dict[str, str]
    ) -> bool:
        """Find canonical splice sites in the input sequence.

        :param hit_sites: canonical splice site at one end
        :param pair_seq: sequence at the other pair end
        :param splice_motif_dict: paired splice sites (same strand or different strand)
        :return: canonical splice sites are paired or not
        :rtype: bool
        """
        paired = False
        for site in hit_sites:
            if site in splice_motif_dict and splice_motif_dict[site] in pair_seq:
                paired = True
                break
        return paired

    if strand_changed:
        splice_motif_dict = {"GT": "CT", "AG": "AC", "CT": "GT", "AC": "AG"}
    else:
        splice_motif_dict = {"GT": "AG", "AG": "GT", "CT": "AC", "AC": "CT"}
    try:
        junc1 = list(cvg[HTSeq.GenomicPosition(chrm1, pos1)])[0]
    except IndexError:
        junc1 = ""
    try:
        junc2 = list(cvg[HTSeq.GenomicPosition(chrm2, pos2)])[0]
    except IndexError:
        junc2 = ""
    # Non-annotated coding exon boundary
    if junc1 not in splice_motif_dict and junc2 not in splice_motif_dict:
        junc_seq1 = genome_fasta[chrm1][pos1 - splice_bin : pos1 + splice_bin].seq
        junc_seq2 = genome_fasta[chrm2][pos2 - splice_bin : pos2 + splice_bin].seq
        _junc1 = canonical_site_finder(junc_seq1)
        _junc2 = canonical_site_finder(junc_seq2)
        if splice_paired_checker(
            _junc1, junc_seq2, splice_motif_dict
        ) or splice_paired_checker(_junc2, junc_seq1, splice_motif_dict):
            return True, 0, 1
        else:
            return (False, 0, 0) if motif_required else (True, 0, 0)
    # pos1 in annotated coding exon boundary, pos2 not.
    elif junc1 in splice_motif_dict and junc2 not in splice_motif_dict:
        junc_seq = genome_fasta[chrm2][pos2 - splice_bin : pos2 + splice_bin].seq
        if splice_motif_dict[junc1] in junc_seq:
            return True, 2, 1
        else:
            return (False, 2, 0) if motif_required else (True, 2, 0)

    # pos2 in annotated coding exon boundary, pos1 not.
    elif junc1 not in splice_motif_dict and junc2 in splice_motif_dict:
        junc_seq = genome_fasta[chrm1][pos1 - splice_bin : pos1 + splice_bin].seq
        if splice_motif_dict[junc2] in junc_seq:
            return True, 1, 1
        else:
            return (False, 1, 0) if motif_required else (True, 1, 0)

    # pos1 and pos2 both in annotated coding exon boundary
    # junc1 in splice_motif_dict and junc2 in splice_motif_dict
    else:
        if splice_motif_dict[junc1] == junc2:
            return True, 3, 1
        else:
            return (False, 3, 0) if motif_required else (True, 3, 0)


def cigar_validity(cigar_str: str) -> str:
    """Merge the first two OR last two 'same' operations in the CIGAR string generated by BLAT.

    :param cigar_str: BLAT generated cigarstring from 'softclipped_seq2SA_tag'
    :type cigar_str: str
    :return: valid cigarstring
    :rtype: str

    ..note ::
        assert cigar_validity('45S50S100M1S') == '95S100M1S'
    """
    cigartuple = list(map(list, re.findall(r"(\d+)(\w)", cigar_str)))
    # first two operations are the same
    if cigartuple[0][1] == cigartuple[1][1]:
        cigartuple[1][0] = str(
            int(cigartuple[0][0]) + int(cigartuple[1][0])  # type: ignore
        )
        del cigartuple[0]

    # last two operations are the same
    elif cigartuple[-1][1] == cigartuple[-2][1]:
        cigartuple[-2][0] = str(
            int(cigartuple[-1][0]) + int(cigartuple[-2][0])  # type: ignore
        )
        del cigartuple[-1]

    valid_cigar = ""
    for len_str, op_str in cigartuple:
        valid_cigar = valid_cigar + len_str + op_str  # type: ignore
    return valid_cigar


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
        __chrm_sa, __pos_sa, __strand_sa, __cigar_sa_partial, __nm_sa = blat.psl2sam(
            top_hsp, in_seq_len
        )
        if read_strand == __strand_sa:
            # same strand: different reads mode
            # MS(1) ~ SM(2) or SM(2) ~ MS(1)
            if read_mode == 1:
                __cigar_sa = "{}S{}".format(
                    read_length - in_seq_len, __cigar_sa_partial
                )  # SM
            else:
                __cigar_sa = "{1}{0}S".format(
                    read_length - in_seq_len, __cigar_sa_partial
                )  # MS
        else:
            # opposite strand: same reads mode
            # MS(1) ~ MS(1) or SM(2) ~ SM(2)
            if read_mode == 1:
                __cigar_sa = "{1}{0}S".format(
                    read_length - in_seq_len, __cigar_sa_partial
                )  # MS
            else:
                __cigar_sa = "{}S{}".format(
                    read_length - in_seq_len, __cigar_sa_partial
                )  # SM
        valid_cigar_sa = cigar_validity(__cigar_sa)
        if __mapq < mapq_cutoff and int(__nm_sa) < max_allowed_nm:
            chimeric_aln_str = "{},{},{},{},{},{};".format(
                __chrm_sa, __pos_sa, __strand_sa, valid_cigar_sa, __mapq, __nm_sa
            )
        else:
            chimeric_aln_str = ""

    return chimeric_aln_str


def strand_mode_checker(strand1: str, strand2: str, mode1: int, mode2: int) -> bool:
    """Check if the two strands are compatible with the two modes."""
    flag = False
    if (strand1 == strand2 and mode1 != mode2) or (
        strand1 != strand2 and mode1 == mode2
    ):
        flag = True
    return flag


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
    flag = False
    if mode == 2:
        if read.lt_soft_len < event_size + bp_region_seq_len:
            flag = True
    else:
        if read.rt_soft_len < event_size + bp_region_seq_len:
            flag = True
    return flag


def obtain_bp_region_seq(read, mode, bp_region_seq_len) -> str:
    """Obtain breakpoint region sequence from read.

    :param bp_region_seq_len: the length of the breakpoint region sequence
    :param read:  the chimeirc read
    :param mode: mode for the chimeirc read
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
            bp_region_seq = read_seq[
                read.lt_soft_len : read.lt_soft_len - bp_region_seq_len
            ]
        elif mode == 1:  # MS
            bp_region_seq = read_seq[
                bp_region_seq_len - read.rt_soft_len : -read.rt_soft_len
            ]
        bp_region_seq = "-" + bp_region_seq
    else:
        bp_region_seq = ""
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

        if bp_region_seq_len > 0:
            query_offset = read_lt.reference_match_size + read_rt.reference_match_size
        else:
            query_offset = (
                read_lt.reference_match_size
                + read_rt.reference_match_size
                + bp_region_seq_len
            )

        lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
        rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)

        evt_size = query_offset - target_offset

        logger.trace(f"{evt_size=}, {query_offset=}")
        if evt_size <= 0:  # deletion
            del_start = read_rt.ref_end
            del_end = del_start + abs(evt_size)
            _, _anno, _can = splicing_confirmation(
                lt_chrm,
                del_start,
                lt_chrm,
                del_end,
                splice_bin,
                genome_fasta,
                cvg,
                False,
                motif_required,
            )
            _genes = gene_annotation(lt_chrm, del_start, lt_chrm, del_end, gene_iv)
            # 1 => 2
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
                chrm_end,
                junc_end,
                splice_bin,
                genome_fasta,
                cvg,
                False,
                motif_required,
            )
            _genes = gene_annotation(
                chrm_start, junc_start, chrm_end, junc_end, gene_iv
            )
            if _nls:
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
                    chrm_end,
                    junc_end,
                    splice_bin,
                    genome_fasta,
                    cvg,
                    False,
                    motif_required,
                )
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                # 2 => 1
                if _nls:
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
):
    """Handler for same chrom and different strand."""
    logger.trace("same_chrom_diff_strand_handler takes over the task.")
    # lt_mode must be equal to rt_mode
    if lt_mode != rt_mode:
        raise SystemExit from ModesNotEqualError

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

    if ra_bp == sa_bp:  # inverted duplication (IDUP)
        chrm_start = lt_chrm
        junc_start = ra_bp
        chrm_end = lt_chrm
        junc_end = ra_bp
        _nls, _anno, _can = splicing_confirmation(
            chrm_start,
            junc_start,
            chrm_end,
            junc_end,
            splice_bin,
            genome_fasta,
            cvg,
            True,
            motif_required,
        )
        strands = (read_lt.strand, read_rt.strand)
        lt_start_end_exons = (read_lt.ref_start, read_lt.ref_end, lt_exons)
        rt_start_end_exons = (read_rt.ref_start, read_rt.ref_end, rt_exons)
        lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
        rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
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
                tuple([*strands]),
                [*_genes],
            )
        else:
            return noreturn
    else:
        chrm_start = lt_chrm
        junc_start = min(ra_bp, sa_bp)
        chrm_end = lt_chrm
        junc_end = junc_start + abs(ra_bp - sa_bp)
        _nls, _anno, _can = splicing_confirmation(
            chrm_start,
            junc_start,
            chrm_end,
            junc_end,
            splice_bin,
            genome_fasta,
            cvg,
            True,
            motif_required,
        )

        if junc_start == ra_bp:
            strands = (read_lt.strand, read_rt.strand)
            lt_start_end_exons = (read_lt.ref_start, read_lt.ref_end, lt_exons)
            rt_start_end_exons = (read_rt.ref_start, read_rt.ref_end, rt_exons)
            lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
            rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
        elif junc_start == sa_bp:
            strands = (read_rt.strand, read_lt.strand)
            lt_start_end_exons = (read_rt.ref_start, read_rt.ref_end, rt_exons)
            rt_start_end_exons = (read_lt.ref_start, read_lt.ref_end, lt_exons)
            lt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
            rt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
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
                tuple([*strands]),
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
    lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
    rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
    _nls, _anno, _can = splicing_confirmation(
        chrm_start,
        junc_start,
        chrm_end,
        junc_end,
        splice_bin,
        genome_fasta,
        cvg,
        False,
        motif_required,
    )
    _genes = gene_annotation(chrm_start, junc_start, chrm_end, junc_end, gene_iv)
    if _nls:
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
    else:
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
):
    """Diff chrom different strand handler."""
    logger.trace("diff_chrom_diff_strand_handler takes over the task.")
    # lt_mode must be equal to rt_mode
    if lt_mode != rt_mode:
        raise SystemExit from ModesNotEqualError

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
    lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
    rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
    _nls, _anno, _can = splicing_confirmation(
        chrm_start,
        junc_start,
        chrm_end,
        junc_end,
        splice_bin,
        genome_fasta,
        cvg,
        True,
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
    else:
        return noreturn
