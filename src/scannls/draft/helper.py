import re
from collections import defaultdict

__funcs__ = {
    "extract_splice_sites",
    "gene_annotation",
    "splicing_confirmation",
    "update_breakpoints",
}

from typing import Tuple, List

import HTSeq  # type: ignore


def extract_splice_sites(in_file: str, bin_size: int) -> Tuple:
    """Extract splice sites and gene regions from input GTF file
    :param in_file: gene annotation file (GTF file)
    :param bin_size: bin size to search splice site
    :type in_file: str
    :type bin_size: int
    :return: annotated splice sites (HTSeq.GenomicArrayOfSets) and annotated gene regions (HTSeq.GenomicArrayOfSets)
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
        exonList = trx_to_exon[trx_id]
        exonList.sort(key=lambda x: x.start)
        exon_num = len(exonList)
        first_exon = exonList[0]
        last_exon = exonList[-1]
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
            for _exon in exonList[1:-1]:
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
) -> Tuple:
    """obtain gene annotations for breakpoints
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
    return (gene1, gene2)


def splicing_confirmation(
    chrm1,
    pos1,
    chrm2,
    pos2,
    splice_bin,
    genome_fasta,
    cvg,
    strand_changed,
    motif_required=True,
) -> tuple:
    """Judge whether the breakpoints are NLS events or not
    if motif_required is ON: it will only report NLS events with 'canonical splice sites';
    otherwise: it will report NLS events whatever the splice sites they used
    :param chrm1: chromosome for breakpoint1
    :param chrm2: chromosome for breakpoint2
    :param pos1: position for breakpoint1
    :param pos2: position for breakpoint2
    :param splice_bin: bin size for splice sites searching
    :param genome_fasta: reference genome (pyfaidx.Fasta object)
    :param cvg: splice site annotations (HTSeq.GenomicArrayOfSets)
    :param strand_changed: whether breakpoint1 and breakpoint2 use the same strand or not
    :param motif_required: canonical splice sites required; if True: considering canonical splice sites only; else: considering canonical and noncanonical splice sites both
    :type chrm1: str
    :type chrm2: str
    :type pos1: int
    :type pos2: int
    :type splice_bin: int
    :type genome_fasta: pyfaidx.Fasta
    :type cvg: HTSeq.GenomicArrayOfSets
    :type strand_changed: bool
    :type motif_required: bool
    :return: report/not report, overlapping boundary in bits, canonical splice site/noncanonical splice site
    :rtype: tuple

    .. note::
        Possible current_output scenarios
        * True,  3(11), 1 => reported, both breakpoints overlap with annotated coding exons boundary, using canonical splice motif
        * True,  2(10), 0 => reported, one breakpoint overlap with annotated coding exons boundary, using noncanonical splice motif
        * True,  1(01), 0 => reported, one breakpoint overlap with annotated coding exons boundary, using noncanonical splice motif
        * False, 0(00), 0 => not reported, none breakpoint overlap with annotated coding exons boundary, using noncanonical splice motif
    """

    def canonical_site_finder(in_seq: str) -> List:
        """find canonical splice sites in the input sequence
        :param in_seq: input sequence (usually sequence nearby the breakpoints)
        :type in_seq: str
        :return: a list of canonical splice sites in the input sequence, it can be a empty list
        :rtype: list
        """
        candidate_sites = ["GT", "AG", "CT", "AC"]
        matches = (i in in_seq for i in candidate_sites)
        hit_sites = [j for i, j in zip(matches, candidate_sites) if i]
        return hit_sites

    def splice_paired_checker(
        hit_sites: List, pair_seq: str, splice_motif_dict
    ) -> bool:
        """find canonical splice sites in the input sequence
        :param hit_sites: canonical splice site at one end
        :param pair_seq: sequence at the other pair end
        :param splice_motif_dict: paired splice sites (same strand or different strand)
        :type hit_sites: list
        :type pair_seq: str
        :type splice_motif_dict: dict
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
        if splice_paired_checker(_junc1, junc_seq2, splice_motif_dict):
            if motif_required:
                return True, 0, 1
            else:
                return True, 0, 1
        elif splice_paired_checker(_junc2, junc_seq1, splice_motif_dict):
            if motif_required:
                return True, 0, 1
            else:
                return True, 0, 1
        else:
            if motif_required:
                return False, 0, 0
            else:
                return True, 0, 0
    # pos1 in annotated coding exon boundary, pos2 not.
    elif junc1 in splice_motif_dict and junc2 not in splice_motif_dict:
        junc_seq = genome_fasta[chrm2][pos2 - splice_bin : pos2 + splice_bin].seq
        if splice_motif_dict[junc1] in junc_seq:
            if motif_required:
                return True, 2, 1
            else:
                return True, 2, 1
        else:
            if motif_required:
                return False, 2, 0
            else:
                return True, 2, 0

    # pos2 in annotated coding exon boundary, pos1 not.
    elif junc1 not in splice_motif_dict and junc2 in splice_motif_dict:
        junc_seq = genome_fasta[chrm1][pos1 - splice_bin : pos1 + splice_bin].seq
        if splice_motif_dict[junc2] in junc_seq:
            if motif_required:
                return True, 1, 1
            else:
                return True, 1, 1
        else:
            if motif_required:
                return False, 1, 0
            else:
                return True, 1, 0

    # pos1 and pos2 both in annotated coding exon boundary
    # junc1 in splice_motif_dict and junc2 in splice_motif_dict
    else:
        if splice_motif_dict[junc1] == junc2:
            if motif_required:
                return True, 3, 1
            else:
                return True, 3, 1
        else:
            if motif_required:
                return False, 3, 0
            else:
                return True, 3, 0


def update_breakpoints(
    bp1_chrm,
    bp1_pos,
    bp2_chrm,
    bp2_pos,
    bp1_strand,
    bp2_strand,
    bp1_mode,
    bp2_mode,
    splice_bin,
    genome_fasta,
) -> tuple:
    """
    Update the breakpoints of NLS events with canonical splice sites
    keep NLS events with noncanonical splice sites unchanged (GT-AG;GC-AG;AT-AC)
    :param bp1_chrm: chromosome for breakpoint1
    :param bp1_pos: position for breakpoint1
    :param bp2_chrm: chromosome for breakpoint2
    :param bp2_pos: position for breakpoint2
    :param bp1_strand: strand for breakpoint1
    :param bp2_strand: strand for breakpoint2
    :param bp1_mode: mode for breakpoint1
    :param bp2_mode: mode for breakpoint2
    :param splice_bin: bin size for splice site searching
    :param genome_fasta: reference genome (pyfaidx.Fasta)
    :type bp1_chrm: str
    :type bp1_pos: int
    :type bp2_chrm: str
    :type bp2_pos: int
    :type bp1_strand: str (-/+)
    :type bp2_strand: str (-/+)
    :type bp1_mode: int (1/2)
    :type bp2_mode: int (1/2)
    :type splice_bin: int
    :type genome_fasta: pyfaidx.Fasta object
    :return: updated position for breakpoint1 and updated position for breakpoint2
    :rtype: tuple
    ..note :
        mode moving rules:
        * mode (SM) [2]: breakpoint move to right
        * mode (MS) [1]: breakpoint move to left
    """

    def splice_site_search(in_str, s_site, splice_bin, for_acceptor=True):
        """search for 's_site' in 'in_str';
        searching for acceptor site (searching from left to right [==>> ...XXX])
        searching for donor site (searching from right to left [XXX... <<==])
        :param in_str: nearby sequence of breakpoints
        :param s_site: splice site
        :param splice_bin: bin size for splice site searching
        :param for_acceptor: searching for acceptor site [True] OR donor site[False]
        :type in_str: str
        :type s_site: str
        :type splice_bin: int
        :type for_acceptor: bool
        :return: optimal splice site position in 'in_str'; not found(-1)
        :rtype: int
        ..note ::
             For donor, splice position will be ZXXXXGTYYY
                                                |||||^
             For acceptor, splice position will be ZXXXAGTYYY
                                                        ^||||
             assert splice_site_search('AGXXAGTXPX', 'AG', 5, True) == 5
             assert splice_site_search('AGXXAGTXPX', 'GT', 5, False) == 4
        """
        shift_positions = []
        if for_acceptor:  # ==>>
            for i in range(len(in_str) - 1):
                _motif = in_str[i : i + 2]
                if _motif == s_site:
                    shift_positions.append(i)
        else:  # <<==
            for i in range(len(in_str) - 1):
                if i == 0:
                    _motif = in_str[-i - 2 :]
                else:
                    _motif = in_str[-i - 2 : -i]
                if _motif == s_site:
                    shift_positions.append(i)
        # No found canonical splice site
        if len(shift_positions) == 0:
            return -1
        else:
            # select position closest to the center point of the 'in_str'
            ordered_shift_positions = sorted(
                shift_positions, key=lambda k: abs(k - (splice_bin - 1))
            )
            return ordered_shift_positions[0] + 1

    def obtain_bps_shift_len(bp1_dict, bp2_dict, bp1_is_upstream=True) -> tuple:
        """obtain shift length for breakpoint1 and breakpoint2 according to putative canonical splice
        site positions in 'bp1_dict' and 'bp2_dict'
        In order to find the correct breakpoint pair, canonical splice site matches is needed (GT-AG, GC-AG, AT-AC)
        :param bp1_dict: breakpoint1 splice site positions
        :param bp2_dict: breakpoint2 splice site positions
        :param bp1_is_upstream: wheather sequence at breakpoint1 is the upstream segment of the transcript
        :type bp1_dict: dict
        :type bp2_dict: dict
        :type bp1_is_upstream: bool
        :return: shift length for breakpoint1, shift length for breakpoint2
        :rtype: tuple
        """
        if bp1_is_upstream:
            if (
                "GT" in bp1_dict
                and "AG" in bp2_dict
                and bp1_dict["GT"] != -1
                and bp2_dict["AG"] != -1
            ):
                return bp1_dict["GT"], bp2_dict["AG"]
            elif (
                "GC" in bp1_dict
                and "AG" in bp2_dict
                and bp1_dict["GC"] != -1
                and bp2_dict["AG"] != -1
            ):
                return bp1_dict["GC"], bp2_dict["AG"]
            elif (
                "AT" in bp1_dict
                and "AC" in bp2_dict
                and bp1_dict["AT"] != -1
                and bp2_dict["AC"] != -1
            ):
                return bp1_dict["AT"], bp2_dict["AC"]
            else:
                # No splice site found
                return True, True
        else:
            if (
                "GT" in bp2_dict
                and "AG" in bp1_dict
                and bp2_dict["GT"] != -1
                and bp1_dict["AG"] != -1
            ):
                return bp1_dict["AG"], bp2_dict["GT"]
            elif (
                "GC" in bp2_dict
                and "AG" in bp1_dict
                and bp2_dict["GC"] != -1
                and bp1_dict["AG"] != -1
            ):
                return bp1_dict["AG"], bp2_dict["GC"]
            elif (
                "AT" in bp2_dict
                and "AC" in bp1_dict
                and bp2_dict["AT"] != -1
                and bp1_dict["AC"] != -1
            ):
                return bp1_dict["AC"], bp2_dict["AT"]
            else:
                # No splice site found
                return None, None

    bp1_pos_dict = {}
    if bp1_mode == 2:
        if bp1_strand == "+":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].seq
            # acceptor site: AG/AC
            bp1_pos_dict["AG"] = splice_site_search(boundary_seq1, "AG", splice_bin)
            bp1_pos_dict["AC"] = splice_site_search(boundary_seq1, "AC", splice_bin)
        elif bp1_strand == "-":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].reverse.complement.seq
            # donor site: GT/GC/AT
            bp1_pos_dict["GT"] = splice_site_search(
                boundary_seq1, "GT", splice_bin, False
            )
            bp1_pos_dict["GC"] = splice_site_search(
                boundary_seq1, "GC", splice_bin, False
            )
            bp1_pos_dict["AT"] = splice_site_search(
                boundary_seq1, "AT", splice_bin, False
            )
    elif bp1_mode == 1:
        if bp1_strand == "+":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].seq
            # donor site: GT/GC/AT
            bp1_pos_dict["GT"] = splice_site_search(
                boundary_seq1, "GT", splice_bin, False
            )
            bp1_pos_dict["GC"] = splice_site_search(
                boundary_seq1, "GC", splice_bin, False
            )
            bp1_pos_dict["AT"] = splice_site_search(
                boundary_seq1, "AT", splice_bin, False
            )
        elif bp1_strand == "-":
            boundary_seq1 = genome_fasta[bp1_chrm][
                bp1_pos - splice_bin : bp1_pos + splice_bin
            ].reverse.complement.seq
            # acceptor site: AG/AC
            bp1_pos_dict["AG"] = splice_site_search(boundary_seq1, "AG", splice_bin)
            bp1_pos_dict["AC"] = splice_site_search(boundary_seq1, "AC", splice_bin)

    bp2_pos_dict = {}
    if bp2_mode == 2:
        if bp2_strand == "+":
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].seq
            # acceptor site: AG/AC
            bp2_pos_dict["AG"] = splice_site_search(boundary_seq2, "AG", splice_bin)
            bp2_pos_dict["AC"] = splice_site_search(boundary_seq2, "AC", splice_bin)
        elif bp2_strand == "-":
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].reverse.complement.seq
            # donor site: GT/GC/AT
            bp2_pos_dict["GT"] = splice_site_search(
                boundary_seq2, "GT", splice_bin, False
            )
            bp2_pos_dict["GC"] = splice_site_search(
                boundary_seq2, "GC", splice_bin, False
            )
            bp2_pos_dict["AT"] = splice_site_search(
                boundary_seq2, "AT", splice_bin, False
            )
    elif bp2_mode == 1:
        if bp2_strand == "+":
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].seq
            # donor site: GT/GC/AT
            bp2_pos_dict["GT"] = splice_site_search(
                boundary_seq2, "GT", splice_bin, False
            )
            bp2_pos_dict["GC"] = splice_site_search(
                boundary_seq2, "GC", splice_bin, False
            )
            bp2_pos_dict["AT"] = splice_site_search(
                boundary_seq2, "AT", splice_bin, False
            )
        elif bp2_strand == "-":
            boundary_seq2 = genome_fasta[bp2_chrm][
                bp2_pos - splice_bin : bp2_pos + splice_bin
            ].reverse.complement.seq
            # acceptor site: AG/AC
            bp2_pos_dict["AG"] = splice_site_search(boundary_seq2, "AG", splice_bin)
            bp2_pos_dict["AC"] = splice_site_search(boundary_seq2, "AC", splice_bin)

    shift1 = 0
    shift2 = 0
    bp1_is_upstream = transcript_upstream_part_determiner(
        bp1_strand, bp2_strand, bp1_mode, bp2_mode
    )
    shift1, shift2 = obtain_bps_shift_len(bp1_pos_dict, bp2_pos_dict, bp1_is_upstream)

    # No canonical splice sites found, so change the annotation to noncanonical splice site
    if shift1 == None and shift2 == None:
        return 0, 0

    new_bp1_pos = bp1_pos
    new_bp2_pos = bp2_pos
    # [SM:2] breakpoint move to right
    # [MS:1] breakpoint move to left
    if bp1_mode == 2:
        new_bp1_pos = bp1_pos + (shift1 - splice_bin + 2) - 1
    elif bp1_mode == 1:
        new_bp1_pos = bp1_pos - (shift1 - splice_bin + 2) + 1

    if bp2_mode == 2:
        new_bp2_pos = bp2_pos + (shift2 - splice_bin + 2) - 1
    elif bp2_mode == 1:
        new_bp2_pos = bp2_pos - (shift2 - splice_bin + 2) + 1

    # print('new_bp1: ',new_bp1_pos, 'new_bp2: ',new_bp2_pos)
    return new_bp1_pos, new_bp2_pos


def transcript_upstream_part_determiner(strand1, strand2, mode1, mode2) -> bool:
    """Determine the transcript upstream part using strand and mode information
    :param strand1: strand for breakpoint1
    :param strand2: strand for breakpoint2
    :param mode1: mode for breakpoint1
    :param mode2: mode for breakpoint2
    :type strand1: str(-/+)
    :type strand2: str(-/+)
    :type mode1: int(1/2)
    :type mode2: int(1/2)
    :return: wheather breakpoint1 is the upstream segment of the transcript or not
    :rtype: bool
    """
    bp1_is_upstream = None
    if strand1 == "+" and strand2 == "-":
        if mode1 == 1 and mode2 == 1:
            is_bp1_upstream = True
        elif mode1 == 2 and mode2 == 2:
            is_bp1_upstream = False
    elif strand1 == "-" and strand2 == "+":
        if mode1 == 1 and mode2 == 1:
            is_bp1_upstream = False
        elif mode1 == 2 and mode2 == 2:
            is_bp1_upstream = True
    elif strand1 == "+" and strand2 == "+":
        if mode1 == 1 and mode2 == 2:
            is_bp1_upstream = True
        elif mode1 == 2 and mode2 == 1:
            is_bp1_upstream = False
    elif strand1 == "-" and strand2 == "-":
        if mode1 == 1 and mode2 == 2:
            is_bp1_upstream = False
        elif mode1 == 2 and mode2 == 1:
            is_bp1_upstream = True
    return is_bp1_upstream


def cigar_validity(cigar_str: str) -> str:
    """merge the first two OR last two 'same' operations in the CIGAR string generated by BLAT
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
        cigartuple[1][0] = str(int(cigartuple[0][0]) + int(cigartuple[1][0]))  # type: ignore
        del cigartuple[0]

    # last two operations are the same
    elif cigartuple[-1][1] == cigartuple[-2][1]:
        cigartuple[-2][0] = str(int(cigartuple[-1][0]) + int(cigartuple[-2][0]))  # type: ignore
        del cigartuple[-1]

    valid_cigar = ""
    for len_str, op_str in cigartuple:
        valid_cigar = valid_cigar + len_str + op_str
    return valid_cigar


def blat2chimeric_alignment(
    in_seq,
    read_length,
    read_strand,
    read_mode,
    blat,
    mapq_cutoff,
    max_allowed_nm,
    blat_ident_pct_cutoff=0.95,
) -> str:
    """
    create chimeric alignments from the alignments which has a long softclipped segment but without SA tag

    :param blat:
    :param in_seq: softclipped segment of the aligned read
    :type in_seq: str
    :param read_length: the length of the aligned read
    :type read_length: int
    :param read_strand: the strand of the aligned read (-/+)
    :type read_strand: str
    :param read_mode: mode of the aligned read (1/2)
    :type read_mode: int
    :param mapq_cutoff: MAPQ cutoff
    :type mapq_cutoff: int
    :param max_allowed_nm: mismatches cutoff used for discarding supplementary alignments
    :type max_allowed_nm: int
    :param blat_ident_pct_cutoff: BLAT HSP identity cutoff
    :type blat_ident_pct_cutoff: float
    :return: putative supplementary alignment of the alignment which is ready for put in the SA tag
    :rtype: str
    """

    chimeric_aln_str = ""
    in_seq_len = len(in_seq)

    top_hsp, __mapq = blat.fetch_mapq(in_seq, blat_ident_pct_cutoff)

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
                __cigar_sa = "{0}S{1}".format(
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
                __cigar_sa = "{0}S{1}".format(
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
