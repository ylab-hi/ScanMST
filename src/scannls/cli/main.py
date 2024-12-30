"""Module contains the main functon of the draft scannls."""

import copy
import inspect
import math
import re
from itertools import chain
from pathlib import Path

import HTSeq
import pyfaidx
import pysam
from Bio import SearchIO
from pyfaidx import Fasta, FastaNotFoundError

from scannls.base import (
    CircRNAFilter,
    Event,
    ExonFilter,
    MappingMode,
    MyLogger,
    ParallelWorker,
    RTSwitchingFilter,
    detect_read_read_connections_from_cigar,
    reverse_complement,
)
from scannls.base.helper import (
    blat2chimeric_alignment,
    extract_splice_sites,
    get_transcriptome_length,
    insertion2chimeric_alignment,
    obtain_variants_stats,
    strand_mode_checker,
)
from scannls.base.nls_inference import infer_nls_from_connected_reads
from scannls.graph import NLPath
from scannls.mtype import LoggerType
from scannls.utils import (
    cigar_validity,
    cigarstring2cigartuples,
    get_longest_insertion_sequence,
    get_softclip_length,
)


class BamScanner:
    """BcamScanner scan the bam file and output the result to a file."""

    def __init__(
        self,
        input_bam,
        mapq_cutoff,
        ref_genome,
        gtf,
        splice_in,
        logger,
        motif_required,
        max_allowed_nm,
        min_soft_seg_len,
        long_indel_length,
        substitutions_num,
        substitutions_fraction,
        indels_fraction,
        aligner,
    ) -> None:
        """Initialize the class."""
        self.in_bam_path = input_bam
        self.in_bam = pysam.AlignmentFile(input_bam, "rb")

        self.bam_chrom_info = set()
        self.mapq_cutoff = mapq_cutoff
        self.ref_genome = ref_genome.expanduser() if "~" in str(ref_genome) else ref_genome
        self.gtf = gtf.expanduser() if "~" in str(gtf) else gtf

        self.splice_bin = splice_in
        self.logger = logger
        self.motif_required = motif_required
        self.max_allowed_nm = max_allowed_nm
        self.min_soft_seg_len = min_soft_seg_len

        self.pat_left_s = re.compile(r"^(\d+)S")
        self.pat_right_s = re.compile(r"(\d+)S$")
        self.header = self._get_bam_header()
        self.total_length = 0

        self.long_indel_length = long_indel_length
        self.substitutions_num = substitutions_num
        self.substitutions_fraction = substitutions_fraction
        self.indels_fraction = indels_fraction

        self.aligner = aligner
        self.representative_alignments_new_cigar = {}
        self.representative_alignments_new_record = {}

    def _check_bam_sort(self, header) -> bool:
        """Check if the bam file is sorted."""
        try:
            return header["HD"]["SO"] == "coordinate"
        except KeyError:
            msg = f"Bam file {self.in_bam} is not sorted"
            raise RuntimeError(msg) from KeyError

    def _count_chrom_info(self, read):
        """Count the chrom and the chrom start and the chrom end."""
        if read.reference_name not in self.bam_chrom_info:
            self.bam_chrom_info.add(read.reference_name)

    def _get_bam_header(self):
        """Get bam header."""
        header = self.in_bam.header.as_dict()  # type: ignore
        self._check_bam_sort(header)
        return header

    def iter_bam(self) -> None:
        """Iterate the bam file."""
        # supplementary alignment cigarstring extraction
        # key: read.query_name + left S + right S
        # For minimap2, "-Y" need to be used, use soft clipping for supplementary alignments
        # "--cs" need to be used, cs tag store information about SNVs and DELs
        self.logger.info("Iter bam file and Extracting supplementary alignments")

        for read in self.in_bam.fetch():
            self._count_chrom_info(read)
            self.total_length += read.query_length

            if read.is_supplementary:
                if read.cigarstring is None:
                    raise ValueError
                read_strand = "-" if read.is_reverse else "+"
                read_query_seq = read.query_sequence
                read_mapq = read.mapping_quality

                left_mat = self.pat_left_s.search(read.cigarstring)
                right_mat = self.pat_right_s.search(read.cigarstring)

                lt_soft_len = int(left_mat.group(1)) if left_mat else 0
                rt_soft_len = int(right_mat.group(1)) if right_mat else 0

                nm = read.get_tag("NM")
                cs_tag = read.get_tag("cs")
                num_of_subs, subs_fraction, ins_fraction, del_fraction = obtain_variants_stats(
                    cs_tag,
                    self.long_indel_length,
                )

                if (
                    not (num_of_subs > self.substitutions_num and subs_fraction > self.substitutions_fraction)
                    and ins_fraction <= self.indels_fraction
                    and del_fraction <= self.indels_fraction
                ):
                    self.representative_alignments_new_cigar[f"{read.query_name}\t{lt_soft_len}\t{rt_soft_len}"] = read.cigarstring
                else:
                    # realignment using BLAT
                    is_passed_qc = False
                    if self.aligner:
                        read_matched_seq = _get_read_matched_sequence(read_query_seq, lt_soft_len, rt_soft_len, read_strand)
                        out_blat = self.aligner.query(in_seq=read_matched_seq)
                        blat_result = None
                        try:
                            blat_result = SearchIO.read(out_blat, "blat-psl", pslx=True)
                        except ValueError as e:
                            self.logger.warning(f"Error reading BLAT result: {e}, read_name: {read.query_name}")
                            blat_result = False

                        if blat_result:
                            hit, top_hsp, mapq_aligner = self.aligner._query_insertion(
                                blat_result,
                                read_matched_seq,
                                threshold_identity=0.9,
                                top=1,
                            )
                            if hit == 1:
                                (
                                    chrom_aligner,
                                    position_aligner,
                                    strand_aligner,
                                    cigar_aligner,
                                    nm_aligner,
                                ) = self.aligner.psl2sam(
                                    top_hsp,
                                    in_seq_len=len(read_matched_seq),
                                )
                                (
                                    substitution_num_aligner,
                                    subs_fraction_aligner,
                                    ins_fraction_aligner,
                                    del_fraction_aligner,
                                ) = self.aligner.obtain_variants_stats(top_hsp, in_seq_len=len(read_matched_seq))
                                self.logger.trace(
                                    f"{substitution_num_aligner=}, {subs_fraction_aligner=}, {ins_fraction_aligner=}, {del_fraction_aligner=}"
                                )

                                if strand_aligner == read_strand:
                                    updated_cigar = cigar_validity(f"{lt_soft_len}S{cigar_aligner}{rt_soft_len}S")
                                else:
                                    updated_cigar = cigar_validity(f"{rt_soft_len}S{cigar_aligner}{lt_soft_len}S")

                                self.logger.trace(f"{cigar_aligner=},{updated_cigar=}")
                                if (
                                    substitution_num_aligner < num_of_subs
                                    and subs_fraction_aligner < subs_fraction
                                    and ins_fraction_aligner <= ins_fraction
                                    and del_fraction_aligner <= del_fraction
                                ):
                                    # use original mapq as realignment mapq temporarily
                                    new_record = f"{chrom_aligner},{position_aligner+1},{strand_aligner},{updated_cigar},{read_mapq},{nm_aligner}"
                                    self.representative_alignments_new_record[f"{read.query_name}\t{lt_soft_len}\t{rt_soft_len}"] = new_record
                                    is_passed_qc = True

                    if not is_passed_qc:
                        self.representative_alignments_new_cigar[f"{read.query_name}\t{lt_soft_len}\t{rt_soft_len}"] = read.cigarstring


def _get_read_matched_sequence(read_query_seq, lt_soft_len, rt_soft_len, read_strand) -> str:
    """Obtain the original read sequence matched part in cigar (M+I)."""
    read_length = len(read_query_seq)
    read_matched_part_in_bam = read_query_seq[lt_soft_len : read_length - rt_soft_len]
    if read_strand == "+":
        read_matched_seq = read_matched_part_in_bam
    elif read_strand == "-":
        read_matched_seq = reverse_complement(read_matched_part_in_bam)
    return read_matched_seq


def _get_genome_fasta(ref_genome):
    """Get the genome fasta file."""
    try:
        return Fasta(str(ref_genome), sequence_always_upper=True)
    except FastaNotFoundError:
        msg = f"Reference File {ref_genome} is Not Found!"
        raise SystemExit(
            msg,
        ) from FastaNotFoundError


def _get_cvg_gene_iv(gtf, splice_bin):
    """Get the gene coverage interval.

    :param gtf: gtf file
    :param splice_bin: splice_bin file
    """
    try:
        return extract_splice_sites(str(gtf), splice_bin)
    except OSError:
        msg = f"Reading GTF file {gtf} error!"
        raise SystemExit(msg) from OSError


def update_position_event_list(event_list: list[Event]) -> list[Event]:
    """Update co-linear exons start/end positions.

    pre_read intersected
                read
      [c1] ---- [c2]
                [c2] ---- [c3]
                          [c3] ---- [c4]
            ^         ^         ^
           evt1      evt2      evt3

     evt1 an evt2: No pre_read_info, so [c1] of evt1 keep unchanged.
                   Update [c2]'s end position of evt1, using [c2]'s junnction-side position of evt2
                   Keep [c2] of evt1 as pre_read_info, Add `evt1` in the list.

     evt2 an evt3: Update [c2] of evt2 using pre_read_info.
                   Update [c3]'s end position of evt2, using [c3]'s junnction-side position of evt3
                   Keep [c3] of evt2 as pre_read_info, Add `evt2` in the list.

     For the last event:
                   Update [c3] of evt3 using pre_read_info, Add `evt3` in the list.
    """

    def is_start_match(event: Event, check_bp1=False) -> bool:
        """Check if event breakpoint equals to event's ref_start"""
        bp1_position = int(event.bp1.split(":")[1])
        bp2_position = int(event.bp2.split(":")[1])
        if check_bp1:
            if bp1_position == event.read1_ref_start:
                return True
            if bp1_position == event.read1_ref_end:
                return False
            return None
        if bp2_position == event.read2_ref_start:
            return True
        if bp2_position == event.read2_ref_end:
            return False
        return None

    hop_number = len(event_list)
    if hop_number < 2:
        return event_list
    pre_read_info = None
    updated_event_list = []
    for idx, (pre_evt, next_evt) in enumerate(zip(event_list[:], event_list[1:]), 1):
        # if is_read_reversed is True, pre_read will be read2, intersected read will be read1
        # otherwise pre_read will be read1, intersected read will be read2
        if pre_evt.is_read_reversed:
            if pre_read_info:
                pre_evt.read2_ref_start, pre_evt.read2_ref_end, pre_evt.read2_exons = pre_read_info
            # intersected read is read1 for previous event
            is_start_match_for_pre_evt = is_start_match(pre_evt, check_bp1=True)
            # if is_read_reversed is False, intersected read will be read1
            # otherwise intersected read will be read2
            if next_evt.is_read_reversed:
                # intersected read is read2
                is_start_match_for_next_evt = is_start_match(next_evt)
                if is_start_match_for_pre_evt and not is_start_match_for_next_evt:
                    pre_evt.read1_ref_end = next_evt.read2_ref_end
                    pre_evt.read1_exons.last.end = next_evt.read2_ref_end
                elif not is_start_match_for_pre_evt and is_start_match_for_next_evt:
                    pre_evt.read1_ref_start = next_evt.read2_ref_start
                    pre_evt.read1_exons.first.start = next_evt.read2_ref_start
            else:
                # intersected read is read1 for next event
                is_start_match_for_next_evt = is_start_match(next_evt, check_bp1=True)
                if is_start_match_for_pre_evt and not is_start_match_for_next_evt:
                    pre_evt.read1_ref_end = next_evt.read1_ref_end
                    pre_evt.read1_exons.last.end = next_evt.read1_ref_end
                elif not is_start_match_for_pre_evt and is_start_match_for_next_evt:
                    pre_evt.read1_ref_start = next_evt.read1_ref_start
                    pre_evt.read1_exons.first.start = next_evt.read1_ref_start
            updated_event_list.append(pre_evt)
            pre_read_info = (
                pre_evt.read1_ref_start,
                pre_evt.read1_ref_end,
                pre_evt.read1_exons,
            )
            if idx == hop_number - 1:
                if next_evt.is_read_reversed:
                    (
                        next_evt.read2_ref_start,
                        next_evt.read2_ref_end,
                        next_evt.read2_exons,
                    ) = pre_read_info
                else:
                    (
                        next_evt.read1_ref_start,
                        next_evt.read1_ref_end,
                        next_evt.read1_exons,
                    ) = pre_read_info
                updated_event_list.append(next_evt)
        else:
            if pre_read_info:
                pre_evt.read1_ref_start, pre_evt.read1_ref_end, pre_evt.read1_exons = pre_read_info

            is_start_match_for_pre_evt = is_start_match(pre_evt)
            if next_evt.is_read_reversed:
                is_start_match_for_next_evt = is_start_match(next_evt)
                if is_start_match_for_pre_evt and not is_start_match_for_next_evt:
                    pre_evt.read2_ref_end = next_evt.read2_ref_end
                    pre_evt.read2_exons.last.end = next_evt.read2_ref_end
                elif not is_start_match_for_pre_evt and is_start_match_for_next_evt:
                    pre_evt.read2_ref_start = next_evt.read2_ref_start
                    pre_evt.read2_exons.first.start = next_evt.read2_ref_start
            else:
                is_start_match_for_next_evt = is_start_match(next_evt, check_bp1=True)
                if is_start_match_for_pre_evt and not is_start_match_for_next_evt:
                    pre_evt.read2_ref_end = next_evt.read1_ref_end
                    pre_evt.read2_exons.last.end = next_evt.read1_ref_end
                elif not is_start_match_for_pre_evt and is_start_match_for_next_evt:
                    pre_evt.read2_ref_start = next_evt.read1_ref_start
                    pre_evt.read2_exons.first.start = next_evt.read1_ref_start
            updated_event_list.append(pre_evt)
            pre_read_info = (
                pre_evt.read2_ref_start,
                pre_evt.read2_ref_end,
                pre_evt.read2_exons,
            )
            if idx == hop_number - 1:
                if next_evt.is_read_reversed:
                    (
                        next_evt.read2_ref_start,
                        next_evt.read2_ref_end,
                        next_evt.read2_exons,
                    ) = pre_read_info
                else:
                    (
                        next_evt.read1_ref_start,
                        next_evt.read1_ref_end,
                        next_evt.read1_exons,
                    ) = pre_read_info
                updated_event_list.append(next_evt)
    return updated_event_list


def detect_sv_from_cigar(
    *,
    read: pysam.AlignedSegment,
    mapq_cutoff: int,
    max_allowed_nm: int,
    splice_bin: int,
    genome_fasta: pyfaidx.Fasta,
    cvg: HTSeq.GenomicArrayOfSets,
    gene_iv: HTSeq.GenomicArrayOfSets,
    motif_required: bool,
    blat_ident_pct_cutoff: float,
    aligner,
    logger: LoggerType,
):
    """Detect SV from cigar string.

    :param logger: logger for logging
    :param read: A read from pysam.AlignedSegment
    :param mapq_cutoff: MAPQ cutoff
    :param max_allowed_nm: NM cutoff
    :param splice_bin: a small bin for splice site searching
    :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
    :param cvg: annotated splice sites (HTSeq.GenomicArrayOfSets) of reference gene
           annotation (GTF file)
    :param gene_iv: annotated gene region (HTSeq.GenomicArrayOfSets) of reference
           gene annotation (GTF file)
    :param motif_required: considering canonical splice sites only OR considering both canonical
           and noncanonical splice sites
    :return: event groups in a list, every group is also a list
    :rtype: list (list of lists)
    """

    if ret := detect_read_read_connections_from_cigar(
        read=read,
        mapq_cutoff=mapq_cutoff,
        max_allowed_nm=max_allowed_nm,
        aligner=aligner,
        blat_ident_pct_cutoff=blat_ident_pct_cutoff,
        logger=logger,
    ):
        (read_chains, reads_pair_mode_dict, num_added_reads) = ret

        event_list: list[Event] = []
        # every chain is a group of connected reads
        # every chain may have a list of events
        for lt, rt in zip(read_chains[:], read_chains[1:]):
            if (lt, rt) in reads_pair_mode_dict:
                lt_mode, rt_mode = reads_pair_mode_dict[(lt, rt)]
            elif (rt, lt) in reads_pair_mode_dict:
                rt_mode, lt_mode = reads_pair_mode_dict[(rt, lt)]
            else:
                logger.warning(f"{lt=}, {rt=} are not in {reads_pair_mode_dict}")
                raise ValueError

            if not strand_mode_checker(lt.strand, rt.strand, lt_mode, rt_mode):
                logger.warning(
                    f"{lt.strand=}, {rt.strand=}, {lt_mode=}, {rt_mode=}",
                )

            original_event_info = infer_nls_from_connected_reads(
                read_lt=lt,
                read_rt=rt,
                lt_mode=lt_mode,
                rt_mode=rt_mode,
                splice_bin=splice_bin,
                genome_fasta=genome_fasta,
                cvg=cvg,
                gene_iv=gene_iv,
                motif_required=motif_required,
            )
            if original_event_info is not None:
                logger.trace(f"{original_event_info=}")
                event = Event(original_event_info)
                event_list.append(event)
                logger.trace(str(event))
            else:  # temporary solution
                logger.warning(f"Event Type is NA {original_event_info=}")

        updated_event_list = update_position_event_list(event_list)

        return updated_event_list, read_chains, num_added_reads

    return None


def _scan_bam_helper(
    identified_key,
    lock,
    *,
    running_mode,
    blat_two_bit,
    blat_port,
    tmp_dir,
    blat_info,
    aligner,
    in_bam_path,
    ref_genome,
    gtf,
    mapq_cutoff,
    representative_alignments_new_cigar,
    representative_alignments_new_record,
    max_allowed_nm,
    min_soft_seg_len,
    blat_ident_pct_cutoff,
    splice_bin,
    motif_required,
    long_indel_length,
    substitutions_num,
    substitutions_fraction,
    indels_fraction,
    circular_rna,
    exon_filter,
    rt_switching_filter_len,
    prune_threshold,
    max_allowed_ins,
):
    """Scan BAM file and write output to file."""
    from loguru import logger

    # Set exon boundary size internally
    boundary_size = 10
    genome_fasta = _get_genome_fasta(ref_genome)
    cvg, gene_iv = _get_cvg_gene_iv(gtf, splice_bin)
    exon_filter = ExonFilter(gtf, boundary_size)
    rt_switching_filter = RTSwitchingFilter(rt_switching_filter_len)
    in_bam_io_object = pysam.AlignmentFile(in_bam_path, "rb")

    if running_mode == "parallel":
        logger = MyLogger(identified_key, logger)
        chrom_bam_io_object = in_bam_io_object.fetch(contig=identified_key)
    else:
        chrom_bam_io_object = chain.from_iterable(
            [in_bam_io_object.fetch(contig=key) for key in identified_key],
        )

    logger.trace(f"{identified_key=} start")

    if aligner is not None:
        aligner.lock = lock

    nls_src_forms_list = []
    # read name: sequence at transcriptional direction
    read_name_to_seq_dict = {}

    pat_left_s = re.compile(r"^(\d+)S")
    pat_right_s = re.compile(r"(\d+)S$")

    # Circular RNA filter
    circ_rna_filter = CircRNAFilter(gtf, boundary_size, prune_threshold)
    logger.trace(f"{representative_alignments_new_cigar=}")
    logger.trace(f"{representative_alignments_new_record=}")
    # update SA tags and iterate the BAM file
    for read in chrom_bam_io_object:
        if (
            read.mapping_quality >= mapq_cutoff
            and not read.is_secondary
            and not read.has_tag("XA")
            and not read.is_unmapped
            and not read.is_supplementary
        ):
            # update SA tag of representative alignments (START)
            if read.has_tag("SA"):
                logger.trace(
                    f"Pre-checking: {read.query_name=} has SA; supplementary read: " f"{read.is_supplementary}",
                )

                updated_chimeric_alns = []
                chimeric_alns = read.get_tag("SA")[:-1].split(";")

                # one representative alignment could have multiple corresponding
                # supplementary alignments
                for _aln in chimeric_alns:
                    (
                        chr_sa,
                        pos_sa,
                        strand_sa,
                        __cigar_sa,
                        mapq_sa,
                        nm_sa,
                    ) = _aln.split(",")

                    left_mat = pat_left_s.search(__cigar_sa)
                    right_mat = pat_right_s.search(__cigar_sa)

                    l_s_len = int(left_mat.group(1)) if left_mat else 0
                    r_s_len = int(right_mat.group(1)) if right_mat else 0

                    tgt_key = f"{read.query_name}\t{l_s_len}\t{r_s_len}"
                    # After realignment, the segment with the opposite strand
                    # will have opposite way of using soft-clipping
                    alt_tgt_key = f"{read.query_name}\t{r_s_len}\t{l_s_len}"

                    # realignment has the same strand as orignal one
                    if tgt_key in representative_alignments_new_record and alt_tgt_key not in representative_alignments_new_record:
                        updated_record_realignment = representative_alignments_new_record[tgt_key]
                        (
                            chrm_realign,
                            pos_realign,
                            strand_realign,
                            cigar_realign,
                            mapq_realign,
                            nm_realign,
                        ) = updated_record_realignment.split(",")
                        if int(nm_realign) <= max_allowed_nm:
                            updated_chimeric_alns.append(updated_record_realignment)
                    # realignment has the opposite strand as orignal one
                    elif alt_tgt_key in representative_alignments_new_record and tgt_key not in representative_alignments_new_record:
                        updated_record_realignment = representative_alignments_new_record[alt_tgt_key]
                        (
                            chrm_realign,
                            pos_realign,
                            strand_realign,
                            cigar_realign,
                            mapq_realign,
                            nm_realign,
                        ) = updated_record_realignment.split(",")
                        if strand_sa != strand_realign and int(nm_realign) <= max_allowed_nm:
                            updated_chimeric_alns.append(updated_record_realignment)
                    # realignment not found use the orignal one
                    elif tgt_key in representative_alignments_new_cigar:
                        updated_cigar = representative_alignments_new_cigar[tgt_key]
                        # discard supplementary alignments with too long edit distance
                        # supplementary alignments with lower MAPQ is allowed
                        if not (int(nm_sa) > max_allowed_nm):
                            updated_chimeric_alns.append(
                                f"{chr_sa},{pos_sa},{strand_sa},{updated_cigar},{mapq_sa},{nm_sa}",
                            )

                if len(updated_chimeric_alns) == 0 | len(updated_chimeric_alns) != len(chimeric_alns):
                    read.set_tag("SA", None)
                else:
                    read.set_tag("SA", "{};".format(";".join(updated_chimeric_alns)))

                # remove SA tags of representative alignments with too much mismatches
                # update SA tag of representative alignments (END)

            # Detect novel chimeric alignments for reads with long softclipped segment
            # but without SA tags using BLAT
            elif not read.has_tag("SA"):
                read_strand = "-" if read.is_reverse else "+"
                read_ori_nm = read.get_tag("NM")
                read_length = int(read.query_length)
                ins_ref_pos, ins_seq, ins_len = get_longest_insertion_sequence(read)

                ret = get_softclip_length(read, mode=MappingMode.Type0)
                if ret is not None and ret[1] and len(ret[1]) >= min_soft_seg_len:
                    soft_seq_ori = reverse_complement(ret[1]) if read.is_reverse else ret[1]
                    read_mode = ret[-1]
                    logger.trace("Funcion blat2chimeric_alignment works on it.")
                    chimeric_aln_str = blat2chimeric_alignment(
                        soft_seq_ori,
                        read_length,
                        read_strand,
                        read_mode,
                        aligner,
                        mapq_cutoff,
                        max_allowed_nm,
                        blat_ident_pct_cutoff,
                    )

                    if chimeric_aln_str:
                        logger.trace(f"auxiliary alignment[2] is effective here. reads_name:{read.query_name} query_sequence:{soft_seq_ori}")

                        logger.trace(
                            f"Pre-checking: {read.query_name=} "
                            f"does not has SA, after BLAT [softclipped segment] (length={len(soft_seq_ori)}bp), it "
                            f"has one SA tag ",
                        )
                        read.set_tag("SA", chimeric_aln_str)

                # _anno:annotated exon boundary (0/1/2); _can: canonical_or_not(1/0);
                # Detect novel chimeric alignments for reads with long insertion (I)
                # but without SA tags using BLAT
                elif ins_ref_pos > 0:
                    logger.trace("Funcion insertion2chimeric_alignment works on it.")
                    (
                        primary_aln_cigarstring,
                        chimeric_aln_str,
                    ) = insertion2chimeric_alignment(
                        read,
                        ins_ref_pos,
                        ins_seq,
                        read_length,
                        read_strand,
                        mapq_cutoff,
                        max_allowed_nm,
                        aligner,
                        blat_ident_pct_cutoff,
                    )

                    if primary_aln_cigarstring:
                        logger.trace(f"auxiliary alignment[3] is effective here. reads_name:{read.query_name} query_sequence:{ins_seq}")
                        logger.trace(
                            f"Pre-checking: {read.query_name=} "
                            f"does not has SA, after BLAT [long insertion] (length={len(ins_seq)}bp), it has one SA tag",
                        )

                        read.cigarstring = primary_aln_cigarstring
                        read.cigartuples = cigarstring2cigartuples(
                            primary_aln_cigarstring,
                        )
                        read.reference_start = ins_ref_pos
                        read.set_tag("NM", read_ori_nm - ins_len)
                        read.set_tag("SA", chimeric_aln_str)

            # select reads with SA tags (original or newly-added), ignore supplementary alignment
            if read.has_tag("SA"):
                logger.trace(
                    f"{read.query_name=} has SA; supplementary read: {read.is_supplementary}",
                )

                nm = read.get_tag("NM")

                read_name = read.query_name
                read_sequence = reverse_complement(read.query_sequence) if read.is_reverse else read.query_sequence

                if read.cigarstring is None:
                    msg = f"{read}'s cigarstring is None"
                    raise ValueError(msg)
                num_of_subs, subs_fraction, ins_fraction, del_fraction = obtain_variants_stats(
                    read.get_tag("cs"),
                    long_indel_length,
                )

                if int(nm) <= max_allowed_nm:
                    if ret := detect_sv_from_cigar(
                        read=read,
                        mapq_cutoff=mapq_cutoff,
                        max_allowed_nm=max_allowed_nm,
                        splice_bin=splice_bin,
                        genome_fasta=genome_fasta,
                        cvg=cvg,
                        gene_iv=gene_iv,
                        motif_required=motif_required,
                        aligner=aligner,
                        blat_ident_pct_cutoff=blat_ident_pct_cutoff,
                        logger=logger,  # type: ignore
                    ):
                        event_lists, read_chains, num_added_reads = ret
                        logger.trace(f"Unordered {event_lists=}, {read_chains=}")
                    else:
                        event_lists, read_chains, num_added_reads = [], [], 0

                    nls_event_list = []
                    for event in event_lists:
                        if event.sv_type in {
                            "TDUP",
                            "INV",
                            "TRA",
                            "DEL",
                            "IDUP",
                        }:
                            if exon_filter:
                                if not exon_filter.is_breakpoints_in_same_exon(
                                    event,
                                ) and not rt_switching_filter.is_from_rt_switching(
                                    event,
                                ):
                                    nls_event_list.append(event)
                                elif rt_switching_filter.is_from_rt_switching(event):
                                    logger.trace(
                                        f"{event} is filtered out owing to RTSwitchingFilter",
                                    )
                                else:
                                    logger.trace(
                                        f"{event} is filtered out owing to ExonFilter",
                                    )
                            elif not rt_switching_filter.is_from_rt_switching(event):
                                nls_event_list.append(event)
                            else:
                                logger.trace(
                                    f"{event} is filtered out owing to RTSwitchingFilter",
                                )

                    # num of alignment segments should be equal to the number of hops + 1
                    # after exon, RT switching and other filtering, the condition may be not satisfied.
                    if len(nls_event_list) > 0 and (len(nls_event_list) == len(read.get_tag("SA")[:-1].split(";")) + num_added_reads):
                        logger.debug(f"{nls_event_list=}")
                        nlpath = NLPath.new(
                            events=nls_event_list,
                            read_chains=read_chains,
                            splice_bin=splice_bin,
                            genome_fasta=genome_fasta,
                            cvg=cvg,
                            gene_iv=gene_iv,
                            motif_required=motif_required,
                            aligner=aligner,
                        )

                        read_name_to_seq_dict[read_name] = read_sequence

                        nlpath.squeeze()
                        nlpath.setup_breakpoints()

                        if (
                            not nlpath.is_all_type_del()
                            and not nlpath.is_forming_circle(prune_threshold)
                            and nlpath.is_maximum_novel_insertion_length_valid(max_allowed_ins)
                            and nlpath.is_minimum_node_length_larger_than_threshold(boundary_size)
                        ):
                            if circular_rna == "remove":
                                if not circ_rna_filter.is_circrna(nlpath):
                                    nls_src_forms_list.append(nlpath)
                                else:
                                    logger.trace(
                                        f"{nlpath} is filtered out owing to CircRNAFilter",
                                    )
                            elif circular_rna == "extract":
                                if circ_rna_filter.is_circrna(nlpath):
                                    nls_src_forms_list.append(nlpath)
                                    logger.trace(f"extracted circular RNA: {nlpath=}")
                            else:
                                nls_src_forms_list.append(nlpath)
                                logger.trace(f"{nlpath=}")

                        # debug purposes only, remove it later
                        elif nlpath.is_forming_circle(prune_threshold):
                            logger.warning(
                                f"{nlpath} is filtered out owing to forming circle",
                            )

                else:
                    logger.trace(
                        f"{read.query_name=} does not pass the num of mismatches(edit distance) cutoff." f"{nm=}",
                    )
    logger.debug(f"Total nlpaths: {nls_src_forms_list}")
    logger.complete()
    in_bam_io_object.close()
    return nls_src_forms_list, read_name_to_seq_dict


def scanbam_run(
    blat_two_bit,
    blat_port,
    tmp_dir,
    blat_info,
    in_bam_path,
    mapq_cutoff,
    ref_genome,
    gtf,
    splice_bin,
    aligner,
    logger,
    motif_required,
    parallel,
    max_allowed_nm,
    min_soft_seg_len,
    blat_ident_pct_cutoff,
    long_indel_length,
    substitutions_num,
    substitutions_fraction,
    indels_fraction,
    species,
    circular_rna,
    exon_filter,
    rt_switching_filter_len,
    prune_threshold,
    max_allowed_ins,
):
    """Main function to run scanbam."""
    bam_scanner = BamScanner(
        input_bam=Path(in_bam_path),
        mapq_cutoff=mapq_cutoff,
        ref_genome=Path(ref_genome),
        gtf=Path(gtf),
        splice_in=splice_bin,
        logger=logger,
        motif_required=motif_required,
        max_allowed_nm=max_allowed_nm,
        min_soft_seg_len=min_soft_seg_len,
        long_indel_length=long_indel_length,
        substitutions_num=substitutions_num,
        substitutions_fraction=substitutions_fraction,
        indels_fraction=indels_fraction,
        aligner=aligner,
    )
    # iterate over all read of the bam file
    bam_scanner.iter_bam()

    representative_alignments_new_cigar = bam_scanner.representative_alignments_new_cigar
    representative_alignments_new_record = bam_scanner.representative_alignments_new_record

    logger.info(f"{representative_alignments_new_cigar=}, {representative_alignments_new_record=}")

    avg_cov = math.ceil(bam_scanner.total_length / get_transcriptome_length(species))

    num_chimeric_reads = len(representative_alignments_new_cigar)
    logger.info(
        f"species: {species}, Reads coverage: {avg_cov:.2f}, Number of chimeric reads: {num_chimeric_reads}",
    )
    # get the chromosome name we want to scan

    contigs = [contig for contig in bam_scanner.bam_chrom_info if "_" not in contig and "M" not in contig]

    logger.info(f" Processing {contigs=}")
    running_mode = "normal" if parallel == 1 else "parallel"
    # get current local namespace
    self_local_namespace = copy.copy(locals())
    # get the keyword arguments for the _scan_bam_helper function
    keyword_parameters_dict = {
        key: self_local_namespace[key] for key, value in inspect.signature(_scan_bam_helper).parameters.items() if value.kind.name == "KEYWORD_ONLY"
    }

    intact_series_list = []
    intact_read_query_name_to_sequence_dict = {}

    if parallel == 1:
        intact_series_list, intact_read_query_name_to_sequence_dict = _scan_bam_helper(contigs, None, **keyword_parameters_dict)
    else:
        parallel_worker = ParallelWorker(_scan_bam_helper, logger, parallel)
        result = parallel_worker.run(*contigs, **keyword_parameters_dict)
        for contig in contigs:
            contig_series_list, contig_read_query_name_to_sequence_dict = result[contig]
            intact_series_list.extend(contig_series_list)
            intact_read_query_name_to_sequence_dict.update(contig_read_query_name_to_sequence_dict)

    bam_scanner.in_bam.close()
    return (
        intact_series_list,
        intact_read_query_name_to_sequence_dict,
        bam_scanner.header,
        avg_cov,
    )
