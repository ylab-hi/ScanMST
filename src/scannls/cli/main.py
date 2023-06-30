"""Module contains the main function of the draft scannls."""
import copy
import inspect
import math
import re
from itertools import chain
from pathlib import Path
from typing import Any

import HTSeq
import pyfaidx
import pysam
from pyfaidx import Fasta
from pyfaidx import FastaNotFoundError

from .. import Blat
from .. import cigarstring2cigartuples
from .. import detect_read_read_connections_from_cigar
from .. import Event
from .. import get_longest_insertion_sequence
from .. import get_softclip_length
from .. import MyLogger
from .. import ParallelWorker
from .. import reverse_complement
from ..base.filters import CircRNAFilter
from ..base.filters import ExonFilter
from ..base.filters import RTSwitchingFilter
from ..base.type import LoggerType
from ..graph.basicGraph import NLPath
from .helper import blat2chimeric_alignment
from .helper import extract_splice_sites
from .helper import get_transcriptome_length
from .helper import insertion2chimeric_alignment
from .helper import obtain_variants_stats
from .helper import strand_mode_checker
from .nls_inference import infer_nls_from_connected_reads


class BamScanner:
    """BcamScanner scan the bam file and output the result to a file."""

    def __init__(
        self,
        input_bam,
        mapq_cutoff,
        ref_genome,
        gtf,
        splice_in,
        blat,
        logger,
        motif_required,
        max_allowed_nm,
        min_soft_seg_len,
        blat_ident_pct_cutoff,
        long_indel_length,
        substitutions_num,
        substitutions_fraction,
        indels_fraction,
    ) -> None:
        """Initialize the class."""
        self.in_bam_path = input_bam
        self.in_bam = pysam.AlignmentFile(input_bam, "rb")

        self.bam_chrom_info = {}

        self.mapq_cutoff = mapq_cutoff

        self.ref_genome = (
            ref_genome.expanduser() if "~" in str(ref_genome) else ref_genome
        )

        self.gtf = gtf.expanduser() if "~" in str(gtf) else gtf

        self.splice_bin = splice_in
        self.blat = blat
        self.logger = logger
        self.motif_required = motif_required
        self.max_allowed_nm = max_allowed_nm
        self.min_soft_seg_len = min_soft_seg_len
        self.blat_ident_pct_cutoff = blat_ident_pct_cutoff

        self.pat_left_s = re.compile(r"^(\d+)S")
        self.pat_right_s = re.compile(r"(\d+)S$")
        self.header = self._get_bam_header()
        self.total_length = 0

        self.long_indel_length = long_indel_length
        self.substitutions_num = substitutions_num
        self.substitutions_fraction = substitutions_fraction
        self.indels_fraction = indels_fraction
        self.representative_alignments_new_cigar = {}

    def _check_bam_sort(self, header) -> bool:
        """Check if the bam file is sorted."""
        try:
            return header["HD"]["SO"] == "coordinate"
        except KeyError:
            raise RuntimeError(f"Bam file {self.in_bam} is not sorted") from KeyError

    def _count_chrom_info(self, read):
        """Count the chrom and the chrom start and the chrom end."""
        if read.reference_name in self.bam_chrom_info:
            if read.reference_end > self.bam_chrom_info[read.reference_name][1]:
                self.bam_chrom_info[read.reference_name][1] = read.reference_end
        else:
            self.bam_chrom_info[read.reference_name] = [
                read.reference_start,
                read.reference_end,
            ]

    def _get_bam_header(self):
        """Get bam header."""
        header = self.in_bam.header.as_dict()
        self._check_bam_sort(header)
        return header

    def iter_bam(self):
        """Iterate the bam file."""
        # supplementary alignment cigarstring extraction
        # key: read.query_name + left S + right S
        # For minimap2, "-Y" need to be used, use soft clipping for supplementary alignments
        # "--MD" need to be used, MD tag store information about SNVs and DELs
        self.logger.info("Iter bam file and Extracting supplementary alignments")

        for read in self.in_bam.fetch():
            self._count_chrom_info(read)
            self.total_length += read.query_length
            if read.is_supplementary:
                sup_aln_cigar = read.cigarstring
                left_mat = self.pat_left_s.search(sup_aln_cigar)
                right_mat = self.pat_right_s.search(sup_aln_cigar)

                l_s_len = left_mat.group(1) if left_mat else ""
                r_s_len = right_mat.group(1) if right_mat else ""

                nm = read.get_tag("NM")
                md_tag = read.get_tag("MD")
                num_of_subs, ins_fraction, del_fraction = obtain_variants_stats(
                    sup_aln_cigar, md_tag, self.long_indel_length
                )

                subs_fraction = 0 if nm == 0 else num_of_subs / nm
                if (
                    not (
                        num_of_subs > self.substitutions_num
                        and subs_fraction > self.substitutions_fraction
                    )
                    and ins_fraction <= self.indels_fraction
                    and del_fraction <= self.indels_fraction
                ):
                    self.representative_alignments_new_cigar[
                        f"{read.qname}\t{l_s_len}\t{r_s_len}"
                    ] = sup_aln_cigar
                else:
                    self.logger.trace(
                        f"{read.query_name=} does not pass the substitutions/indel cutoff. "
                        f"{nm=}, {num_of_subs=}, {subs_fraction=}, {ins_fraction=}, {del_fraction=}"
                    )
        return self.representative_alignments_new_cigar


def _get_genome_fasta(ref_genome):
    """Get the genome fasta file."""
    try:
        return Fasta(str(ref_genome), sequence_always_upper=True)
    except FastaNotFoundError:
        raise SystemExit(
            f"Reference File {ref_genome} is Not Found!"
        ) from FastaNotFoundError


def _get_cvg_gene_iv(gtf, splice_bin):
    """Get the gene coverage interval.

    :param gtf: gtf file
    :param splice_bin: splice_bin file
    """
    try:
        return extract_splice_sites(str(gtf), splice_bin)
    except OSError:
        raise SystemExit(f"Reading GTF file {gtf} error!") from OSError


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
    blat: Blat,
    logger: LoggerType,
) -> Any:
    """Detect SV from cigar string.

    :param logger: logger for logging
    :param blat: `class.Blat`
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
    (
        read_chains,
        reads_pair_mode_dict,
        num_added_reads,
    ) = detect_read_read_connections_from_cigar(
        read=read,
        mapq_cutoff=mapq_cutoff,
        max_allowed_nm=max_allowed_nm,
        blat=blat,
        logger=logger,
    )

    event_list: list[Event] = []

    if read_chains:
        # every chain is a group of connected reads
        # every chain may have a list of events
        for _lt, _rt in zip(read_chains[::1], read_chains[1::1]):
            if (_lt, _rt) in reads_pair_mode_dict:
                _lt_mode, _rt_mode = reads_pair_mode_dict[(_lt, _rt)]

            elif (_rt, _lt) in reads_pair_mode_dict:
                _rt_mode, _lt_mode = reads_pair_mode_dict[(_rt, _lt)]

            if not strand_mode_checker(_lt.strand, _rt.strand, _lt_mode, _rt_mode):
                logger.warning(
                    f"{_lt.strand=}, {_rt.strand=}, {_lt_mode=}, {_rt_mode=}"
                )

            event = Event(
                infer_nls_from_connected_reads(
                    read_lt=_lt,
                    read_rt=_rt,
                    lt_mode=_lt_mode,
                    rt_mode=_rt_mode,
                    splice_bin=splice_bin,
                    genome_fasta=genome_fasta,
                    cvg=cvg,
                    gene_iv=gene_iv,
                    motif_required=motif_required,
                )
            )

            if not event.is_type_na():
                event_list.append(event)
                logger.trace(str(event))
            else:  # temporary solution
                logger.warning(f"Event Type is NA {event=}")

    return event_list, read_chains, num_added_reads


def _scan_bam_helper(
    identified_key,
    lock,
    *,
    running_mode,
    two_bit,
    port,
    tmp_dir,
    blat_info,
    in_bam_path,
    ref_genome,
    gtf,
    mapq_cutoff,
    representative_alignments_new_cigar,
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
):
    """Scan BAM file and write output to file."""
    from loguru import logger

    genome_fasta = _get_genome_fasta(ref_genome)
    cvg, gene_iv = _get_cvg_gene_iv(gtf, splice_bin)
    exon_filter = ExonFilter(gtf, 10, logger)
    rt_switching_filter = RTSwitchingFilter(rt_switching_filter_len, logger)
    in_bam_io_object = pysam.AlignmentFile(in_bam_path, "rb")

    if running_mode == "parallel":
        logger = MyLogger(identified_key, logger)
        chrom_bam_io_object = in_bam_io_object.fetch(contig=identified_key)
    else:
        chrom_bam_io_object = chain.from_iterable(
            [in_bam_io_object.fetch(contig=key) for key in identified_key]
        )

    logger.trace(f"{identified_key=} start")

    blat_log_file, blat_is_start_server = blat_info
    blat = Blat(two_bit, port, tmp_dir, blat_log_file, blat_is_start_server, lock)

    nls_src_forms_list = []

    pat_left_s = re.compile(r"^(\d+)S")
    pat_right_s = re.compile(r"(\d+)S$")

    # Circular RNA filter
    circ_rna_filter = CircRNAFilter(gtf, 10, logger)
    # update SA tags and iterate the BAM file
    for read in chrom_bam_io_object:
        if (
            read.mapq >= mapq_cutoff
            and not read.is_secondary
            and not read.has_tag("XA")
            and not read.is_unmapped
            and not read.is_supplementary
        ):
            # update SA tag of representative alignments (START)
            if read.has_tag("SA"):
                logger.trace(
                    f"Pre-checking: {read.query_name= } has SA; supplementary read: "
                    f"{read.is_supplementary}"
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

                    l_s_len = left_mat.group(1) if left_mat else ""
                    r_s_len = right_mat.group(1) if right_mat else ""

                    tgt_key = f"{read.qname}\t{l_s_len}\t{r_s_len}"

                    if tgt_key in representative_alignments_new_cigar:
                        updated_cigar = representative_alignments_new_cigar[tgt_key]
                        # discard supplementary alignments with too many mismatches
                        # supplementary alignments with lower MAPQ is allowed
                        if not (int(nm_sa) > max_allowed_nm):
                            updated_chimeric_alns.append(
                                f"{chr_sa},{pos_sa},{strand_sa},{updated_cigar},{mapq_sa},{nm_sa}"
                            )

                if (
                    len(updated_chimeric_alns)
                    == 0 | len(updated_chimeric_alns)
                    != len(chimeric_alns)
                ):
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
                _, _soft_seq, _, read_mode = get_softclip_length(read, mode=0)
                ins_ref_pos, ins_seq, ins_len = get_longest_insertion_sequence(read)

                soft_seq_ori = (
                    reverse_complement(_soft_seq) if read.is_reverse else _soft_seq
                )

                if (
                    read_mode in {1, 2}
                    and soft_seq_ori
                    and len(soft_seq_ori) >= min_soft_seg_len
                ):
                    chimeric_aln_str = blat2chimeric_alignment(
                        soft_seq_ori,
                        read_length,
                        read_strand,
                        read_mode,
                        blat,
                        mapq_cutoff,
                        max_allowed_nm,
                        blat_ident_pct_cutoff,
                    )

                    if chimeric_aln_str:
                        logger.trace(
                            f"Pre-checking: {read.query_name= } "
                            f"does not has SA, after BLAT [softclipped segment] (length={len(soft_seq_ori)}bp), it "
                            f"has one SA tag "
                        )
                        read.set_tag("SA", chimeric_aln_str)

                # _anno:annotated exon boundary (0/1/2); _can: canonical_or_not(1/0);
                # Detect novel chimeric alignments for reads with long insertion (I)
                # but without SA tags using BLAT
                elif ins_ref_pos > 0:
                    (
                        primary_aln_cigarstring,
                        chimeric_aln_str,
                    ) = insertion2chimeric_alignment(
                        read,
                        ins_ref_pos,
                        ins_seq,
                        read_length,
                        read_strand,
                        max_allowed_nm,
                        blat,
                        blat_ident_pct_cutoff,
                    )

                    if primary_aln_cigarstring:
                        logger.trace(
                            f"Pre-checking: {read.query_name= } "
                            f"does not has SA, after BLAT [long insertion] (length={len(ins_seq)}bp), it has one SA tag"
                        )

                        read.cigarstring = primary_aln_cigarstring
                        read.cigartuples = cigarstring2cigartuples(
                            primary_aln_cigarstring
                        )
                        read.reference_start = ins_ref_pos
                        read.set_tag("NM", read_ori_nm - ins_len)
                        read.set_tag("SA", chimeric_aln_str)

            # select reads with SA tags (original or newly-added), ignore supplementary alignment
            if read.has_tag("SA"):
                logger.trace(
                    f"{read.query_name= } has SA; supplementary read: {read.is_supplementary}"
                )

                nm = read.get_tag("NM")

                num_of_subs, ins_fraction, del_fraction = obtain_variants_stats(
                    read.cigarstring, read.get_tag("MD"), long_indel_length
                )

                subs_fraction = 0 if nm == 0 else num_of_subs / nm

                if (
                    not (
                        num_of_subs > substitutions_num
                        and subs_fraction > substitutions_fraction
                    )
                    and ins_fraction <= indels_fraction
                    and del_fraction <= indels_fraction
                ):
                    event_lists, read_chains, num_added_reads = detect_sv_from_cigar(
                        read=read,
                        mapq_cutoff=mapq_cutoff,
                        max_allowed_nm=max_allowed_nm,
                        splice_bin=splice_bin,
                        genome_fasta=genome_fasta,
                        cvg=cvg,
                        gene_iv=gene_iv,
                        motif_required=motif_required,
                        blat=blat,
                        logger=logger,
                    )

                    logger.trace(f"{read_chains=}")

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
                                    event
                                ) and not rt_switching_filter.is_from_rt_switching(
                                    event
                                ):
                                    nls_event_list.append(event)
                            else:
                                if not rt_switching_filter.is_from_rt_switching(event):
                                    nls_event_list.append(event)

                    if nls_event_list:
                        logger.debug(f"{nls_event_list=}")
                        nlpath = NLPath.new(
                            events=nls_event_list,
                            read_chains=read_chains,
                            splice_bin=splice_bin,
                            genome_fasta=genome_fasta,
                            cvg=cvg,
                            gene_iv=gene_iv,
                            motif_required=motif_required,
                            blat=blat,
                        )
                        if (
                            not nlpath.is_all_type_del()
                            and nlpath.is_minimum_node_length_larger_than_threshold()
                        ):
                            if circular_rna == "remove":
                                if not circ_rna_filter.is_circRNA(nlpath):
                                    nls_src_forms_list.append(nlpath)
                                    logger.trace(f"{nlpath=}")
                            elif circular_rna == "extract":
                                if circ_rna_filter.is_circRNA(nlpath):
                                    nls_src_forms_list.append(nlpath)
                                    logger.trace(f"extracted circular RNA: {nlpath=}")
                            else:
                                nls_src_forms_list.append(nlpath)
                                logger.trace(f"{nlpath=}")
                else:
                    logger.trace(
                        f"{read.query_name= } does not pass the substitutions/indel cutoff. "
                        f"{nm=}, {num_of_subs=}, {ins_fraction=}, {del_fraction=}"
                    )
    logger.debug(f"Total Series: {nls_src_forms_list}")
    logger.complete()
    in_bam_io_object.close()
    return nls_src_forms_list


def scanbam_run(
    two_bit,
    port,
    tmp_dir,
    blat_info,
    in_bam_path,
    mapq_cutoff,
    ref_genome,
    gtf,
    splice_bin,
    blat,
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
):
    """Main function to run scanbam."""
    bam_scanner = BamScanner(
        input_bam=Path(in_bam_path),
        mapq_cutoff=mapq_cutoff,
        ref_genome=Path(ref_genome),
        gtf=Path(gtf),
        splice_in=splice_bin,
        blat=blat,
        logger=logger,
        motif_required=motif_required,
        max_allowed_nm=max_allowed_nm,
        min_soft_seg_len=min_soft_seg_len,
        blat_ident_pct_cutoff=blat_ident_pct_cutoff,
        long_indel_length=long_indel_length,
        substitutions_num=substitutions_num,
        substitutions_fraction=substitutions_fraction,
        indels_fraction=indels_fraction,
    )
    # iterate over all read of the bam file
    representative_alignments_new_cigar = bam_scanner.iter_bam()

    avg_cov = math.ceil(bam_scanner.total_length / get_transcriptome_length(species))

    num_chimeric_reads = len(representative_alignments_new_cigar)
    logger.info(
        f"species: {species}, Reads coverage: {avg_cov:.2f}, Number of chimeric reads: {num_chimeric_reads}"
    )
    # get the chromosome name we want to scan

    contigs = [
        contig
        for contig in bam_scanner.bam_chrom_info
        if "_" not in contig and "M" not in contig
    ]

    logger.info(f" Processing {contigs=}")
    # get running mode
    running_mode = "normal" if parallel == 1 else "parallel"

    # get current local namespace
    self_local_namespace = copy.copy(locals())
    # get the keyword arguments for the _scan_bam_helper function
    keyword_parameters_dict = {
        key: self_local_namespace[key]
        for key, value in inspect.signature(_scan_bam_helper).parameters.items()
        if value.kind.name == "KEYWORD_ONLY"
    }

    intact_series_list = []

    if parallel == 1:
        intact_series_list = _scan_bam_helper(contigs, None, **keyword_parameters_dict)

    else:
        parallel_worker = ParallelWorker(_scan_bam_helper, logger, parallel)
        result = parallel_worker.run(*contigs, **keyword_parameters_dict)

        for contig in contigs:
            contig_series_list = result[contig]
            intact_series_list.extend(contig_series_list)

    bam_scanner.in_bam.close()
    return intact_series_list, bam_scanner.header, avg_cov
