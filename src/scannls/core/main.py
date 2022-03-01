# !/usr/bin/env python
"""This module contains the main function of the draft scannls."""
import copy
import inspect
import re
from pathlib import Path
from typing import Any
from typing import List

import HTSeq  # type: ignore
import pyfaidx  # type: ignore
import pysam  # type: ignore
from pyfaidx import Fasta
from pyfaidx import FastaNotFoundError

from .. import Blat
from .. import detect_read_read_connections_from_cigar
from .. import Event
from .. import get_softclip_length
from .. import MyLogger
from .. import ParallelWorker
from .. import reverse_complement
from .. import Series
from .._class.type import LoggerType
from .helper import blat2chimeric_alignment
from .helper import extract_splice_sites
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
    ):
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

        self.representative_alignments_new_cigar = {}
        self.supplementary_alignments_md_tag = {}

    def _check_bam_sort(self, header) -> bool:
        """Check if the bam file is sorted."""
        try:
            return header["HD"]["SO"] == "coordinate"
        except KeyError:
            raise RuntimeError(f"Bam file {self.in_bam} is not sorted") from None

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
        self.logger.info("Iter bam file and Extracting supplementary alignments")
        try:
            for read in self.in_bam.fetch():
                self._count_chrom_info(read)
                self.total_length += read.query_length
                if read.is_supplementary:
                    sup_aln_cigar = read.cigarstring
                    left_mat = self.pat_left_s.search(sup_aln_cigar)
                    right_mat = self.pat_right_s.search(sup_aln_cigar)
                    if left_mat:
                        l_s_len = left_mat.group(1)
                    else:
                        l_s_len = ""
                    if right_mat:
                        r_s_len = right_mat.group(1)
                    else:
                        r_s_len = ""
                    self.representative_alignments_new_cigar[
                        f"{read.qname}\t{l_s_len}\t{r_s_len}"
                    ] = sup_aln_cigar

                    self.supplementary_alignments_md_tag[
                        f"{read.qname}\t{l_s_len}\t{r_s_len}"
                    ] = read.get_tag("MD")
        except ValueError:
            raise SystemExit("BAM index file is not found!") from None
        else:
            return (
                self.representative_alignments_new_cigar,
                self.supplementary_alignments_md_tag,
            )


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
        raise SystemExit(f"Reading GTF file {gtf} error!") from None


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

    event_list: List[Event] = []
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
                    logger=logger,
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
):
    """Scan BAM file and write output to file."""
    from loguru import logger

    if running_mode == "parallel":
        logger = MyLogger(identified_key, logger)

    logger.trace(f"{identified_key=} start")

    genome_fasta = _get_genome_fasta(ref_genome)
    cvg, gene_iv = _get_cvg_gene_iv(gtf, splice_bin)
    in_bam_io_object = pysam.AlignmentFile(in_bam_path, "rb")
    chrom_bam_io_object = in_bam_io_object.fetch(contig=identified_key)

    blat_log_file, blat_is_start_server = blat_info
    blat = Blat(two_bit, logger, port, tmp_dir, blat_log_file, blat_is_start_server)

    nls_src_forms_list = []

    pat_left_s = re.compile(r"^(\d+)S")
    pat_right_s = re.compile(r"(\d+)S$")

    # update SA tags and iterate the BAM file
    for read in chrom_bam_io_object:
        if (
            read.mapq >= mapq_cutoff
            and not read.is_secondary
            and not read.has_tag("XA")
            and not read.is_unmapped
            and not read.is_supplementary
        ):
            chimeric_alns_num = 2
            # update SA tag of representative alignments (START)
            if read.has_tag("SA"):
                logger.trace(
                    f"Pre-checking: {read.query_name= } has SA; supplementary read: "
                    f"{read.is_supplementary}"
                )
                updated_chimeric_alns = []
                chimeric_alns = read.get_tag("SA")[:-1].split(";")
                chimeric_alns_num = len(chimeric_alns) + 1
                # one representative alignment could have multiple corresponding
                # supplementary alignments
                for _aln in chimeric_alns:
                    (
                        __chr_sa,
                        __pos_sa,
                        __strand_sa,
                        __cigar_sa,
                        __mapq_sa,
                        __nm_sa,
                    ) = _aln.split(",")
                    left_mat = pat_left_s.search(__cigar_sa)
                    right_mat = pat_right_s.search(__cigar_sa)

                    l_s_len = left_mat.group(1) if left_mat else ""
                    r_s_len = right_mat.group(1) if right_mat else ""

                    tgt_key = f"{read.qname}\t{l_s_len}\t{r_s_len}"
                    if tgt_key in representative_alignments_new_cigar:
                        __updated_cigar = representative_alignments_new_cigar[tgt_key]
                        # discard supplementary alignments with too many mismatches or lower MAPQ
                        if not (
                            int(__nm_sa) > max_allowed_nm
                            or int(__mapq_sa) < mapq_cutoff
                        ):
                            updated_chimeric_alns.append(
                                "{},{},{},{},{},{}".format(
                                    __chr_sa,
                                    __pos_sa,
                                    __strand_sa,
                                    __updated_cigar,
                                    __mapq_sa,
                                    __nm_sa,
                                )
                            )
                if len(updated_chimeric_alns) == 0:
                    read.set_tag("SA", None)
                else:
                    read.set_tag("SA", "{};".format(";".join(updated_chimeric_alns)))
                    after_set_sa_chimeric_alns = read.get_tag("SA")[:-1].split(";")
                    after_set_sa_chimeric_alns_num = len(after_set_sa_chimeric_alns) + 1
                # remove SA tags of representative alignments with too much mismatches
                # update SA tag of representative alignments (END)

            # Detect novel chimeric alignments for reads with long softclipped segment
            # but without SA tags using BLAT
            elif not read.has_tag("SA"):
                chimeric_alns_num = 1
                read_strand = "-" if read.is_reverse else "+"
                read_length = int(read.query_length)
                _, _soft_seq, _, read_mode = get_softclip_length(read, mode=0)

                if read.is_reverse:
                    soft_seq_ori = reverse_complement(_soft_seq)
                else:
                    soft_seq_ori = _soft_seq

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
                            f"does not has SA, after BLAT it has one SA tag"
                        )
                        read.set_tag("SA", chimeric_aln_str)
                        after_set_sa_chimeric_alns_num = 1

                        # _anno:annotated exon boundary (0/1/2); _can: canonical_or_not(1/0);

            # select reads with SA tags (original or newly-added), ignore supplementary alignment
            if (
                read.has_tag("SA")
                and chimeric_alns_num == after_set_sa_chimeric_alns_num
            ):

                logger.trace(
                    f"{read.query_name= } has SA; supplementary read: {read.is_supplementary}"
                )
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
                    if event.sv_type in {"TDUP", "INV", "TRA", "DEL", "IDUP"}:
                        nls_event_list.append(event)

                chimeric_alns_num += num_added_reads
                if nls_event_list and (len(nls_event_list) + 1 == chimeric_alns_num):
                    series = Series(blat=blat, logger=logger)
                    logger.debug(f"{nls_event_list=}")
                    series.init(
                        nls_event_list,
                        read_chains,
                        splice_bin,
                        genome_fasta,
                        cvg,
                        gene_iv,
                        motif_required,
                    )
                    series.disable_blat_logger()
                    if not series.is_all_type_del():
                        nls_src_forms_list.append(series)
                        logger.trace(f"{series=}")

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
    )
    # iterate over all read of the bam file
    (
        representative_alignments_new_cigar,
        supplementary_alignments_md_tag,
    ) = bam_scanner.iter_bam()
    avg_cov = bam_scanner.total_length / 150000000
    num_chimeric_reads = len(representative_alignments_new_cigar)
    logger.info(
        f"Reads coverage: {avg_cov:.4f}, Number of chimeric reads: {num_chimeric_reads}"
    )
    # get the chromosome name we want to scan
    filter_chrom_list = [f"chr{i}" for i in range(1, 23)]
    filter_chrom_list.extend(["chrX", "chrY"])

    contigs = [
        contig
        for contig in bam_scanner.bam_chrom_info.keys()
        if contig in filter_chrom_list
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

        intact_series_list = _scan_bam_helper(None, **keyword_parameters_dict)

    else:

        parallel_worker = ParallelWorker(_scan_bam_helper, logger, parallel)
        result = parallel_worker.run(*contigs, **keyword_parameters_dict)

        for contig in contigs:
            contig_series_list = result[contig]
            intact_series_list.extend(contig_series_list)

    return intact_series_list, bam_scanner.in_bam
