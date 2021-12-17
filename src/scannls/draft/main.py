# !/usr/bin/env python
# -*- coding: utf-8 -*-
# ===========================================================
import copy
import inspect
import re
import subprocess
import time
from collections import defaultdict
from pathlib import Path
from typing import Any
from typing import List
from typing import Tuple

import pysam  # type: ignore
from pyfaidx import Fasta  # type: ignore
from pyfaidx import FastaNotFoundError  # type: ignore

from .._class.basicClass import Series  # type: ignore
from .._class.blat import Blat  # type: ignore
from .._class.myLogger import MyLogger  # type: ignore
from .._class.parallel import ParallelWorker  # type: ignore
from .._class.readConnector import detect_read_read_connections_from_cigar  # type: ignore
from ..utils import get_softclip_length  # type: ignore
from ..utils import reverse_complement  # type: ignore
from ..utils import write_series_to_file  # type: ignore
from .helper import blat2chimeric_alignment  # type: ignore
from .helper import extract_splice_sites  # type: ignore
from .nls_inference import infer_nls_from_connected_reads  # type: ignore


class BamScanner:
    """"""

    def __init__(
        self,
        input_bam,
        mapq_cutoff,
        output,
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

        self.in_bam_path = input_bam
        self.in_bam = pysam.AlignmentFile(input_bam, "rb")

        self.bam_chrom_info = {}

        self.output = output
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

        self.pat_left_S = re.compile(r"^(\d+)S")
        self.pat_right_S = re.compile(r"(\d+)S$")
        self.header = self._get_bam_header()

        self.representative_alignments_new_cigar = {}

    def _check_bam_sort(self, header):
        """
        check if the bam file is sorted
        """
        try:
            return True if header["HD"]["SO"] == "coordinate" else False
        except KeyError:
            self.logger.error(f"Bam file {self.in_bam} is not sorted")
            raise SystemExit

    def _count_chrom_info(self, read):
        """
        count the chrom and the chrom start and the chrom end
        """
        if read.reference_name in self.bam_chrom_info:
            if read.reference_end > self.bam_chrom_info[read.reference_name][1]:
                self.bam_chrom_info[read.reference_name][1] = read.reference_end
        else:
            self.bam_chrom_info[read.reference_name] = [
                read.reference_start,
                read.reference_end,
            ]

    def _get_bam_header(self):

        header = self.in_bam.header.as_dict()
        self._check_bam_sort(header)
        return header

    def iter_bam(self):
        # supplementary alignment cigarstring extraction
        # key: read.query_name + left S + right S
        # For minimap2, "-Y" need to be used, use soft clipping for supplementary alignments
        self.logger.info("Iter bam file and Extracting supplementary alignments")
        try:
            for read in self.in_bam.fetch():
                self._count_chrom_info(read)
                if read.is_supplementary:
                    sup_aln_cigar = read.cigarstring
                    left_mat = self.pat_left_S.search(sup_aln_cigar)
                    right_mat = self.pat_right_S.search(sup_aln_cigar)
                    if left_mat:
                        l_S_len = left_mat.group(1)
                    else:
                        l_S_len = ""
                    if right_mat:
                        r_S_len = right_mat.group(1)
                    else:
                        r_S_len = ""
                    self.representative_alignments_new_cigar[
                        "{}\t{}\t{}".format(read.qname, l_S_len, r_S_len)
                    ] = sup_aln_cigar
        except ValueError as e:
            self.logger.error(
                f"BAM index file is not found in supplementary alignments! {e}"
            )
            raise SystemExit
        else:
            return self.representative_alignments_new_cigar


def _get_genome_fasta(ref_genome, logger):
    """
    get the genome fasta file
    """
    try:
        return Fasta(str(ref_genome), sequence_always_upper=True)
    except FastaNotFoundError:
        logger.error(f"Cannot find the reference genome {ref_genome}")
        raise SystemExit


def _get_cvg_gene_iv(gtf, splice_bin, logger):
    try:
        cvg, gene_iv = extract_splice_sites(str(gtf), splice_bin)
        logger.success(f"{gtf} loaded successfully")
    except IOError as e:
        logger.error(f"read GTF file {gtf} error!", e)
        raise SystemExit
    else:
        return cvg, gene_iv


def detect_sv_from_cigar(
    *,
    read,
    mapq_cutoff,
    max_allowed_nm,
    splice_bin,
    genome_fasta,
    cvg,
    gene_iv,
    motif_required,
    blat,
    logger,
) -> Any:
    """
    :param logger: logger for logging
    :param blat: `class.Blat`
    :param read: A read from pysam.AlignedSegment
    :param mapq_cutoff: MAPQ cutoff
    :param max_allowed_nm: NM cutoff
    :param splice_bin: a small bin for splice site searching
    :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
    :param cvg: annotated splice sites (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param gene_iv: annotated gene region (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param motif_required: considering canonical splice sites only OR considering both canonical and noncanonical splice sites
    :type read: pysam.AlignedSegment
    :type mapq_cutoff: int
    :type max_allowed_nm: int
    :type splice_bin: int
    :type genome_fasta: pyfaidx.Fasta
    :type cvg: HTSeq.GenomicArrayOfSets
    :type gene_iv: HTSeq.GenomicArrayOfSets
    :type motif_required: bool
    :return: event groups in a list, every group is also a list
    :rtype: list (list of lists)
    """
    (
        read_to_read_chains,
        reads_pair_mode_dict,
    ) = detect_read_read_connections_from_cigar(
        read=read,
        mapq_cutoff=mapq_cutoff,
        max_allowed_nm=max_allowed_nm,
        blat=blat,
        logger=logger,
    )

    read_to_read_chains = [read_to_read_chains]

    event_list = []  # type: List[Tuple[Any,...]]
    if read_to_read_chains:
        # every chain is a group of connected reads
        # every chain may have a list of events
        for chain in read_to_read_chains:
            event_list = []
            for _lt, _rt in zip(chain[::1], chain[1::1]):
                # print(_lt, _rt)
                if (_lt, _rt) in reads_pair_mode_dict:
                    _lt_mode, _rt_mode = reads_pair_mode_dict[(_lt, _rt)]
                elif (_rt, _lt) in reads_pair_mode_dict:
                    _rt_mode, _lt_mode = reads_pair_mode_dict[(_rt, _lt)]
                (
                    nls_type,
                    _anno,
                    _canonical,
                    positions,
                    lt_info,
                    rt_info,
                    bp_seqs,
                    strands,
                    genes,
                ) = infer_nls_from_connected_reads(
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

                if nls_type != "NA":
                    event_list.append(
                        (
                            nls_type,
                            _anno,
                            _canonical,
                            positions,
                            lt_info,
                            rt_info,
                            bp_seqs,
                            strands,
                            genes,
                        )
                    )
                else:  # temporary solution
                    logger.warning(f"{nls_type=}")
    return event_list, read_to_read_chains[0]


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
    output,
    header,
    mapq_cutoff,
    representative_alignments_new_cigar,
    max_allowed_nm,
    min_soft_seg_len,
    blat_ident_pct_cutoff,
    splice_bin,
    motif_required,
):
    from loguru import logger

    if running_mode == "parallel":
        logger = MyLogger(identified_key, logger)

    logger.info(f"{identified_key= } start")

    output = Path(output)

    genome_fasta = _get_genome_fasta(ref_genome, logger)
    cvg, gene_iv = _get_cvg_gene_iv(gtf, splice_bin, logger)
    in_bam_io_object = pysam.AlignmentFile(in_bam_path, "rb")
    chrom_bam_io_object = in_bam_io_object.fetch(contig=identified_key)

    blat_log_file, blat_is_start_server = blat_info
    blat = Blat(two_bit, logger, port, tmp_dir, blat_log_file, blat_is_start_server)

    temp_id = int(time.time_ns())
    current_output = output.parent.joinpath(f"{identified_key}_{temp_id}_{output.name}")
    output_bam = pysam.AlignmentFile(f"{current_output}", "wb", header=header)

    nls_src_forms_list = []
    candidate_ins_dict = defaultdict(int)
    candidate_ao_dict = defaultdict(int)

    pat_left_S = re.compile(r"^(\d+)S")
    pat_right_S = re.compile(r"(\d+)S$")

    # update SA tags and iterate the BAM file
    for read in chrom_bam_io_object:
        if (
            read.mapq >= mapq_cutoff
            and not read.is_secondary
            and not read.has_tag("XA")
            and not read.is_unmapped
            and not read.is_supplementary
        ):
            chrom = read.reference_name
            chimeric_alns_num = 2
            # update SA tag of representative alignments (START)
            if read.has_tag("SA"):
                logger.trace(
                    f"Pre-checking: {read.query_name= } has SA; supplementary read: {read.is_supplementary}"
                )
                updated_chimeric_alns = []
                chimeric_alns = read.get_tag("SA")[:-1].split(";")
                chimeric_alns_num = len(chimeric_alns) + 1
                # one representative alignment could have multiple corresponding supplementary alignments
                for _aln in chimeric_alns:
                    (
                        __chr_sa,
                        __pos_sa,
                        __strand_sa,
                        __cigar_sa,
                        __mapq_sa,
                        __nm_sa,
                    ) = _aln.split(",")
                    left_mat = pat_left_S.search(__cigar_sa)
                    right_mat = pat_right_S.search(__cigar_sa)

                    l_S_len = left_mat.group(1) if left_mat else ""
                    r_S_len = right_mat.group(1) if right_mat else ""

                    tgt_key = "{}\t{}\t{}".format(read.qname, l_S_len, r_S_len)
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

            # Detect novel chimeric alignments for reads with long softclipped segment but without SA tags using BLAT
            elif not read.has_tag("SA"):
                chimeric_alns_num = 1
                read_strand = "-" if read.is_reverse else "+"
                read_length = int(read.query_length)
                _, _soft_seq, _, read_mode = get_softclip_length(read)

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
                            f"Pre-checking: {read.query_name= } does not has SA, after BLAT it has one SA tag"
                        )
                        read.set_tag("SA", chimeric_aln_str)
                        after_set_sa_chimeric_alns_num = 1

                        # _anno:annotated exon boundary (0/1/2); _can: canonical_or_not(1/0);
            # newpos=[pos,size/pos2_of_translocation, rep_aln_mode, sup_aln_mode]

            # select reads with SA tags (original or newly-added), ignore supplementary alignment
            if (
                read.has_tag("SA")
                and chimeric_alns_num == after_set_sa_chimeric_alns_num
            ):

                logger.trace(
                    f"{read.query_name= } has SA; supplementary read: {read.is_supplementary}"
                )
                event_lists, read_chains = detect_sv_from_cigar(
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
                sv_tag_list = []
                ot_tag_list = []
                nls_event_list = []
                for event in event_lists:
                    (
                        _type,
                        _anno,
                        _canonical,
                        _positions,
                        read1_info,
                        read2_info,
                        bp_seqs,
                        strands,
                        genes,
                    ) = event
                    _bp1, _bp2, _mode1, _mode2 = _positions
                    _strand1, _strand2 = strands
                    _gene1, _gene2 = genes
                    if _type in {"TDUP", "INV", "TRA", "DEL"}:
                        _chrm1, _pos1 = _bp1.split(":")
                        _chrm2, _pos2 = _bp2.split(":")
                        # SV tag uses SA tag corrdinate system (start with 1)
                        # So, position should always add 1
                        sv_tag_list.append(
                            f"{_type},{_anno}|{_canonical},{_chrm1}:{int(_pos1) + 1},{_chrm2}:{int(_pos2) + 1},{_mode1}{_mode2},{_strand1}{_strand2},{_gene1}|{_gene2};"
                        )
                        nls_event_list.append(event)

                        event_key = f"{_type}\t{_canonical}\t{_chrm1}:{int(_pos1) + 1}\t{_chrm2}:{int(_pos2) + 1}\t{_strand1}{_strand2}"
                        reversed_event_key = f"{_type}\t{_canonical}\t{_chrm2}:{int(_pos2) + 1}\t{_chrm1}:{int(_pos1) + 1}\t{_strand2}{_strand1}"
                        if event_key in candidate_ao_dict:
                            candidate_ao_dict[event_key] += 1
                        elif reversed_event_key in candidate_ao_dict:
                            candidate_ao_dict[reversed_event_key] += 1

                    elif _type in {"INS"}:
                        _end_pos = int(_bp1) + int(_bp2)
                        ot_tag_list.append(
                            f"{_type},{_anno}|{_canonical},{_bp1},{_end_pos},{_mode1}{_mode2},{_strand1}{_strand2},{_gene1}|{_gene2};"
                        )

                        candidate_ins_dict[
                            f"{_type}\t{_anno}\t{_canonical}\t{chrom}:{_bp1}\t{chrom}:{_end_pos}\t{_strand1}{_strand2}"
                        ] += 1

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
                    nls_src_forms_list.append(series)
                    logger.debug(f"{series=}")

                if sv_tag_list:
                    read.set_tag("SV", "".join(sv_tag_list))
                if ot_tag_list:
                    read.set_tag("OT", "".join(ot_tag_list))

        output_bam.write(read)

    output_bam.close()

    subprocess.check_call("samtools index {}".format(current_output), shell=True)
    logger.debug(f"{nls_src_forms_list=}")
    logger.complete()
    return current_output, nls_src_forms_list


def scanbam_run(
    two_bit,
    port,
    tmp_dir,
    blat_info,
    in_bam_path,
    mapq_cutoff,
    output,
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
    bam_scanner = BamScanner(
        input_bam=Path(in_bam_path),
        mapq_cutoff=mapq_cutoff,
        output=Path(output),
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
    representative_alignments_new_cigar = bam_scanner.iter_bam()
    # get the chromosome name we want to scan
    filter_chrom_list = [f"chr{i}" for i in range(1, 23)]
    filter_chrom_list.extend(["chrX", "chrY"])

    contigs = [
        contig
        for contig in bam_scanner.bam_chrom_info.keys()
        if contig in filter_chrom_list
    ]

    logger.info(f" Processing {contigs=}")
    # get the header of the bam file in order to write the new bam file
    header = bam_scanner.header
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

    # create a temporary directory for storing temporary files of bam
    output = Path(output)
    temp_id = int(time.time_ns())
    temp_dirname = output.parent.joinpath(f"temp_{temp_id}")
    temp_dirname.mkdir()
    temp_output = temp_dirname.joinpath(output.name)

    keyword_parameters_dict["output"] = temp_output

    intact_series_list = []
    temp_bamfiles = []

    if parallel == 1:

        for contig in contigs:
            result = _scan_bam_helper(contig, **keyword_parameters_dict)
            contig_output, contig_series_list = result
            temp_bamfiles.append(contig_output)
            intact_series_list.extend(contig_series_list)

    else:

        parallel_worker = ParallelWorker(_scan_bam_helper, logger, parallel)
        result = parallel_worker.run(*contigs, **keyword_parameters_dict)

        for contig in contigs:
            contig_output, contig_series_list = result[contig]
            temp_bamfiles.append(contig_output)
            intact_series_list.extend(contig_series_list)

    merge_cmd = f"samtools merge -f {output} {temp_dirname}/*.bam"
    subprocess.check_call(merge_cmd, shell=True)

    write_series_to_file(
        f"{output.parent.joinpath(output.stem)}_series.txt", intact_series_list
    )
    return intact_series_list
