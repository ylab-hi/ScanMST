import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Any
from typing import List
from typing import Tuple

import HTSeq  # type: ignore
import pyfaidx  # type: ignore
from loguru._logger import Logger

from .._class.exception import SeqNotFoundError  # type: ignore
from .helper import gene_annotation  # type: ignore
from .helper import splicing_confirmation  # type: ignore


@dataclass
class AlignerResult:
    """ Class to store the result of aligner. """

    seq1: str
    seq2: str
    start1: int
    start2: int
    end1: int
    end2: int
    score: float
    n_gaps: float
    n_mismatches: float


class Aligner(object):
    """Aligner class for aligning two sequences by semi-global algorithm

    :param seqa: the sequence to be aligned
    :param seqb:  another sequence to be aligned

    .. note::
        This class is wrapper of `gapmis` based on C implementation so that
        we can simply call it from Python.
    """

    def __init__(self, seqa: str, seqb: str):
        self.seqa = seqa
        self.seqb = seqb
        self.cmd = "gapmis -a {seqa} -b {seqb} -o {out}".format

    def __repr__(self):
        return f"Aligner(seqa={self.seqa}, seqb={self.seqb})"

    def run(self) -> AlignerResult:
        """Run the aligner with temp file and return the result

        .. note::
           1. Using gapmis to align two sequences
           2. Using `module::subprocess` to run the shell command of gapmis
           3. All files are stored in temp directory and will be deleted after aligning
           4. Return the result of aligning including the aligned sequence, start, end,
               score, number of gaps, number of mismatches for two sequences
        """
        with tempfile.TemporaryDirectory() as tmpdirname:
            tempfile_seq1 = os.path.join(tmpdirname, "seq1.fa")
            tempfile_seq2 = os.path.join(tmpdirname, "seq2.fa")
            with open(tempfile_seq1, "w") as f1, open(tempfile_seq2, "w") as f2:
                f1.write(f">seq1\n{self.seqa}\n")
                f2.write(f">seq2\n{self.seqb}\n")
            tempfile_name = os.path.join(tmpdirname, "tempfile.txt")
            subprocess.check_call(
                self.cmd(seqa=tempfile_seq1, seqb=tempfile_seq2, out=tempfile_name),
                shell=True,
            )
            align_result = self.parse_gapmis_result(tempfile_name)
        return align_result

    def parse_gapmis_result(self, result_file: str) -> AlignerResult:
        """Parse the result of gapmis

        .. note::
            1. Return the result of aligning including the aligned sequence, start, end,
                score, number of gaps, number of mismatches for two sequences
            2. start and end coordinates are 0-based and end is not included, same as Python

        """
        seqa_coords: List[Tuple] = []
        seqb_coords: List[Tuple] = []
        with open(result_file, "r") as f:
            for line in [line.strip() for line in f if not line.startswith("#")]:
                if line.startswith("seq1"):
                    # (1, 50)
                    seqa_coords.append(
                        (int(line.split()[1]) - 1, int(line.split()[-1]))
                    )
                elif line.startswith("seq2"):
                    seqb_coords.append(
                        (int(line.split()[1]) - 1, int(line.split()[-1]))
                    )
                elif line.startswith("Alignment"):
                    score = float(line.split()[-1])
                elif line.startswith("Number"):
                    mismatches = float(line.split()[-1])
                elif line.startswith("Length"):
                    gaps = float(line.split()[-1])
        return AlignerResult(
            seq1=self.seqa[seqa_coords[0][0] : seqa_coords[-1][-1]],
            seq2=self.seqb[seqb_coords[0][0] : seqb_coords[-1][-1]],
            start1=seqa_coords[0][0],
            start2=seqb_coords[0][0],
            end1=seqa_coords[-1][-1],
            end2=seqb_coords[-1][-1],
            score=score,
            n_gaps=gaps,
            n_mismatches=mismatches,
        )


def short_tdup_or_not(
    chrm, ra_mode, sa_start, sa_end, ins_seq_in_read, fastafile, logger
) -> bool:
    """judge the ins_seq_in_read is a TDUP (TDUP size < reads length)
    OR novel sequence insertion using reference sequence inferred
    from chimeric alignment start position and indel_size from 'query_offset - target_offset'

    :param chrm: the chromosome
    :param ra_mode: representative alignment mode
    :param sa_start: supplementary alignment reference start position
    :param sa_end: supplementary alignment reference end positions
    :param ins_seq_in_read: putative insertion sequence from the read
    :param fastafile: pyfaidx.Fasta object of reference genome (FASTA file)
    :type chrm: str
    :type ra_mode: int
    :type sa_start: int
    :type sa_end: int
    :type ins_seq_in_read: str
    :type fastafile: pyfaidx.Fasta object
    :returns: True if it is a short TDUP
    :rtype: bool
    """
    indel_size = len(ins_seq_in_read)

    ref_seq = (
        fastafile[chrm][sa_start - 10 : sa_start + indel_size].seq
        if ra_mode == 1
        else fastafile[chrm][sa_end - indel_size : sa_end + 10].seq
    )
    if not ref_seq or not ins_seq_in_read:
        logger.error("Gapmis: Sequence not found for semi-global alignment")
        raise SeqNotFoundError
    logger.trace(f"Aligner is working, {ref_seq=}, {ins_seq_in_read=}")
    aligner = Aligner(ref_seq, ins_seq_in_read)
    alignment_result = aligner.run()
    search_seq = alignment_result.seq1

    search_start, search_end = alignment_result.start1, alignment_result.end1
    aln_len = search_end - search_start

    total_mismatches = len(search_seq) - aln_len + alignment_result.n_mismatches
    return True if total_mismatches <= 3 else False


def infer_nls_from_connected_reads(
    read_lt,
    read_rt,
    lt_mode: int,
    rt_mode: int,
    splice_bin: int,
    genome_fasta: pyfaidx.Fasta,
    cvg: HTSeq.GenomicArrayOfSets,
    gene_iv: HTSeq.GenomicArrayOfSets,
    motif_required: bool,
    logger: Logger,
) -> Any:
    """
    :param logger:
    :param read_lt: Read 1
    :param read_rt: Read 2
    :param lt_mode: mode of Read 1
    :param rt_mode: mode of Read 2
    :param splice_bin: a small bin for splice site searching
    :param genome_fasta: pyfaidx.Fasta object of reference genome (FASTA file)
    :param cvg: annotated splice sites (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param gene_iv: annotated gene region (HTSeq.GenomicArrayOfSets) of reference gene annotation (GTF file)
    :param motif_required: considering canonical splice sites only OR considering both canonical and noncanonical splice sites
    :return: putative event from reads-pair

    .. note:: putative event

    examples: * 'NA', 0, 0, (), (), (), (), (), [] * 'TDUP', annotation, canonical/noncanonical, ('bp_chrm1:bp_pos1',
    'bp_chrm2:bp_pos2', bp_mode1, bp_mode2), (bp_read1_ref_start, bp_read1_ref_end, bp_read1_exons),
    (bp_read2_ref_start, bp_read2_ref_end, bp_read2_exons), (lt_bp_seq, rt_bp_seq), (strand1, strand2), [gene1, gene2]

       annotation explanation:
       3(11) => both breakpoints overlap with coding exons boundary
       2(10) => one breakpoint overlap with coding exons boundary
       1(01) => one breakpoint overlap with coding exons boundary
       0(00) => none breakpoint overlap with coding exons boundary
    """

    def obtain_ins_seq_from_softclipped_part_read(read, mode, indel_size) -> str:
        """
        :param read: a chimeric read
        :param mode: mode for the chimeric read
        :param indel_size: indel size inferred from 'query_offset - target_offset'
        :type read : Read
        :type mode: int
        :type indel_size: int
        :return: putative insertion sequence from the read
        :rtype: str
        """
        read_seq = read.query_sequence

        ins_seq_in_read = (
            read_seq[: read.lt_soft_len][-indel_size:]
            if mode == 2
            else read_seq[-read.rt_soft_len :][:indel_size]
        )

        return ins_seq_in_read

    def obtain_bp_region_seq(read, mode, bp_region_seq_len) -> Any:
        """
        :param bp_region_seq_len: the length of the breakpoint region sequence
        :param read:  the chimeirc read
        :param mode: mode for the chimeirc read
        :type mode: int
        :return: (putative insertion/microhomology sequence from the read; + means insertion, - means microhomology, mode)
        :rtype: tuple
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

    NAN = "NA", 0, 0, (), (), (), (), (), []  # type: ignore

    if lt_mode == 3 or rt_mode == 3:
        return NAN
    lt_chrm, lt_strand, lt_start, lt_end, lt_cigartuples, lt_cigarstring = (
        read_lt.chrom,
        read_lt.strand,
        read_lt.ref_start,
        read_lt.ref_end,
        read_lt.cigartuples,
        read_lt.cigarstring,
    )
    rt_chrm, rt_strand, rt_start, rt_end, rt_cigartuples, rt_cigarstring = (
        read_rt.chrom,
        read_rt.strand,
        read_rt.ref_start,
        read_rt.ref_end,
        read_rt.cigartuples,
        read_rt.cigarstring,
    )

    lt_exons, lt_introns = read_lt.get_exons_and_introns()
    rt_exons, rt_introns = read_rt.get_exons_and_introns()

    logger.trace(f"{lt_chrm=} {rt_chrm=}")
    logger.trace(f"{lt_strand=} {rt_strand=}")
    logger.trace(f"{lt_mode=} {rt_mode=}")

    if lt_chrm == rt_chrm:
        if lt_strand == rt_strand:  # deletion, insertion, duplication
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
                    query_offset = (
                        read_lt.reference_match_size + read_rt.reference_match_size
                    )
                else:
                    query_offset = (
                        read_lt.reference_match_size
                        + read_rt.reference_match_size
                        + bp_region_seq_len
                    )

                lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
                rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)

                evt_size = query_offset - target_offset
                if evt_size == 0:  # micro-inversion
                    return NAN
                elif evt_size < 0:  # deletion
                    del_start = read_rt.ref_start + read_rt.reference_match_size
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
                    _genes = gene_annotation(
                        lt_chrm, del_start, lt_chrm, del_end, gene_iv
                    )
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
                        (rt_start, rt_end, rt_exons),
                        (lt_start, lt_end, lt_exons),
                        (rt_bp_seq, lt_bp_seq),
                        (rt_strand, lt_strand),
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
                            (lt_start, lt_end, lt_exons),
                            (rt_start, rt_end, rt_exons),
                            (lt_bp_seq, rt_bp_seq),
                            (lt_strand, rt_strand),
                            [*_genes],
                        )
                    else:
                        return NAN
                else:  # read length > tandem duplication size
                    ins_start = read_lt.ref_start
                    ref_allele = genome_fasta[lt_chrm][ins_start : ins_start + 1].seq
                    ins_seq_in_read = obtain_ins_seq_from_softclipped_part_read(
                        read_lt, lt_mode, evt_size
                    )

                    is_dup = None
                    if short_tdup_or_not(
                        lt_chrm,
                        lt_mode,
                        rt_start,
                        rt_end,
                        ins_seq_in_read,
                        genome_fasta,
                        logger,
                    ):
                        is_dup = True
                    else:
                        is_dup = False
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
                                (lt_start, lt_end, lt_exons),
                                (rt_start, rt_end, rt_exons),
                                (lt_bp_seq, rt_bp_seq),
                                (lt_strand, rt_strand),
                                [*_genes],
                            )
                        else:
                            return NAN
                    else:  # it's a short insertion
                        _genes = gene_annotation(
                            lt_chrm, ins_start, lt_chrm, ins_start, gene_iv
                        )
                        return (
                            "INS",
                            ref_allele,
                            ins_seq_in_read,
                            (ins_start, len(ins_seq_in_read), 2, 1),
                            (lt_start, lt_end, lt_exons),
                            (rt_start, rt_end, rt_exons),
                            (lt_bp_seq, rt_bp_seq),
                            (lt_strand, rt_strand),
                            [*_genes],
                        )
            elif lt_mode == 1 and rt_mode == 2:
                target_start = read_lt.ref_start
                target_end = read_rt.ref_start + read_rt.reference_match_size
                target_offset = target_end - target_start
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.lt_soft_len
                    - read_rt.rt_soft_len
                    - read_lt.read_match_size
                    - read_rt.read_match_size
                )

                logger.trace(f"{bp_region_seq_len=}")

                if bp_region_seq_len > 0:
                    query_offset = (
                        read_lt.reference_match_size + read_rt.reference_match_size
                    )
                else:
                    query_offset = (
                        read_lt.reference_match_size
                        + read_rt.reference_match_size
                        + bp_region_seq_len
                    )

                lt_bp_seq = obtain_bp_region_seq(read_lt, lt_mode, bp_region_seq_len)
                rt_bp_seq = obtain_bp_region_seq(read_rt, rt_mode, bp_region_seq_len)
                evt_size = query_offset - target_offset
                if evt_size == 0:  # micro-inversion
                    return NAN
                elif evt_size < 0:  # deletion
                    del_start = read_lt.ref_start + read_lt.reference_match_size
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
                    _genes = gene_annotation(
                        lt_chrm, del_start, lt_chrm, del_end, gene_iv
                    )
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
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                elif evt_size >= query_offset:
                    chrm_start = rt_chrm
                    junc_start = read_rt.ref_start
                    chrm_end = rt_chrm
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
                                f"{rt_chrm}:{junc_start}",
                                f"{rt_chrm}:{junc_end}",
                                2,
                                1,
                            ),
                            (rt_start, rt_end, rt_exons),
                            (lt_start, lt_end, lt_exons),
                            (rt_bp_seq, lt_bp_seq),
                            (rt_strand, lt_strand),
                            [*_genes],
                        )
                    else:
                        return NAN
                # indel_size < query_offset
                else:
                    ins_start = read_lt.ref_start + read_lt.reference_match_size
                    ref_allele = genome_fasta[lt_chrm][ins_start : ins_start + 1].seq
                    ins_seq_in_read = obtain_ins_seq_from_softclipped_part_read(
                        read_lt, lt_mode, evt_size
                    )

                    is_dup = None

                    if short_tdup_or_not(
                        lt_chrm,
                        lt_mode,
                        rt_start,
                        rt_end,
                        ins_seq_in_read,
                        genome_fasta,
                        logger,
                    ):
                        is_dup = True
                    else:
                        is_dup = False
                    if is_dup:
                        chrm_start = rt_chrm
                        junc_start = rt_start
                        chrm_end = rt_chrm
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
                                    f"{rt_chrm}:{junc_start}",
                                    f"{rt_chrm}:{junc_end}",
                                    2,
                                    1,
                                ),
                                (rt_start, rt_end, rt_exons),
                                (lt_start, lt_end, lt_exons),
                                (rt_bp_seq, lt_bp_seq),
                                (rt_strand, lt_strand),
                                [*_genes],
                            )
                        else:
                            return NAN
                    # it is a short insertion
                    else:
                        _genes = gene_annotation(
                            rt_chrm, ins_start, rt_chrm, ins_start, gene_iv
                        )
                        return (
                            "INS",
                            ref_allele,
                            ins_seq_in_read,
                            (ins_start, len(ins_seq_in_read), 1, 2),
                            (rt_start, rt_end, rt_exons),
                            (lt_start, lt_end, lt_exons),
                            (rt_bp_seq, lt_bp_seq),
                            (rt_strand, lt_strand),
                            [*_genes],
                        )
            else:
                return NAN
        else:  # lt_strand != rt_strand
            if lt_mode == rt_mode == 1:
                ra_bp = read_lt.ref_start + read_lt.reference_match_size
                sa_bp = read_rt.ref_start + read_rt.reference_match_size
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.lt_soft_len
                    - read_rt.lt_soft_len
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
                    strands = (lt_strand, rt_strand)
                    lt_start_end_exons = (lt_start, lt_end, lt_exons)
                    rt_start_end_exons = (rt_start, rt_end, rt_exons)
                    lt_bp_seq = obtain_bp_region_seq(
                        read_lt, lt_mode, bp_region_seq_len
                    )
                    rt_bp_seq = obtain_bp_region_seq(
                        read_rt, rt_mode, bp_region_seq_len
                    )
                    if _nls:
                        return (
                            "IDUP",
                            _anno,
                            _can,
                            (
                                f"{lt_chrm}:{junc_start}",
                                f"{lt_chrm}:{junc_end}",
                                1,
                                1,
                            ),
                            lt_start_end_exons,
                            rt_start_end_exons,
                            (lt_bp_seq, rt_bp_seq),
                            tuple([*strands]),
                            [*_genes],
                        )
                    else:
                        return NAN
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
                        strands = (lt_strand, rt_strand)
                        lt_start_end_exons = (lt_start, lt_end, lt_exons)
                        rt_start_end_exons = (rt_start, rt_end, rt_exons)
                        lt_bp_seq = obtain_bp_region_seq(
                            read_lt, lt_mode, bp_region_seq_len
                        )
                        rt_bp_seq = obtain_bp_region_seq(
                            read_rt, rt_mode, bp_region_seq_len
                        )
                    elif junc_start == sa_bp:
                        strands = (rt_strand, lt_strand)
                        lt_start_end_exons = (rt_start, rt_end, rt_exons)
                        rt_start_end_exons = (lt_start, lt_end, lt_exons)
                        lt_bp_seq = obtain_bp_region_seq(
                            read_rt, rt_mode, bp_region_seq_len
                        )
                        rt_bp_seq = obtain_bp_region_seq(
                            read_lt, lt_mode, bp_region_seq_len
                        )
                    _genes = gene_annotation(
                        chrm_start, junc_start, chrm_end, junc_end, gene_iv
                    )
                    if _nls:
                        # check whether the chimeric alignments uses canonical splice sites or not (60% fraction by default)
                        if read_lt.splice_site_checker(
                            genome_fasta
                        ) and read_rt.splice_site_checker(genome_fasta):
                            return (
                                "INV",
                                _anno,
                                _can,
                                (
                                    f"{lt_chrm}:{junc_start}",
                                    f"{lt_chrm}:{junc_end}",
                                    1,
                                    1,
                                ),
                                lt_start_end_exons,
                                rt_start_end_exons,
                                (lt_bp_seq, rt_bp_seq),
                                tuple([*strands]),
                                [*_genes],
                            )
                        else:
                            return NAN
                    else:
                        return NAN
            elif lt_mode == rt_mode == 2:  # inversion
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
                    strands = (lt_strand, rt_strand)
                    lt_start_end_exons = (lt_start, lt_end, lt_exons)
                    rt_start_end_exons = (rt_start, rt_end, rt_exons)
                    lt_bp_seq = obtain_bp_region_seq(
                        read_lt, lt_mode, bp_region_seq_len
                    )
                    rt_bp_seq = obtain_bp_region_seq(
                        read_rt, rt_mode, bp_region_seq_len
                    )
                    if _nls:
                        return (
                            "IDUP",
                            _anno,
                            _can,
                            (
                                f"{lt_chrm}:{junc_start}",
                                f"{lt_chrm}:{junc_end}",
                                2,
                                2,
                            ),
                            lt_start_end_exons,
                            rt_start_end_exons,
                            (lt_bp_seq, rt_bp_seq),
                            tuple([*strands]),
                            [*_genes],
                        )
                    else:
                        return NAN
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
                        strands = (lt_strand, rt_strand)
                        lt_start_end_exons = (lt_start, lt_end, lt_exons)
                        rt_start_end_exons = (rt_start, rt_end, rt_exons)
                        lt_bp_seq = obtain_bp_region_seq(
                            read_lt, lt_mode, bp_region_seq_len
                        )
                        rt_bp_seq = obtain_bp_region_seq(
                            read_rt, rt_mode, bp_region_seq_len
                        )
                    elif junc_start == sa_bp:
                        strands = (rt_strand, lt_strand)
                        lt_start_end_exons = (rt_start, rt_end, rt_exons)
                        rt_start_end_exons = (lt_start, lt_end, lt_exons)
                        lt_bp_seq = obtain_bp_region_seq(
                            read_rt, rt_mode, bp_region_seq_len
                        )
                        rt_bp_seq = obtain_bp_region_seq(
                            read_lt, lt_mode, bp_region_seq_len
                        )
                    _genes = gene_annotation(
                        chrm_start, junc_start, chrm_end, junc_end, gene_iv
                    )
                    if _nls:
                        # check whether the chimeric alignments uses canonical splice sites or not (60% fraction by default)
                        if read_lt.splice_site_checker(
                            genome_fasta
                        ) and read_rt.splice_site_checker(genome_fasta):
                            return (
                                "INV",
                                _anno,
                                _can,
                                (
                                    f"{lt_chrm}:{junc_start}",
                                    f"{lt_chrm}:{junc_end}",
                                    2,
                                    2,
                                ),
                                lt_start_end_exons,
                                rt_start_end_exons,
                                (lt_bp_seq, rt_bp_seq),
                                tuple([*strands]),
                                [*_genes],
                            )
                        else:
                            return NAN
                    else:
                        return NAN
            else:
                return NAN
    else:  # lt_chrm != rt_chrm
        if lt_strand == rt_strand:
            if lt_mode == 1 and rt_mode == 2:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start + read_lt.reference_match_size
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.lt_soft_len
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
                    False,
                    motif_required,
                )
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                if _nls:
                    return (
                        "TRA",
                        _anno,
                        _can,
                        (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 1, 2),
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                else:
                    return NAN
            elif lt_mode == 2 and rt_mode == 1:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start
                chrm_end = rt_chrm
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
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                if _nls:
                    return (
                        "TRA",
                        _anno,
                        _can,
                        (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 2, 1),
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                else:
                    return NAN
            else:
                return NAN
        else:  # lt_strand != rt_strand
            if lt_mode == rt_mode == 1:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start + read_lt.reference_match_size
                chrm_end = rt_chrm
                junc_end = read_rt.ref_start + read_rt.reference_match_size
                bp_region_seq_len = (
                    read_lt.query_length
                    - read_lt.lt_soft_len
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
                    True,
                    motif_required,
                )
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                if _nls:
                    return (
                        "TRA",
                        _anno,
                        _can,
                        (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 1, 1),
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                else:
                    return NAN
            elif lt_mode == rt_mode == 2:
                chrm_start = lt_chrm
                junc_start = read_lt.ref_start
                chrm_end = rt_chrm
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
                _genes = gene_annotation(
                    chrm_start, junc_start, chrm_end, junc_end, gene_iv
                )
                if _nls:
                    return (
                        "TRA",
                        _anno,
                        _can,
                        (f"{lt_chrm}:{junc_start}", f"{rt_chrm}:{junc_end}", 2, 2),
                        (lt_start, lt_end, lt_exons),
                        (rt_start, rt_end, rt_exons),
                        (lt_bp_seq, rt_bp_seq),
                        (lt_strand, rt_strand),
                        [*_genes],
                    )
                else:
                    return NAN
            else:
                return NAN
