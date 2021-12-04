__funcs__ = {"detect_read_read_connections_from_cigar"}

from typing import Any

from loguru._logger import Logger
from pysam import AlignedSegment  # type: ignore
from ..classes import ReadsConnecter, Blat, Read  # type: ignore
from ..utils import reverse_complement  # type: ignore


def detect_read_read_connections_from_cigar(
    read: AlignedSegment, mapq_cutoff: int, blat: Blat, logger: Logger
) -> Any:
    """Detecting read-read connections with chimeric alignments CIGAR string

    :param logger:
    :param blat:
    :param mapq_cutoff: MAPQ cutoff
    :type read: pysam.AlignedSegment object
    :type mapq_cutoff: int
    :return: Read-to-Read chain (a list of lists), a dictionary of Read-pair(Read1, Read2) => mode-of-Read1, mode-of-Read2
    :rtype: tuple
    .. note::
        Read-to-Read chain scenarios
        * [[Read1, Read2, Read3]]
        * [[Read1, Read2, Read3],[Read4,Read5]]

        Dictionary of Read-pair scenarios
        * (Read1, Read2) => mode-of-Read1, mode-of-Read2
        * (Read2, Read1) => mode-of-Read2, mode-of-Read1

    .. important::
        If no 'SA' tag is found in this read, read-to-read chain and the read-pair => mode dictionary will become empty.

    #return: NLS_type(TDUP/INV), exon_boundary(0/1/2/3), canonical_or_not (1/0), [position, size, rep_aln_mode, sup_aln_mode], [++]
    #        TRA, canonical_or_not (1/0), [position, sup_position, rep_aln_mode, sup_aln_mode], [+-]
    #        e.g., INV,1,43947377,181934993,1,1,++
    #              TRA,1,160289623,chr17:17189212,1,1,+-
    """

    def format_sa_tag(in_str: str) -> Any:
        """
        To keep read.reference_start and start position of SA alignment consistent, start position of SA alignment need to substract 1
        :param in_str: string of supplementary read item in the SA tag
        :type in_str: str
        :return: chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa
        :rtype: tuple
        .. note::
             pos_sa, mapq_sa and nm_sa are integral variables now.
        """
        chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = in_str.split(",")
        pos_sa = int(pos_sa) - 1  # type: ignore
        mapq_sa = int(mapq_sa)  # type: ignore
        nm_sa = int(nm_sa)  # type: ignore
        return chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa

    def obtain_sa_query_seq_from_ra(
        query_seq_ra: str, strand_ra: str, strand_sa: str
    ) -> str:
        """a helper function to define query_seq for the supplementary alignment
        :param query_seq_ra: query sequence of representative alignment
        :param strand_ra: direction of representative read (-|+)
        :type strand_ra: str
        :param strand_sa: direction of supplementary read (-|+)
        :type strand_sa: str
        :return: query sequence of supplementary alignment
        :rtype: str
        """
        if strand_ra == strand_sa:
            return query_seq_ra
        else:
            return reverse_complement(query_seq_ra)

    if read.has_tag("SV"):
        return [], {}

    if read.is_supplementary:
        return [], {}

    # if no 'SA' tag was found, read-to-read chain will be empty
    try:
        chimeric_aln = read.get_tag("SA")[:-1].split(";")
    except KeyError:
        return [], {}

    # chimeric alignments for a chimeric read
    # a chimeric read can have multiple chimeric alignments
    chimeric_aln_list = []

    chrm_ra = read.reference_name
    pos_ra = read.reference_start
    if read.is_reverse:
        strand_ra = "-"
    else:
        strand_ra = "+"
    cigar_ra = read.cigarstring
    mapq_ra = read.mapping_quality
    nm_ra = read.get_tag("NM")
    seq_ra = read.query_sequence

    if mapq_ra > mapq_cutoff:
        chimeric_aln_list.append(
            Read.init(chrm_ra, pos_ra, strand_ra, cigar_ra, mapq_ra, nm_ra, seq_ra)
        )

    for sa_string in chimeric_aln:
        chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = format_sa_tag(sa_string)
        seq_sa = obtain_sa_query_seq_from_ra(seq_ra, strand_ra, strand_sa)
        if mapq_sa > mapq_cutoff:
            chimeric_aln_list.append(
                Read.init(chrm_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa, seq_sa)
            )

    if len(chimeric_aln_list) < 1 + len(chimeric_aln):
        return [], {}
    else:

        read_connecter = ReadsConnecter(
            aln_list=chimeric_aln_list, blat=blat, logger=logger
        )
        flag = read_connecter.run()
        if flag:
            logger.debug(f"reads chain: {read_connecter.reads_chain}")
            logger.debug(f"reads pair mode: {read_connecter.read_pair_mode_dict}")
            return (
                read_connecter.reads_chain,
                read_connecter.read_pair_mode_dict,
            )
        else:
            return [], {}
