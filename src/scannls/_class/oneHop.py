#!/usr/bin/env python
# ===============================================================================
import datetime
import secrets
from dataclasses import dataclass
from functools import singledispatchmethod
from pathlib import Path
from typing import Any
from typing import IO
from typing import List
from typing import Optional

from pyfaidx import Fasta  # type: ignore
from pyfaidx import FastaNotFoundError

from .. import __version__
from ..type import LoggerType
from .basicClass import reverse_complement
from .exception import NumberOfHopIsNotValidError
from .writer import Writer


@dataclass
class MetaExon:
    """Class to store the information of one meta-exon."""

    chrom: str
    locus: str
    strand: str
    exons: list
    p5_pos: int
    p3_pos: int
    wt_seq: str
    mt_seq: str
    nls_type: Optional[str] = None


class OneHop:
    """Generate one-hop NLS transcript (GTF and FASTA)."""

    def __init__(
        self,
        chrom_to_genes: dict,
        gene_to_trx: dict,
        trx_to_exons: dict,
        trx_to_introns: dict,
        gene_to_intergenic: dict,
        reference: str,
        logger: LoggerType,
    ) -> None:
        """Initialize OneHop.

        :param chrom_to_genes: chrom => list((gene_name, strand))
        :param gene_to_trx: gene => set(transcript_id))
        :param trx_to_exons: transcript_id => [exon:HTSeq.GenomicInterval]
        :param trx_to_introns: transcript_id => [intron:HTSeq.GenomicInterval]
        :param gene_to_intergenic: gene => {"upstream": (start, end), "downstream": (start, end)}
        :param reference: reference fasta file
        :param logger: logger
        """
        self.chrom_to_genes = chrom_to_genes
        self.available_chroms = list(chrom_to_genes)
        self.gene_to_trx = gene_to_trx
        self.gene_to_intergenic = gene_to_intergenic
        self.trx_to_exons = trx_to_exons
        self.trx_to_introns = trx_to_introns
        self.reference = Path(reference)
        if not self.reference.exists():
            raise FastaNotFoundError
        self.reference_io = Fasta(reference, sequence_always_upper=True)
        self.logger = logger

    def _metaexon(
        self,
        locus: str,
        strand: str,
        locus_type: str,
        direction: str,
        nls_type: str = None,
    ) -> MetaExon:
        """Generate meta-exon according to locus type (exonic, intronic, intergenic).

        ..note: direction options: upstream OR downstream
              * For exonic, upstream: upstream exons of selected annotated transcript
                downstream: downstream exons of selected annotated transcript
              * For intergenic, upstream: upstream intergenic region of selected gene
                downstream: downstream intergenic region of selected gene
        """
        _trxs = self.gene_to_trx[locus]
        _tgt_trx = secrets.choice(_trxs)
        _chrom = self.trx_to_exons[_tgt_trx][0].chrom
        _num_exon = len(self.trx_to_exons[_tgt_trx])
        _num_intron = len(self.trx_to_introns[_tgt_trx])
        wt_seq = ""
        if locus_type == "exonic" or locus_type == "intronic":
            if strand == "+":
                for exon in self.trx_to_exons[_tgt_trx]:
                    wt_seq += self.reference_io[exon.chrom][exon.start : exon.end].seq
            else:
                for exon in self.trx_to_exons[_tgt_trx][::-1]:
                    wt_seq += self.reference_io[exon.chrom][
                        exon.start : exon.end
                    ].reverse.complement.seq
        # intergenic/intronic wt_seq == ''

        r_metaexon = []
        mt_seq = ""
        if locus_type == "exonic":
            _upstream_choice_exon_num = round(_num_exon / 2)
            _downstream_choice_exon_num = _num_exon - _upstream_choice_exon_num
            if direction == "upstream":
                if strand == "+":
                    tgt_exon_index = secrets.choice(range(_upstream_choice_exon_num))
                    # list of GenomicInterval
                    _metaexon = self.trx_to_exons[_tgt_trx][: tgt_exon_index + 1]
                    for exon in _metaexon:
                        mt_seq += self.reference_io[exon.chrom][
                            exon.start : exon.end
                        ].seq
                    p5_pos, p3_pos = _metaexon[0].start, _metaexon[-1].end
                else:
                    tgt_exon_index = secrets.choice(
                        range(_downstream_choice_exon_num, _num_exon)
                    )
                    _metaexon = self.trx_to_exons[_tgt_trx][tgt_exon_index:]
                    for exon in _metaexon[::-1]:
                        mt_seq += self.reference_io[exon.chrom][
                            exon.start : exon.end
                        ].reverse.complement.seq
                    p5_pos, p3_pos = _metaexon[-1].end, _metaexon[0].start
            # downstream direction
            else:
                if strand == "+":
                    tgt_exon_index = secrets.choice(
                        range(_upstream_choice_exon_num, _num_exon)
                    )
                    _metaexon = self.trx_to_exons[_tgt_trx][tgt_exon_index:]
                    for exon in _metaexon:
                        mt_seq += self.reference_io[exon.chrom][
                            exon.start : exon.end
                        ].seq
                    p5_pos, p3_pos = _metaexon[0].start, _metaexon[-1].end
                else:
                    tgt_exon_index = secrets.choice(range(_downstream_choice_exon_num))
                    _metaexon = self.trx_to_exons[_tgt_trx][: tgt_exon_index + 1]
                    for exon in _metaexon[::-1]:
                        mt_seq += self.reference_io[exon.chrom][
                            exon.start : exon.end
                        ].reverse.complement.seq
                    p5_pos, p3_pos = _metaexon[-1].end, _metaexon[0].start
            r_metaexon = _metaexon
        elif locus_type == "intronic":
            tgt_intron_index = secrets.choice(range(_num_intron))
            # GenomicInterval
            _metaexon = self.trx_to_introns[_tgt_trx][tgt_intron_index]
            if strand == "+":
                mt_seq += self.reference_io[_metaexon.chrom][
                    _metaexon.start : _metaexon.end
                ].seq
                p5_pos, p3_pos = _metaexon.start, _metaexon.end
            else:
                mt_seq += self.reference_io[_metaexon.chrom][
                    _metaexon.start : _metaexon.end
                ].reverse.complement.seq
                p5_pos, p3_pos = _metaexon.end, _metaexon.start
            r_metaexon = [_metaexon]

        elif locus_type == "intergenic":
            # GenomicInterval
            _metaexon = self.gene_to_intergenic[locus][direction]
            if strand == "+":
                mt_seq = self.reference_io[_metaexon.chrom][
                    _metaexon.start : _metaexon.end
                ].seq
                p5_pos, p3_pos = _metaexon.start, _metaexon.end
            else:
                mt_seq = self.reference_io[_metaexon.chrom][
                    _metaexon.start : _metaexon.end
                ].reverse.complement.seq
                p5_pos, p3_pos = _metaexon.end, _metaexon.start
            r_metaexon = [_metaexon]

        return MetaExon(
            chrom=_chrom,
            locus=locus,
            strand=strand,
            exons=r_metaexon,
            p5_pos=p5_pos,
            p3_pos=p3_pos,
            wt_seq=wt_seq,
            mt_seq=mt_seq,
            nls_type=nls_type,
        )

    @staticmethod
    def reverse_metaexon(input_metaexon: MetaExon) -> MetaExon:
        """Reverse metaexon for IDUP."""
        strand = "+" if input_metaexon.strand == "-" else "-"
        return MetaExon(
            chrom=input_metaexon.chrom,
            locus=input_metaexon.locus,
            strand=strand,
            exons=input_metaexon.exons,
            p5_pos=input_metaexon.p3_pos,
            p3_pos=input_metaexon.p5_pos,
            wt_seq=input_metaexon.wt_seq,
            mt_seq=reverse_complement(input_metaexon.mt_seq),
            nls_type=input_metaexon.nls_type,
        )

    def _tdup_hopper(self, current_metaexons) -> None:
        """TDUP hopper."""
        if current_metaexons == []:
            _select_chrom = secrets.choice(self.available_chroms)
            while True:
                _select_gene1, _select_strand1 = secrets.choice(
                    self.chrom_to_genes[_select_chrom]
                )
                _index1 = self.chrom_to_genes[_select_chrom].index(
                    (_select_gene1, _select_strand1)
                )
                _metaexon1 = self._metaexon(
                    locus=_select_gene1,
                    strand=_select_strand1,
                    locus_type="exonic",
                    direction="upstream",
                    nls_type="TDUP",
                )

                _shift = secrets.choice(range(2, 5))
                _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
                if _select_strand1 == "+":
                    _index2 = _index1 - _shift
                else:
                    _index2 = _index1 + _shift
                _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                    _index2
                ]
                if _select_strand1 == _select_strand2:
                    _metaexon2 = self._metaexon(
                        locus=_select_gene2,
                        strand=_select_strand2,
                        locus_type=_locus_type,
                        direction="downstream",
                    )
                    current_metaexons.append(_metaexon1)
                    current_metaexons.append(_metaexon2)
                    break
        else:
            current_metaexons[-1].nls_type = "TDUP"
            last_metaexon = current_metaexons[-1]
            _select_chrom = last_metaexon.chrom
            _select_gene1, _select_strand1 = last_metaexon.locus, last_metaexon.strand
            _index1 = self.chrom_to_genes[_select_chrom].index(
                (_select_gene1, _select_strand1)
            )

            while True:
                _shift = secrets.choice(range(2, 12))
                _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
                if _select_strand1 == "+":
                    _index2 = _index1 - _shift
                else:
                    _index2 = _index1 + _shift
                _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                    _index2
                ]
                if _select_strand1 == _select_strand2:
                    _metaexon2 = self._metaexon(
                        locus=_select_gene2,
                        strand=_select_strand2,
                        locus_type=_locus_type,
                        direction="downstream",
                    )
                    current_metaexons.append(_metaexon2)
                    break
            return current_metaexons

    def _idup_hopper(self, current_metaexons) -> None:
        """IDUP hopper."""
        if current_metaexons == []:
            _select_chrom = secrets.choice(self.available_chroms)
            _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
            _select_gene1, _select_strand1 = secrets.choice(
                self.chrom_to_genes[_select_chrom]
            )
            _metaexon1 = self._metaexon(
                locus=_select_gene1,
                strand=_select_strand1,
                locus_type=_locus_type,
                direction="upstream",
                nls_type="IDUP",
            )
            _metaexon2 = OneHop.reverse_metaexon(_metaexon1)
            current_metaexons.append(_metaexon1)
            current_metaexons.append(_metaexon2)
        else:
            current_metaexons[-1].nls_type = "IDUP"
            last_metaexon = current_metaexons[-1]
            _metaexon2 = OneHop.reverse_metaexon(last_metaexon)
            current_metaexons.append(_metaexon2)

    def _inv_hopper(self, current_metaexons) -> None:
        """INV hopper."""
        if current_metaexons == []:
            _select_chrom = secrets.choice(self.available_chroms)
            while True:
                _select_gene1, _select_strand1 = secrets.choice(
                    self.chrom_to_genes[_select_chrom]
                )
                _index1 = self.chrom_to_genes[_select_chrom].index(
                    (_select_gene1, _select_strand1)
                )
                _metaexon1 = self._metaexon(
                    locus=_select_gene1,
                    strand=_select_strand1,
                    locus_type="exonic",
                    direction="upstream",
                    nls_type="INV",
                )

                _shift = secrets.choice(
                    list(map(lambda x: -x, range(2, 5))) + list(range(2, 5))
                )
                _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
                _index2 = _index1 + _shift
                _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                    _index2
                ]
                if _select_strand1 != _select_strand2:
                    _metaexon2 = self._metaexon(
                        locus=_select_gene2,
                        strand=_select_strand2,
                        locus_type=_locus_type,
                        direction="downstream",
                    )
                    current_metaexons.append(_metaexon1)
                    current_metaexons.append(_metaexon2)
                    break
        else:
            current_metaexons[-1].nls_type = "INV"
            last_metaexon = current_metaexons[-1]
            _select_chrom = last_metaexon.chrom
            _select_gene1, _select_strand1 = last_metaexon.locus, last_metaexon.strand
            _index1 = self.chrom_to_genes[_select_chrom].index(
                (_select_gene1, _select_strand1)
            )

            while True:
                _shift = secrets.choice(
                    list(map(lambda x: -x, range(2, 12))) + list(range(2, 12))
                )
                _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
                _index2 = _index1 + _shift
                _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                    _index2
                ]
                if _select_strand1 != _select_strand2:
                    _metaexon2 = self._metaexon(
                        locus=_select_gene2,
                        strand=_select_strand2,
                        locus_type=_locus_type,
                        direction="downstream",
                    )
                    current_metaexons.append(_metaexon2)
                    break
            return current_metaexons

    def _tra_hopper(self, current_metaexons) -> None:
        """TRA hopper."""
        if current_metaexons == []:
            _select_chrom1, _select_chrom2 = secrets.SystemRandom().sample(
                self.available_chroms, k=2
            )
            _select_gene1, _select_strand1 = secrets.choice(
                self.chrom_to_genes[_select_chrom1]
            )
            _select_gene2, _select_strand2 = secrets.choice(
                self.chrom_to_genes[_select_chrom2]
            )
            _metaexon1 = self._metaexon(
                locus=_select_gene1,
                strand=_select_strand1,
                locus_type="exonic",
                direction="upstream",
                nls_type="TRA",
            )

            _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
            _metaexon2 = self._metaexon(
                locus=_select_gene2,
                strand=_select_strand2,
                locus_type=_locus_type,
                direction="downstream",
            )
            current_metaexons.append(_metaexon1)
            current_metaexons.append(_metaexon2)
        else:
            current_metaexons[-1].nls_type = "TRA"
            last_metaexon = current_metaexons[-1]
            _select_chrom1 = last_metaexon.chrom
            _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
            while True:
                _select_chrom2 = secrets.choice(self.available_chroms)
                if _select_chrom1 != _select_chrom2:
                    _select_gene2, _select_strand2 = secrets.choice(
                        self.chrom_to_genes[_select_chrom2]
                    )
                    _metaexon2 = self._metaexon(
                        locus=_select_gene2,
                        strand=_select_strand2,
                        locus_type=_locus_type,
                        direction="downstream",
                    )
                    current_metaexons.append(_metaexon2)
                    break
            return current_metaexons

    def _del_hopper(self, current_metaexons) -> None:
        """DEL hopper."""
        if current_metaexons == []:
            _select_chrom = secrets.choice(self.available_chroms)
            while True:
                _select_gene1, _select_strand1 = secrets.choice(
                    self.chrom_to_genes[_select_chrom]
                )
                _index1 = self.chrom_to_genes[_select_chrom].index(
                    (_select_gene1, _select_strand1)
                )
                _metaexon1 = self._metaexon(
                    locus=_select_gene1,
                    strand=_select_strand1,
                    locus_type="exonic",
                    direction="upstream",
                    nls_type="DEL",
                )

                _shift = secrets.choice(range(2, 5))
                _locus_type = "exonic"
                if _select_strand1 == "+":
                    _index2 = _index1 + _shift
                else:
                    _index2 = _index1 - _shift
                _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                    _index2
                ]
                if _select_strand1 == _select_strand2:
                    _metaexon2 = self._metaexon(
                        locus=_select_gene2,
                        strand=_select_strand2,
                        locus_type=_locus_type,
                        direction="downstream",
                    )
                    current_metaexons.append(_metaexon1)
                    current_metaexons.append(_metaexon2)
                    break
        else:
            current_metaexons[-1].nls_type = "DEL"
            last_metaexon = current_metaexons[-1]
            _select_chrom = last_metaexon.chrom
            _select_gene1, _select_strand1 = last_metaexon.locus, last_metaexon.strand
            _index1 = self.chrom_to_genes[_select_chrom].index(
                (_select_gene1, _select_strand1)
            )

            while True:
                _shift = secrets.choice(range(2, 12))
                _locus_type = "exonic"
                if _select_strand1 == "+":
                    _index2 = _index1 + _shift
                else:
                    _index2 = _index1 - _shift
                _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                    _index2
                ]
                if _select_strand1 == _select_strand2:
                    _metaexon2 = self._metaexon(
                        locus=_select_gene2,
                        strand=_select_strand2,
                        locus_type=_locus_type,
                        direction="downstream",
                    )
                    current_metaexons.append(_metaexon2)
                    break
            return current_metaexons

    def _microhomology_checker(self) -> bool:
        """Check Microhomology."""
        pass

    def one_hop_generator(self, hop_type: str) -> List[MetaExon]:
        """Generate one-hop NLS event.

        .. note:: candidate_hop_types = ["TDUP", "IDUP", "INV", "TRA"]
        """
        total_metaexons: List[MetaExon] = []
        _select_type = hop_type
        if _select_type == "TDUP":
            self._tdup_hopper(total_metaexons)
        elif _select_type == "IDUP":
            self._idup_hopper(total_metaexons)
        elif _select_type == "INV":
            self._inv_hopper(total_metaexons)
        elif _select_type == "TRA":
            self._tra_hopper(total_metaexons)
        return total_metaexons

    def multi_hop_generator(self, num_of_hops: int) -> list:
        """Generate multi-hop NLS event.

        .. note::
              All the hop types are 'DEL' is not allowed.
        """
        candidate_hop_types = ["TDUP", "IDUP", "INV", "TRA", "DEL"]
        while True:
            total_metaexons: List[MetaExon] = []
            hops_type_list = []
            for _hop_idx in range(num_of_hops):
                _select_type = secrets.choice(candidate_hop_types)
                hops_type_list.append(_select_type)
                if _select_type == "TDUP":
                    self._tdup_hopper(total_metaexons)
                elif _select_type == "IDUP":
                    self._idup_hopper(total_metaexons)
                elif _select_type == "INV":
                    self._inv_hopper(total_metaexons)
                elif _select_type == "TRA":
                    self._tra_hopper(total_metaexons)
                elif _select_type == "DEL":
                    self._del_hopper(total_metaexons)
            if hops_type_list.count("DEL") < num_of_hops:
                break
        return total_metaexons

    def transcripts_generator(
        self, num_of_hops: int, num_of_transcripts: int, hop_type: str
    ) -> dict:
        """Generate transcripts."""
        transcripts_dict = {}
        if num_of_hops < 1:
            raise SystemExit from NumberOfHopIsNotValidError
        elif num_of_hops == 1:  # user must provide hop_type
            for trx_idx in range(num_of_transcripts):
                transcripts_dict[f"nls_{trx_idx}"] = self.one_hop_generator(hop_type)
        else:  # hop number > 1
            for trx_idx in range(num_of_transcripts):
                transcripts_dict[f"nls_{trx_idx}"] = self.multi_hop_generator(
                    num_of_hops
                )
        return transcripts_dict


class SimVCFWriter(Writer):
    """Writer for VCF files for Simulated data.

    .. note::

        1. CHROM: The name of the sequence (typically a chromosome) on which the variation
            is being called. This sequence is usually known as 'the reference sequence',
            i.e. the sequence against which the given sample varies.
        2. POS: The 1-based position of the variation on the given sequence.
        3. ID: The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown
            a ".". Multiple identifiers should be separated by semi-colons without white-space.
        4. REF:The reference base (or bases in the case of an indel) at the given position
            on the given reference sequence.
        5. ALT: The list of alternative alleles at this position.
        6. QUAL: A quality score associated with the inference of the given alleles.
        7. FILTER: A flag indicating which of a given set of filters the variation has
            failed or PASS if all the filters were passed successfully.
        8. INFO: An extensible list of key-value pairs (fields) describing the variation.
            See below for some common fields. Multiple fields are separated by semicolons
            with optional values in the format: <key>=<data>[,data].
        9. FORMAT: An (optional) extensible list of fields for describing the samples.
            See below for some common fields.
        10. SAMPLE: For each (optional) sample described in the file,
            values are given for the fields listed in FORMAT
    """

    num_fields = 10

    reserved_info = {
        "SVMETHOD": "String",
        "SVTYPE": "String",
        "STRAND1": "String",
        "STRAND2": "String",
        "SVLEN": "Integer",
        "CHR2": "String",
        "END": "Integer",
    }
    reserved_format = {"GT": "String"}
    reserved_alt = ["INS", "DEL", "TDUP", "IDUP", "INV", "TRA"]

    description = {
        "CANONICAL": "Canonical splice site",
        "NONCANONICAL": "Noncanonical splice site",
        "BOUNDARY": "The coding exon boundary type of event, BOTH, LEFT, RIGHT, NEITHER.",
        "DP": "Total read depth at the breakpoint for insertion",
        "DP1": "Total read depth at the breakpoint1",
        "DP2": "Total read depth at the breakpoint2",
        "SR": "The number of support reads for the breakpoints",
        "AO": "Alternate allele observations, "
        "with partial observations recorded fractionally",
        "AF": "Estimated allele frequency in the range (0,1], "
        "representing the ratio of reads showing the alternative allele to all reads",
        "PSO": "Estimated Percent splice-out in the range (0,1], "
        "representing the percentage of NLS transcripts",
        "SVTYPE": "The type of event, INS, DEL, TDUP, IDUP, INV, TRA.",
        "SVLEN": "Difference in length between REF and ALT alleles",
        "CHR2": "Chromosome for END coordinate in case of a translocation",
        "END": "2nd position of the structural variant",
        "GENE": "Overlapped coding gene for insertion",
        "GENE1": "Overlapped coding gene for breakpoint1",
        "GENE2": "Overlapped coding gene for breakpoint2",
        "TRANSCRIPT_ID": "Transcript ID",
        "SVMETHOD": "Type of approach used to detect SV",
        "STRAND": "Strand for insertion",
        "STRAND1": "Strand for breakpoint1",
        "STRAND2": "Strand for breakpoint2",
        "MODE1": "Mode for softclipped reads at breakpoint1",
        "MODE2": "Mode for softclipped reads at breakpoint2",
        "GT": "Genotype",
        "INS": "Insertion",
        "DEL": "Deletion",
        "TDUP": "Tandem duplication",
        "IDUP": "Inverted duplication",
        "INV": "Inversion",
        "TRA": "Translocation",
    }

    def __init__(
        self,
        file_path: str,
        output_prefix: str,
        logger: LoggerType,
    ) -> None:
        """Initialize SimVCFWriter object."""
        super().__init__(file_path, logger)
        self.id = 1
        self.sample_name = output_prefix

    @property
    def is_opened(self) -> bool:
        """Check if file is opened."""
        return self.io is not None and not self.io.closed

    def formatter(self, fields: List[str], delimiter: str = "\t") -> str:
        """Formatter for writing data."""
        if fields is None or len(fields) != SimVCFWriter.num_fields:
            self.logger.warning(
                f"{self.__class__.__name__}: Number of fields is not equal to 10."
            )
            raise SystemExit
        return delimiter.join(fields) + "\n"

    def open(self, mode: str = "w") -> IO:
        """Open file."""
        if self.is_opened:
            self.logger.warning(f"{self.__class__.__name__}: File is already opened.")
        with open(self.file_path, mode) as self.io:
            pass
        return self.io

    def close(self) -> None:
        """Close file."""
        if self.is_opened:
            self.io.close()  # type: ignore
            self.io = None

    def write_line(self, line: str) -> None:
        """Write line to file."""
        if self.is_opened:
            self.io.write(line)
        else:
            self.logger.warning(f"{self.__class__.__name__}: File is not opened.")

    def write_header(self) -> None:
        """Write header to VCF file."""
        self.write_line(self.header)

    @singledispatchmethod
    def write_data(self, data_object: Any) -> None:
        """Write data to file.

        :param: data_object: Data to write to file.
        """

    @write_data.register
    def _(self, data_object: List[MetaExon]) -> None:
        """Write metaexons to VCF file.

        :param data_object: metaexons to write to file.
        """
        if len(data_object.nodes) == 0:
            self.logger.warning(
                f"{self.__class__.__name__}: No nodes to write to VCF file."
            )
        self.logger.trace(f"{self.__class__.__name__}: Writing metaexons to VCF file.")
        for hop_vcf_feature in get_vcf_features_from_metaexons(data_object, self.id):
            self.write_line(self.formatter(hop_vcf_feature))
        self.id += 1

    @property
    def header(self) -> str:
        """VCF header provides metadata describing the body of the file."""
        # Metadata parsers/constants

        date = datetime.datetime.today().strftime("%Y%m%d")
        source = f"ScanNLS v{__version__}"

        header_lines = [
            "##fileformat=VCFv4.3",
            f"##fileDate={date}",
            f"##source={source}",
        ]

        for _id in SimVCFWriter.reserved_info:
            _number = 0 if SimVCFWriter.reserved_info[_id] == "Flag" else 1
            header_lines.append(
                f"##INFO=<ID={_id},Number={_number},Type={SimVCFWriter.reserved_info[_id]},"
                f'Description="{SimVCFWriter.description[_id]}">'
            )

        for _id in SimVCFWriter.reserved_format:
            header_lines.append(
                f"##FORMAT=<ID={_id},Number=1,Type={SimVCFWriter.reserved_format[_id]},"
                f'Description="{SimVCFWriter.description[_id]}">'
            )

        for _id in SimVCFWriter.reserved_alt:
            header_lines.append(
                f'##ALT=<ID={_id},Description="{SimVCFWriter.description[_id]}">'
            )
        header_lines.append(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}"
        )

        return "\n".join(header_lines) + "\n"


def get_vcf_features_from_metaexons(
    series: List[MetaExon],
    series_id: int,
) -> List[List[str]]:
    """Obtain hop vcf features from one list of metaexons."""
    series_hops_features = []

    for event_id, current_node in enumerate(series[:-1], 1):
        next_node = series[event_id]

        _chrom1 = current_node.chrom
        _chrom2 = next_node.chrom

        _strand1, _pos1 = current_node.strand, current_node.p3_pos
        _strand2, _pos2 = next_node.strand, next_node.p5_pos

        sv_distance = abs(_pos1 - _pos2) if current_node.nls_type != "TRA" else 0

        series_hops_features.append(
            [
                _chrom1,
                f"{int(_pos1) + 1}",
                f"HOP_{event_id}",
                ".",
                f"<{current_node.nls_type}>",
                ".",
                ".",
                (
                    f"SVTYPE={current_node.nls_type};"
                    f"CHR2={_chrom2};END={int(_pos2) + 1};"
                    f"SVLEN={sv_distance};"
                    f"STRAND1={_strand1};STRAND2={_strand2};"
                    f"TRANSCRIPT_ID={series_id};SVMETHOD=ScanNLS_Simulator"
                ),
                "GT",
                "0/1",
            ]
        )
        # current_node insertion_seq: #TODO
    return series_hops_features
