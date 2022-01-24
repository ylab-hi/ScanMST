#!/usr/bin/env python
# ===============================================================================
import secrets
from dataclasses import dataclass
from pathlib import Path
from typing import List
from typing import Optional

from pyfaidx import Fasta  # type: ignore
from pyfaidx import FastaNotFoundError

from ..type import LoggerType
from .basicClass import reverse_complement
from .exception import NumberOfHopIsNotValidError


@dataclass
class MetaExon:
    """Class to store the information of one meta-exon."""

    chrom: str = None
    locus: str = None
    locus_strand: str = None
    strand: str = None
    exons: list = None
    p5_pos: int = None
    p3_pos: int = None
    wt_seq: str = None
    mt_seq: str = None
    nls_type: Optional[str] = None

    def __repr__(self) -> str:
        """Get a string representation of a MetaExon."""
        exons_repr = "|".join([f"{i.start}-{i.end}" for i in self.exons]) if self.exons else "None"  # type: ignore
        return (
            (
                f"MetaExon({self.chrom}:{self.exons[0].start}-{self.exons[-1].end}:{self.strand}, "
                f"{exons_repr}, {self.nls_type}, 5p:{self.p5_pos}, "
                f"3p:{self.p3_pos}, "
                f"locus={self.locus}:{self.locus_strand}) "
            )
            if self.exons
            else "None"
        )


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

    def gen_metaexon(
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
        _trxs = list(self.gene_to_trx[locus])
        _tgt_trx = secrets.choice(_trxs)
        _chrom = self.trx_to_exons[_tgt_trx][0].chrom
        _num_exon = len(self.trx_to_exons[_tgt_trx])
        # We only consider transcript with >= 3 exons
        if _num_exon <= 2:
            return MetaExon()

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
        # intergenic wt_seq == ''

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
            _num_intron = len(self.trx_to_introns[_tgt_trx])
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
            locus_strand=strand,
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
        print(f"{input_metaexon.mt_seq=}")
        return MetaExon(
            chrom=input_metaexon.chrom,
            locus=input_metaexon.locus,
            locus_strand=input_metaexon.locus_strand,
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
            gene_num_on_select_chrom = len(self.chrom_to_genes[_select_chrom])
            while True:
                _select_gene1, _select_strand1 = secrets.choice(
                    self.chrom_to_genes[_select_chrom]
                )
                _index1 = self.chrom_to_genes[_select_chrom].index(
                    (_select_gene1, _select_strand1)
                )
                _metaexon1 = self.gen_metaexon(
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
                if 0 <= _index2 < gene_num_on_select_chrom:
                    _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                        _index2
                    ]
                    if _metaexon1.chrom and _select_strand1 == _select_strand2:
                        _metaexon2 = self.gen_metaexon(
                            locus=_select_gene2,
                            strand=_select_strand2,
                            locus_type=_locus_type,
                            direction="downstream",
                        )
                        self.logger.trace(f"TDUP: {_metaexon2=}")
                        if _metaexon2.chrom:
                            current_metaexons.append(_metaexon1)
                            current_metaexons.append(_metaexon2)
                            break
        else:
            current_metaexons[-1].nls_type = "TDUP"
            last_metaexon = current_metaexons[-1]
            _select_chrom = last_metaexon.chrom
            _select_gene1, _select_strand1 = (
                last_metaexon.locus,
                last_metaexon.locus_strand,
            )
            _index1 = self.chrom_to_genes[_select_chrom].index(
                (_select_gene1, _select_strand1)
            )
            gene_num_on_select_chrom = len(self.chrom_to_genes[_select_chrom])

            while True:
                _shift = secrets.choice(range(2, 12))
                _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
                if _select_strand1 == "+":
                    _index2 = _index1 - _shift
                else:
                    _index2 = _index1 + _shift
                if 0 <= _index2 < gene_num_on_select_chrom:
                    _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                        _index2
                    ]
                    if _select_strand1 == _select_strand2:
                        _metaexon2 = self.gen_metaexon(
                            locus=_select_gene2,
                            strand=_select_strand2,
                            locus_type=_locus_type,
                            direction="downstream",
                        )
                        self.logger.trace(f"TDUP: {_metaexon2=}")
                        if _metaexon2.chrom:
                            current_metaexons.append(_metaexon2)
                            break
            return current_metaexons

    def _idup_hopper(self, current_metaexons) -> None:
        """IDUP hopper."""
        if current_metaexons == []:
            _select_chrom = secrets.choice(self.available_chroms)
            while True:
                _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
                _select_gene1, _select_strand1 = secrets.choice(
                    self.chrom_to_genes[_select_chrom]
                )
                _metaexon1 = self.gen_metaexon(
                    locus=_select_gene1,
                    strand=_select_strand1,
                    locus_type=_locus_type,
                    direction="upstream",
                    nls_type="IDUP",
                )
                if _metaexon1.chrom:
                    self.logger.trace(f"IDUP metaexon1: {_metaexon1=}")
                    _metaexon2 = OneHop.reverse_metaexon(_metaexon1)
                    self.logger.trace(f"IDUP: {_metaexon2=}")
                    if _metaexon2.chrom:
                        current_metaexons.append(_metaexon1)
                        current_metaexons.append(_metaexon2)
                        break
        else:
            current_metaexons[-1].nls_type = "IDUP"
            last_metaexon = current_metaexons[-1]
            _metaexon2 = OneHop.reverse_metaexon(last_metaexon)
            self.logger.trace(f"IDUP: {_metaexon2=}")
            if _metaexon2.chrom:
                current_metaexons.append(_metaexon2)

    def _inv_hopper(self, current_metaexons) -> None:
        """INV hopper."""
        if current_metaexons == []:
            _select_chrom = secrets.choice(self.available_chroms)
            gene_num_on_select_chrom = len(self.chrom_to_genes[_select_chrom])
            while True:
                _select_gene1, _select_strand1 = secrets.choice(
                    self.chrom_to_genes[_select_chrom]
                )
                _index1 = self.chrom_to_genes[_select_chrom].index(
                    (_select_gene1, _select_strand1)
                )
                _metaexon1 = self.gen_metaexon(
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
                if 0 <= _index2 < gene_num_on_select_chrom:
                    _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                        _index2
                    ]
                    if _metaexon1.chrom and _select_strand1 != _select_strand2:
                        _metaexon2 = self.gen_metaexon(
                            locus=_select_gene2,
                            strand=_select_strand2,
                            locus_type=_locus_type,
                            direction="downstream",
                        )

                        self.logger.trace(f"INV: {_metaexon2=}")
                        if _metaexon2.chrom:
                            current_metaexons.append(_metaexon1)
                            current_metaexons.append(_metaexon2)
                            break
        else:
            current_metaexons[-1].nls_type = "INV"
            last_metaexon = current_metaexons[-1]
            _select_chrom = last_metaexon.chrom
            gene_num_on_select_chrom = len(self.chrom_to_genes[_select_chrom])
            _select_gene1, _select_strand1 = (
                last_metaexon.locus,
                last_metaexon.locus_strand,
            )
            _index1 = self.chrom_to_genes[_select_chrom].index(
                (_select_gene1, _select_strand1)
            )

            while True:
                _shift = secrets.choice(
                    list(map(lambda x: -x, range(2, 12))) + list(range(2, 12))
                )
                _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
                _index2 = _index1 + _shift
                if 0 <= _index2 < gene_num_on_select_chrom:
                    _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                        _index2
                    ]
                    if _select_strand1 != _select_strand2:
                        _metaexon2 = self.gen_metaexon(
                            locus=_select_gene2,
                            strand=_select_strand2,
                            locus_type=_locus_type,
                            direction="downstream",
                        )
                        self.logger.trace(f"INV: {_metaexon2=}")
                        if _metaexon2.chrom:
                            current_metaexons.append(_metaexon2)
                            break
            return current_metaexons

    def _tra_hopper(self, current_metaexons) -> None:
        """TRA hopper."""
        if current_metaexons == []:
            _select_chrom1, _select_chrom2 = secrets.SystemRandom().sample(
                self.available_chroms, k=2
            )
            while True:
                _select_gene1, _select_strand1 = secrets.choice(
                    self.chrom_to_genes[_select_chrom1]
                )
                _select_gene2, _select_strand2 = secrets.choice(
                    self.chrom_to_genes[_select_chrom2]
                )
                _metaexon1 = self.gen_metaexon(
                    locus=_select_gene1,
                    strand=_select_strand1,
                    locus_type="exonic",
                    direction="upstream",
                    nls_type="TRA",
                )

                _locus_type = secrets.choice(["exonic", "intronic", "intergenic"])
                _metaexon2 = self.gen_metaexon(
                    locus=_select_gene2,
                    strand=_select_strand2,
                    locus_type=_locus_type,
                    direction="downstream",
                )
                self.logger.trace(f"TRA: {_metaexon2=}")
                if _metaexon1.chrom and _metaexon2.chrom:
                    current_metaexons.append(_metaexon1)
                    current_metaexons.append(_metaexon2)
                    break
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
                    _metaexon2 = self.gen_metaexon(
                        locus=_select_gene2,
                        strand=_select_strand2,
                        locus_type=_locus_type,
                        direction="downstream",
                    )
                    self.logger.trace(f"TRA: {_metaexon2=}")
                    if _metaexon2.chrom:
                        current_metaexons.append(_metaexon2)
                        break
            return current_metaexons

    def _del_hopper(self, current_metaexons) -> None:
        """DEL hopper."""
        if current_metaexons == []:
            _select_chrom = secrets.choice(self.available_chroms)
            gene_num_on_select_chrom = len(self.chrom_to_genes[_select_chrom])
            while True:
                _select_gene1, _select_strand1 = secrets.choice(
                    self.chrom_to_genes[_select_chrom]
                )
                _index1 = self.chrom_to_genes[_select_chrom].index(
                    (_select_gene1, _select_strand1)
                )
                _metaexon1 = self.gen_metaexon(
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
                if 0 <= _index2 < gene_num_on_select_chrom:
                    _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                        _index2
                    ]
                    if _metaexon1.chrom and _select_strand1 == _select_strand2:
                        _metaexon2 = self.gen_metaexon(
                            locus=_select_gene2,
                            strand=_select_strand2,
                            locus_type=_locus_type,
                            direction="downstream",
                        )
                        self.logger.trace(f"DEL: {_metaexon2=}")
                        if _metaexon2.chrom:
                            current_metaexons.append(_metaexon1)
                            current_metaexons.append(_metaexon2)
                            break
        else:
            current_metaexons[-1].nls_type = "DEL"
            last_metaexon = current_metaexons[-1]
            _select_chrom = last_metaexon.chrom
            gene_num_on_select_chrom = len(self.chrom_to_genes[_select_chrom])
            _select_gene1, _select_strand1 = (
                last_metaexon.locus,
                last_metaexon.locus_strand,
            )
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
                if 0 <= _index2 < gene_num_on_select_chrom:
                    _select_gene2, _select_strand2 = self.chrom_to_genes[_select_chrom][
                        _index2
                    ]
                    if _select_strand1 == _select_strand2:
                        _metaexon2 = self.gen_metaexon(
                            locus=_select_gene2,
                            strand=_select_strand2,
                            locus_type=_locus_type,
                            direction="downstream",
                        )
                        self.logger.trace(f"DEL: {_metaexon2=}")
                        if _metaexon2.chrom:
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
        repr_metaexons = [repr(i) for i in total_metaexons]
        self.logger.trace(f"{'; '.join(repr_metaexons)}")
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
        repr_metaexons = [repr(i) for i in total_metaexons]
        self.logger.trace(f"{';'.join(repr_metaexons)}")
        return total_metaexons

    def transcripts_generator(
        self, num_of_hops: int, num_of_transcripts: int, hop_type: str = None
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
