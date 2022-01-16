#!/usr/bin/env python
# ===============================================================================
from pathlib import Path

from pyfaidx import Fasta  # type: ignore
from pyfaidx import FastaNotFoundError

from ..type import LoggerType


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
        self.gene_to_trx = gene_to_trx
        self.gene_to_intergenic = gene_to_intergenic
        self.trx_to_exons = trx_to_exons
        self.trx_to_introns = trx_to_introns
        self.reference = Path(reference)
        if not self.reference.exists():
            raise FastaNotFoundError
        self.reference_io = Fasta(reference, sequence_always_upper=True)

    def _megaexon(self, locus: str) -> None:
        """Generate Megaexon according to locus type (exon, intron, intergenic)."""
        pass

    def _tdup_hopper(self, prev_locus) -> None:
        pass

    def _idup_hopper(self, prev_locus) -> None:
        pass

    def _inv_hopper(self, prev_locus) -> None:
        pass

    def _tra_hopper(self, prev_locus) -> None:
        pass

    def _microhomology_checker(self) -> bool:
        pass
