#!/usr/bin/env python
"""Generate NLS transcripts."""
# ===============================================================================
from .. import GTFReader
from .. import Intergenic
from .. import OneHop
from .. import SimVCFWriter
from ..type import LoggerType
from .helper import to_wt_and_mt_fasta


def prepare_intergenic_gtf(input_gtf: str, output_gtf: str, logger: LoggerType) -> str:
    """Prepare the intergenic GTF file using annotation GTF file."""
    intergenic = Intergenic(input_gtf=input_gtf, output_gtf=output_gtf, logger=logger)
    outfile = intergenic.run()
    return outfile


def single_hop_generator(
    annotation_gtf: str,
    reference: str,
    sv_type: str,
    num_of_transcripts: int,
    output_prefix: str,
    logger: LoggerType,
):
    """Generate One-hop NLS events."""
    gtf = GTFReader(input_gtf=annotation_gtf, logger=logger)
    gtf.parser()

    # initialize OneHop class using GTF information
    one_hop = OneHop(
        chrom_to_genes=gtf.chrom_to_genes,
        gene_to_trx=gtf.gene_to_trx,
        trx_to_exons=gtf.trx_to_exons,
        trx_to_introns=gtf.trx_to_introns,
        gene_to_intergenic=gtf.gene_to_intergenic,
        reference=reference,
        logger=logger,
    )

    nls_dict = one_hop.transcripts_generator(
        num_of_hops=1, num_of_transcripts=num_of_transcripts, hop_type=sv_type
    )
    to_wt_and_mt_fasta(input_dict=nls_dict, output_prefix=output_prefix)

    vcf_writer = SimVCFWriter(f"{output_prefix}.vcf", output_prefix, logger)
    with vcf_writer.open() as _:
        vcf_writer.write_header()
        for trx_idx in nls_dict:
            vcf_writer.write_data(nls_dict[trx_idx])


def multi_hop_generator(
    annotation_gtf: str,
    reference: str,
    num_of_hops: int,
    num_of_transcripts: int,
    output_prefix: str,
    logger: LoggerType,
):
    """Generate Multi-hop NLS events."""
    gtf = GTFReader(input_gtf=annotation_gtf, logger=logger)
    gtf.parser()

    # initialize OneHop class using GTF information
    one_hop = OneHop(
        chrom_to_genes=gtf.chrom_to_genes,
        gene_to_trx=gtf.gene_to_trx,
        trx_to_exons=gtf.trx_to_exons,
        trx_to_introns=gtf.trx_to_introns,
        gene_to_intergenic=gtf.gene_to_intergenic,
        reference=reference,
        logger=logger,
    )

    nls_dict = one_hop.transcripts_generator(
        num_of_hops=num_of_hops, num_of_transcripts=num_of_transcripts
    )

    to_wt_and_mt_fasta(input_dict=nls_dict, output_prefix=output_prefix)

    vcf_writer = SimVCFWriter(f"{output_prefix}.vcf", output_prefix, logger)

    with vcf_writer.open() as _:
        vcf_writer.write_header()
        for trx_idx in nls_dict:
            vcf_writer.write_data(nls_dict[trx_idx])
