#!/usr/bin/env python
"""Generate NLS transcripts."""
from .GTFReader import GTFReader
from .helper import to_wt_and_mt_fasta
from .intergenicGTF import Intergenic
from .oneHop import OneHop
from .writer import GTFWriter
from .writer import SimVCFWriter
from scannls import LoggerType
from scannls import Writers


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
    shift: int,
    max_len: int,
    logger: LoggerType,
):
    """Generate One-hop NLS events."""
    gtf = GTFReader(input_gtf=annotation_gtf, logger=logger)
    gtf.parser()

    if gtf.chrom_to_genes is None:
        raise SystemExit from ValueError("chrom_to_genes is None")
    # initialize OneHop class using GTF information
    one_hop = OneHop(
        chrom_to_genes=gtf.chrom_to_genes,
        gene_to_trx=gtf.gene_to_trx,
        trx_to_exons=gtf.trx_to_exons,
        trx_to_introns=gtf.trx_to_introns,
        gene_to_intergenic=gtf.gene_to_intergenic,
        reference=reference,
        shift=shift,
        max_length=max_len,
        logger=logger,
    )

    nls_dict = one_hop.transcripts_generator(
        num_of_hops=1, num_of_transcripts=num_of_transcripts, hop_type=sv_type
    )
    to_wt_and_mt_fasta(input_dict=nls_dict, output_prefix=output_prefix)

    vcf_writer = SimVCFWriter(f"{output_prefix}.vcf", output_prefix, logger)
    gtf_writer = GTFWriter(f"{output_prefix}.gtf", logger)

    writers = Writers((gtf_writer, vcf_writer))

    with writers.open() as _:
        for trx_idx in nls_dict:
            writers.write_series(nls_dict[trx_idx])


def multi_hop_generator(
    annotation_gtf: str,
    reference: str,
    num_of_hops: int,
    num_of_transcripts: int,
    output_prefix: str,
    shift: int,
    max_len: int,
    logger: LoggerType,
):
    """Generate Multi-hop NLS events."""
    gtf = GTFReader(input_gtf=annotation_gtf, logger=logger)
    gtf.parser()

    if gtf.chrom_to_genes is None:
        raise SystemExit from ValueError("chrom_to_genes is None")
    # initialize OneHop class using GTF information
    one_hop = OneHop(
        chrom_to_genes=gtf.chrom_to_genes,
        gene_to_trx=gtf.gene_to_trx,
        trx_to_exons=gtf.trx_to_exons,
        trx_to_introns=gtf.trx_to_introns,
        gene_to_intergenic=gtf.gene_to_intergenic,
        reference=reference,
        shift=shift,
        max_length=max_len,
        logger=logger,
    )

    nls_dict = one_hop.transcripts_generator(
        num_of_hops=num_of_hops, num_of_transcripts=num_of_transcripts
    )

    to_wt_and_mt_fasta(input_dict=nls_dict, output_prefix=output_prefix)

    gtf_writer = GTFWriter(f"{output_prefix}.gtf", logger)

    writers = Writers(gtf_writer)

    with writers.open() as _:
        for trx_idx in nls_dict:
            writers.write_series(nls_dict[trx_idx])

    vcf_writer = SimVCFWriter(f"{output_prefix}.vcf", output_prefix, logger)
    vcf_writer.open()
    vcf_writer.write_data(nls_dict)
    vcf_writer.close()
