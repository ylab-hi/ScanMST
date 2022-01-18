#!/usr/bin/env python
# ===============================================================================
import operator
import os
from collections import defaultdict

import HTSeq
from loguru._logger import Logger  # type: ignore

from .intergenicGTF import Intergenic


class GTFReader:
    """GTFReader."""

    # hg38
    chrm_size = {
        "chr1": 248956422,
        "chr2": 242193529,
        "chr3": 198295559,
        "chr4": 190214555,
        "chr5": 181538259,
        "chr6": 170805979,
        "chr7": 159345973,
        "chr8": 145138636,
        "chr9": 138394717,
        "chr10": 133797422,
        "chr11": 135086622,
        "chr12": 133275309,
        "chr13": 114364328,
        "chr14": 107043718,
        "chr15": 101991189,
        "chr16": 90338345,
        "chr17": 83257441,
        "chr18": 80373285,
        "chr19": 58617616,
        "chr20": 64444167,
        "chr22": 50818468,
        "chr21": 46709983,
        "chrX": 156040895,
        "chrY": 57227415,
    }

    def __init__(self, input_gtf: str, logger: Logger, input_gtf_source: str) -> None:
        """Init."""
        self.gtf = input_gtf
        self.logger = logger
        self.source = input_gtf_source
        self.chrm_to_genes = None
        self.gene_to_trx = None
        self.trx_to_exon = None
        self.trx_to_intron = None
        self.gene_to_intergenic = None

    @staticmethod
    def _gtf_parser(in_gtf):
        """Gene => exons in genomicinterval."""
        protein_coding = {"protein_coding"}
        trx_to_exon = defaultdict(list)
        gene_to_trx = defaultdict(set)
        gtf_file = HTSeq.GFF_Reader(in_gtf)
        gene_positions = []
        gene_positions_dict = {}
        available_chroms = set()
        for feature in gtf_file:
            biotype = feature.attr["gene_type"]
            gene_name = feature.attr["gene_name"]
            if biotype in protein_coding:
                if feature.type == "exon":
                    trx_id = feature.attr["transcript_id"]
                    trx_to_exon[trx_id].append(feature.iv)
                    gene_to_trx[gene_name].add(trx_id)
                if feature.type == "gene":
                    gene_positions.append(
                        {
                            "chrom": feature.iv.chrom,
                            "pos": int(feature.iv.start),
                            "strand": feature.iv.strand,
                            "gene_name": gene_name,
                        }
                    )
                    gene_positions_dict[gene_name] = (
                        int(feature.iv.start),
                        int(feature.iv.end),
                    )

                    available_chroms.add(feature.iv.chrom)
        if "chrM" in available_chroms:
            available_chroms.remove("chrM")
        if "MT" in available_chroms:
            available_chroms.remove("MT")
        # remove chrX and chrY due to the Human Pseudoautosomal Region
        if "chrX" in available_chroms:
            available_chroms.remove("chrX")
        if "chrY" in available_chroms:
            available_chroms.remove("chrY")

        gene_positions.sort(key=operator.itemgetter("chrm"))
        gene_positions.sort(key=operator.itemgetter("pos"))

        # chrom => (gene_name, strand)
        chrm_to_ordered_genes = defaultdict(list)
        for item in gene_positions:
            chrom = item["chrom"]
            gene_name = item["gene_name"]
            strand = item["strand"]
            chrm_to_ordered_genes[chrom].append((gene_name, strand))

        sorted_trx_to_exon = {}
        for trx_id in trx_to_exon:
            tmp_exons = trx_to_exon[trx_id]
            tmp_exons.sort(key=lambda x: x.start)
            sorted_trx_to_exon[trx_id] = tmp_exons
        trx_to_exon = None
        return chrm_to_ordered_genes, gene_to_trx, sorted_trx_to_exon

    @staticmethod
    def _obtain_trx_to_intron(trx_to_exon):
        """Obtain transcript to introns dictionary."""
        trx_to_intron = defaultdict(list)
        for trx_id in trx_to_exon:
            strand = trx_to_exon[trx_id][0].strand
            chrm = trx_to_exon[trx_id][0].chrom
            tmp_list = []
            for i in trx_to_exon[trx_id]:
                tmp_list.append(i.start)
                tmp_list.append(i.end)
            tmp_list.pop(0)
            tmp_list.pop(-1)
            if len(tmp_list) >= 2:
                for j, k in zip(tmp_list[0::2], tmp_list[1::2]):
                    trx_to_intron[trx_id].append(
                        HTSeq.GenomicInterval(chrm, j, k, strand)
                    )
        return trx_to_intron

    @staticmethod
    def gene_to_intergenic_parser(input_gtf):
        """Parse gene to upstream/downstream intergenic GTF file."""
        protein_coding = {"protein_coding"}
        gene_to_intergenic = {}
        gtf_file = HTSeq.GFF_Reader(input_gtf)
        for feature in gtf_file:
            biotype = feature.attr["gene_type"]
            gene_name = feature.attr["gene_name"]
            chrom = feature.iv.chrom
            if biotype in protein_coding and feature.type == "gene":
                upstream_region = feature.attr["upstream_intergenic"]
                downstream_region = feature.attr["downstream_intergenic"]
                _up_start, _up_end = upstream_region.split("-")
                _down_start, _down_end = downstream_region.split("-")
                gene_to_intergenic[gene_name] = {
                    "upstream": HTSeq.GenomicInterval(
                        chrom, int(_up_start), int(_up_end), "."
                    ),
                    "downstream": HTSeq.GenomicInterval(
                        chrom, int(_down_start), int(_down_end), "."
                    ),
                }
        return gene_to_intergenic

    def parser(self):
        """Parse the annotation GTF file to generate serval useful dictionaries."""
        self.chrm_to_genes, self.gene_to_trx, self.trx_to_exon = GTFReader._gtf_parser(
            self.gtf
        )
        self.trx_to_intron = GTFReader._obtain_trx_to_intron(self.trx_to_exon)
        gtf_name = os.path.splitext(os.path.basename(self.gtf))[0]
        intergenic = Intergenic(
            input_gtf=self.gtf,
            output_gtf=f"{gtf_name}.intergenic.gtf",
            logger=self.logger,
        )
        intergenic_gtf = intergenic.run()
        self.gene_to_intergenic = GTFReader.gene_to_intergenic_parser(intergenic_gtf)
