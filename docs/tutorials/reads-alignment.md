# Long reads RNA-seq Alignment

Learn how to perform read alignment and prepare BAM files for use with ScanMST.

!!! info "Learning Objectives"
By the end of this tutorial, you will be able to:

    - Necesary files preparation
    - Alignment file (BAM) files generation

    **Prerequisites**:

    - SAMtools installed
    - Minimap2 installed
    - Basic command-line experience

    **Time**: Approximately 30 minutes to several hours, depending on the size of the FASTQ file.

## Necessary files preparation

```bash
# Download reference genome FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Download reference annotation GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
gunzip gencode.v48.annotation.gtf.gz

# Convert GTF/GFF to BED12
paftools.js gff2bed gencode.v48.annotation.gtf > annotation.bed12
```

## PacBio Iso-Seq Data

```bash
minimap2 -t 16 -R "@RG\tID:sample\tSM:hs\tLB:ga\tPL:PACBIO" -Y --cs -ax splice:hq -uf --secondary=no --junc-bed annotation.bed12 hg38.fa FLNC_reads.fastq | samtools sort -@ 8 -O BAM -o sample.bam - && samtools index sample.bam sample.bai
```

## ONT direct RNA Data

```bash
minimap2 -t 16 -R "@RG\tID:sample\tSM:hs\tLB:ga\tPL:NANOPORE" -Y --cs -ax splice -uf -k14 --junc-bed annotation.bed12 hg38.fa directRNA.fastq | samtools sort -@ 8 -O BAM -o sample.bam - && samtools index sample.bam sample.bai
```

## ONT direct cDNA Data

The raw FASTQ file from Nanopore caller should not used, you have to prepossess with [Pychopper](https://github.com/epi2me-labs/pychopper)

```bash
minimap2 -t 16 -R "@RG\tID:sample\tSM:hs\tLB:ga\tPL:NANOPORE" -Y --cs -ax splice -uf -k14 --junc-bed annotation.bed12 hg38.fa directcDNA_pychopper.fastq | samtools sort -@ 8 -O BAM -o sample.bam - && samtools index sample.bam sample.bai
```

### Expected Output

1. **BAM**: BAM file with SA and cs tags.

### Check if the BAM file with SA and cs tags

=== "Check SA tag"

    ```bash
    samtools view isoseq_data.bam|grep "SA:Z"|less
    ```

=== "Check cs tag"

    ```bash
    samtools view isoseq_data.bam|grep "cs:Z"|less
    ```

## Summary

You've learned how to:

- ✅ Long read RNA-seq alignment for ScanMST.

!!! success "Long read alignment files are Ready!"
The alignment files (BAM format) are now ready for MST calling with ScanMST!
