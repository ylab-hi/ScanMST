# Running ScanMST

Learn how to perform multi-segment transcript identificaiton with ScanMST.

!!! info "Learning Objectives"

    By the end of this tutorial, you will be able to:

    - Prepare the necessary reference files
    - Run ScanMST on your own data

    **Prerequisites**:

    - ScanMST installed
    - BAM file prepared, see the [Reads Alignment Tutorial](reads-alignment.md).
    - Basic command-line experience

    **Time**: Approximately 30 minutes to several hours, depending on the size of the BAM file

!!! tip "New to ScanMST?"

    Try the [Quick Start](../getting-started/quick-start.md) first — it runs a complete analysis on bundled example data in about a minute, with no genome download, and walks through how to read each output file.

## Necessary files preparation

```bash
# Download reference genome FASTA
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Download reference genome 2bit file
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

# Download reference annotation GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
gunzip gencode.v48.annotation.gtf.gz
```

## Begin to run ScanMST

```bash
scanmst --thread 1 --input sample.bam --output sample --ref hg38.fa --gtf gencode.v48.annotation.gtf --ncan --aligner blat --blat-2bit hg38.2bit --blat-port 88890 --log-level debug --graph --refine &> sample.log
```

### Expected Output

1. **GTF**: Store transcript segments (detailed exons) and segment links.
2. **VCF**: Aggregated segment links by positions
3. **FASTA**: Consensus transcript sequences
4. **JSON**: Transcript segment graph

## Summary

You've learned how to:

- ✅ Prepare necessary files for ScanMST
- ✅ Run ScanMST to identify MSTs.

!!! success "Identified MSTs Data Ready!"
The MSTs identified are now ready for high-quality downstream analysis!
