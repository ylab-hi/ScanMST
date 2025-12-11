# Quick Start

Get started with ScanMST in 5 minutes! This tutorial will guide you through your first multi-segment transcript detection.

!!! info "What you'll learn"
    - How to run ScanMST on BAM files
    - Understanding ScanMST output format
    - Verifying your results

    **Time**: ~5 minutes

## Prerequisites

- ScanMST installed ([Installation Guide](installation.md))
- Basic command-line experience
- A BAM file to analyze (we'll provide sample data)
- A reference genome FASTA file
- A gene annotation GTF file

!!! tip "ScanMST needs specific tags in the BAM file"
For comprehensive guide of reads alignment, see the [Reads Alignment Tutorial](../tutorials/reads-alignment.md).

## Step 1: Get Sample Data

ScanMST includes test data in the repository. If you installed from source:

```bash
# Sample data is already available
ls tests/data/isoseq_test.bam
```

If you installed via pip, download the sample data:

```bash
# Download sample BAM file with index
wget https://github.com/ylab-hi/ScanMST/raw/main/tests/data/isoseq_test.bam
wget https://github.com/ylab-hi/ScanMST/raw/main/tests/data/isoseq_test.bam.bai

# Or using curl
curl -L -o isoseq_test.bam https://github.com/ylab-hi/ScanMST/raw/main/tests/data/isoseq_test.bam
curl -L -o isoseq_test.bam.bai https://github.com/ylab-hi/ScanMST/raw/main/tests/data/isoseq_test.bam.bai

# Verify files downloaded correctly
ls -lh isoseq_test.bam*
```

!!! tip "About the Sample Data"
    The sample file `isoseq_test.bam` contains 10 chimeric reads sequenced using PacBio Iso-seq protocol.

## Step 2: Your First Run

Run ScanMST on the sample data:

=== "CPU Mode"

    ```bash
    scanmst --input isoseq_test.bam --ref ref_test.fasta --gtf anno_test.gtf --output test
    ```

    **Expected output**:
    ```text

    ```

    MSTs are saved to: `test.vcf, test.gtf, test.fasta`

## Step 3: Understand the Output

ScanMST creates a bunch of files (VCF, GTF, FASTA) per sample:

```bash
# View predictions from first batch
head -10 test.gtf
```

**Output format** (GTF file):

```text

.	scanmst	transcript	.	.	.	.	.	sr "1"; osr "1"; transcript_id "TSP4137634832"; gene_id "TSG0000000001"; extend "False";
chr10	scanmst	exon	100267571	100267671	.	-	.	exon_id "001"; segment_id "TSN3076341880"; ptc "1"; ptf "1.0"; transcript_id "TSP4137634832"; gene_id "TSG0000000001";
chr10	scanmst	exon	100261979	100262063	.	-	.	exon_id "002"; segment_id "TSN3076341880"; ptc "1"; ptf "1.0"; transcript_id "TSP4137634832"; gene_id "TSG0000000001";
...
```

**Labels**:

- **segment_id**:
- **transcript_id**:
- **gene_id**:


## Checkpoint: Verify Your Identifications Worked

✅ **Success indicators**:

- [ ] VCF, GTF and FASTA files created
- [ ] Files are not empty

!!! success "Congratulations!"
    You've successfully run your first ScanMST identification! :tada:

## Troubleshooting

Encountered an issue? Check our [Troubleshooting Guide](troubleshooting.md) for common problems and solutions.

!!! question "Need Help?"
    - :material-github: [Open an issue](https://github.com/ylab-hi/ScanMST/issues)
    - :material-chat: [GitHub Discussions](https://github.com/ylab-hi/ScanMST/discussions)
