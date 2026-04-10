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
ls tests/data/data_test.bam
```

If you installed via pip, download the sample data:

```bash
# Download sample BAM file with index
wget https://github.com/ylab-hi/ScanMST/raw/main/tests/data/data_test.bam
wget https://github.com/ylab-hi/ScanMST/raw/main/tests/data/data_test.bam.bai

# Or using curl
curl -L -o data_test.bam https://github.com/ylab-hi/ScanMST/raw/main/tests/data/data_test.bam
curl -L -o data_test.bam.bai https://github.com/ylab-hi/ScanMST/raw/main/tests/data/data_test.bam.bai

# Verify files downloaded correctly
ls -lh data_test.bam*
```

!!! tip "About the Sample Data"
    The sample file `data_test.bam` contains 10 chimeric reads sequenced using PacBio Iso-seq protocol.

## Step 2: Your First Run

Run ScanMST on the sample data:


```bash
scanmst --input data_test.bam --ref hg38.fa --gtf anno_test.gtf --output data_test
```

**Expected output**:

1. **GTF file (data_test.gtf)**: Store transcript segments (detailed exons) and segment links.
2. **VCF file (data_test.vcf)**: Aggregated segment links by positions
3. **FASTA file (data_test.fasta)**: Consensus transcript sequences
4. **JSON file (data_test_TSG0000000001_cy.json)**:  Transcript segment graph


## Step 3: Understand the Output

ScanMST creates a bunch of files (including VCF, GTF, FASTA files) per sample:


=== "GTF file"

    ```bash
    # View first 8 lines in GTF file.
    head -n 8 data_test.gtf
    ```

    **Output format** (GTF file):

    ```text
    .	scanmst	transcript	.	.	.	.	.	sr "4"; osr "4"; transcript_id "TSP2127603018"; gene_id "TSG0000000001"; extend "False";
    chr6	scanmst	exon	157823215	157823446	.	+	.	exon_id "001"; segment_id "TSN5024142628"; ptc "3"; ptf "1.0"; transcript_id "TSP2127603018"; gene_id "TSG0000000001";
    chr6	scanmst	exon	157867547	157867633	.	+	.	exon_id "002"; segment_id "TSN5024142628"; ptc "3"; ptf "1.0"; transcript_id "TSP2127603018"; gene_id "TSG0000000001";
    chr6	scanmst	exon	157873102	157873176	.	+	.	exon_id "003"; segment_id "TSN5024142628"; ptc "3"; ptf "1.0"; transcript_id "TSP2127603018"; gene_id "TSG0000000001";
    chr6	scanmst	exon	157875051	157875176	.	+	.	exon_id "004"; segment_id "TSN5024142628"; ptc "3"; ptf "1.0"; transcript_id "TSP2127603018"; gene_id "TSG0000000001";
    chr10	scanmst	exon	94820496	94820637	.	+	.	exon_id "001"; segment_id "TSN2095345182"; ptc "2"; ptf "0.6666666666666666"; transcript_id "TSP2127603018"; gene_id "TSG0000000001";
    chr10	scanmst	exon	93636995	93637062	.	+	.	exon_id "001"; segment_id "TSN7383915218"; ptc "2"; ptf "0.6666666666666666"; transcript_id "TSP2127603018"; gene_id "TSG0000000001";
    chr10	scanmst	exon	89136383	89139408	.	-	.	exon_id "001"; segment_id "TSN1003188885"; ptc "1"; ptf "0.3333333333333333"; transcript_id "TSP2127603018"; gene_id "TSG0000000001";
    ```

    **Attributes**:

    - **segment_id**: Identifier of an individual transcript segment
    - **transcript_id**: Identifier of a multi-segment transcript (MST)
    - **gene_id**: Identifier grouping MSTs that share transcript segments

=== "VCF file"

    ```bash
    # View the first transcript segment link.
    sed -n '492,493p' isoseq_data.vcf
    ```

    **Output format** (GTF file):

    ```text
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	isoseq_data
    chr6	157875177	1	.	<ITPL>	.	.	CANONICAL;LINKTYPE=ITPL;SR=7;OSR=7;CHR2=chr10;SVEND=94820496;DP1=10;DP2=7;PSI=0.452;SVLEN=0;GENE1=SNX9;GENE2=CYP2C19;SVTYPE=TRA;SEGMENT1=TSN5024142628,TSN5024142628;SEGMENT2=TSN2095345182,TSN2095345182;STRAND1=+;STRAND2=+;MODE1=MS;MODE2=SM;HOMSEQ=AGG;INSSEQ=.;TRANSCRIPT_ID=TSP2127603018,TSP1746195317;GENE_ID=TSG0000000001;SR_ID=m64135_220622_211525/166592832/ccs|m64135_220622_211525/81921303/ccs|m64135_220622_211525/2623225/ccs|m64135_220622_211525/157876553/ccs|m64135_220622_211525/6686285/ccs|m64135_220622_211525/26215277/ccs|m64135_220622_211525/114754580/ccs,m64135_220622_211525/166592832/ccs|m64135_220622_211525/81921303/ccs|m64135_220622_211525/2623225/ccs|m64135_220622_211525/157876553/ccs|m64135_220622_211525/6686285/ccs|m64135_220622_211525/26215277/ccs|m64135_220622_211525/114754580/ccs;SVMETHOD=ScanMST	GT	0/1
    ...
    ```

    **Fields**:

    - **LINKTYPE**: The type of link, ICRL, ICTL, ITPL, ITTL
    - **SR**: Supporting reads number
    - **SEGMENT1**: ID for source transcript segment
    - **SEGMENT2**: ID for target transcript segment

=== "FASTA file"

    ```bash
    # View predictions from first batch
    head -n 2 data_test.fasta
    ```

    **Output format** (FASTA file):

    ```text
    >TSP2127603018 520|142|68|3026
    GAGTAGCCGAGCGCCCAGCGGCTGGGCCTGAGCGTCGAGACTCGGGGCCGAGGCGGAGGAGCGGCCGCCGCGCCGGGGCCCAGCCGGAGCCGCCGCCCTCGCCC
    ...
    ...
    ```


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
