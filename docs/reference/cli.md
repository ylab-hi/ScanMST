# CLI Commands Reference

Complete reference for all ScanMLST command-line interface parameters.

## Overview

______________________________________________________________________

```console
❯ scanmst -h

usage: scanmst [-h] [--version] --input INPUT --ref REF --gtf GTF --output OUTPUT [--output-seq {consensus,reference,both}] [--sr SUPPORT_READS]
               [--splice-bin SPLICE_BIN] [--mapq MAPQ] [--log-level {info,debug,trace,warning}] [--parallel PARALLEL] [--aligner {blat,}]
               [--blat-identity IDENT_CUTOFF] [--blat-2bit BLAT_TWO_BIT] [--blat-nclosed] [--blat-nsleep] [--blat-port BLAT_PORT] [--species {human,mouse}]
               [--circular-rna-filter {remove,keep,extract}] [--off-exon-filter] [--rt-switching-filter RT_SWITCHING_FILTER_LEN] [--ncan] [--graph] [--refine]
               [--nbound] [--max-allowed-nm MAX_ALLOWED_NM] [--max-allowed-ins MAX_ALLOWED_INS] [--min-required-ins MIN_REQUIRED_INS] [--long-indel-length LONG_INDEL_LENGTH]
               [--indel-fraction INDEL_FRACTION] [--prune-threshold PRUNE_THRESHOLD] [--soft-len SOFT_LEN]
               [--mismatch MISMATCH] [--min-soft-seg-len MIN_SOFT_SEG_LEN] [--alignment-fraction ALIGNMENT_FRACTION]
               [--substitution-fraction SUBSTITUTIONS_FRACTION] [--ignore-circle] [--rescue-sr]

scanmst 🚀 Multi-segment transcript (MST) identification using transcriptomic long reads data

options:
  -h, --help                              show this help message and exit
  --version                               show program's version number and exit
  --input INPUT                           Input alignment BAM file, which must contain both cs and SA tags.
  --ref REF                               Reference genome in FASTA format (with fai index)
  --gtf GTF                               Gene annotations in GTF format
  --output OUTPUT                         Output file prefix
  --output-seq {consensus,reference,both}
                                          Output sequence type (default: consensus)
  --sr SUPPORT_READS                      The minimum number of supporting reads required for calling MST. (default: 1)
  --splice-bin SPLICE_BIN                 Bin size for searching canonical splice sites. (default: 5)
  --mapq MAPQ                             Minimum MAPQ of reads required for calling MST. (default: 20)
  --log-level {info,debug,trace,warning}  Set log level (default: warning)
  --thread THREAD                         Set the thread number (default: 1)
  --aligner {blat,}                       Aligner used for additional realignment to recover missing chimeric alignments. (default: None)
  --blat-identity IDENT_CUTOFF            BLAT identity cutoff (default: 0.9)
  --blat-2bit BLAT_TWO_BIT                Reference genome in 2bit format for BLAT aligner
  --blat-nclosed                          Close BLAT server when the job is complete (default: True)
  --blat-nsleep                           Whether to sleep randomly before starting BLAT server (default: True)
  --blat-port BLAT_PORT                   Port for BLAT server (default: 88888)
  --species {human,mouse}                 Name of the species for the reference genome (default: human)
  --circular-rna-filter {remove,keep,extract}
                                          The way of dealing with putative circular RNAs (default: remove)
  --off-exon-filter                       Turn on exon filter (default: True)
  --rt-switching-filter RT_SWITCHING_FILTER_LEN
                                          Set the length threshold for RT switching filter. (default length: 10)
  --ncan                                  Considering non-canonical splice sites (default: False)
  --graph                                 Whether to output graph (default: False)
  --refine                                Whether to refine the graph (default: False)
  --nbound                                Whether to add maximum increment limit using average reads depth when rescuing SR (default: True)
  --max-allowed-nm MAX_ALLOWED_NM         Maximum allowed edit distance (NM tag). (default: 100)
  --max-allowed-ins MAX_ALLOWED_INS       Maximum allowed micro-insertion length (default: 50)
  --min-required-ins MIN_REQUIRED_INS     Minimum required insertion length in read to infer chimeric alignment (default: 100)
  --long-indel-length LONG_INDEL_LENGTH   Length cutoff for defining long indels in reads. (default: 10)
  --indel-fraction INDEL_FRACTION         Maximum allowed fraction of long indels in the reads. (default: 0.001)
  --prune-threshold PRUNE_THRESHOLD       Length threshold for pruning the transcript segment graph (default: 10)
  --soft-len SOFT_LEN                     Minimum length of soft-clipped portion to be rescued (default: 5)
  --mismatch MISMATCH                     Maximum number of mismatched bases allowed in a rescued segment (default: 3)
  --min-soft-seg-len MIN_SOFT_SEG_LEN     Minimum length of soft-clipped portion required to trigger BLAT alignment. (default: 200)
  --alignment-fraction ALIGNMENT_FRACTION
                                          Minimum fraction of the sequence that must align in Smith-Waterman local alignment. (default: 0.8)
  --substitution-fraction SUBSTITUTIONS_FRACTION
                                          Maximum allowed fraction of substitutions in the reads (default: 0.05)
  --ignore-circle                         Whether to export result when the transcript segment graph contains a circle (default: False)
  --rescue-sr                             Whether to rescue SR for segment links (default: False)

```

### `--version`

Display ScanMST version information.

```bash
scanmst --version
```

**Output:**

```text
scanmst 0.1.7
```

### `--help`, `-h`

Display help information for all commands.

```bash
scanmst --help
```

## Essential Arguments

### `--input`

Input the alignment BAM file, which must have cs and SA tags in it.

### `--ref`

Reference genome FASTA file.

### `--gtf`

Gene annotation GTF file

### `--output`

Output file prefix

## Optional Arguments

### `--output-seq`

Output sequence type

**Default:** `consensus`

 -  `consensus`: consensus sequences derived from reads
 -  `reference`: sequences based on the reference genome
 -  `both`: consensus and reference sequences

### `--sr`

**Default:** `1`

The minimum number of supporting reads (SRs) required to report an MST.

It is defined at the transcript level as the minimum SR across all its segment links, each segment link has a corresponding SR and must meet this minimum threshold.

### `--splice-bin`

**Default:** `5`

Bin size for searching canonical splice sites.

### `--aligner`

**Default:** `None`

Aligner used for additional realignment to recover missing chimeric alignments.

______________________________________________________________________

## See Also

- [Quick Start Tutorial](../getting-started/quick-start.md)
- [Running ScanMST](../tutorials/running-scanmst.md)
- [Reads alignment Tutorial](../tutorials/reads-alignment.md)
