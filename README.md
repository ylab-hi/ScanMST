# <img src="https://raw.githubusercontent.com/ylab-hi/ScanMST/main/images/logo.png" alt="ScanMST logo" height="100"/> [![social](https://img.shields.io/github/stars/ylab-hi/ScanMST?style=social)](https://github.com/ylab-hi/ScanMST/stargazers)

[![pypi](https://img.shields.io/pypi/v/scanmst.svg?style=for-the-badge)](https://pypi.org/project/scanmst/)
[![conda](https://img.shields.io/conda/vn/bioconda/scanmst?style=for-the-badge)](https://anaconda.org/channels/bioconda/packages/scanmst/overview)
<!-- [![publication](https://img.shields.io/badge/published%20in-Nature-green.svg?style=for-the-badge)][paper]
[paper]: https://www.nature.com/articles/d41586-023-03067-6
-->

## What is ScanMST?

ScanMST is a powerful tool for detecting Multi-segment transcripts (MSTs) with long reads and transcript segment graphs.
It is able to detect and classify the multi-segment transcripts with the following four forms of transcript segment links: ICRL, ICTL, ITPL, and ITTL (see the figure below).

<div align="center">
<img src="https://raw.githubusercontent.com/ylab-hi/ScanMST/main/images/segment_links.png" alt="Modeling segment connectivity" width="60%">
</div>


## 🧬 BLAT Aligner (Automatic Setup)
ScanMST utilizes BLAT (BLAST-like alignment tool) for auxiliary alignments.

**You do not need to install BLAT manually.** When you run ScanMST with the `--aligner blat` option, the tool will automatically detect your operating system (Linux or macOS/Darwin) and chip architecture (Intel or Apple Silicon). It will then download the appropriate executables (gfServer, gfClient, and faToTwoBit) from the UCSC Genome Browser servers into the installation directory.

## 🚀 **Getting Started**

The first step in starting your journey with `ScanMST` is to install the tool.
To do this, there are two options shown below:

- **PyPI**

```bash
pip install scanmst
```

- **CONDA** via [Bioconda](https://bioconda.github.io/)

```bash
conda install scanmst
```

Congratulations! You've successfully installed `ScanMST` on your local machine.
If you have some issues, please check the [document](https://ylab-hi.github.io/ScanMST/) first before opening an issue.

### 🤖 **Using ScanMST**

```console
❯ scanmst -h
usage: scanmst [-h] [--version] --input INPUT --ref REF --gtf GTF --output OUTPUT [--output-seq {consensus,reference,both}] [--sr SUPPORT_READS]
               [--splice-bin SPLICE_BIN] [--mapq MAPQ] [--log-level {info,debug,trace,warning}] [--thread THREAD] [--aligner {blat,}]
               [--blat-identity IDENT_CUTOFF] [--blat-2bit BLAT_TWO_BIT] [--blat-nclosed] [--blat-nsleep] [--blat-port BLAT_PORT] [--species {human,mouse}]
               [--circular-rna-filter {remove,keep,extract}] [--off-exon-filter] [--rt-switching-filter RT_SWITCHING_FILTER_LEN] [--ncan] [--graph]
               [--refine] [--refine-threshold REFINE_THRESHOLD] [--prune-threshold PRUNE_THRESHOLD] [--max-allowed-nm MAX_ALLOWED_NM]
               [--max-allowed-ins MAX_ALLOWED_INS] [--min-required-ins MIN_REQUIRED_INS] [--min-soft-seg-len MIN_SOFT_SEG_LEN]
               [--long-indel-length LONG_INDEL_LENGTH] [--indel-fraction INDEL_FRACTION] [--substitution-fraction SUBSTITUTIONS_FRACTION] [--rescue-sr]
               [--soft-len SOFT_LEN] [--mismatch MISMATCH] [--alignment-fraction ALIGNMENT_FRACTION] [--nbound] [--ignore-circle]

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
  --graph                                 Whether to output transcript segment graph. (default: False)
  --refine                                Whether to refine the transcript segment graph after construction. (default: False)
  --refine-threshold REFINE_THRESHOLD     Threshold for merging nodes during refinement (default: 3)
  --prune-threshold PRUNE_THRESHOLD       Length threshold for pruning the transcript segment graph (default: 10)
  --max-allowed-nm MAX_ALLOWED_NM         Maximum allowed edit distance (NM tag). (default: 100)
  --max-allowed-ins MAX_ALLOWED_INS       Maximum allowed micro-insertion length (default: 50)
  --min-required-ins MIN_REQUIRED_INS     Minimum required insertion length in read to infer chimeric alignment (default: 100)
  --min-soft-seg-len MIN_SOFT_SEG_LEN     Minimum length of soft-clipped portion required to trigger BLAT alignment. (default: 200)
  --long-indel-length LONG_INDEL_LENGTH   Length cutoff for defining long indels in reads. (default: 10)
  --indel-fraction INDEL_FRACTION         Maximum allowed fraction of long indels in the reads. (default: 0.001)
  --substitution-fraction SUBSTITUTIONS_FRACTION
                                          Maximum allowed fraction of substitutions in the reads (default: 0.05)
  --rescue-sr                             Whether to rescue SR for segment links (default: False)
  --soft-len SOFT_LEN                     Minimum length of soft-clipped portion to be rescued (default: 5)
  --mismatch MISMATCH                     Maximum number of mismatched bases allowed in a rescued segment (default: 3)
  --alignment-fraction ALIGNMENT_FRACTION
                                          Minimum fraction of the sequence that must align in Smith-Waterman local alignment. (default: 0.8)
  --nbound                                Whether to add maximum increment limit using average reads depth when rescuing SR (default: True)
  --ignore-circle                         Whether to export result when the transcript segment graph contains a circle (default: False)
```

Please refer to the [document](https://ylab-hi.github.io/ScanMST/) for details and more examples.

## Contributing

Contributions are very welcome. To learn more, see the [Contributor Guide].

## 🪪 **License**

ScanMST is free software available under the GNU General Public License v3.0 (GPLv3). You are free to modify and redistribute this software under the terms of this license.

### ⚠️ External Dependency Licensing (BLAT)
While ScanMST itself is GPL-licensed, it utilizes the **BLAT** aligner for specific functionality. BLAT is **not** distributed with ScanMST; it is downloaded automatically from the University of California, Santa Cruz (UCSC) servers upon first use.

### BLAT License Terms:

- Academic/Non-Profit: Free for use.

- Commercial: A license is required from [Kent Informatics](https://kentinformatics.com/).

By using this software with the `--aligner blat` option, you acknowledge that you are responsible for adhering to the [UCSC](https://genome.ucsc.edu/license/) and [Kent Informatics](https://kentinformatics.com/) license terms regarding the use of BLAT executables.

## 🤝 **Contact**

If you experience any problems or have suggestions, please create an issue or a pull request.

## Credits


[hypermodern python cookiecutter]: https://github.com/cjolowicz/cookiecutter-hypermodern-python
[file an issue]: https://github.com/ylab-hi/ScanMST/issues
[pip]: https://pip.pypa.io/
[contributor guide]: CONTRIBUTING.md
[command-line reference]: https://ylab-hi.github.io/ScanMST/
