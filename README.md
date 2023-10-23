# ScanNLS: A powerful tool for detecting non-co-linear transcripts with long reads and nonlinear splicing graphs

[![pypi](https://img.shields.io/pypi/v/scannls.svg?style=for-the-badge)][pypi]
![conda](https://img.shields.io/conda/vn/bioconda/scannls?style=for-the-badge)
[![publication](https://img.shields.io/badge/published%20in-Nature-green.svg?style=for-the-badge)][paper]

[paper]: https://www.nature.com/articles/d41586-023-03067-6

## What is ScanNLS?

ScanNLS is a non-co-linear transcript caller for third-generation sequencing reads.
It is able to detect and classify the non-co-linear transcripts with the following four forms of nonlinear splicing junctions: reversed junctions, inverted junctions, translocated junctions (same strand), and translocated junctions (different strand) (see the figure below).

## 🚀 **Getting Started**

The first step in starting your journey with `ScanNLS` is to install the tool.
To do this, there are two options shown below:

- **PyPI**

```bash
pip install scannls
```

- **CONDA** via [Bioconda](https://bioconda.github.io/)

```bash
conda install scannls
```

Congratulations! You've successfully installed `ScanNLS` on your local machine.
If you have some issues, please check the [document](https://scannls.readthedocs.io/en/latest/installation.html) first before opening an issue.

### 🤖 **Using ScanNLS**

```console
❯ scannls -h

 sage: scannls [-h] [--version] --input INPUT --ref REF --gtf GTF --output OUTPUT [--sr SUPPORT_READS] [--splice-bin SPLICE_BIN] [--mapq MAPQ]
               [--log-level {info,debug,trace}] [--parallel PARALLEL] [--2bit TWO_BIT] [--species {human,mouse}]
               [--circular-rna-filter {remove,keep,extract}] [--off-exon-filter] [--rt-switching-filter RT_SWITCHING_FILTER_LEN] [--ncan] [--nclosed]
               [--nsleep] [--graph] [--nbound] [--port PORT] [--max-allowed-nm MAX_ALLOWED_NM] [--identity IDENT_CUTOFF]
               [--long-indel-length LONG_INDEL_LENGTH] [--substitution-num SUBSTITUTIONS_NUM] [--indel-fraction INDEL_FRACTION]
               [--prune-threshold PRUNE_THRESHOLD] [--soft-len SOFT_LEN] [--mismatch MISMATCH] [--min-soft-seg-len MIN_SOFT_SEG_LEN]
               [--alignment-fraction ALIGNMENT_FRACTION] [--substitution-fraction SUBSTITUTIONS_FRACTION]

scannls 🚀 Nonlinear splicing (NLS) events identification using transcriptomic long reads data

options:
  -h, --help                              show this help message and exit
  --version                               show program's version number and exit
  --input INPUT                           input BAM file
  --ref REF                               reference genome in FASTA format (with fai index)
  --gtf GTF                               gene annotations in GTF format
  --output OUTPUT                         output prefix
  --sr SUPPORT_READS                      minimum number of support reads for reporting NLS (default: 1)
  --splice-bin SPLICE_BIN                 splice site bin size (default: 5)
  --mapq MAPQ                             minimum MAPQ of reads for calling NLS (default: 15)
  --log-level {info,debug,trace}          set log level (default: info)
  --parallel PARALLEL                     set working mode in processor (default: 1)
  --2bit TWO_BIT                          reference genome in 2bit format
  --species {human,mouse}                 species name for reference genome (default: human)
  --circular-rna-filter {remove,keep,extract}
                                          The way of dealing with circular RNAs (default: remove)
  --off-exon-filter                       Turn on exon filter (default: True)
  --rt-switching-filter RT_SWITCHING_FILTER_LEN
                                          Set RT switching filter (default length: 10)
  --ncan                                  considering Non canonical spliced sites (default: False)
  --nclosed                               close BLAT server when job has done (default: True)
  --nsleep                                if sleep randomly before starting BLAT server (default: True)
  --graph                                 if output graph (default: False)
  --nbound                                if add maximum increment limit using average reads depth when rescuing sr (default: True)
  --port PORT                             port for BLAT server (default: 88888)
  --max-allowed-nm MAX_ALLOWED_NM         maximum allowed NM to keep AS tag (default: 60)
  --identity IDENT_CUTOFF                 BLAT identity cutoff (default: 0.99)
  --long-indel-length LONG_INDEL_LENGTH   the length cutoff of defining long indel in the reads (default: 5)
  --substitution-num SUBSTITUTIONS_NUM    the allowed maximum substitution number in the reads (default: 20)
  --indel-fraction INDEL_FRACTION         the allowed maximum long indel fraction in the reads (default: 0.1)
  --prune-threshold PRUNE_THRESHOLD       splice graph pruning length threshold (default: 3)
  --soft-len SOFT_LEN                     minimum softclipped segment length to be rescued (default: 5)
  --mismatch MISMATCH                     maximum allowed mismatch bases of rescued segment (default: 3)
  --min-soft-seg-len MIN_SOFT_SEG_LEN     minimum softclipped segment length to trigger BLAT alignment (default: 200)
  --alignment-fraction ALIGNMENT_FRACTION
                                          minimal fraction of aligned part for smith waterman local alignment (default: 0.8)
  --substitution-fraction SUBSTITUTIONS_FRACTION
                                          the allowed maximum substitution fraction in the reads (default: 0.1)
```

Please see the [document](https://scannls.readthedocs.io/en/latest/) for details and more examples.

## 📎 **Citation**

Feel free to read and cite our paper in [Nature](https://www.nature.com/articles/d41586-023-03067-6).

## Contributing

Contributions are very welcome. To learn more, see the [Contributor Guide].

## 🪪 **License**

The project is licensed under the GNU General Public License.

## 🤝 **Contact**

If you experience any problems or have suggestions please create an issue or a pull request.

## Credits

[mit license]: https://opensource.org/licenses/mit
[pypi]: https://pypi.org/
[hypermodern python cookiecutter]: https://github.com/cjolowicz/cookiecutter-hypermodern-python
[file an issue]: https://github.com/ylab-hi/ScanNLS/issues
[pip]: https://pip.pypa.io/
[contributor guide]: CONTRIBUTING.md
[command-line reference]: https://ScanNLS.readthedocs.io/en/latest/usage.html
