# ScanMST

<div class="hero" markdown>

## ScanMST: Detection of Multi-Segment Transcripts from Long Reads

A powerful tool for detecting Multi-segment transcripts (MSTs) with long reads and transcript segment graphs.

[Get Started](getting-started/quick-start.md){ .md-button .md-button--primary }
[View on GitHub](https://github.com/ylab-hi/ScanMST){ .md-button }

</div>

---

## :material-star: Key Features

<div class="feature-grid" markdown>

<div class="feature-item" markdown>
### :material-speedometer: High Accuracy
Robust detection of multi-segment transcripts at base resolution.
</div>

<div class="feature-item" markdown>
### :material-console-line: Easy to Use
A streamlined CLI with sensible defaults—start running analyses.
</div>

<div class="feature-item" markdown>
### :material-flash: Visualization-Ready JSON Output
Produce structured JSON files optimized for transcript segment graph visualization.
</div>

</div>

---

## Quick Start

Get up and running with ScanMST in under 5 minutes:

```bash
# Install ScanMST
pip install scanmst

# Run ScanMST on your long-read rna-seq data
scanmst --input input_data.bam --ref ref.fasta --gtf annotation.gtf --output output_name
```

Ready to dive in? Check out our [Quick Start Guide](getting-started/quick-start.md).

---

## What is ScanMST?

ScanMST is a multi-segment transcript caller for third-generation sequencing reads. It is able to detect and classify the multi-segment transcripts with the following four forms of transcript segment links: ICRL, ICTL, ITPL, and ITTL (see the figure below).

<div align="center">
<img src="./images/segment_links.png" alt="Modeling segment connectivity" width="60%">
</div>
______________________________________________________________________

## Citation

If you use ScanMST in your research, please cite:

```bibtex
@software{scanmst2026,
  title={Structural variation reshapes the transcriptome by generating pervasive multi-segment transcripts},
  author={Wang, Ting-You, Li, Yangyang, Liu, Qi and Yang, Rendong},
  year={2026},
  url={https://github.com/ylab-hi/ScanMST}
}
```

---

## License

ScanMST is licensed under the GNU General Public License 3.0. See [License](about/license.md) for details.
