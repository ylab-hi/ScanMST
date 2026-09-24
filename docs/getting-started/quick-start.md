# Quick Start

Run ScanMST end to end on the bundled example data, and learn to read its output.

!!! info "What you'll learn"

    - How to run ScanMST on a real BAM file
    - What the GTF, VCF, FASTA and graph JSON outputs contain
    - How to confirm your run produced the expected result

    **Time**: ~5 minutes

## Prerequisites

- ScanMST installed ([Installation Guide](installation.md))
- Basic command-line experience

Nothing else is needed. The reads, the reference sequence and the gene annotation for this example all ship with the repository, so you do not have to download a genome to try it.

## Step 1: Get the example data

=== "From a clone"

    ```bash
    git clone https://github.com/ylab-hi/ScanMST.git
    cd ScanMST
    ls example/
    ```

=== "Download individually"

    ```bash
    BASE=https://github.com/ylab-hi/ScanMST/raw/main/example
    mkdir -p example && cd example

    wget $BASE/example.bam
    wget $BASE/example.bam.bai
    wget $BASE/example_ref.windows.fasta
    wget $BASE/example_ref.gtf
    wget $BASE/make_reference.py
    ```

!!! tip "About the example data"

    `example.bam` contains 10 chimeric PacBio Iso-Seq CCS reads from the VCaP prostate cancer cell line (25 alignment records, 15 of them supplementary), aligned to hg38 across chr12, chr16 and chr17.

    Those reads encode two multi-segment transcripts: one joining **ZDHHC7** on chr16 to chr17, and one joining **USP10** and **ZDHHC7** on chr16 to **ABCB9** on chr12. Five reads support each.

## Step 2: Build the example reference

```bash
python example/make_reference.py
```

```text
chr12: 122,961,323 bp (33,331 bp of real sequence from 122,927,993)
chr16: 85,012,535 bp (313,536 bp of real sequence from 84,699,000)
chr17: 75,787,455 bp (5,349 bp of real sequence from 75,782,107)

Wrote example/example_ref.fasta (288 MB)
```

!!! note "Why this step exists"

    ScanMST looks up reference sequence by absolute genomic coordinate, so each contig in the FASTA has to start at position 1 — even though this example only touches a few small windows of hg38. Padding those windows out produces a 288 MB file, too large to keep in git, so the repository ships only the 359 KB of real sequence (`example_ref.windows.fasta`) and this script regenerates the padding locally.

    `example/example_ref.fasta` is disposable — delete it when you are done and rerun the script whenever you need it back.

## Step 3: Run ScanMST

```bash
mkdir -p out
scanmst --input example/example.bam --output out/example \
        --ref example/example_ref.fasta --gtf example/example_ref.gtf \
        --ncan --graph --refine --prune-threshold 20
```

The output directory must already exist, which is what `mkdir -p out` is for. The run takes about a second and writes four things:

- `out/example.gtf` — transcript segments and their exons
- `out/example.vcf` — segment links, aggregated by position
- `out/example.fasta` — the reconstructed transcript sequences
- `out/graph_example/` — one transcript segment graph per gene, as JSON

!!! tip "Where the graph JSON goes"

    `--graph` writes to `graph_<input BAM name>/` next to your output prefix, so `--input example/example.bam --output out/example` produces `out/graph_example/`.

## Step 4: Understand the output

ScanMST uses four kinds of identifier throughout its output:

- **TSP** — a multi-segment transcript (MST)
- **TSN** — a single transcript segment, a node in the graph
- **TSE** — a link between two segments, an edge in the graph
- **TSG** — a gene, grouping MSTs that share segments

=== "GTF file"

    ```bash
    head -n 5 out/example.gtf
    ```

    ```text
    .	scanmst	transcript	.	.	.	.	.	sr "5"; osr "5"; transcript_id "TSP1011239472"; gene_id "TSG0000000001"; extend "False";
    chr16	scanmst	exon	85011286	85011535	.	-	.	exon_id "001"; segment_id "TSN7236494217"; ptc "1"; ptf "1.0"; transcript_id "TSP1011239472"; gene_id "TSG0000000001";
    chr16	scanmst	exon	84995922	84996007	.	-	.	exon_id "002"; segment_id "TSN7236494217"; ptc "1"; ptf "1.0"; transcript_id "TSP1011239472"; gene_id "TSG0000000001";
    chr17	scanmst	exon	75783107	75786454	.	-	.	exon_id "001"; segment_id "TSN4224182903"; ptc "1"; ptf "1.0"; transcript_id "TSP1011239472"; gene_id "TSG0000000001";
    .	scanmst	transcript	.	.	.	.	.	sr "5"; osr "5"; transcript_id "TSP1451884094"; gene_id "TSG0000000002"; extend "False";
    ```

    Each MST begins with a `transcript` line. Because an MST spans multiple loci, it has no single chromosome or coordinate, so those columns are `.` and the line carries only counts and identifiers. The `exon` lines that follow give the real coordinates, grouped by `segment_id`.

    The first transcript above has two segments: `TSN7236494217` (two exons on chr16) joined to `TSN4224182903` (one exon on chr17).

    **Attributes**:

    - **sr** — reads supporting this transcript, and **osr** the count before any rescue step
    - **segment_id** — which transcript segment the exon belongs to
    - **ptc** / **ptf** — how many paths through the graph use this segment, as a count and a fraction

=== "VCF file"

    Every segment link is one VCF record:

    ```bash
    grep -v '^##' out/example.vcf | cut -f1,2,3,5
    ```

    ```text
    #CHROM	POS	ID	ALT
    chr16	84995922	1	<ITPL>
    chr16	84700112	1	<ICTL>
    chr16	84990304	2	<ITPL>
    ```

    The ALT allele names the kind of link, and the file defines all four in its own header:

    - **ICRL** — Intra-Chromosomal Reverse Link
    - **ICTL** — Intra-Chromosomal Trans-strand Link
    - **ITPL** — InTer-chromosomal Parallel Link
    - **ITTL** — InTer-chromosomal Trans-strand Link

    The detail lives in the INFO column. Here is the first record — one long line in the file, wrapped here for readability, with the `SR_ID` read list shortened:

    ```text
    CANONICAL;LINKTYPE=ITPL;SR=5;OSR=5;CHR2=chr17;SVEND=75786455;DP1=5;DP2=5;SVLEN=0;
    GENE1=ZDHHC7;GENE2=UNK;SVTYPE=TRA;SEGMENT1=TSN7236494217;SEGMENT2=TSN4224182903;
    STRAND1=-;STRAND2=-;MODE1=SM;MODE2=MS;HOMSEQ=CT;INSSEQ=.;
    TRANSCRIPT_ID=TSP1011239472;GENE_ID=TSG0000000001;SR_ID=...;SVMETHOD=ScanMST
    ```

    So this link joins `TSN7236494217` to `TSN4224182903`, from chr16:84995922 to chr17:75786455, across a `ZDHHC7`–`UNK` gene pair, supported by 5 reads, at a canonical splice site with a 2 bp `CT` microhomology.

    **Fields**:

    - **LINKTYPE** — one of the four link types above; **SVTYPE** is the matching structural variant class (`TDUP`, `INV` or `TRA`)
    - **SR** / **OSR** — supporting reads, after and before rescue processing
    - **SEGMENT1** / **SEGMENT2** — source and target transcript segment
    - **CHR2** / **SVEND** — the other end of the link
    - **HOMSEQ** / **INSSEQ** — microhomology or microinsertion at the breakpoint
    - **SR_ID** — the names of the supporting reads

    Every INFO and ALT field is described in the VCF's own `##INFO` and `##ALT` header lines, so `grep '^##' out/example.vcf` is a complete reference.

=== "FASTA file"

    ```bash
    grep '>' out/example.fasta
    ```

    ```text
    >TSP1011239472 336|3348
    >TSP1451884094 112|332|3306
    ```

    The sequence is the full reconstructed transcript. The pipe-separated numbers after the ID are the lengths of each segment in order, so they sum to the sequence length and tell you where one segment ends and the next begins — 336 + 3348 = 3684 bases for the first transcript, and three segments totalling 3750 for the second.

=== "Graph JSON"

    ```bash
    ls out/graph_example/
    ```

    ```text
    example_TSG0000000001_cy.json  example_TSG0000000002_cy.json
    ```

    One file per gene (`TSG`), in [Cytoscape.js](https://js.cytoscape.org/) format, so it can be loaded straight into a graph viewer. Nodes are transcript segments (`TSN`) carrying their coordinates, strand, exons and supporting reads; edges are segment links (`TSE`) carrying the breakpoints, supporting read count and any microhomology. The `data` block at the top lists the paths through the graph, which is where the `TSP` transcripts come from.

## Step 5: Verify your run

The repository ships the expected output, so you can check your run exactly:

```bash
diff out/example.gtf example/expected_output/example.gtf
diff out/example.fasta example/expected_output/example.fasta
```

Both should print nothing. The VCF matches too, except for its `##fileDate` header, which is regenerated on every run:

```bash
diff <(grep -v '^##fileDate' out/example.vcf) \
     <(grep -v '^##fileDate' example/expected_output/example.vcf)
```

✅ **Success indicators**:

- [ ] `out/` contains `example.gtf`, `example.vcf` and `example.fasta`
- [ ] `out/graph_example/` contains two JSON files
- [ ] The GTF reports two transcripts, and the VCF three segment links
- [ ] All three diffs above are clean

!!! success "Congratulations!"

    You've run your first multi-segment transcript identification. :tada:

## Next steps

Ready to use your own data? Two things change: you need a full reference genome and annotation instead of the windowed example ones, and your BAM has to carry the `SA` and `cs` tags ScanMST relies on.

- [Reads Alignment](../tutorials/reads-alignment.md) — how to align reads so the BAM is compatible. The example BAM here was produced with `minimap2 -Y --cs -ax splice:hq -uf --secondary=no --junc-bed`, recorded in the `##reference` header of its VCF.
- [Running ScanMST](../tutorials/running-scanmst.md) — a full run against hg38 and GENCODE
- [CLI Reference](../reference/cli.md) — every option, including the `--ncan`, `--graph`, `--refine` and `--prune-threshold` flags used above

## Troubleshooting

Encountered an issue? Check our [Troubleshooting Guide](troubleshooting.md) for common problems and solutions.

!!! question "Need Help?"

    - :material-github: [Open an issue](https://github.com/ylab-hi/ScanMST/issues)
    - :material-chat: [GitHub Discussions](https://github.com/ylab-hi/ScanMST/discussions)
