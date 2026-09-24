# !/usr/bin/env python
"""End-to-end tests for the CLI.
"""
import json
import re
from pathlib import Path

import pytest
from scanmst.cli import DefaultOptions, cli

DATA_DIR = Path(__file__).parent.parent / "data"

WINDOW_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")


def _read_fasta(path):
    seqs, name, chunks = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks)
            name, chunks = line[1:], []
        else:
            chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


@pytest.fixture(scope="module")
def padded_ref(tmp_path_factory):
    """Reconstruct an N-padded reference FASTA (preserving absolute hg38
    coordinates, since pyfaidx.Fasta is always indexed by absolute position)
    from the small windows-only fixture, so the multi-hundred-MB padded file
    never has to be checked into the repo."""
    windows = _read_fasta(DATA_DIR / "example_ref.windows.fasta")
    dest = tmp_path_factory.mktemp("ref") / "example_ref.fasta"
    with open(dest, "w") as out:
        for header, seq in windows.items():
            match = WINDOW_RE.match(header)
            chrom, start = match["chrom"], int(match["start"])
            padded = ("N" * (start - 1)) + seq
            out.write(f">{chrom}\n")
            for i in range(0, len(padded), 60):
                out.write(padded[i : i + 60] + "\n")
    return dest


@pytest.mark.parametrize("thread", [1, 8])
def test_cli_example(tmp_path, padded_ref, thread):
    """Run the CLI against example.bam and compare against the checked-in expected output."""
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    op = DefaultOptions(
        input=str(DATA_DIR / "example.bam"),
        output=str(out_dir / "example"),
        ref=str(padded_ref),
        gtf=str(DATA_DIR / "example_ref.gtf"),
        blat_two_bit=str(DATA_DIR / "dummy.2bit"),
        noncanonical=True,
        log="trace",
        substitutions_fraction=0.05,
        indel_fraction=0.001,
        prune_threshold=20,
        circular_rna="remove",
        graph=True,
        refine=True,
        thread=thread,
    )
    cli.cli(op)

    # FASTA / GTF: fully deterministic, compare byte-for-byte
    assert (out_dir / "example.fasta").read_text() == (DATA_DIR / "example.fasta").read_text()
    assert (out_dir / "example.gtf").read_text() == (DATA_DIR / "example.gtf").read_text()

    # VCF: skip ##-meta lines (##fileDate is regenerated every run); compare the rest
    def data_lines(path):
        return [line for line in path.read_text().splitlines() if not line.startswith("##")]

    assert data_lines(out_dir / "example.vcf") == data_lines(DATA_DIR / "example.vcf")

    # --graph output: graph_{input_bam_stem}/ dir, one *_cy.json per TSG cluster
    out_graph_dir = out_dir / "graph_example"
    expect_graph_dir = DATA_DIR / "graph_example"
    out_files = sorted(p.name for p in out_graph_dir.glob("*_cy.json"))
    expect_files = sorted(p.name for p in expect_graph_dir.glob("*_cy.json"))
    assert out_files == expect_files
    for name in expect_files:
        out_json = json.loads((out_graph_dir / name).read_text())
        expect_json = json.loads((expect_graph_dir / name).read_text())
        assert out_json == expect_json
