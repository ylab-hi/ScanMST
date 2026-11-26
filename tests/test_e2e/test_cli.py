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


def _run_example(out_dir, padded_ref, **overrides):
    op = DefaultOptions(
        input=str(DATA_DIR / "example.bam"),
        output=str(out_dir / "example"),
        ref=str(padded_ref),
        gtf=str(DATA_DIR / "example_ref.gtf"),
        blat_two_bit=str(DATA_DIR / "dummy.2bit"),
        noncanonical=True,
        substitutions_fraction=0.05,
        indel_fraction=0.001,
        prune_threshold=20,
        circular_rna="remove",
        refine=True,
        **overrides,
    )
    cli.cli(op)
    return [line for line in (out_dir / "example.vcf").read_text().splitlines() if not line.startswith("#")]


def test_cli_blacklist_without_overlap_keeps_output(tmp_path, padded_ref):
    """A blacklist that overlaps no junction must not change the result."""
    bed = tmp_path / "blacklist.bed"
    bed.write_text("chr16\t1\t100\nchrUn\t1\t1000000\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    records = _run_example(out_dir, padded_ref, blacklist=str(bed))
    expected = [line for line in (DATA_DIR / "example.vcf").read_text().splitlines() if not line.startswith("#")]
    assert records == expected


@pytest.mark.parametrize("thread", [1, 8])
def test_cli_blacklist_removes_overlapping_junction(tmp_path, padded_ref, thread):
    """A junction with a breakpoint in a blacklisted region is excluded from the output."""
    # chr16:84995922 (VCF, 1-based) is breakpoint 1 of a chr16-chr17 ITPL in example.vcf
    bed = tmp_path / "blacklist.bed"
    bed.write_text("chr16\t84995900\t84995950\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    records = _run_example(out_dir, padded_ref, blacklist=str(bed), thread=thread)
    expected = [line for line in (DATA_DIR / "example.vcf").read_text().splitlines() if not line.startswith("#")]
    assert any(line.split("\t")[:2] == ["chr16", "84995922"] for line in expected)
    assert not any(line.split("\t")[:2] == ["chr16", "84995922"] for line in records)
    assert "CHR2=chr17" not in "\n".join(records)
    assert 0 < len(records) < len(expected)


def test_cli_missing_blacklist_file_exits(tmp_path, padded_ref):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    with pytest.raises(SystemExit, match="blacklist file"):
        _run_example(out_dir, padded_ref, blacklist=str(tmp_path / "missing.bed"))


def test_cli_output_seq_both_and_rescue_sr(tmp_path, padded_ref):
    """Verify --output-seq both and --rescue-sr produce both FASTA variants and populated SR/OSR records."""
    out_dir = tmp_path / "out_both"
    out_dir.mkdir()
    op = DefaultOptions(
        input=str(DATA_DIR / "example.bam"),
        output=str(out_dir / "example"),
        ref=str(padded_ref),
        gtf=str(DATA_DIR / "example_ref.gtf"),
        blat_two_bit=str(DATA_DIR / "dummy.2bit"),
        output_sequence_choice="both",
        noncanonical=True,
        substitutions_fraction=0.05,
        indel_fraction=0.001,
        prune_threshold=20,
        circular_rna="remove",
        refine=True,
        rescue_sr=True,
    )
    cli.cli(op)

    assert (out_dir / "example.consensus.fasta").read_text() == (DATA_DIR / "example.fasta").read_text()
    assert (out_dir / "example.reference.fasta").read_text().startswith(">TSP")
    vcf_records = [line for line in (out_dir / "example.vcf").read_text().splitlines() if not line.startswith("#")]
    assert len(vcf_records) > 0
    assert all("SR=" in line and "OSR=" in line for line in vcf_records)
