#!/usr/bin/env python3
"""Rebuild the example reference FASTA from the windowed sequences.

ScanMST looks up reference sequence by absolute genomic coordinate, so the
FASTA needs every contig padded from position 1 even though the example only
uses a few small windows of hg38. The padded file is ~288 MB, which is too
large to keep in git, so only the real sequence windows are committed
(example_ref.windows.fasta) and the padding is regenerated here.
"""

import re
from pathlib import Path

WINDOW_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
LINE_WIDTH = 60

HERE = Path(__file__).parent
SOURCE = HERE / "example_ref.windows.fasta"
DEST = HERE / "example_ref.fasta"


def read_windows(path):
    """Read a FASTA whose headers are `chrom:start-end` regions."""
    windows, header, chunks = [], None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                windows.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        else:
            chunks.append(line)
    if header is not None:
        windows.append((header, "".join(chunks)))
    return windows


def main():
    with open(DEST, "w", encoding="utf-8") as out:
        for header, sequence in read_windows(SOURCE):
            match = WINDOW_RE.match(header)
            if match is None:
                msg = f"unexpected FASTA header {header!r}, expected 'chrom:start-end'"
                raise ValueError(msg)
            chrom, start = match["chrom"], int(match["start"])
            padded = ("N" * (start - 1)) + sequence
            out.write(f">{chrom}\n")
            for i in range(0, len(padded), LINE_WIDTH):
                out.write(padded[i : i + LINE_WIDTH] + "\n")
            print(f"{chrom}: {len(padded):,} bp ({len(sequence):,} bp of real sequence from {start:,})")
    print(f"\nWrote {DEST} ({DEST.stat().st_size / 1e6:.0f} MB)")


if __name__ == "__main__":
    main()
