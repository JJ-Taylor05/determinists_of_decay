#!/usr/bin/env python3
"""
sum_of_bits.py
==============
Calculate per-sequence "sum of bits" scores from a FIMO best_site.narrowPeak file.

Method
------
Implements the "best match score per motif" approach described in:
    Stormo GD & Fields DS (1998). Specificity, free energy and information
    content in protein-DNA interactions. Trends Biochem Sci 23(3):109-13.
    DOI: 10.1089/cmb.1998.5.211

For each sequence:
  1. The best_site file already contains ONE hit per (sequence, motif) pair —
     the highest-scoring match in that sequence for that motif.
  2. The per-motif bit score (column 7) is summed across all motifs that
     produced a hit in that sequence.

The output is a tab-separated file with one row per sequence, sorted
descending by sum-of-bits score.

Input format (narrowPeak, 10 columns, tab-separated)
-----------------------------------------------------
Col 1  sequence ID
Col 4  motif name
Col 7  bit score (float — THIS is the score we sum)

Usage
-----
    python sum_of_bits.py --input best_site.narrowPeak --output sum_of_bits_scores.tsv
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Sum-of-bits score per sequence from FIMO best_site.narrowPeak",
    )
    p.add_argument("-i", "--input", required=True,
                   help="FIMO best_site.narrowPeak file")
    p.add_argument("-o", "--output", required=True,
                   help="Output TSV file (per-sequence sum-of-bits scores)")
    return p.parse_args()


def read_best_site(path):
    """
    Parse the narrowPeak file.

    Returns
    -------
    scores : dict[seq_id -> dict[motif_name -> float]]
        Best-site bit score for each (sequence, motif) pair.
    """
    scores = defaultdict(dict)

    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")

            if len(cols) < 7:
                print(f"  [WARNING] Line {lineno}: expected >=7 columns, "
                      f"got {len(cols)} — skipping", file=sys.stderr)
                continue

            seq_id = cols[0]
            motif  = cols[3]

            try:
                bit_score = float(cols[6])
            except ValueError:
                print(f"  [WARNING] Line {lineno}: cannot parse score '{cols[6]}' "
                      f"— skipping", file=sys.stderr)
                continue

            # The best_site file should already have one hit per (seq, motif),
            # but if duplicates exist, keep the maximum score.
            existing = scores[seq_id].get(motif)
            if existing is None or bit_score > existing:
                scores[seq_id][motif] = bit_score

    return scores


def compute_sum_of_bits(scores):
    """Sum bit scores across all motifs for each sequence."""
    return {
        seq_id: sum(motif_scores.values())
        for seq_id, motif_scores in scores.items()
    }


def write_output(out_path, sum_scores, scores):
    """Write per-sequence TSV: seq_id, n_motifs_hit, sum_of_bits_score."""
    rows = sorted(
        [(seq_id, len(scores[seq_id]), s) for seq_id, s in sum_scores.items()],
        key=lambda r: r[2],
        reverse=True,
    )
    with open(out_path, "w") as fh:
        fh.write("sequence_id\tn_motifs_hit\tsum_of_bits_score\n")
        for seq_id, n_motifs, s in rows:
            fh.write(f"{seq_id}\t{n_motifs}\t{s:.6f}\n")


def main():
    args = parse_args()
    input_path  = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        sys.exit(f"ERROR: Input file not found: {input_path}")

    print(f"Reading: {input_path}")
    scores = read_best_site(input_path)
    sum_scores = compute_sum_of_bits(scores)

    n_motifs = len({m for motif_scores in scores.values() for m in motif_scores})
    print(f"  {len(scores):,} unique sequences")
    print(f"  {n_motifs:,} unique motifs")

    write_output(output_path, sum_scores, scores)
    print(f"Done. Wrote: {output_path}")


if __name__ == "__main__":
    main()
