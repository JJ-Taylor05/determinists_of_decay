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
    p.add_argument("-s", "--streme", required=True,
                   help="STREME output file (used to look up nsites per motif)")
    p.add_argument("--min-sites", type=int, default=278,
                   help="Minimum nsites a motif must have to be included "
                        "in the sum-of-bits calculation (default: 278, ~20%% of sequences)")
    return p.parse_args()


def parse_streme(path):
    """
    Parse a STREME output file and return the number of sites per motif.

    Returns
    -------
    nsites : dict[motif_name -> int]
        e.g. {"1-AUAUUUAUUU": 428, "2-UGUACAUAU": 606, ...}

    How it works
    ------------
    The STREME file has a consistent two-line block for each motif:
        MOTIF 1-AUAUUUAUUU STREME-1
        letter-probability matrix: alength= 4 w= 10 nsites= 428 E= 1.2e-003
    We read line-by-line, and whenever we see a MOTIF line we store the
    motif name, then on the very next line we extract the nsites= value.
    """
    nsites = {}
    pending_motif = None  # motif name waiting for its nsites line

    with open(path) as fh:
        for line in fh:
            line = line.strip()

            if line.startswith("MOTIF "):
                # "MOTIF 1-AUAUUUAUUU STREME-1"  → take the second token
                parts = line.split()
                pending_motif = parts[1]

            elif pending_motif and line.startswith("letter-probability matrix:"):
                # "letter-probability matrix: alength= 4 w= 10 nsites= 428 E= ..."
                # Note: nsites= and the number are separate tokens (space after =)
                tokens = line.split()
                for i, token in enumerate(tokens):
                    if token == "nsites=" and i + 1 < len(tokens):
                        nsites[pending_motif] = int(tokens[i + 1])
                        break
                    # Also handle the case where nsites=428 (no space) appears as one token
                    if token.startswith("nsites=") and len(token) > 7:
                        nsites[pending_motif] = int(token.split("=")[1])
                        break
                pending_motif = None  # reset; done with this motif's header

    return nsites


def read_best_site(path, nsites_map, min_sites):
    """
    Parse the narrowPeak file, filtering out motifs below the nsites threshold.

    Parameters
    ----------
    path : Path
        FIMO best_site.narrowPeak file.
    nsites_map : dict[motif_name -> int]
        nsites for each motif, parsed from the STREME file.
    min_sites : int
        Motifs with nsites < min_sites are excluded from scoring.

    Returns
    -------
    scores : dict[seq_id -> dict[motif_name -> float]]
        Best-site bit score for each (sequence, motif) pair.
    """
    scores = defaultdict(dict)
    rejected_motifs = set()   # motifs that failed the threshold (reported once)
    unknown_motifs  = set()   # motifs not found in the STREME file at all

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

            # --- nsites threshold check ---
            if motif not in nsites_map:
                # Motif name from narrowPeak not found in STREME file at all
                if motif not in unknown_motifs:
                    print(f"  [WARNING] Motif '{motif}' not found in STREME file "
                          f"— excluding from scoring", file=sys.stderr)
                    unknown_motifs.add(motif)
                continue

            if nsites_map[motif] < min_sites:
                # Motif exists but doesn't meet the site-count threshold
                rejected_motifs.add(motif)
                continue
            # --- end threshold check ---

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

    if rejected_motifs:
        print(f"  Motifs excluded (nsites < {min_sites}): "
              f"{', '.join(sorted(rejected_motifs))}", file=sys.stderr)

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
    streme_path = Path(args.streme)

    if not input_path.exists():
        sys.exit(f"ERROR: Input file not found: {input_path}")
    if not streme_path.exists():
        sys.exit(f"ERROR: STREME file not found: {streme_path}")

    print(f"Reading STREME file: {streme_path}")
    nsites_map = parse_streme(streme_path)
    print(f"  {len(nsites_map):,} motifs found in STREME file")
    print(f"  nsites threshold: >= {args.min_sites}")

    print(f"Reading: {input_path}")
    scores = read_best_site(input_path, nsites_map, args.min_sites)
    sum_scores = compute_sum_of_bits(scores)

    n_motifs = len({m for motif_scores in scores.values() for m in motif_scores})
    print(f"  {len(scores):,} unique sequences")
    print(f"  {n_motifs:,} unique motifs (after threshold filter)")

    write_output(output_path, sum_scores, scores)
    print(f"Done. Wrote: {output_path}")


if __name__ == "__main__":
    main()
