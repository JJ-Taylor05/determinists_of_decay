#!/usr/bin/env python3
"""
Calculate a differential "sum of bits" score per sequence from two FIMO
best_site.narrowPeak files: one containing stable motif hits and one
containing unstable motif hits.

Scoring logic: sum_of_bits_combined_score = stable_sum_of_bits − unstable_sum_of_bits

Sequences are rewarded for containing stable motifs and penalised for
containing unstable motifs.  Sequences present in only one file receive
a contribution of 0.0 from the absent file.

The output is a tab-separated file with one row per sequence, sorted
descending by score.

Usage
    python sum_of_bits_differential.py \\
        --stable   stable_best_site.narrowPeak \\
        --unstable unstable_best_site.narrowPeak \\
        --output   sum_of_bits_combined_scores.tsv
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# 1. ARGUMENT PARSING
# ---------------------------------------------------------------------------
# Accept two input files instead of one.
# --stable   → hits for motifs associated with stable behaviour  (rewarded)
# --unstable → hits for motifs associated with unstable behaviour (penalised)

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Differential sum-of-bits score per sequence from two FIMO "
            "best_site.narrowPeak files (stable vs unstable motifs)."
        ),
    )
    p.add_argument("-s", "--stable", required=True,
                   help="FIMO best_site.narrowPeak file — stable motifs")
    p.add_argument("-u", "--unstable", required=True,
                   help="FIMO best_site.narrowPeak file — unstable motifs")
    p.add_argument("-o", "--output", required=True,
                   help="Output TSV file (per-sequence differential scores)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# 2. FILE PARSING  (unchanged from original)
# ---------------------------------------------------------------------------
# This function reads one narrowPeak file and returns a nested dictionary:
#
#   { sequence_id: { motif_name: best_bit_score, ... }, ... }
#
# If a (sequence, motif) pair appears more than once we keep the maximum
# score — the best_site file should already guarantee uniqueness, but this
# guards against edge cases.

def read_best_site(path):
    """
    Parse a narrowPeak file.

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

            existing = scores[seq_id].get(motif)
            if existing is None or bit_score > existing:
                scores[seq_id][motif] = bit_score

    return scores


# ---------------------------------------------------------------------------
# 3. SUM OF BITS  (unchanged from original)
# ---------------------------------------------------------------------------
# Collapses the nested dict into a flat {seq_id: total_score} dict by
# summing across all motifs that hit that sequence.

def compute_sum_of_bits(scores):
    """Sum bit scores across all motifs for each sequence."""
    return {
        seq_id: sum(motif_scores.values())
        for seq_id, motif_scores in scores.items()
    }


# ---------------------------------------------------------------------------
# 4. OUTPUT  (extended to include component scores)
# ---------------------------------------------------------------------------
# We write four columns so you can see exactly where each sequence's
# differential score comes from:
#
#   sequence_id | n_stable_motifs | stable_score | n_unstable_motifs |
#   unstable_score | sum_of_bits_combined_score
#
# Rows are sorted descending by sum_of_bits_combined_score (best sequences first).

def write_output(out_path, all_seq_ids,
                 stable_scores, stable_raw,
                 unstable_scores, unstable_raw):
    """
    Write per-sequence TSV with stable, unstable, and differential scores.

    Parameters
    ----------
    all_seq_ids     : set of every sequence ID seen across both files
    stable_scores   : dict[seq_id -> float]  — summed stable bits
    stable_raw      : dict[seq_id -> dict]   — raw stable hits (for motif count)
    unstable_scores : dict[seq_id -> float]  — summed unstable bits
    unstable_raw    : dict[seq_id -> dict]   — raw unstable hits (for motif count)
    """
    rows = []
    for seq_id in all_seq_ids:
        s_score  = stable_scores.get(seq_id, 0.0)
        u_score  = unstable_scores.get(seq_id, 0.0)
        diff     = s_score - u_score
        n_stable   = len(stable_raw.get(seq_id, {}))
        n_unstable = len(unstable_raw.get(seq_id, {}))
        rows.append((seq_id, n_stable, s_score, n_unstable, u_score, diff))

    # Sort descending by differential score — highest-scoring sequences first
    rows.sort(key=lambda r: r[5], reverse=True)

    with open(out_path, "w") as fh:
        fh.write(
            "sequence_id\t"
            "n_stable_motifs_hit\tstable_sum_of_bits\t"
            "n_unstable_motifs_hit\tunstable_sum_of_bits\t"
            "sum_of_bits_combined_score\n"
        )
        for seq_id, n_s, s_sc, n_u, u_sc, diff in rows:
            fh.write(
                f"{seq_id}\t{n_s}\t{s_sc:.6f}\t"
                f"{n_u}\t{u_sc:.6f}\t{diff:.6f}\n"
            )


# ---------------------------------------------------------------------------
# 5. MAIN  (orchestrates everything)
# ---------------------------------------------------------------------------
# Reads both files → computes sums → takes the union of all sequence IDs →
# subtracts unstable from stable → writes output.

def main():
    args = parse_args()
    stable_path   = Path(args.stable)
    unstable_path = Path(args.unstable)
    output_path   = Path(args.output)

    for p in (stable_path, unstable_path):
        if not p.exists():
            sys.exit(f"ERROR: Input file not found: {p}")

    # --- Read both files ---
    print(f"Reading stable:   {stable_path}")
    stable_raw = read_best_site(stable_path)

    print(f"Reading unstable: {unstable_path}")
    unstable_raw = read_best_site(unstable_path)

    # --- Compute per-file sums ---
    stable_scores   = compute_sum_of_bits(stable_raw)
    unstable_scores = compute_sum_of_bits(unstable_raw)

    # --- Union of all sequence IDs across both files ---
    # Sequences absent from one file simply score 0.0 for that component.
    all_seq_ids = set(stable_scores) | set(unstable_scores)

    # --- Summary stats ---
    n_stable_motifs   = len({m for md in stable_raw.values()   for m in md})
    n_unstable_motifs = len({m for md in unstable_raw.values() for m in md})
    print(f"  Stable file:   {len(stable_raw):,} sequences, "
          f"{n_stable_motifs:,} unique motifs")
    print(f"  Unstable file: {len(unstable_raw):,} sequences, "
          f"{n_unstable_motifs:,} unique motifs")
    print(f"  Total unique sequences (union): {len(all_seq_ids):,}")

    # --- Write output ---
    write_output(output_path, all_seq_ids,
                 stable_scores, stable_raw,
                 unstable_scores, unstable_raw)
    print(f"Done. Wrote: {output_path}")


if __name__ == "__main__":
    main()
