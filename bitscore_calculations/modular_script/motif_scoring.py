"""
motif_scoring.py
================
Core reading and scoring functions shared across all sum-of-bits scripts.

Functions
---------
read_best_site(path, allowed_motifs=None)
    Parse a FIMO best_site.narrowPeak file into a nested dict.
    If ``allowed_motifs`` is provided (a set of motif name strings),
    only motifs whose name is in that set are retained.

compute_sum_of_bits(scores)
    Collapse a nested {seq_id: {motif: score}} dict into {seq_id: total}.

get_all_seq_ids(stable_raw, unstable_raw)
    Return the union of sequence IDs across both files.
"""

import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# FILE PARSING
# ---------------------------------------------------------------------------
# Reads one narrowPeak file and returns a nested dictionary:
#
#   { sequence_id: { motif_name: best_bit_score, ... }, ... }
#
# The optional `allowed_motifs` argument is the hook all filter modules use:
# they each produce a set of motif names that pass their threshold, and pass
# that set in here so only approved motifs contribute to the score.
#
# Column layout (0-based):
#   0  seq_id
#   1  start
#   2  end
#   3  motif_name
#   4  (unused score field)
#   5  strand
#   6  bit_score  ← this is what we accumulate

def read_best_site(path, allowed_motifs=None):
    """
    Parse a narrowPeak file.

    Parameters
    ----------
    path : str or Path
        Path to the FIMO best_site.narrowPeak file.
    allowed_motifs : set[str] or None
        If provided, only motifs whose name (col 3) is in this set are kept.
        Pass None to keep all motifs (no filtering).

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
                print(
                    f"  [WARNING] Line {lineno}: expected >=7 columns, "
                    f"got {len(cols)} — skipping",
                    file=sys.stderr,
                )
                continue

            seq_id = cols[0]
            motif  = cols[3]

            # ── Filter gate ──────────────────────────────────────────────
            # If a filter module has produced an allowed-motif set, skip any
            # motif not in it. If no filter was applied, allowed_motifs is
            # None and every motif passes through.
            if allowed_motifs is not None and motif not in allowed_motifs:
                continue

            try:
                bit_score = float(cols[6])
            except ValueError:
                print(
                    f"  [WARNING] Line {lineno}: cannot parse score "
                    f"'{cols[6]}' — skipping",
                    file=sys.stderr,
                )
                continue

            # Keep the highest score if the same (seq, motif) appears twice
            existing = scores[seq_id].get(motif)
            if existing is None or bit_score > existing:
                scores[seq_id][motif] = bit_score

    return scores


# ---------------------------------------------------------------------------
# SUM OF BITS
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
# UNION OF SEQUENCE IDs
# ---------------------------------------------------------------------------
# Sequences absent from one file simply score 0.0 for that component.

def get_all_seq_ids(stable_raw, unstable_raw):
    """Return the union of sequence IDs seen across both narrowPeak files."""
    return set(stable_raw) | set(unstable_raw)
