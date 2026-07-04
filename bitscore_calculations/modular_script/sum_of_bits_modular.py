#!/usr/bin/env python3
"""
sum_of_bits_combined.py
=======================
Calculate a differential "sum of bits" score per sequence from two FIMO
best_site.narrowPeak files (stable vs unstable motifs).

    score = stable_sum_of_bits − unstable_sum_of_bits

This script is now a thin orchestrator. All reading, scoring, and filtering
logic lives in the imported modules below — mix and match filters as needed.

Available filters (import from filters/)
-----------------------------------------
  streme_filters.filter_by_evalue(xml_path, threshold)
  streme_filters.filter_by_pvalue(xml_path, threshold)
  streme_filters.filter_by_num_sites(xml_path, min_sites)
  solo_motif.filter_by_motif_name(motif_name)
  region_filter.partition_hits_by_region(narrowpeak, fasta, allowed_motifs)

Combining filters
-----------------
Because each filter returns a set of motif names, you can intersect them
with the & operator before passing to read_best_site():

    allowed = filter_by_evalue(xml, 0.05) & filter_by_num_sites(xml, 10)
    stable_raw = read_best_site(stable_path, allowed_motifs=allowed)

Usage (no filters — reproduces original behaviour)
--------------------------------------------------
    python sum_of_bits_combined.py \\
        --stable   stable_best_site.narrowPeak \\
        --unstable unstable_best_site.narrowPeak \\
        --output   sum_of_bits_combined_scores.tsv

Usage (with E-value filter)
---------------------------
    python sum_of_bits_combined.py \\
        --stable   stable_best_site.narrowPeak \\
        --unstable unstable_best_site.narrowPeak \\
        --output   sum_of_bits_combined_scores.tsv \\
        --stable-streme   stable_streme_out/streme.xml \\
        --unstable-streme unstable_streme_out/streme.xml \\
        --evalue 0.05

Usage (with E-value + minimum sites)
-------------------------------------
    python sum_of_bits_combined.py \\
        ... \\
        --stable-streme   stable_streme_out/streme.xml \\
        --unstable-streme unstable_streme_out/streme.xml \\
        --evalue 0.05 --min-sites 10

Usage (solo motif)
------------------
    python sum_of_bits_combined.py \\
        ... \\
        --motif "10-GCUGGAGCUGGGRU"
"""

import argparse
import sys
from pathlib import Path

from motif_scoring import read_best_site, compute_sum_of_bits, get_all_seq_ids
from filters.streme_filters import (
    filter_by_evalue,
    filter_by_pvalue,
    filter_by_num_sites,
)
from filters.solo_motif import filter_by_motif_name


# ---------------------------------------------------------------------------
# 1. ARGUMENT PARSING
# ---------------------------------------------------------------------------
# --stable / --unstable / --output are the original required arguments.
# All filter arguments are optional; omitting them replicates original behaviour.

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Differential sum-of-bits score per sequence from two FIMO "
            "best_site.narrowPeak files (stable vs unstable motifs)."
        ),
    )
    # Core inputs (unchanged)
    p.add_argument("-s", "--stable",   required=True,
                   help="FIMO best_site.narrowPeak — stable motifs")
    p.add_argument("-u", "--unstable", required=True,
                   help="FIMO best_site.narrowPeak — unstable motifs")
    p.add_argument("-o", "--output",   required=True,
                   help="Output TSV file (per-sequence differential scores)")

    # Optional filter arguments
    p.add_argument("--stable-streme", default=None,
                   help="Path to streme.xml from the stable STREME run "
                        "(required for --evalue, --pvalue, --min-sites)")
    p.add_argument("--unstable-streme", default=None,
                   help="Path to streme.xml from the unstable STREME run "
                        "(required for --evalue, --pvalue, --min-sites)")
    p.add_argument("--evalue",    type=float, default=None,
                   help="Retain only motifs with E-value <= this threshold")
    p.add_argument("--pvalue",    type=float, default=None,
                   help="Retain only motifs with P-value <= this threshold")
    p.add_argument("--min-sites", type=int,   default=None,
                   help="Retain only motifs found in >= this many positive sequences")
    p.add_argument("--motif",     default=None,
                   help="Score a single named motif only (e.g. '10-GCUGGAGCUGGGRU')")

    return p.parse_args()


# ---------------------------------------------------------------------------
# 2. BUILD ALLOWED-MOTIF SETS FROM WHICHEVER FILTERS WERE REQUESTED
# ---------------------------------------------------------------------------
# Because stable and unstable motifs come from separate STREME runs, filters
# are applied independently to each: stable motifs are checked against
# --stable-streme, unstable motifs against --unstable-streme.
#
# This function returns a (stable_allowed, unstable_allowed) tuple. Each
# element is either a set of motif names that passed all filters for that
# file, or None (meaning no filtering — keep all motifs from that file).
#
# The --motif (solo motif) filter is the exception: it produces a single
# one-element set that is applied identically to both files, since it
# selects by motif name rather than by statistical threshold.
#
# If multiple stat filters are active, their results are intersected (&)
# within each file so only motifs passing ALL filters for that file are kept.

def _build_allowed_for_one(xml_path, label, args):
    """
    Build the allowed-motifs set for one STREME file (stable or unstable).

    Parameters
    ----------
    xml_path : str or None   — path to this file's streme.xml
    label    : str           — "stable" or "unstable" (for printed messages)
    args     : Namespace     — parsed CLI args

    Returns None if no stat filters were requested, otherwise a set[str].
    The --motif filter is handled separately in build_allowed_motifs().
    """
    allowed = None

    if args.evalue is not None:
        print(f"  Applying E-value filter (<= {args.evalue}) to {label} motifs:")
        ev_set  = filter_by_evalue(xml_path, args.evalue)
        allowed = ev_set if allowed is None else allowed & ev_set

    if args.pvalue is not None:
        print(f"  Applying P-value filter (<= {args.pvalue}) to {label} motifs:")
        pv_set  = filter_by_pvalue(xml_path, args.pvalue)
        allowed = pv_set if allowed is None else allowed & pv_set

    if args.min_sites is not None:
        print(f"  Applying site-count filter (>= {args.min_sites}) to {label} motifs:")
        sc_set  = filter_by_num_sites(xml_path, args.min_sites)
        allowed = sc_set if allowed is None else allowed & sc_set

    return allowed


def build_allowed_motifs(args):
    """
    Construct per-file allowed-motifs sets from CLI arguments.

    Returns
    -------
    stable_allowed, unstable_allowed : set[str] or None
        None means no filtering (all motifs from that file are kept).
    """
    streme_filters_requested = any([
        args.evalue is not None,
        args.pvalue is not None,
        args.min_sites is not None,
    ])

    # Both STREME xml paths are required if any stat filter is active.
    # They are independent: missing one is always an error, even if the other
    # is provided, because each file's motifs must be checked against their
    # own STREME run.
    if streme_filters_requested:
        missing = []
        if args.stable_streme is None:
            missing.append("--stable-streme")
        if args.unstable_streme is None:
            missing.append("--unstable-streme")
        if missing:
            sys.exit(
                f"ERROR: {' and '.join(missing)} "
                f"{'is' if len(missing) == 1 else 'are'} required when using "
                "--evalue, --pvalue, or --min-sites."
            )

    # Build per-file stat-filter allowlists (None if no stat filters active)
    stable_allowed   = _build_allowed_for_one(args.stable_streme,   "stable",   args)
    unstable_allowed = _build_allowed_for_one(args.unstable_streme, "unstable", args)

    # --motif selects a single motif by name and is applied to both files
    # identically (it doesn't depend on STREME statistics).
    if args.motif is not None:
        print("Applying solo-motif filter to both files:")
        mo_set           = filter_by_motif_name(args.motif)
        stable_allowed   = mo_set if stable_allowed   is None else stable_allowed   & mo_set
        unstable_allowed = mo_set if unstable_allowed is None else unstable_allowed & mo_set

    # Summary
    if stable_allowed is not None or unstable_allowed is not None:
        s_n = len(stable_allowed)   if stable_allowed   is not None else "all"
        u_n = len(unstable_allowed) if unstable_allowed is not None else "all"
        print(f"  → {s_n} stable motif(s) and {u_n} unstable motif(s) will be scored.")
    else:
        print("No filters applied — all motifs included.")

    return stable_allowed, unstable_allowed


# ---------------------------------------------------------------------------
# 3. OUTPUT
# ---------------------------------------------------------------------------
# Writes six-column TSV: stable and unstable component scores + differential.
# Sorted descending by differential score (best sequences first).

def write_output(out_path, all_seq_ids,
                 stable_scores, stable_raw,
                 unstable_scores, unstable_raw):
    rows = []
    for seq_id in all_seq_ids:
        s_score  = stable_scores.get(seq_id, 0.0)
        u_score  = unstable_scores.get(seq_id, 0.0)
        diff     = s_score - u_score
        n_stable   = len(stable_raw.get(seq_id, {}))
        n_unstable = len(unstable_raw.get(seq_id, {}))
        rows.append((seq_id, n_stable, s_score, n_unstable, u_score, diff))

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
# 4. MAIN
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    for p in (args.stable, args.unstable):
        if not Path(p).exists():
            sys.exit(f"ERROR: Input file not found: {p}")

    # Build per-file motif allowlists (each is a set or None)
    stable_allowed, unstable_allowed = build_allowed_motifs(args)

    # Read both narrowPeak files, each filtered by its own allowlist
    print(f"Reading stable:   {args.stable}")
    stable_raw   = read_best_site(args.stable,   allowed_motifs=stable_allowed)

    print(f"Reading unstable: {args.unstable}")
    unstable_raw = read_best_site(args.unstable, allowed_motifs=unstable_allowed)

    # Sum scores and take the union of sequence IDs
    stable_scores   = compute_sum_of_bits(stable_raw)
    unstable_scores = compute_sum_of_bits(unstable_raw)
    all_seq_ids     = get_all_seq_ids(stable_raw, unstable_raw)

    # Summary
    n_stable_motifs   = len({m for md in stable_raw.values()   for m in md})
    n_unstable_motifs = len({m for md in unstable_raw.values() for m in md})
    print(f"  Stable file:   {len(stable_raw):,} sequences, "
          f"{n_stable_motifs:,} unique motifs")
    print(f"  Unstable file: {len(unstable_raw):,} sequences, "
          f"{n_unstable_motifs:,} unique motifs")
    print(f"  Total unique sequences (union): {len(all_seq_ids):,}")

    write_output(
        args.output, all_seq_ids,
        stable_scores, stable_raw,
        unstable_scores, unstable_raw,
    )
    print(f"Done. Wrote: {args.output}")


if __name__ == "__main__":
    main()