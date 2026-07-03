#!/usr/bin/env python3
"""
sum_of_bits_regional.py
=======================
Calculate per-region differential sum-of-bits scores using ORF boundaries
derived from a FASTA file to separate motif hits into 5'UTR, CDS, and 3'UTR.

Produces three output TSV files (one per region), each with the same column
structure as sum_of_bits_combined.py.

Can be combined with streme.xml filters so only statistically significant
motifs contribute to regional scores.

Usage (no streme filters)
--------------------------
    python sum_of_bits_regional.py \\
        --stable   stable_best_site.narrowPeak \\
        --unstable unstable_best_site.narrowPeak \\
        --fasta    transcripts.fa \\
        --output   regional_scores

    # Produces:
    #   regional_scores_utr5.tsv
    #   regional_scores_cds.tsv
    #   regional_scores_utr3.tsv

Usage (with E-value filter)
---------------------------
    python sum_of_bits_regional.py \\
        --stable   stable_best_site.narrowPeak \\
        --unstable unstable_best_site.narrowPeak \\
        --fasta    transcripts.fa \\
        --output   regional_scores \\
        --streme-xml streme_out/streme.xml \\
        --evalue 0.05
"""

import argparse
import sys
from pathlib import Path

from motif_scoring import compute_sum_of_bits
from filters.streme_filters import (
    filter_by_evalue,
    filter_by_pvalue,
    filter_by_num_sites,
)
from filters.solo_motif import filter_by_motif_name
from filters.region_filter import (
    partition_hits_by_region,
    read_best_site_for_region,
)


# ---------------------------------------------------------------------------
# ARGUMENT PARSING
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Per-region differential sum-of-bits from two narrowPeak files. "
            "Produces three output TSVs: 5'UTR, CDS, and 3'UTR."
        ),
    )
    p.add_argument("-s", "--stable",   required=True,
                   help="FIMO best_site.narrowPeak — stable motifs")
    p.add_argument("-u", "--unstable", required=True,
                   help="FIMO best_site.narrowPeak — unstable motifs")
    p.add_argument("-f", "--fasta",    required=True,
                   help="FASTA file of transcript sequences (for ORF finding)")
    p.add_argument("-o", "--output",   required=True,
                   help="Output prefix (e.g. 'results' → results/results_utr5.tsv etc.)")

    p.add_argument("--streme-xml", default=None,
                   help="Path to streme.xml (required for --evalue, --pvalue, --min-sites)")
    p.add_argument("--evalue",    type=float, default=None)
    p.add_argument("--pvalue",    type=float, default=None)
    p.add_argument("--min-sites", type=int,   default=None)
    p.add_argument("--motif",     default=None,
                   help="Score a single named motif only")

    return p.parse_args()


# ---------------------------------------------------------------------------
# BUILD MOTIF ALLOWLIST  (same logic as sum_of_bits_combined.py)
# ---------------------------------------------------------------------------

def build_allowed_motifs(args):
    streme_requested = any([
        args.evalue is not None,
        args.pvalue is not None,
        args.min_sites is not None,
    ])
    if streme_requested and args.streme_xml is None:
        sys.exit("ERROR: --streme-xml is required when using --evalue, --pvalue, or --min-sites.")

    allowed = None
    if args.evalue is not None:
        print(f"Applying E-value filter (<= {args.evalue}):")
        s = filter_by_evalue(args.streme_xml, args.evalue)
        allowed = s if allowed is None else allowed & s
    if args.pvalue is not None:
        print(f"Applying P-value filter (<= {args.pvalue}):")
        s = filter_by_pvalue(args.streme_xml, args.pvalue)
        allowed = s if allowed is None else allowed & s
    if args.min_sites is not None:
        print(f"Applying site-count filter (>= {args.min_sites}):")
        s = filter_by_num_sites(args.streme_xml, args.min_sites)
        allowed = s if allowed is None else allowed & s
    if args.motif is not None:
        print("Applying solo-motif filter:")
        s = filter_by_motif_name(args.motif)
        allowed = s if allowed is None else allowed & s
    if allowed is not None:
        print(f"  → {len(allowed)} motif(s) will be used for scoring.")
    else:
        print("No filters applied — all motifs included.")
    return allowed


# ---------------------------------------------------------------------------
# OUTPUT  (same format as sum_of_bits_combined.py, labelled by region)
# ---------------------------------------------------------------------------

def write_output(out_path, region_label, all_seq_ids,
                 stable_scores, stable_raw,
                 unstable_scores, unstable_raw):
    rows = []
    for seq_id in all_seq_ids:
        s_score    = stable_scores.get(seq_id, 0.0)
        u_score    = unstable_scores.get(seq_id, 0.0)
        diff       = s_score - u_score
        n_stable   = len(stable_raw.get(seq_id, {}))
        n_unstable = len(unstable_raw.get(seq_id, {}))
        rows.append((seq_id, n_stable, s_score, n_unstable, u_score, diff))

    rows.sort(key=lambda r: r[5], reverse=True)

    with open(out_path, "w") as fh:
        fh.write(
            f"# Region: {region_label}\n"
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
    print(f"  Wrote {len(rows):,} rows → {out_path}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    for p in (args.stable, args.unstable, args.fasta):
        if not Path(p).exists():
            sys.exit(f"ERROR: Input file not found: {p}")

    # Create the output directory (named after the prefix) if it doesn't exist.
    # parents=True allows nested paths like "results/run1"; exist_ok=True means
    # re-running the script won't raise an error if the directory is already there.
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {out_dir}/")

    allowed_motifs = build_allowed_motifs(args)

    # Partition hits from BOTH files into regions.
    # partition_hits_by_region returns sets of (seq_id, motif) tuples —
    # we call it separately for stable and unstable so each gets its own
    # region assignment (the same hit appears in at most one region per file).

    print(f"\nPartitioning stable hits by region:")
    s_utr5, s_cds, s_utr3 = partition_hits_by_region(
        args.stable, args.fasta, allowed_motifs=allowed_motifs
    )

    print(f"\nPartitioning unstable hits by region:")
    u_utr5, u_cds, u_utr3 = partition_hits_by_region(
        args.unstable, args.fasta, allowed_motifs=allowed_motifs
    )

    # For each region: read scores for hits assigned to that region,
    # compute sums, take union of seq IDs, write output.
    regions = [
        ("5'UTR", "utr5", s_utr5, u_utr5),
        ("CDS",   "cds",  s_cds,  u_cds),
        ("3'UTR", "utr3", s_utr3, u_utr3),
    ]

    print()
    for label, suffix, s_hits, u_hits in regions:
        print(f"Scoring {label}:")
        s_raw    = read_best_site_for_region(args.stable,   s_hits)
        u_raw    = read_best_site_for_region(args.unstable, u_hits)
        s_scores = compute_sum_of_bits(s_raw)
        u_scores = compute_sum_of_bits(u_raw)
        all_ids  = set(s_scores) | set(u_scores)

        out_path = out_dir / f"{out_dir.name}_{suffix}.tsv"
        write_output(out_path, label, all_ids,
                     s_scores, s_raw, u_scores, u_raw)

    print("\nDone.")


if __name__ == "__main__":
    main()
