#!/usr/bin/env python3
"""
sum_of_bits_regional.py
=======================
Calculate per-region differential sum-of-bits scores using ORF boundaries
derived from a FASTA file to separate motif hits into 5'UTR, CDS, and 3'UTR.

Produces three output TSV files (one per region), each with the same column
structure as sum_of_bits_combined.py. Outputs are organised by region, and 
any subsequent runs of this script will add their files to the appropriate
output directory. Also produces a CSV containing the boundaries of each 
region identified in each transcript. This will be upserted on subsequent 
runs of the script. 

Can be combined with streme.xml filters so only statistically significant
motifs contribute to regional scores.

Usage (no streme filters)
--------------------------
    python sum_of_bits_regional.py \\
        --stable   stable_best_site.narrowPeak \\
        --unstable unstable_best_site.narrowPeak \\
        --fasta    transcripts.fa \\
        --output   q1 \\
        --outdir   regional_summed_scores

    # Produces:
    #   regional_summed_scores/utr5/q1_utr5.tsv
    #   regional_summed_scores/cds/q1_cds.tsv
    #   regional_summed_scores/utr3/q1_utr3.tsv
    #   regional_summed_scores/boundaries/boundaries.csv   (shared, upserted)

    # Re-running with --output q2 will produce:
    #   regional_summed_scores/utr5/q2_utr5.tsv
    #   regional_summed_scores/cds/q2_cds.tsv
    #   regional_summed_scores/utr3/q2_utr3.tsv
    #   The boundaries.csv file will be updated with any new transcripts found in q2.

Usage (with E-value filter)
---------------------------
    python sum_of_bits_regional.py \\
        --stable   stable_best_site.narrowPeak \\
        --unstable unstable_best_site.narrowPeak \\
        --fasta    transcripts.fa \\
        --output   q1 \\
        --outdir  regional_summed_scores \\
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
    compute_region_boundaries,
    write_region_boundaries_csv,
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
    p.add_argument("-o", "--outdir",   required=True,
                   help=(
                        "Label for this run, used as the filename prefix "
                        "e.g. 'q1' → q1_utr5.tsv etc. "
                        "Distinguish repeated runs of the script by using a different " 
                        "label each time (1 per quantile)"
                   ))
    p.add_argument("--outdir", default=".",
                     help=(
                        "Parent directory for regional subdirectories (utr5/, cds/, utr3/) " 
                        "and boundaries/), created if it doesn't exist and then reused on "
                        "subsequent runs. Default: current working directory."
                     ))
    p.add_argument("--streme-xml", default=None,
                   help="Path to streme.xml (required for --evalue, --pvalue, --min-sites)")
    p.add_argument("--evalue",    type=float, default=None)
    p.add_argument("--pvalue",    type=float, default=None)
    p.add_argument("--min-sites", type=int,   default=None)
    p.add_argument("--motif",     default=None,
                   help="Score a single named motif only")

    return p.parse_args()


# ---------------------------------------------------------------------------
# BUILD MOTIF ALLOWLIST 
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
# OUTPUT 
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

    # Create regional output directories (plus boundaries/) if it doesn't exist.
    # These are created on initial run of this script (parents=True, exist_ok=True)  
    # and are reused on subsequent runs. The boundaries/ directory is shared across 
    # all runs of the script.
    out_dir = Path(args.outdir)
    subdirs = {
        name: outdir / name
        for name in ("utr5", "cds", "utr3", "boundaries")
    }
    for name, d in subdirs.items():
        d.mkdir(parents=True, exist_ok=True)
    print(f"Output directories (created if missing, reused if present):")
    for name, d in subdirs.items():
        print(f"  {d}/")


    allowed_motifs = build_allowed_motifs(args)

    # Compute ORF-derived region boundaries once from the FASTA, then reuse
    # for both the stable and unstable partitioning (and for the boundaries CSV).
    print(f"\nComputing region boundaries from FASTA:")
    region_map = compute_region_boundaries(args.fasta)

    boundaries_csv_path = subdirs["boundaries"] / "boundaries.csv"
    write_region_boundaries_csv(boundaries_csv_path, region_map)

    # Partition hits from BOTH files into regions.
    # partition_hits_by_region returns sets of (seq_id, motif) tuples —
    # we call it separately for stable and unstable so each gets its own
    # region assignment (the same hit appears in at most one region per file).

    print(f"\nPartitioning stable hits by region:")
    s_utr5, s_cds, s_utr3 = partition_hits_by_region(
        args.stable, args.fasta, allowed_motifs=allowed_motifs, region_map=region_map,
    )

    print(f"\nPartitioning unstable hits by region:")
    u_utr5, u_cds, u_utr3 = partition_hits_by_region(
        args.unstable, args.fasta, allowed_motifs=allowed_motifs, region_map=region_map,
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

        out_path = subdirs[suffix] / f"{args.output}_{suffix}.tsv"
        write_output(out_path, label, all_ids,
                     s_scores, s_raw, u_scores, u_raw)

    print("\nDone.")


if __name__ == "__main__":
    main()
