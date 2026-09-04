#!/usr/bin/env python3
"""
map_rnaprot_peaks_to_metatranscript.py

Joins RNAProt `predict --mode 2` peak calls (peak_regions.tsv) with
transcript UTR5/CDS/UTR3 boundary annotations (boundaries.csv) to produce
a "hits" TSV in the exact format expected by
metatranscript_stable_vs_unstable.R's prepare_hits() function:

    start   end   seq_len   orf_start   orf_end

Usage:
    python map_rnaprot_peaks_to_metatranscript.py \
        --peaks stable_train_predict_out/peak_regions.tsv \
        --boundaries boundaries.csv \
        --out plot_data/stable_motif_mapping_hits.tsv

Run once for the stable set and once for the unstable set, pointing
--out at the two paths your R script already reads:
    plot_data/stable_motif_mapping_hits.tsv
    plot_data/unstable_motif_mapping_hits.tsv
"""
import argparse
import os
import sys

import pandas as pd


def find_column(columns, exact_names, contains_fallback):
    """Return the first matching column name: try exact matches first,
    then fall back to a substring search."""
    for name in exact_names:
        if name in columns:
            return name
    for c in columns:
        if contains_fallback in c.lower():
            return c
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--peaks", required=True, help="RNAProt peak_regions.tsv file")
    ap.add_argument(
        "--boundaries", required=True, help="boundaries.csv with UTR5/CDS/UTR3 coords"
    )
    ap.add_argument(
        "--out",
        required=True,
        help="Output hits TSV (start, end, seq_len, orf_start, orf_end)",
    )
    ap.add_argument(
        "--peak-id-col",
        default=None,
        help="Column name in the peaks file holding the transcript/sequence ID "
        "(auto-detected if not given)",
    )
    ap.add_argument(
        "--start-col",
        default=None,
        help="Column name in the peaks file holding the peak start position "
        "(auto-detected if not given; prefers 'peak_region_s' over 'window_s')",
    )
    ap.add_argument(
        "--end-col",
        default=None,
        help="Column name in the peaks file holding the peak end position "
        "(auto-detected if not given; prefers 'peak_region_e' over 'window_e')",
    )
    args = ap.parse_args()

    peaks = pd.read_csv(args.peaks, sep="\t")
    bounds = pd.read_csv(args.boundaries)

    if "sequence_id" not in bounds.columns:
        sys.exit(
            f"Expected a 'sequence_id' column in {args.boundaries}, "
            f"found: {list(bounds.columns)}"
        )

    # --- Identify the transcript ID column in the peaks file ---
    id_col = args.peak_id_col
    if id_col is None:
        id_col = find_column(
            peaks.columns,
            exact_names=["reference_id", "ref_id", "seq_id", "site_id", "chr", "id"],
            contains_fallback="id",
        )
        if id_col is None:
            sys.exit(
                f"Could not auto-detect the transcript ID column in {args.peaks}.\n"
                f"Found columns: {list(peaks.columns)}\n"
                f"Re-run with --peak-id-col <name> to specify it explicitly."
            )
        print(f"Using '{id_col}' as the transcript ID column in the peaks file.")

    # --- Identify start/end columns ---
    # RNAProt's peak_regions.tsv uses 'peak_region_s'/'peak_region_e' (the
    # called peak itself) and 'window_s'/'window_e' (the wider scan window).
    # Prefer the peak region bounds for metatranscript mapping, since that's
    # the actual called structural element, not the surrounding window.
    start_col = args.start_col or find_column(
        peaks.columns,
        exact_names=["peak_region_s", "start", "window_s"],
        contains_fallback="start",
    )
    end_col = args.end_col or find_column(
        peaks.columns,
        exact_names=["peak_region_e", "end", "window_e"],
        contains_fallback="end",
    )
    if start_col is None or end_col is None:
        sys.exit(
            f"Could not find start/end columns in {args.peaks}.\n"
            f"Found columns: {list(peaks.columns)}\n"
            f"Re-run with --start-col/--end-col to specify explicitly."
        )
    print(f"Using '{start_col}'/'{end_col}' as peak start/end columns.")

    hits = peaks[[id_col, start_col, end_col]].rename(
        columns={id_col: "sequence_id", start_col: "start", end_col: "end"}
    )

    # --- Merge with boundaries ---
    merged = hits.merge(bounds, on="sequence_id", how="left", validate="many_to_one")

    missing = merged["seq_len"].isna().sum()
    if missing:
        print(
            f"Warning: {missing} peak(s) had no matching entry in "
            f"{args.boundaries} and will be dropped.",
            file=sys.stderr,
        )
        merged = merged.dropna(subset=["seq_len"])

    # Transcripts with no called ORF (has_orf == False) get NA orf_start/end,
    # which triggers the R script's is.na(orf_start) fallback (hit_centre / seq_len).
    has_orf = merged["has_orf"].astype(str).str.lower().isin(["true", "1", "yes"])
    merged["orf_start"] = merged["cds_start"].where(has_orf)
    merged["orf_end"] = merged["cds_end"].where(has_orf)

    out = merged[["start", "end", "seq_len", "orf_start", "orf_end"]]

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(out)} hits to {args.out}")


if __name__ == "__main__":
    main()
