#!/usr/bin/env python3
"""
Analyze the amplified dinucleotide-shuffle negative control for discriminative
STREME motif discovery.

Null design: shuffled foreground vs REAL (unshuffled) background, matching
the real search's direction each time (e.g. shuffled-stable vs real-unstable
for the stable-direction null). Only the foreground class's dinucleotide
composition is preserved from the real data each iteration; the background
class is exactly the real training data STREME sees in the actual search.

For each direction (stable, unstable), produces:
  1. A "global null" summary across N null iterations: how many motifs
     passed STREME's E-value threshold and what the best E-value was in
     each null iteration, with your real run's values placed against that
     empirical distribution.
  2. A per-real-motif empirical p-value: the fraction of the N null
     iterations in which Tomtom found a significant match to that motif
     among that iteration's null motifs, then BH-FDR corrected across the
     direction's motif set.

Run this after all 01_null_array.sbatch array tasks have finished (handled
automatically by submit_all.sh via --dependency=afterany).

Usage:
    python 02_analyze_null.py --project-dir /path/to/streme_null_pipeline \
                               --data-dir /path/to/data --n-iter 100
"""
import argparse
import os
import re
from collections import defaultdict

import pandas as pd
from statsmodels.stats.multitest import multipletests

MOTIF_HEADER_RE = re.compile(r"^MOTIF\s+(\S+)")
EVALUE_RE = re.compile(r"E=\s*([\d.eE+-]+)")
TOMTOM_EVALUE_THRESHOLD = 0.5


def parse_streme_motifs(streme_txt_path):
    """Return list of (motif_id, evalue) from a streme.txt file."""
    motifs = []
    if not os.path.exists(streme_txt_path):
        return motifs
    with open(streme_txt_path) as fh:
        current_id = None
        for line in fh:
            m = MOTIF_HEADER_RE.match(line)
            if m:
                current_id = m.group(1)
                continue
            if line.startswith("letter-probability matrix"):
                em = EVALUE_RE.search(line)
                evalue = float(em.group(1)) if em else None
                motifs.append((current_id, evalue))
    return motifs


def summarize_real_run(streme_txt_path, label):
    motifs = parse_streme_motifs(streme_txt_path)
    print(f"\n=== Real {label} run: {streme_txt_path} ===")
    print(f"  {len(motifs)} motifs reported")
    for mid, ev in motifs:
        print(f"    {mid}\tE={ev:.2e}")
    return motifs


def build_global_null(null_dir, direction, n_iter):
    """Returns DataFrame: iteration, n_motifs, best_evalue (NaN if missing)."""
    rows = []
    for i in range(1, n_iter + 1):
        streme_path = os.path.join(
            null_dir, f"iter_{i}", f"streme_{direction}_null", "streme.txt"
        )
        motifs = parse_streme_motifs(streme_path)
        best_e = min((e for _, e in motifs if e is not None), default=float("nan"))
        rows.append({"iteration": i, "n_motifs": len(motifs), "best_evalue": best_e})
    df = pd.DataFrame(rows)
    n_missing = df["n_motifs"].isna().sum() if df.empty else (df["best_evalue"].isna() & (df["n_motifs"] == 0)).sum()
    if n_missing:
        print(f"  NOTE: {n_missing}/{n_iter} {direction} null iterations had no streme.txt "
              f"(job failed, timed out, or STREME found zero motifs - check logs/)")
    return df


def empirical_percentile(real_value, null_values, higher_is_more_extreme):
    """Fraction of null values at least as extreme as the real value."""
    null_values = [v for v in null_values if pd.notna(v)]
    if not null_values or pd.isna(real_value):
        return float("nan")
    if higher_is_more_extreme:
        n_extreme = sum(1 for v in null_values if v >= real_value)
    else:
        n_extreme = sum(1 for v in null_values if v <= real_value)
    return n_extreme / len(null_values)


def get_significant_query_hits(tomtom_tsv_path):
    """
    Return the set of real (query) motif IDs that had >=1 significant hit
    in this single iteration's Tomtom run.

    Tomtom itself was called with -evalue -thresh 0.5, so tomtom.tsv already
    only contains rows with E-value <= 0.5 - this filter re-applies the same
    threshold on the "E-value" column explicitly (rather than relying only
    on what Tomtom already wrote) so the criterion stays correct even if the
    Tomtom invocation's threshold is changed later.
    """
    if not os.path.exists(tomtom_tsv_path):
        return set()
    df = pd.read_csv(tomtom_tsv_path, sep="\t", comment="#")
    if df.empty or "Query_ID" not in df.columns:
        return set()
    df = df.dropna(subset=["Query_ID"])
    sig = df[df["E-value"] < TOMTOM_EVALUE_THRESHOLD]
    return set(sig["Query_ID"].unique())


def compute_recurrence(project_dir, direction, real_motif_ids, n_iter):
    """
    For each real motif ID, count how many of the n_iter null iterations
    produced a significant Tomtom match to it (exact per-iteration mapping,
    since Tomtom was run separately per iteration against that iteration's
    null motifs only).
    """
    counts = defaultdict(int)
    n_iters_with_tomtom_output = 0
    for i in range(1, n_iter + 1):
        tsv_path = os.path.join(
            project_dir, "null_runs", f"iter_{i}", f"tomtom_{direction}", "tomtom.tsv"
        )
        if os.path.exists(tsv_path):
            n_iters_with_tomtom_output += 1
        hits = get_significant_query_hits(tsv_path)
        for mid in hits:
            counts[mid] += 1
    print(f"  {direction}: Tomtom output found for {n_iters_with_tomtom_output}/{n_iter} iterations")
    return {mid: counts.get(mid, 0) for mid in real_motif_ids}, n_iters_with_tomtom_output


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-dir", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--n-iter", type=int, default=100)
    args = ap.parse_args()

    null_dir = os.path.join(args.project_dir, "null_runs")
    out_dir = os.path.join(args.project_dir, "results")
    os.makedirs(out_dir, exist_ok=True)

    all_results = []

    for direction in ["stable", "unstable"]:
        real_path = os.path.join(args.data_dir, f"streme_{direction}.txt")
        real_motifs = summarize_real_run(real_path, direction)
        real_motif_ids = [mid for mid, _ in real_motifs]

        # --- Global null: motif counts and best E-value per iteration ---
        null_df = build_global_null(null_dir, direction, args.n_iter)
        null_df.to_csv(os.path.join(out_dir, f"global_null_{direction}.csv"), index=False)

        real_n_motifs = len(real_motifs)
        real_best_e = min((e for _, e in real_motifs if e is not None), default=float("nan"))

        pct_count = empirical_percentile(real_n_motifs, null_df["n_motifs"], higher_is_more_extreme=True)
        pct_evalue = empirical_percentile(real_best_e, null_df["best_evalue"], higher_is_more_extreme=False)

        print(f"\n--- Global null summary: {direction} ---")
        print(f"  Real run: {real_n_motifs} motifs, best E={real_best_e:.2e}")
        print(f"  Null (n={len(null_df)}): mean motifs={null_df['n_motifs'].mean():.2f} "
              f"(sd={null_df['n_motifs'].std():.2f}), mean best-E={null_df['best_evalue'].mean():.2e}")
        print(f"  Empirical p (motif count as extreme as real): {pct_count:.4f}")
        print(f"  Empirical p (best E-value as extreme as real): {pct_evalue:.4f}")

        # --- Per-motif recurrence via per-iteration Tomtom ---
        recurrence, n_valid_iters = compute_recurrence(
            args.project_dir, direction, real_motif_ids, args.n_iter
        )

        for mid, evalue in real_motifs:
            n_hits = recurrence[mid]
            denom = n_valid_iters if n_valid_iters > 0 else args.n_iter
            emp_p = n_hits / denom
            all_results.append({
                "direction": direction,
                "motif_id": mid,
                "real_evalue": evalue,
                "n_null_iterations_with_match": n_hits,
                "n_iterations_evaluated": denom,
                "empirical_p_recurrence": emp_p,
            })

    results_df = pd.DataFrame(all_results)

    # BH correction within each direction's motif family separately
    results_df["bh_qvalue"] = float("nan")
    for direction in results_df["direction"].unique():
        mask = results_df["direction"] == direction
        pvals = results_df.loc[mask, "empirical_p_recurrence"].values
        _, qvals, _, _ = multipletests(pvals, method="fdr_bh")
        results_df.loc[mask, "bh_qvalue"] = qvals

    results_df["robust_call"] = (
        (results_df["bh_qvalue"] < 0.05) & (results_df["empirical_p_recurrence"] < 0.05)
    ).map({True: "PASS (not explained by shuffled null)", False: "FAIL (recurs in null / not significant after FDR)"})

    out_csv = os.path.join(out_dir, "motif_robustness_summary.csv")
    results_df.to_csv(out_csv, index=False)

    print(f"\n=== Final per-motif table written to {out_csv} ===")
    print(results_df.to_string(index=False))
    print(
        "\nNOTE: empirical p-value resolution is limited to 1/n_iter (e.g. 1/100 = 0.01 "
        "minimum non-zero value) - a motif with 0 null matches gets p=0, which BH treats "
        "as maximally significant but really just means 'not once in n_iter tries'. "
        "Consider this when n_iter is small.\n"
        "Combine this table with your SEA test-set enrichment result: a motif should "
        "ideally be SEA-significant in the held-out test set AND get a PASS here to be "
        "reported as a robust discriminative motif."
    )


if __name__ == "__main__":
    main()
