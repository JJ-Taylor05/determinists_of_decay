#!/usr/bin/env python3
"""
merge_scores.py
===============
Merges FIMO sum-of-bits scores with mRNA degradation data by mapping
transcript IDs (ENST) to gene IDs (ENSG) via the Ensembl REST API.

Steps
-----
1. Load transcript IDs from the sum-of-bits TSV.
2. POST them to the Ensembl REST API in batches of 1000 to retrieve
   the parent gene ID for each transcript.
3. Load the half-life CSV (Gene ID, Gene Name, Degradation value).
4. Join on Gene ID and write the merged output.

Output columns
--------------
gene_id | gene_name | transcript_id | degradation_value | sum_of_bits_score

Usage
-----
    python merge_scores.py \\
        --bits   sum_of_bits_scores.tsv \\
        --halflife human_halflife_data_sorted.csv \\
        --output merged_scores.tsv
"""

import argparse
import csv
import json
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id"
BATCH_SIZE = 1000
RETRY_LIMIT = 3
RETRY_BACKOFF = 5  # seconds


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Merge sum-of-bits scores with degradation data via Ensembl REST API"
    )
    p.add_argument("--bits", required=True,
                   help="Sum-of-bits TSV (sequence_id, n_motifs_hit, sum_of_bits_score)")
    p.add_argument("--halflife", required=True,
                   help="Half-life CSV (no header: gene_id, gene_name, degradation_value)")
    p.add_argument("--output", required=True,
                   help="Output TSV file")
    return p.parse_args()


# ---------------------------------------------------------------------------
# File loading
# ---------------------------------------------------------------------------

def load_bits(path):
    """
    Returns dict: transcript_id (with version) -> sum_of_bits_score
    e.g. {'ENST00000252999.7': 545.931, ...}
    """
    bits = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            bits[row["sequence_id"]] = float(row["sum_of_bits_score"])
    print(f"  Loaded {len(bits):,} transcript scores from {path}")
    return bits


def load_halflife(path):
    """
    Returns dict: gene_id (no version) -> (gene_name, degradation_value)
    e.g. {'ENSG00000110848': ('CD69', -17.18)}
    """
    halflife = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        for row in reader:
            if len(row) < 3:
                continue
            gene_id, gene_name, deg_val = row[0].strip(), row[1].strip(), row[2].strip()
            try:
                halflife[gene_id] = (gene_name, float(deg_val))
            except ValueError:
                print(f"  [WARNING] Could not parse degradation value '{deg_val}' "
                      f"for {gene_id} — skipping", file=sys.stderr)
    print(f"  Loaded {len(halflife):,} gene entries from {path}")
    return halflife


# ---------------------------------------------------------------------------
# Ensembl REST API
# ---------------------------------------------------------------------------

def ensembl_lookup_batch(transcript_ids):
    """
    POST a batch of transcript IDs to the Ensembl REST API.
    Returns dict: transcript_id -> gene_id (both without version numbers)
    """
    payload = json.dumps({"ids": transcript_ids}).encode("utf-8")
    req = urllib.request.Request(
        ENSEMBL_LOOKUP_URL,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept":       "application/json",
        },
        method="POST",
    )

    for attempt in range(1, RETRY_LIMIT + 1):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            result = {}
            for tid, info in data.items():
                if info and "Parent" in info:
                    # Strip version suffixes from both IDs
                    clean_tid  = tid.split(".")[0]
                    clean_gene = info["Parent"].split(".")[0]
                    result[clean_tid] = clean_gene
            return result
        except urllib.error.HTTPError as e:
            # 429 = rate limited; back off and retry
            if e.code == 429:
                wait = int(e.headers.get("Retry-After", RETRY_BACKOFF))
                print(f"  Rate limited — waiting {wait}s (attempt {attempt}/{RETRY_LIMIT})")
                time.sleep(wait)
            else:
                print(f"  [WARNING] HTTP {e.code} from Ensembl API: {e.reason}",
                      file=sys.stderr)
                return {}
        except Exception as e:
            print(f"  [WARNING] API error (attempt {attempt}/{RETRY_LIMIT}): {e}",
                  file=sys.stderr)
            if attempt < RETRY_LIMIT:
                time.sleep(RETRY_BACKOFF)

    print("  [ERROR] Ensembl API failed after all retries.", file=sys.stderr)
    return {}


def map_transcripts_to_genes(transcript_ids):
    """
    Look up gene IDs for all transcript IDs in batches.
    Transcript IDs are sent without version suffixes.
    Returns dict: transcript_id (no version) -> gene_id (no version)
    """
    # Strip versions for the API query
    stripped = [tid.split(".")[0] for tid in transcript_ids]
    total = len(stripped)
    mapping = {}

    print(f"  Querying Ensembl REST API for {total:,} transcripts "
          f"in batches of {BATCH_SIZE}...")

    for i in range(0, total, BATCH_SIZE):
        batch = stripped[i : i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        n_batches = (total + BATCH_SIZE - 1) // BATCH_SIZE
        print(f"    Batch {batch_num}/{n_batches} ({len(batch)} IDs)...")
        result = ensembl_lookup_batch(batch)
        mapping.update(result)
        # Be polite to the API between batches
        if i + BATCH_SIZE < total:
            time.sleep(0.5)

    print(f"  Mapped {len(mapping):,} / {total:,} transcripts to gene IDs")
    return mapping


# ---------------------------------------------------------------------------
# Merging
# ---------------------------------------------------------------------------

def merge(bits, halflife, transcript_to_gene):
    """
    Join on gene ID. Returns list of row dicts.
    Transcripts that couldn't be mapped or whose gene has no half-life
    entry are reported but excluded from output.
    """
    rows = []
    no_gene_map = 0
    no_halflife = 0

    for transcript_id, bit_score in bits.items():
        stripped_tid = transcript_id.split(".")[0]
        gene_id = transcript_to_gene.get(stripped_tid)

        if gene_id is None:
            no_gene_map += 1
            continue

        if gene_id not in halflife:
            no_halflife += 1
            continue

        gene_name, deg_val = halflife[gene_id]
        rows.append({
            "gene_id":           gene_id,
            "gene_name":         gene_name,
            "transcript_id":     transcript_id,
            "degradation_value": deg_val,
            "sum_of_bits_score": bit_score,
        })

    if no_gene_map:
        print(f"  [INFO] {no_gene_map:,} transcripts had no Ensembl gene mapping")
    if no_halflife:
        print(f"  [INFO] {no_halflife:,} transcripts mapped to a gene not in the half-life data")

    return rows


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_output(path, rows):
    fieldnames = ["gene_id", "gene_name", "transcript_id",
                  "degradation_value", "sum_of_bits_score"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows):,} rows → {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    for f in [args.bits, args.halflife]:
        if not Path(f).exists():
            sys.exit(f"ERROR: File not found: {f}")

    print("Loading input files...")
    bits     = load_bits(args.bits)
    halflife = load_halflife(args.halflife)

    print("\nMapping transcripts to genes via Ensembl REST API...")
    transcript_to_gene = map_transcripts_to_genes(list(bits.keys()))

    print("\nMerging...")
    rows = merge(bits, halflife, transcript_to_gene)

    print("\nWriting output...")
    write_output(args.output, rows)
    print("\nDone.")


if __name__ == "__main__":
    main()
