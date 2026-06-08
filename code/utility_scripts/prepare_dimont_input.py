#!/usr/bin/env python3
"""
prepare_dimont_input.py

Converts a plain-header FASTA (with Ensembl Transcript IDs) into the
annotated FASTA format required by Dimont, using half-life values from
human_halflife_data_sorted.csv as the signal.

Usage:
    python prepare_dimont_input.py \\
        --fasta stable_train.fa \\
        --out dimont_stable_input.fa \\
        [--signal-mode stability|instability]

Signal modes:
    stability    — higher half-life → higher signal (use for stability motifs)
    instability  — lower half-life → higher signal (use for instability motifs)

Requirements: pip install requests biopython
"""

import argparse
import csv
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Config (mirrors get_mrnas_minimise_paralogs_interpro.py)
# ---------------------------------------------------------------------------
ENSEMBL_URL  = "https://rest.ensembl.org"
JSON_HDR     = {"Content-Type": "application/json"}
MAX_WORKERS  = 20
MAX_RETRIES  = 3
RETRY_DELAY  = 2
HALFLIFE_CSV = "human_halflife_data_sorted.csv"

# ---------------------------------------------------------------------------
# HTTP helper (reused from original script)
# ---------------------------------------------------------------------------
_local = threading.local()

def get_session():
    if not hasattr(_local, "session"):
        _local.session = requests.Session()
    return _local.session

def api_get(url, headers=None, params=None):
    session = get_session()
    for attempt in range(MAX_RETRIES):
        try:
            r = session.get(url, headers=headers, params=params, timeout=30)
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", RETRY_DELAY * (attempt + 1)))
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r
        except requests.RequestException:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                raise

# ---------------------------------------------------------------------------
# ENST → ENSG mapping (reverse of get_canonical_transcript in original script)
# ---------------------------------------------------------------------------
def get_gene_id(enst_id):
    """
    Return the ENSG gene ID for a given ENST transcript ID.
    The Ensembl /lookup/id endpoint returns 'Parent' for transcripts.
    """
    data = api_get(
        f"{ENSEMBL_URL}/lookup/id/{enst_id}",
        JSON_HDR
    ).json()
    gene_id = data.get("Parent", "")
    if not gene_id:
        raise ValueError(f"no Parent gene ID found for {enst_id}")
    return gene_id

# ---------------------------------------------------------------------------
# Load half-life CSV  (ENSG_ID, gene_name, half_life_value — no header)
# ---------------------------------------------------------------------------
def load_halflife(csv_path):
    halflife = {}
    with open(csv_path) as f:
        for row in csv.reader(f):
            if row and len(row) >= 3:
                ensg = row[0].strip()
                try:
                    halflife[ensg] = float(row[2].strip())
                except ValueError:
                    pass  # skip malformed rows
    return halflife

# ---------------------------------------------------------------------------
# Anchor position — utr_bias only
# ---------------------------------------------------------------------------
def anchor_utr_bias(seq_len):
    """
    80% of sequence length — mild bias toward 3'UTR region.
    Pair with sd=500 in Dimont to allow spread across CDS and 3'UTR.
    """
    return int(seq_len * 0.80)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--fasta", required=True, help="Input FASTA (ENST headers)")
    parser.add_argument("--out",   required=True, help="Output annotated FASTA for Dimont")
    parser.add_argument("--signal-mode", choices=["stability", "instability"],
                        default="stability",
                        help="Direction of signal weighting (default: stability)")
    args = parser.parse_args()

    # ── 1. Load half-life values ─────────────────────────────────────────────
    print("Loading half-life values...", file=sys.stderr)
    halflife = load_halflife(HALFLIFE_CSV)
    print(f"  Loaded {len(halflife)} ENSG entries", file=sys.stderr)

    # ── 2. Parse FASTA ───────────────────────────────────────────────────────
    print("Parsing FASTA...", file=sys.stderr)
    records = list(SeqIO.parse(args.fasta, "fasta"))
    print(f"  Found {len(records)} sequences", file=sys.stderr)

    # Strip version numbers from transcript IDs (e.g. ENST00000236147.6 → ENST00000236147)
    enst_ids = [r.id.split(".")[0] for r in records]

    # ── 3. Map ENST → ENSG concurrently ─────────────────────────────────────
    print(f"Mapping {len(enst_ids)} transcript IDs to gene IDs via Ensembl...",
          file=sys.stderr)

    enst_to_ensg = {}
    failed_ids   = []

    def _lookup(enst):
        return enst, get_gene_id(enst)

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(_lookup, enst): enst for enst in enst_ids}
        done = 0
        for future in as_completed(futures):
            done += 1
            if done % 100 == 0 or done == len(enst_ids):
                print(f"  Mapped {done}/{len(enst_ids)}...", file=sys.stderr)
            try:
                enst, ensg = future.result()
                enst_to_ensg[enst] = ensg
            except Exception as e:
                enst = futures[future]
                failed_ids.append(enst)
                print(f"  WARNING: could not map {enst} — {e}", file=sys.stderr)

    print(f"  Successfully mapped: {len(enst_to_ensg)}", file=sys.stderr)
    print(f"  Failed to map:       {len(failed_ids)}", file=sys.stderr)

    # ── 4. Compute signal shift so all values are positive ───────────────────
    # Dimont uses signal as a confidence/weight — must be positive.
    # We shift the full CSV value range so the global minimum becomes 0.01.
    global_min = min(halflife.values())
    global_max = max(halflife.values())
    shift = abs(global_min) + 0.01 if global_min < 0 else 0.0
    print(f"\nHalf-life value range: [{global_min:.4f}, {global_max:.4f}]",
          file=sys.stderr)
    print(f"Applying shift of {shift:.4f} to ensure positive signal values",
          file=sys.stderr)

    # ── 5. Write reformatted FASTA ───────────────────────────────────────────
    print(f"\nWriting Dimont input to {args.out}...", file=sys.stderr)

    written  = 0
    no_map   = 0
    no_value = 0

    with open(args.out, "w") as out:
        for record, enst_clean in zip(records, enst_ids):

            # ENST → ENSG
            ensg = enst_to_ensg.get(enst_clean)
            if ensg is None:
                no_map += 1
                continue

            # ENSG → half-life value
            raw_value = halflife.get(ensg)
            if raw_value is None:
                no_value += 1
                continue

            # Apply shift then optionally invert for instability mode
            signal = raw_value + shift
            if args.signal_mode == "instability":
                # Invert: the most unstable sequences get the highest signal
                signal = (global_max + shift) - signal + 0.01

            seq      = str(record.seq)
            position = anchor_utr_bias(len(seq))

            # Write annotated FASTA header
            out.write(f"> peak: {position}; signal: {signal:.6f}\n")

            # Write sequence in 60-character lines (standard FASTA)
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")

            written += 1

    # ── 6. Summary ───────────────────────────────────────────────────────────
    print(f"\n========== Summary ==========", file=sys.stderr)
    print(f"  Written:               {written}", file=sys.stderr)
    print(f"  Skipped (no ENST map): {no_map}", file=sys.stderr)
    print(f"  Skipped (no HL value): {no_value}", file=sys.stderr)


if __name__ == "__main__":
    main()
