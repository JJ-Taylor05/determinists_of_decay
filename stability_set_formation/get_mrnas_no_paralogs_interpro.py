#!/usr/bin/env python3
"""
Usage:
    python get_mrnas.py <input.csv> --unstable N --stable N --average N

Downloads canonical mRNA (cDNA) sequences from Ensembl for gene IDs in
column 1 of <input.csv>, which must be sorted least-stable (top) to
most-stable (bottom).

Produces up to five FASTA files:
    unstable_train.fa   — top of CSV,    90% of --unstable target
    unstable_test.fa    — top of CSV,    10% of --unstable target
    stable_train.fa     — bottom of CSV, 90% of --stable target
    stable_test.fa      — bottom of CSV, 10% of --stable target
    average.fa          — middle of CSV, --average target (no split)

Domain filtering:
    For each gene, the canonical protein is mapped to its UniProt accession
    via Ensembl's xrefs endpoint, then Pfam domain annotations are fetched
    from the InterPro REST API. Sequences sharing any Pfam domain with an
    already-accepted sequence in the same set are rejected. Filtering is
    independent per set.

Train/test assignment is random (fixed seed for reproducibility).

Requirements: pip install requests
"""

import argparse
import csv
import random
import sys
import time
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
ENSEMBL_URL  = "https://rest.ensembl.org"
INTERPRO_URL = "https://www.ebi.ac.uk/interpro/api"
JSON_HDR     = {"Content-Type": "application/json"}
FASTA_HDR    = {"Content-Type": "text/x-fasta"}
MAX_WORKERS  = 50   # pure network I/O — high concurrency is safe and fast
MAX_RETRIES  = 3
RETRY_DELAY  = 2    # seconds base delay for retries
LOOKAHEAD    = 200  # how many genes to prefetch ahead of the collector
TRAIN_FRAC   = 0.9
RANDOM_SEED  = 42


# ---------------------------------------------------------------------------
# HTTP helper — one persistent session per thread via thread-local storage
# ---------------------------------------------------------------------------
import threading
_local = threading.local()

def get_session():
    """Return a requests.Session local to the current thread."""
    if not hasattr(_local, "session"):
        _local.session = requests.Session()
    return _local.session

def api_get(url, headers=None, params=None):
    """GET with retry and rate-limit handling, using a per-thread session."""
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
# Ensembl helpers
# ---------------------------------------------------------------------------
def get_canonical_transcript(gene_id):
    """Return the version-stripped canonical transcript ID for a gene."""
    data = api_get(f"{ENSEMBL_URL}/lookup/id/{gene_id}", JSON_HDR,
                   params={"expand": 1}).json()
    canonical = data.get("canonical_transcript", "").split(".")[0]
    if not canonical:
        raise ValueError("no canonical transcript")
    return canonical


def get_uniprot_accession(canonical_transcript_id):
    """
    Map an Ensembl transcript ID to a UniProt accession via the xrefs
    endpoint. Returns the first Swiss-Prot accession found, falling back
    to TrEMBL if none, or raises ValueError if nothing is found.

    Fetches all xrefs without a server-side filter (the external_db wildcard
    is double-encoded by requests and silently returns nothing). Filtering is
    done locally against the known Ensembl dbname strings for UniProt.
    """
    xrefs = api_get(
        f"{ENSEMBL_URL}/xrefs/id/{canonical_transcript_id}",
        JSON_HDR,
        params={"all_levels": 1}
    ).json()

    # Ensembl uses 'dbname' for the database identifier; the values for
    # UniProt entries are 'Uniprot/SWISSPROT' and 'Uniprot/SPTREMBL'
    swissprot = [x["primary_id"] for x in xrefs
                 if "SWISSPROT" in x.get("dbname", "").upper()]
    trembl    = [x["primary_id"] for x in xrefs
                 if "SPTREMBL"  in x.get("dbname", "").upper()
                 or "TREMBL"    in x.get("dbname", "").upper()]

    accession = (swissprot or trembl or [None])[0]
    if not accession:
        raise ValueError(f"no UniProt accession for {canonical_transcript_id}")
    return accession


def get_cdna_sequence(canonical_transcript_id):
    """Fetch the cDNA FASTA string for a transcript."""
    sequence = api_get(
        f"{ENSEMBL_URL}/sequence/id/{canonical_transcript_id}",
        FASTA_HDR,
        params={"type": "cdna"}
    ).text.strip()
    if not sequence:
        raise ValueError(f"no cDNA sequence for {canonical_transcript_id}")
    return sequence


# ---------------------------------------------------------------------------
# InterPro helper
# ---------------------------------------------------------------------------
def get_pfam_domains(uniprot_accession):
    """
    Fetch Pfam domain accessions for a protein from the InterPro API.
    Returns a set of Pfam accession strings (e.g. {'PF00001', 'PF00002'}).
    """
    r = api_get(
        f"{INTERPRO_URL}/entry/pfam/protein/uniprot/{uniprot_accession}/",
        params={"page_size": 200}
    )
    results = r.json().get("results", [])
    domains = {entry["metadata"]["accession"]
               for entry in results
               if entry.get("metadata", {}).get("accession", "").startswith("PF")}
    if not domains:
        raise ValueError(f"no Pfam domains in InterPro for {uniprot_accession}")
    return domains


# ---------------------------------------------------------------------------
# Combined fetch — all 4 API calls for one gene, runs in a worker thread
# ---------------------------------------------------------------------------
def fetch_gene_data(gene_id):
    """
    Return (gene_id, cdna, domains) or raise an exception with a description
    of what failed. All four API calls are sequential within the thread, but
    many threads run concurrently so overall throughput is high.
    """
    canonical = get_canonical_transcript(gene_id)
    uniprot   = get_uniprot_accession(canonical)
    domains   = get_pfam_domains(uniprot)
    cdna      = get_cdna_sequence(canonical)
    return gene_id, cdna, domains


# ---------------------------------------------------------------------------
# Set collector — unchanged from previous version
# ---------------------------------------------------------------------------
class SetCollector:
    """
    Collects sequences for one logical group, applying an independent
    Pfam domain filter, and optionally splitting into train/test files.
    """

    def __init__(self, name, target, out_train, out_test=None,
                 train_frac=TRAIN_FRAC):
        self.name       = name
        self.target     = target
        self.out_train  = open(out_train, "w")
        self.out_test   = open(out_test, "w") if out_test else None
        self.train_frac = train_frac
        self.seen       = set()
        self.accepted   = 0
        self.rejected   = 0
        self._rng       = random.Random(RANDOM_SEED)

    @property
    def full(self):
        return self.accepted >= self.target

    def try_accept(self, gene_id, cdna, domains):
        """Evaluate one gene. Returns True if accepted."""
        if self.full:
            return False

        overlap = domains & self.seen
        if overlap:
            shared = next(iter(overlap))
            print(f"  REJECTED [{self.name}]: {gene_id} — shares Pfam domain "
                  f"'{shared}' with an accepted sequence", file=sys.stderr)
            self.rejected += 1
            return False

        self.seen |= domains
        self.accepted += 1

        if self.out_test and self._rng.random() > self.train_frac:
            self.out_test.write(cdna + "\n")
            dest = "test"
        else:
            self.out_train.write(cdna + "\n")
            dest = "train"

        print(f"  Accepted [{self.name}/{dest}] [{self.accepted}/{self.target}]: "
              f"{gene_id}", file=sys.stderr)
        return True

    def close(self):
        self.out_train.close()
        if self.out_test:
            self.out_test.close()

    def summary(self):
        print(f"  {self.name}: accepted={self.accepted}  rejected={self.rejected}")


# ---------------------------------------------------------------------------
# Pipelined collect — fetches ahead while evaluating, preserves CSV order
# ---------------------------------------------------------------------------
def collect(gene_ids, collector, executor):
    """
    Stream gene_ids into the collector using a sliding window of in-flight
    futures. Up to LOOKAHEAD fetches run concurrently at any time. Results
    are evaluated in strict CSV order so the domain-overlap check is
    deterministic regardless of which futures complete first.

    Uses the shared persistent executor passed in from main() so there is
    no per-call thread pool overhead.
    """
    gene_ids   = list(gene_ids)
    # OrderedDict preserves submission order for in-order evaluation
    in_flight  = OrderedDict()
    idx        = 0

    def submit_next():
        """Submit the next gene if we have headroom and genes remaining."""
        nonlocal idx
        while idx < len(gene_ids) and len(in_flight) < LOOKAHEAD:
            gid = gene_ids[idx]
            in_flight[gid] = executor.submit(fetch_gene_data, gid)
            idx += 1

    submit_next()

    while in_flight:
        if collector.full:
            # Cancel any pending futures that haven't started yet
            for f in in_flight.values():
                f.cancel()
            break

        # Evaluate the oldest submitted gene first (preserves CSV order)
        gene_id, future = next(iter(in_flight.items()))
        in_flight.pop(gene_id)

        try:
            _, cdna, domains = future.result()
            collector.try_accept(gene_id, cdna, domains)
        except Exception as e:
            print(f"  WARNING: {gene_id} skipped — {e}", file=sys.stderr)

        # Refill the window
        submit_next()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("csv", help="Input CSV (sorted least→most stable)")
    parser.add_argument("--unstable", type=int, default=0, metavar="N",
                        help="Total sequences for unstable set (split 90/10 train/test)")
    parser.add_argument("--stable",   type=int, default=0, metavar="N",
                        help="Total sequences for stable set (split 90/10 train/test)")
    parser.add_argument("--average",  type=int, default=0, metavar="N",
                        help="Sequences for average set (no train/test split)")
    args = parser.parse_args()

    if not any([args.unstable, args.stable, args.average]):
        parser.error("Specify at least one of --unstable, --stable, --average")

    input_path = Path(args.csv)
    with open(input_path) as f:
        all_ids = [row[0] for row in csv.reader(f) if row]

    n = len(all_ids)
    collectors = []

    # Single shared executor for the entire run — no repeated spin-up cost
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:

        # -- Unstable: read from top of CSV
        if args.unstable:
            c = SetCollector("unstable", args.unstable,
                             out_train="unstable_train.fa",
                             out_test="unstable_test.fa")
            collectors.append(c)
            print(f"\n--- Collecting {args.unstable} unstable sequences "
                  f"(top of CSV) ---", file=sys.stderr)
            collect(all_ids, c, executor)

        # -- Stable: read from bottom of CSV
        if args.stable:
            c = SetCollector("stable", args.stable,
                             out_train="stable_train.fa",
                             out_test="stable_test.fa")
            collectors.append(c)
            print(f"\n--- Collecting {args.stable} stable sequences "
                  f"(bottom of CSV) ---", file=sys.stderr)
            collect(reversed(all_ids), c, executor)

        # -- Average: read from middle of CSV
        if args.average:
            c = SetCollector("average", args.average,
                             out_train="average.fa", out_test=None)
            collectors.append(c)
            mid   = n // 2
            half  = (args.average * 3) // 2
            start = max(0, mid - half)
            end   = min(n, mid + half)
            print(f"\n--- Collecting {args.average} average sequences "
                  f"(middle of CSV, rows {start}–{end}) ---", file=sys.stderr)
            collect(all_ids[start:end], c, executor)

    # -- Summary
    print("\n========== Summary ==========")
    for c in collectors:
        c.summary()
        c.close()


if __name__ == "__main__":
    main()
