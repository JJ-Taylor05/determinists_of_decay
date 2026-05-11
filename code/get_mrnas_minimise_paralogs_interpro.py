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

Domain filtering with sequence identity check:
    For each gene, the canonical protein is mapped to its UniProt accession
    via Ensembl's xrefs endpoint, then Pfam domain annotations are fetched
    from the InterPro REST API. If a candidate shares a Pfam domain with an
    already-accepted sequence, a pairwise amino acid sequence identity check
    is performed against all accepted sequences that share that domain. The
    candidate is only rejected if identity exceeds IDENTITY_THRESHOLD (30%)
    with any of those sequences. Filtering is independent per set.

Train/test assignment is random (fixed seed for reproducibility).

Requirements: pip install requests biopython
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
from Bio import Align
from Bio.Align import substitution_matrices

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
ENSEMBL_URL        = "https://rest.ensembl.org"
INTERPRO_URL       = "https://www.ebi.ac.uk/interpro/api"
JSON_HDR           = {"Content-Type": "application/json"}
FASTA_HDR          = {"Content-Type": "text/x-fasta"}
MAX_WORKERS        = 50    # pure network I/O — high concurrency is safe and fast
MAX_RETRIES        = 3
RETRY_DELAY        = 2     # seconds base delay for retries
LOOKAHEAD          = 200   # how many genes to prefetch ahead of the collector
TRAIN_FRAC         = 0.9
RANDOM_SEED        = 42
IDENTITY_THRESHOLD = 0.30  # reject if amino acid identity exceeds this value


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


def get_protein_sequence(canonical_transcript_id):
    """
    Fetch the amino acid (protein) sequence for a transcript from Ensembl.
    Returns the sequence as a plain string (no FASTA header), or raises
    ValueError if none is found.
    """
    fasta = api_get(
        f"{ENSEMBL_URL}/sequence/id/{canonical_transcript_id}",
        FASTA_HDR,
        params={"type": "protein"}
    ).text.strip()
    if not fasta:
        raise ValueError(f"no protein sequence for {canonical_transcript_id}")
    # Strip the FASTA header line(s) and join the sequence lines
    lines = fasta.splitlines()
    seq = "".join(line for line in lines if not line.startswith(">"))
    if not seq:
        raise ValueError(f"empty protein sequence for {canonical_transcript_id}")
    return seq


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
# Sequence identity helper
# ---------------------------------------------------------------------------

# Build a single global aligner so BLOSUM62 is loaded once, not per call.
# Global alignment (Needleman-Wunsch) is appropriate here: we want identity
# over the full length of both sequences, not just a local high-scoring patch.
_aligner = Align.PairwiseAligner()
_aligner.mode              = "global"
_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
_aligner.open_gap_score    = -10   # standard affine gap penalties for proteins
_aligner.extend_gap_score  = -0.5


def sequence_identity(seq_a, seq_b):
    """
    Compute pairwise amino acid sequence identity between two protein
    sequences using a global Needleman-Wunsch alignment with BLOSUM62.

    A length-ratio pre-filter skips the alignment entirely when sequences
    differ so much in length that identity cannot possibly exceed
    IDENTITY_THRESHOLD (the shorter sequence could account for at most
    len_short/len_long of the longer one's residues).

    Identity is defined as:
        (identical aligned residues) / (length of the longer sequence)

    Using the longer sequence as the denominator penalises large length
    differences, consistent with the alignment-length sensitivity noted in
    the twilight-zone literature.

    Returns a float in [0, 1].
    """
    if not seq_a or not seq_b:
        return 0.0

    len_a, len_b = len(seq_a), len(seq_b)
    len_short, len_long = sorted((len_a, len_b))

    # Pre-filter: if even perfect alignment of the shorter sequence cannot
    # reach the threshold, skip the expensive alignment entirely.
    if len_short / len_long <= IDENTITY_THRESHOLD:
        return 0.0

    # Use only the single best alignment to avoid enumerating all ties.
    alignment = next(iter(_aligner.align(seq_a, seq_b)))

    # Count identical (non-gap) residue pairs in the alignment
    aligned_a, aligned_b = alignment.aligned
    identical = sum(
        sum(1 for ra, rb in zip(seq_a[s_a:e_a], seq_b[s_b:e_b]) if ra == rb)
        for (s_a, e_a), (s_b, e_b) in zip(aligned_a, aligned_b)
    )

    return identical / len_long


# ---------------------------------------------------------------------------
# Combined fetch — all API calls for one gene, runs in a worker thread
# ---------------------------------------------------------------------------
def fetch_gene_data(gene_id):
    """
    Return (gene_id, cdna, protein_seq, domains) or raise an exception with
    a description of what failed. All API calls are sequential within the
    thread, but many threads run concurrently so overall throughput is high.
    """
    canonical  = get_canonical_transcript(gene_id)
    uniprot    = get_uniprot_accession(canonical)
    domains    = get_pfam_domains(uniprot)
    cdna       = get_cdna_sequence(canonical)
    protein    = get_protein_sequence(canonical)
    return gene_id, cdna, protein, domains


# ---------------------------------------------------------------------------
# Set collector
# ---------------------------------------------------------------------------
class SetCollector:
    """
    Collects sequences for one logical group, applying an independent
    Pfam domain filter backed by a sequence identity check, and optionally
    splitting into train/test files.

    Filter logic:
        1. If the candidate shares no Pfam domain with any accepted sequence
           → accept immediately.
        2. If there is a domain overlap, compute pairwise amino acid identity
           against every accepted sequence that shares the overlapping domain.
        3. If any comparison exceeds IDENTITY_THRESHOLD → reject (likely
           paralog). Otherwise → accept (domain overlap but sequences are
           sufficiently divergent).
    """

    def __init__(self, name, target, out_train, out_test=None,
                 train_frac=TRAIN_FRAC):
        self.name       = name
        self.target     = target
        self.out_train  = open(out_train, "w")
        self.out_test   = open(out_test, "w") if out_test else None
        self.train_frac = train_frac
        self.accepted   = 0
        self.rejected   = 0
        self._rng       = random.Random(RANDOM_SEED)

        # domain → list of accepted protein sequences that carry that domain
        self._domain_to_seqs: dict[str, list[str]] = {}

    @property
    def full(self):
        return self.accepted >= self.target

    def try_accept(self, gene_id, cdna, protein, domains):
        """
        Evaluate one gene. Returns True if accepted.

        protein  — amino acid sequence string (no FASTA header)
        domains  — set of Pfam accession strings
        """
        if self.full:
            return False

        # Find which of the candidate's domains are already seen
        overlapping_domains = domains & self._domain_to_seqs.keys()

        if overlapping_domains:
            # Collect the unique accepted protein sequences that share at
            # least one domain with the candidate (deduplicated by id())
            seen_ids = set()
            seqs_to_check = []
            for domain in overlapping_domains:
                for seq in self._domain_to_seqs[domain]:
                    if id(seq) not in seen_ids:
                        seen_ids.add(id(seq))
                        seqs_to_check.append(seq)

            # Check pairwise identity against each relevant accepted sequence
            for accepted_seq in seqs_to_check:
                identity = sequence_identity(protein, accepted_seq)
                if identity > IDENTITY_THRESHOLD:
                    shared_domain = next(iter(overlapping_domains))
                    print(
                        f"  REJECTED [{self.name}]: {gene_id} — shares Pfam domain "
                        f"'{shared_domain}' AND sequence identity "
                        f"{identity:.1%} > {IDENTITY_THRESHOLD:.0%} threshold",
                        file=sys.stderr
                    )
                    self.rejected += 1
                    return False

            # Domain overlap but all identity checks passed — log and continue
            shared_domain = next(iter(overlapping_domains))
            print(
                f"  ACCEPTED [{self.name}]: {gene_id} — shared Pfam domain "
                f"'{shared_domain}' but sequence identity below threshold; "
                f"treating as sufficiently divergent",
                file=sys.stderr
            )

        # Register the candidate's domains → protein sequence mapping
        for domain in domains:
            self._domain_to_seqs.setdefault(domain, []).append(protein)

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
    are evaluated in strict CSV order so the domain-overlap and identity
    checks are deterministic regardless of which futures complete first.

    Uses the shared persistent executor passed in from main() so there is
    no per-call thread pool overhead.
    """
    gene_ids  = list(gene_ids)
    in_flight = OrderedDict()
    idx       = 0

    def submit_next():
        nonlocal idx
        while idx < len(gene_ids) and len(in_flight) < LOOKAHEAD:
            gid = gene_ids[idx]
            in_flight[gid] = executor.submit(fetch_gene_data, gid)
            idx += 1

    submit_next()

    while in_flight:
        if collector.full:
            for f in in_flight.values():
                f.cancel()
            break

        gene_id, future = next(iter(in_flight.items()))
        in_flight.pop(gene_id)

        try:
            _, cdna, protein, domains = future.result()
            collector.try_accept(gene_id, cdna, protein, domains)
        except Exception as e:
            print(f"  WARNING: {gene_id} skipped — {e}", file=sys.stderr)

        submit_next()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    global IDENTITY_THRESHOLD
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
    parser.add_argument("--identity-threshold", type=float,
                        default=IDENTITY_THRESHOLD, metavar="F",
                        help=f"Max amino acid sequence identity (0–1) allowed "
                             f"between domain-sharing sequences before rejection "
                             f"(default: {IDENTITY_THRESHOLD})")
    args = parser.parse_args()

    if not any([args.unstable, args.stable, args.average]):
        parser.error("Specify at least one of --unstable, --stable, --average")

    # Allow the threshold to be overridden at runtime
    IDENTITY_THRESHOLD = args.identity_threshold

    input_path = Path(args.csv)
    with open(input_path) as f:
        all_ids = [row[0] for row in csv.reader(f) if row]

    n = len(all_ids)
    collectors = []

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:

        if args.unstable:
            c = SetCollector("unstable", args.unstable,
                             out_train="unstable_train.fa",
                             out_test="unstable_test.fa")
            collectors.append(c)
            print(f"\n--- Collecting {args.unstable} unstable sequences "
                  f"(top of CSV) ---", file=sys.stderr)
            collect(all_ids, c, executor)

        if args.stable:
            c = SetCollector("stable", args.stable,
                             out_train="stable_train.fa",
                             out_test="stable_test.fa")
            collectors.append(c)
            print(f"\n--- Collecting {args.stable} stable sequences "
                  f"(bottom of CSV) ---", file=sys.stderr)
            collect(reversed(all_ids), c, executor)

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

    print("\n========== Summary ==========")
    for c in collectors:
        c.summary()
        c.close()


if __name__ == "__main__":
    main()
