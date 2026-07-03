"""
filters/region_filter.py
========================
Assign each motif hit in a narrowPeak file to a transcript region (5'UTR,
CDS, or 3'UTR) and return per-region hit sets for use in scoring.

The region boundaries are determined by finding the longest ORF (open reading
frame) in each sequence from a provided FASTA file. The ORF finder and
classification logic are adapted directly from visualise_fimo.py so that
region boundaries are defined consistently across the pipeline.

Longest-ORF definition (from visualise_fimo.py)
------------------------------------------------
Scans all three forward reading frames for ATG..stop codons (TAA, TAG, TGA).
The ORF with the greatest nucleotide length is selected. If no ORF is found
in a sequence, the entire sequence is treated as 3'UTR (unclassified CDS and
5'UTR have zero length), matching the behaviour in visualise_fimo.py.

Region boundaries (0-based, half-open intervals)
-------------------------------------------------
  5'UTR : [0,       orf_start)
  CDS   : [orf_start, orf_end)
  3'UTR : [orf_end,   seq_len)

Hit assignment
--------------
A hit is assigned to the region it overlaps *most* (by nucleotide count).
This mirrors the classify() function in visualise_fimo.py. Hits with no
overlap with any defined region are labelled "unclassified" and excluded
from all three output sets.

Main public function
--------------------
partition_hits_by_region(narrowpeak_path, fasta_path, allowed_motifs=None)

Returns
-------
Three sets of (seq_id, motif_name) tuples:
  utr5_hits, cds_hits, utr3_hits

These are passed to read_best_site_for_region() (also in this module) instead
of the standard read_best_site(), because standard filtering works on motif
names alone, whereas region filtering must also check coordinates.

Usage example
-------------
    from filters.region_filter import partition_hits_by_region, read_best_site_for_region
    from motif_scoring import compute_sum_of_bits

    utr5_hits, cds_hits, utr3_hits = partition_hits_by_region(
        "stable_best_site.narrowPeak", "transcripts.fa"
    )
    utr5_raw   = read_best_site_for_region("stable_best_site.narrowPeak", utr5_hits)
    utr5_score = compute_sum_of_bits(utr5_raw)
    # ... repeat for cds and utr3, then write_output() three times
"""

import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# FASTA PARSER  (adapted from visualise_fimo.py)
# ---------------------------------------------------------------------------
# Reads a multi-FASTA file into a dict: { seq_id: sequence_string }
# The sequence ID is everything after ">" up to the first whitespace,
# which is exactly what appears in narrowPeak col 0.

def _parse_fasta(path):
    seqs = {}
    cur_id, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(buf).upper()
                cur_id = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line)
    if cur_id:
        seqs[cur_id] = "".join(buf).upper()
    return seqs


# ---------------------------------------------------------------------------
# LONGEST ORF FINDER  (copied from visualise_fimo.py unchanged)
# ---------------------------------------------------------------------------
# Scans all three forward reading frames.
# Returns (orf_start, orf_end) as 0-based half-open coordinates,
# or (None, None) if no complete ATG..stop ORF is found.

_STOPS = {"TAA", "TAG", "TGA"}

def _longest_orf(seq):
    best = (None, None, 0)
    for frame in range(3):
        start = None
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon == "ATG" and start is None:
                start = i
            elif codon in _STOPS and start is not None:
                end = i + 3
                if end - start > best[2]:
                    best = (start, end, end - start)
                start = None
    return best[0], best[1]


# ---------------------------------------------------------------------------
# REGION BOUNDARIES  (adapted from visualise_fimo.py)
# ---------------------------------------------------------------------------
# Returns a dict of region name → (start, end) 0-based half-open intervals.
# If no ORF was found, UTR5 and CDS are empty; entire sequence is UTR3.

def _region_boundaries(seq_len, orf_s, orf_e):
    if orf_s is None:
        return {"utr5": (0, 0), "cds": (0, 0), "utr3": (0, seq_len)}
    return {
        "utr5": (0, orf_s),
        "cds":  (orf_s, orf_e),
        "utr3": (orf_e, seq_len),
    }


# ---------------------------------------------------------------------------
# HIT CLASSIFIER  (adapted from visualise_fimo.py)
# ---------------------------------------------------------------------------
# Assigns a hit to the region it overlaps most.
# Returns "utr5", "cds", "utr3", or "unclassified".

def _classify(hit_s, hit_e, regions):
    best, best_ov = "unclassified", 0
    for name, (r0, r1) in regions.items():
        ov = max(0, min(hit_e, r1) - max(hit_s, r0))
        if ov > best_ov:
            best, best_ov = name, ov
    return best


# ---------------------------------------------------------------------------
# NARROWPEAK HIT READER  (internal helper)
# ---------------------------------------------------------------------------
# Reads each line of the narrowPeak file and yields a small dict per hit.
# We need coordinates here (not just the bit score) so we can classify hits,
# hence this is a separate reader from read_best_site() in motif_scoring.py.
#
# narrowPeak columns (0-based):
#   0: seq_id   1: start (0-based)   2: end (exclusive)
#   3: motif    6: bit_score

def _iter_hits(narrowpeak_path):
    """Yield one dict per valid hit line in the narrowPeak file."""
    with open(narrowpeak_path) as fh:
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
            try:
                yield {
                    "seq_id":    cols[0],
                    "start":     int(cols[1]),
                    "end":       int(cols[2]),
                    "motif":     cols[3],
                    "bit_score": float(cols[6]),
                }
            except ValueError:
                print(
                    f"  [WARNING] Line {lineno}: could not parse coordinates "
                    "or score — skipping",
                    file=sys.stderr,
                )


# ---------------------------------------------------------------------------
# PUBLIC FUNCTION 1 — partition hits by region
# ---------------------------------------------------------------------------
# This is the main entry point. It:
#   1. Parses the FASTA to find ORF boundaries for each sequence.
#   2. Reads every hit in the narrowPeak file.
#   3. Classifies each hit as utr5 / cds / utr3 / unclassified.
#   4. Returns three sets of (seq_id, motif_name) tuples, one per region.
#
# The optional `allowed_motifs` argument lets you combine this filter with a
# streme_filters allowlist — only hits whose motif is in allowed_motifs AND
# whose coordinates fall in the target region are included.

def partition_hits_by_region(narrowpeak_path, fasta_path, allowed_motifs=None):
    """
    Classify every hit in a narrowPeak file into a transcript region.

    Parameters
    ----------
    narrowpeak_path : str or Path
    fasta_path      : str or Path
    allowed_motifs  : set[str] or None
        Optional pre-filter. If provided, only motifs in this set are
        considered (e.g. the output of filter_by_evalue()).

    Returns
    -------
    utr5_hits, cds_hits, utr3_hits : set of (seq_id, motif_name) tuples
        Each set contains hits belonging to that region.
        Hits whose seq_id is not in the FASTA are skipped with a warning.
        Hits in regions of length 0 (e.g. no 5'UTR before the ORF) are
        labelled "unclassified" and excluded from all three sets.
    """
    # Step 1 — build ORF boundary lookup from FASTA
    print(f"  Reading FASTA for ORF boundaries: {fasta_path}")
    seqs = _parse_fasta(fasta_path)

    # Pre-compute region boundaries for every sequence
    # { seq_id: {"utr5": (s,e), "cds": (s,e), "utr3": (s,e)} }
    region_map = {}
    no_orf = 0
    for seq_id, seq in seqs.items():
        orf_s, orf_e = _longest_orf(seq)
        if orf_s is None:
            no_orf += 1
        region_map[seq_id] = _region_boundaries(len(seq), orf_s, orf_e)

    print(
        f"  {len(seqs):,} sequences in FASTA; "
        f"{no_orf:,} had no detectable ORF (entire sequence treated as 3'UTR)"
    )

    # Step 2 — read hits and classify
    utr5_hits = set()
    cds_hits  = set()
    utr3_hits = set()
    n_skipped_fasta  = 0
    n_skipped_motif  = 0
    n_unclassified   = 0

    for hit in _iter_hits(narrowpeak_path):
        seq_id = hit["seq_id"]
        motif  = hit["motif"]

        # Optional motif pre-filter
        if allowed_motifs is not None and motif not in allowed_motifs:
            n_skipped_motif += 1
            continue

        # Check seq is in FASTA (needed for ORF boundaries)
        if seq_id not in region_map:
            n_skipped_fasta += 1
            continue

        region = _classify(hit["start"], hit["end"], region_map[seq_id])

        if region == "utr5":
            utr5_hits.add((seq_id, motif))
        elif region == "cds":
            cds_hits.add((seq_id, motif))
        elif region == "utr3":
            utr3_hits.add((seq_id, motif))
        else:
            n_unclassified += 1

    if n_skipped_fasta:
        print(
            f"  [WARNING] {n_skipped_fasta:,} hits skipped — seq_id not in FASTA.",
            file=sys.stderr,
        )
    if n_skipped_motif:
        print(f"  {n_skipped_motif:,} hits excluded by motif pre-filter.")
    if n_unclassified:
        print(f"  {n_unclassified:,} hits unclassified (zero-length region or no overlap).")

    print(
        f"  Hits assigned — 5'UTR: {len(utr5_hits):,}  "
        f"CDS: {len(cds_hits):,}  3'UTR: {len(utr3_hits):,}"
    )
    return utr5_hits, cds_hits, utr3_hits


# ---------------------------------------------------------------------------
# PUBLIC FUNCTION 2 — region-aware narrowPeak reader
# ---------------------------------------------------------------------------
# Standard read_best_site() filters by motif name only. For regional scoring
# we need to filter by (seq_id, motif_name) pair, because the same motif can
# hit the same sequence in multiple regions and we want only the hits that
# were assigned to a specific region.
#
# Like read_best_site(), we keep the highest bit score if a (seq, motif) pair
# appears more than once in the allowed set.

def read_best_site_for_region(path, region_hit_set):
    """
    Parse a narrowPeak file, keeping only hits in ``region_hit_set``.

    Parameters
    ----------
    path           : str or Path
    region_hit_set : set of (seq_id, motif_name) tuples
        One of the three sets returned by partition_hits_by_region().

    Returns
    -------
    scores : dict[seq_id -> dict[motif_name -> float]]
        Same structure as read_best_site() in motif_scoring.py.
    """
    scores = defaultdict(dict)

    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 7:
                continue

            seq_id = cols[0]
            motif  = cols[3]

            # Only keep this hit if it was assigned to the target region
            if (seq_id, motif) not in region_hit_set:
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

            existing = scores[seq_id].get(motif)
            if existing is None or bit_score > existing:
                scores[seq_id][motif] = bit_score

    return scores
