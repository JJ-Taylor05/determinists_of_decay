"""
extract_regions_rnaprot.py

Builds RNAProt-format FASTA input files for a stable-vs-unstable mRNA
structure search around three landmarks:
  - the CDS start (start codon)
  - the CDS end (stop codon)
  - the 3' UTR end

Design summary (see conversation for the reasoning behind these choices):
  - Every transcript with a called ORF (has_orf == True) is used -- see
    load_orf_transcripts().
  - Each site gets ONE 100nt window, extracted directly from the transcript
    with NO padding:
      - cds_start: 100nt starting AT the start codon, extending into the
                   CDS. Never touches the 5'UTR, so a short/absent 5'UTR
                   is no longer a problem at all for this site.
      - cds_end:   50nt either side of the stop codon.
      - utr3_end:  the last 100nt of the transcript (entirely upstream --
                   no sequence exists past the annotated 3' end in a
                   mature mRNA).
  - Since there's no padding, a transcript's window at a given site is
    simply dropped if it doesn't fit (e.g. a CDS shorter than 100nt for
    cds_start's downstream bound, or a 3'UTR shorter than 50nt for
    cds_end's downstream bound, or a transcript shorter than 100nt for
    utr3_end's upstream bound). Sites are extracted independently per
    transcript, so failing one site doesn't drop that transcript from the
    other two.
  - A transcript is dropped entirely (all sites) only if the CSV/FASTA
    disagree on its length -- a version/isoform mismatch where the
    coordinates can't be trusted against that sequence at all.
  - SIX output FASTA files total (one per class x site). Unlike PRIESSTESS,
    RNAProt doesn't need pre-built fg/bg pairs baked into separate
    directories per direction: you point `rnaprot gt --in --neg-in` at
    whichever two of these six files you want to compare, and swap them to
    get the other direction. So each class/site combination is written
    exactly once.

Usage:
    python3 extract_regions_rnaprot.py \
        --csv transcript_boundaries.csv \
        --stable stable_train.fa \
        --unstable unstable_train.fa \
        --out-dir rnaprot_inputs
"""

import argparse
import csv
import os


WINDOW_SIZE = 100

# site name (internal) -> short name used in output filenames
SITE_FILENAME = {
    'cds_start': 'start_cds',
    'cds_end':   'end_cds',
    'utr3_end':  'end_utr3',
}


def site_windows(coords):
    """Given a transcript's cds_start/cds_end/utr3_end coordinates, return
    the (start, end) bounds -- in the transcript's own coordinate system --
    of each site's 100nt window."""
    return {
        'cds_start': (coords['cds_start'], coords['cds_start'] + 100),
        'cds_end':   (coords['cds_end'] - 50, coords['cds_end'] + 50),
        'utr3_end':  (coords['utr3_end'] - 100, coords['utr3_end']),
    }


def load_orf_transcripts(csv_path):
    """Read transcript_boundaries.csv and keep every transcript with a
    called ORF. Returns {unversioned_id: {cds_start, cds_end, utr3_end, seq_len}}.

    Not a length filter -- has_orf == True is the only requirement, since
    there's no coordinate to build a window around otherwise. Any window
    that doesn't fit within a given transcript is handled per-site in
    extract_regions(), not filtered out here.
    """
    orf_transcripts = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['has_orf'] != 'True':
                continue
            orf_transcripts[row['sequence_id']] = {
                'cds_start': int(row['cds_start']),
                'cds_end':   int(row['cds_end']),
                'utr3_end':  int(row['utr3_end']),
                'seq_len':   int(row['seq_len']),
            }
    return orf_transcripts


def load_fasta(path):
    """Read a FASTA file into {unversioned_id: full_sequence}."""
    seqs = {}
    current_id = None
    buffer = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(buffer)
                current_id = line[1:].split('.')[0]
                buffer = []
            else:
                buffer.append(line)
        if current_id is not None:
            seqs[current_id] = ''.join(buffer)
    return seqs


def extract_regions(seq, coords):
    """Cut the three 100nt windows (cds_start, cds_end, utr3_end) out of a
    full transcript sequence, with no padding.

    Returns {site_name: 100nt region string or None}. A site is None if its
    window doesn't fit within this specific transcript (extracted
    independently per site -- one site failing doesn't affect the others).

    Returns None (not a dict) if the CSV/FASTA disagree on this
    transcript's length -- a version/isoform mismatch means none of the
    coordinates can be trusted against this sequence at all.
    """
    if coords['seq_len'] != len(seq):
        return None

    regions = {}
    for site_name, (start, end) in site_windows(coords).items():
        if start < 0 or end > len(seq):
            regions[site_name] = None
            continue
        region = seq[start:end]
        assert len(region) == WINDOW_SIZE
        regions[site_name] = region

    return regions


def to_rna(seq):
    """DNA -> RNA alphabet for RNAProt (T -> U, uppercase)."""
    return seq.upper().replace('T', 'U')


def write_fasta(path, records):
    """records: list of (id, seq). Writes standard FASTA with unique headers,
    as required by RNAProt."""
    with open(path, 'w') as f:
        for seq_id, seq in records:
            f.write(f'>{seq_id}\n{seq}\n')


def build_rnaprot_inputs(boundaries_csv, stable_fa, unstable_fa, out_dir):
    """Full pipeline: filter -> load -> extract -> write six FASTA files,
    one per (class, site) combination:
        out_dir/stable_start_cds.fa    out_dir/unstable_start_cds.fa
        out_dir/stable_end_cds.fa      out_dir/unstable_end_cds.fa
        out_dir/stable_end_utr3.fa     out_dir/unstable_end_utr3.fa
    """
    orf_transcripts = load_orf_transcripts(boundaries_csv)
    stable = load_fasta(stable_fa)
    unstable = load_fasta(unstable_fa)

    collected = {'stable': {site: [] for site in ('cds_start', 'cds_end', 'utr3_end')},
                 'unstable': {site: [] for site in ('cds_start', 'cds_end', 'utr3_end')}}
    n_seq_len_mismatch = 0
    n_site_out_of_bounds = {'cds_start': 0, 'cds_end': 0, 'utr3_end': 0}

    for class_name, fasta_dict in [('stable', stable), ('unstable', unstable)]:
        for tid, seq in fasta_dict.items():
            if tid not in orf_transcripts:
                continue
            regions = extract_regions(seq, orf_transcripts[tid])
            if regions is None:
                n_seq_len_mismatch += 1
                continue
            for site, region in regions.items():
                if region is None:
                    n_site_out_of_bounds[site] += 1
                    continue
                collected[class_name][site].append((tid, region))

    print(f"  (dropped {n_seq_len_mismatch} transcripts for CSV/FASTA length mismatch)")
    for site, n in n_site_out_of_bounds.items():
        print(f"  (dropped {n} {site} windows for not fitting within their transcript)")

    os.makedirs(out_dir, exist_ok=True)

    summary = {}
    for class_name in ('stable', 'unstable'):
        for site in ('cds_start', 'cds_end', 'utr3_end'):
            records = collected[class_name][site]
            if not records:
                print(f"  [skip] {class_name}/{site}: no sequences extracted")
                continue

            filename = f"{class_name}_{SITE_FILENAME[site]}.fa"
            out_path = os.path.join(out_dir, filename)
            write_fasta(out_path, [(tid, to_rna(seq)) for tid, seq in records])
            summary[filename] = len(records)

    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', required=True, help='transcript_boundaries.csv')
    parser.add_argument('--stable', required=True, help='stable_train.fa')
    parser.add_argument('--unstable', required=True, help='unstable_train.fa')
    parser.add_argument('--out-dir', required=True,
                         help='output directory for the six RNAProt input FASTA files')
    args = parser.parse_args()

    summary = build_rnaprot_inputs(
        args.csv, args.stable, args.unstable, args.out_dir,
    )

    print(f"Wrote {len(summary)} FASTA files to {args.out_dir}/")
    for filename in sorted(summary):
        print(f"  {filename}: {summary[filename]} sequences")


if __name__ == '__main__':
    main()
