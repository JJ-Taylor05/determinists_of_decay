"""
extract_regions.py

Builds PRIESSTESS-format foreground/background input files for a
stable-vs-unstable mRNA motif search around three landmarks:
  - the CDS start (start codon)
  - the CDS end (stop codon)
  - the 3' UTR end

Design summary (see conversation for the reasoning/citations behind these
choices):
  - Every transcript with a called ORF (has_orf == True) is used -- there's
    no minimum 5'UTR/CDS/3'UTR length filter anymore (see
    load_orf_transcripts()).
  - CDS-start and CDS-end windows are +/-100nt, symmetric around the site.
    The 3'UTR-end window is the last 200nt of the transcript (no downstream
    sequence exists past the annotated 3' end in a mature mRNA). Windowing
    around these three landmarks follows the logic in Geisberg et al. 2014
    (Cell, https://doi.org/10.1016/j.cell.2013.12.026), which mapped
    stabilizing/destabilizing elements to exactly these positions.
  - If a transcript's UTR/CDS is too short to fill a window, the missing
    portion is padded with 'N' at the boundary furthest from the
    biological site, rather than dropping the transcript (see
    extract_regions() docstring for why this is safe -- Sample et al. 2019,
    Nat Biotechnol, doi:10.1038/s41587-019-0164-5; Bailey et al. 2021,
    Bioinformatics, doi:10.1093/bioinformatics/btab203). A transcript is
    only dropped if a window's REAL sequence already contains an 'N' (a
    genuine assembly gap), or if the CSV/FASTA disagree on that
    transcript's length (a version/isoform mismatch -- see extract_regions()).
    Note: if a CDS is shorter than 200nt, the cds_start and cds_end windows
    will overlap each other -- padding doesn't change that.
  - Within each 200nt region, a 100nt window is slid in 10nt steps,
    producing 11 sub-windows per site (33 total across the 3 sites), each
    run through PRIESSTESS separately to localize positional signal.

Usage:
    python3 extract_regions.py \
        --csv transcript_boundaries.csv \
        --stable stable_train.fa \
        --unstable unstable_train.fa \
        --out priesstess_input
"""

import argparse
import csv
import os


SITE_REGION_START = {
    'cds_start': -100,
    'cds_end':   -100,
    'utr3_end':  -200,
}


def load_orf_transcripts(csv_path):
    """Read transcript_boundaries.csv and keep every transcript with a
    called ORF. Returns {unversioned_id: {cds_start, cds_end, utr3_end, seq_len}}.

    This is no longer a length filter -- earlier versions of this function
    dropped transcripts with a short 5'UTR/CDS/3'UTR, but extract_regions()
    now pads any shortfall with 'N' instead of discarding the transcript,
    so that filtering was only throwing away recoverable data. The one
    thing padding can't invent is a set of CDS coordinates to build a
    window around, so has_orf == True is the only requirement left here.
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


def extract_regions(seq, coords, window_size=200):
    """Cut the three 200nt windows (cds_start, cds_end, utr3_end) out of a
    full transcript sequence.

    Previously this returned None (dropping the whole transcript) if a
    window ran off either end of the sequence. Now it instead PADS the
    missing part with 'N' so every transcript can still be used.

    Why this is safe to do here: each window is anchored at a real
    biological landmark (e.g. cds_start sits at relative position 100
    within its window), and any shortfall in available sequence always
    occurs at the far edge of the window -- furthest from that anchor --
    simply because that's where the real transcript sequence runs out.
    Padding is therefore automatically confined to the boundary furthest
    from the site of interest, never the centre. See Sample et al. 2019
    (Nat Biotechnol, doi:10.1038/s41587-019-0164-5) for the same
    boundary-padding logic applied to 5' UTRs of varying length, and
    Bailey et al. 2021 (Bioinformatics, doi:10.1093/bioinformatics/btab203)
    for STREME's explicit handling of ambiguous 'N' characters without
    introducing artifacts.

    A transcript is only dropped now if the REAL (non-padded) part of a
    window already contains an 'N' -- e.g. a genuine sequencing/assembly
    gap in the source data -- since that's a data-quality issue padding
    can't fix.

    Returns {site_name: 200nt region (str, possibly N-padded)}, or None
    if a genuine embedded 'N' was found in the real sequence.
    """
    cds_start = coords['cds_start']
    cds_end = coords['cds_end']
    utr3_end = coords['utr3_end']

    # Sanity check: the boundaries table's seq_len should always equal the
    # actual FASTA sequence length for this ID, since utr3_end (the CDS/UTR
    # coordinates' frame of reference) is defined relative to that same
    # transcript. If they disagree, the FASTA and CSV rows are describing
    # different transcript versions/isoforms that happen to share an
    # unversioned ID -- coordinates from one can't be trusted against
    # sequence from the other, so padding math would be meaningless here.
    # This is rare (~9 in 2,698 in this dataset) but real, so we drop these
    # explicitly rather than let padding silently produce garbage.
    if coords['seq_len'] != len(seq):
        return None

    windows = {
        'cds_start': (cds_start - 100, cds_start + 100),
        'cds_end':   (cds_end - 100, cds_end + 100),
        'utr3_end':  (utr3_end - 200, utr3_end),
    }

    regions = {}
    for site_name, (start, end) in windows.items():
        # How much of the requested window actually falls within the
        # transcript, and how much is missing off each end?
        left_pad = max(0, -start)
        right_pad = max(0, end - len(seq))
        real_start = max(0, start)
        real_end = min(len(seq), end)

        real_seq = seq[real_start:real_end]
        if 'N' in real_seq.upper():
            return None  # genuine data-quality gap, not something padding fixes

        region = ('N' * left_pad) + real_seq + ('N' * right_pad)
        assert len(region) == window_size, (
            f"{site_name}: built a {len(region)}nt region, expected {window_size}nt"
        )
        regions[site_name] = region

    return regions


def to_rna(seq):
    """DNA -> RNA alphabet for PRIESSTESS (T -> U, uppercase)."""
    return seq.upper().replace('T', 'U')


def sliding_windows(region, window=100, step=10):
    """Yield (extraction_relative_start, subsequence) across a region."""
    n_windows = (len(region) - window) // step + 1
    for i in range(n_windows):
        start = i * step
        yield start, region[start:start + window]


def build_priesstess_inputs(boundaries_csv, stable_fa, unstable_fa, out_dir,
                             window=100, step=10, min_fg_seqs=1000):
    """Full pipeline: filter -> load -> extract -> slide -> write fg/bg files.

    min_fg_seqs: PRIESSTESS/STREME needs a minimum number of foreground
    sequences to run motif discovery reliably. Any site/offset combination
    whose fg file would contain fewer than min_fg_seqs sequences is skipped
    entirely (not written), with a warning printed, rather than silently
    handing PRIESSTESS a file that will just fail later.

    Output layout:
        out_dir/<site>/rel_<site_relative_offset>/fg.txt   (stable)
        out_dir/<site>/rel_<site_relative_offset>/bg.txt   (unstable)
    """
    orf_transcripts = load_orf_transcripts(boundaries_csv)
    stable = load_fasta(stable_fa)
    unstable = load_fasta(unstable_fa)

    collected = {'fg': {}, 'bg': {}}

    for class_name, fasta_dict in [('fg', stable), ('bg', unstable)]:
        for tid, seq in fasta_dict.items():
            if tid not in orf_transcripts:
                continue
            regions = extract_regions(seq, orf_transcripts[tid])
            if regions is None:
                continue
            for site, region_seq in regions.items():
                collected[class_name].setdefault(site, {})
                for win_start, subseq in sliding_windows(region_seq, window, step):
                    offset = SITE_REGION_START[site] + win_start
                    collected[class_name][site].setdefault(offset, []).append(subseq)

    summary = {}
    for site in SITE_REGION_START:
        # Union of offsets seen on either side -- a site could in principle
        # be entirely absent from one class (e.g. if every bg transcript
        # happened to fail the seq_len/embedded-N checks in extract_regions
        # for this particular site), so don't assume fg's offsets are a
        # superset of bg's.
        all_offsets = set(collected['fg'].get(site, {})) | set(collected['bg'].get(site, {}))
        for offset in sorted(all_offsets):
            fg_seqs = collected['fg'].get(site, {}).get(offset, [])
            bg_seqs = collected['bg'].get(site, {}).get(offset, [])

            if not bg_seqs:
                print(f"  [skip] {site}/rel_{offset}: no background sequences available")
                continue
            if len(fg_seqs) < min_fg_seqs:
                print(f"  [skip] {site}/rel_{offset}: only {len(fg_seqs)} foreground "
                      f"sequences, below min_fg_seqs={min_fg_seqs}")
                continue

            win_dir = os.path.join(out_dir, site, f"rel_{offset}")
            os.makedirs(win_dir, exist_ok=True)

            with open(os.path.join(win_dir, 'fg.txt'), 'w') as f:
                f.write('\n'.join(to_rna(s) for s in fg_seqs) + '\n')
            with open(os.path.join(win_dir, 'bg.txt'), 'w') as f:
                f.write('\n'.join(to_rna(s) for s in bg_seqs) + '\n')

            summary[f"{site}/rel_{offset}"] = (len(fg_seqs), len(bg_seqs))

    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', required=True, help='transcript_boundaries.csv')
    parser.add_argument('--stable', required=True, help='stable_train.fa')
    parser.add_argument('--unstable', required=True, help='unstable_train.fa')
    parser.add_argument('--out', required=True, help='output directory for PRIESSTESS input files')
    parser.add_argument('--window', type=int, default=100, help='sliding sub-window size (default 100)')
    parser.add_argument('--step', type=int, default=10, help='sliding sub-window step (default 10)')
    parser.add_argument('--min-fg-seqs', type=int, default=1000,
                         help='minimum number of foreground sequences required for '
                              'PRIESSTESS to run; site/window combos below this are '
                              'skipped (default 1000)')
    args = parser.parse_args()

    summary = build_priesstess_inputs(
        args.csv, args.stable, args.unstable, args.out,
        window=args.window, step=args.step, min_fg_seqs=args.min_fg_seqs,
    )

    print(f"Wrote {len(summary)} fg/bg file pairs to {args.out}/")
    for k in sorted(summary):
        fg_n, bg_n = summary[k]
        print(f"  {k}: fg={fg_n}  bg={bg_n}")


if __name__ == '__main__':
    main()
