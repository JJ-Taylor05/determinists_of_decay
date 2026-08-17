"""
test_extract_regions.py -- standalone sanity checks for extract_regions.py.
Not part of the delivered pipeline; just used to verify correctness.
Run with: python3 test_extract_regions.py
"""
import sys
import tempfile
import os

sys.path.insert(0, os.path.dirname(__file__))
from extract_regions import (
    load_orf_transcripts, load_fasta, extract_regions,
    to_rna, build_priesstess_inputs, WINDOW_SIZE,
)

n_pass = 0
n_fail = 0


def check(name, condition):
    global n_pass, n_fail
    if condition:
        n_pass += 1
        print(f"  PASS: {name}")
    else:
        n_fail += 1
        print(f"  FAIL: {name}")


# ---------------------------------------------------------------------------
print("load_fasta")
with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as f:
    f.write(">ENST001.1\nACGT\nACGT\n>ENST002.3\nTTTT\n")
    fa_path = f.name
seqs = load_fasta(fa_path)
check("strips version suffix", 'ENST001' in seqs and 'ENST002' in seqs)
check("joins multi-line sequence", seqs['ENST001'] == 'ACGTACGT')
os.unlink(fa_path)

# ---------------------------------------------------------------------------
print("load_orf_transcripts")
csv_header = "sequence_id,seq_len,has_orf,utr5_start,utr5_end,utr5_len,cds_start,cds_end,cds_len,utr3_start,utr3_end,utr3_len\n"
with tempfile.NamedTemporaryFile('w', suffix='.csv', delete=False) as f:
    f.write(csv_header)
    f.write("ENST001,300,True,0,50,50,50,250,200,250,300,50\n")
    f.write("ENST002,300,False,0,0,0,0,0,0,0,300,300\n")
    csv_path = f.name
orf_t = load_orf_transcripts(csv_path)
check("keeps has_orf==True", 'ENST001' in orf_t)
check("drops has_orf==False", 'ENST002' not in orf_t)
os.unlink(csv_path)

# ---------------------------------------------------------------------------
print("extract_regions - fully in-bounds transcript")
# 5'UTR(50) CDS(200) 3'UTR(150) = 400nt total
seq = 'A' * 50 + 'C' * 200 + 'G' * 150
coords = {'cds_start': 50, 'cds_end': 250, 'utr3_end': 400, 'seq_len': 400}
regions = extract_regions(seq, coords)
check("all 3 sites extracted (none None)", all(v is not None for v in regions.values()))
check("cds_start window is 100nt of CDS, starting AT the start codon",
      regions['cds_start'] == 'C' * 100)
check("cds_start window does NOT include any 5'UTR", 'A' not in regions['cds_start'])
check("cds_end window is 50 C's then 50 G's", regions['cds_end'] == 'C' * 50 + 'G' * 50)
check("utr3_end window is the last 100nt (all G)", regions['utr3_end'] == 'G' * 100)

# ---------------------------------------------------------------------------
print("extract_regions - very short 5'UTR no longer matters for cds_start")
seq2 = 'A' * 2 + 'C' * 200 + 'G' * 150  # only 2nt of 5'UTR
coords2 = {'cds_start': 2, 'cds_end': 202, 'utr3_end': 352, 'seq_len': 352}
regions2 = extract_regions(seq2, coords2)
check("cds_start still fully extracted despite 2nt 5'UTR", regions2['cds_start'] == 'C' * 100)

print("extract_regions - zero-length 5'UTR also fine for cds_start")
seq2b = 'C' * 200 + 'G' * 150
coords2b = {'cds_start': 0, 'cds_end': 200, 'utr3_end': 350, 'seq_len': 350}
regions2b = extract_regions(seq2b, coords2b)
check("cds_start extracted with 0nt 5'UTR", regions2b['cds_start'] == 'C' * 100)

# ---------------------------------------------------------------------------
print("extract_regions - CDS+3'UTR too short for cds_start window, but cds_end/utr3_end still fit")
# Long 5'UTR (200) pushes cds_start near the transcript's own end, leaving only
# CDS(30)+3'UTR(60)=90nt after it -- less than the 100nt cds_start needs.
# cds_end and utr3_end each only need <=100nt of LOCAL context, which is present.
seq3 = 'A' * 200 + 'C' * 30 + 'G' * 60  # 290nt total
coords3 = {'cds_start': 200, 'cds_end': 230, 'utr3_end': 290, 'seq_len': 290}
regions3 = extract_regions(seq3, coords3)
check("cds_start dropped (only 90nt remains after the start codon, need 100)", regions3['cds_start'] is None)
check("cds_end still extracted (independent of cds_start)", regions3['cds_end'] is not None)
check("utr3_end still extracted (independent of cds_start)", regions3['utr3_end'] is not None)

print("extract_regions - 3'UTR shorter than 50nt: cds_end site dropped, others unaffected")
seq4 = 'A' * 100 + 'C' * 100 + 'G' * 20  # only 20nt 3'UTR, need 50 downstream of cds_end
coords4 = {'cds_start': 100, 'cds_end': 200, 'utr3_end': 220, 'seq_len': 220}
regions4 = extract_regions(seq4, coords4)
check("cds_end dropped (window would run past the transcript)", regions4['cds_end'] is None)
check("cds_start still extracted", regions4['cds_start'] is not None)
check("utr3_end still extracted", regions4['utr3_end'] is not None)

print("extract_regions - transcript shorter than 100nt: utr3_end dropped, others may survive")
seq5 = 'A' * 10 + 'C' * 60 + 'G' * 10  # 80nt total
coords5 = {'cds_start': 10, 'cds_end': 70, 'utr3_end': 80, 'seq_len': 80}
regions5 = extract_regions(seq5, coords5)
check("utr3_end dropped (transcript shorter than the 100nt window)", regions5['utr3_end'] is None)

# ---------------------------------------------------------------------------
print("extract_regions - CSV/FASTA length mismatch drops everything")
seq6 = 'A' * 300
coords6 = {'cds_start': 50, 'cds_end': 250, 'utr3_end': 300, 'seq_len': 999}
regions6 = extract_regions(seq6, coords6)
check("entire transcript dropped (returns None, not a dict)", regions6 is None)

print("extract_regions - embedded N in real sequence passes through unchanged")
seq7 = 'A' * 50 + 'C' * 40 + 'N' + 'C' * 59 + 'G' * 150  # 300nt total
coords7 = {'cds_start': 50, 'cds_end': 150, 'utr3_end': 300, 'seq_len': 300}
regions7 = extract_regions(seq7, coords7)
check("N preserved, not dropped", regions7['cds_start'] is not None and 'N' in regions7['cds_start'])

# ---------------------------------------------------------------------------
print("to_rna")
check("T -> U", to_rna("ACGT") == "ACGU")
check("uppercases", to_rna("acgt") == "ACGU")

# ---------------------------------------------------------------------------
print("build_priesstess_inputs - end to end with a tiny synthetic dataset")
tmpdir = tempfile.mkdtemp()
csv_path = os.path.join(tmpdir, 'boundaries.csv')
stable_path = os.path.join(tmpdir, 'stable.fa')
unstable_path = os.path.join(tmpdir, 'unstable.fa')
stable_motifs_dir = os.path.join(tmpdir, 'stable_motifs')
unstable_motifs_dir = os.path.join(tmpdir, 'unstable_motifs')

with open(csv_path, 'w') as f:
    f.write(csv_header)
    for i in range(10):
        f.write(f"ENST{i:03d},400,True,0,50,50,50,250,200,250,400,150\n")

with open(stable_path, 'w') as f:
    for i in range(5):
        f.write(f">ENST{i:03d}.1\n" + ('A' * 50 + 'C' * 200 + 'G' * 150) + "\n")
with open(unstable_path, 'w') as f:
    for i in range(5, 10):
        f.write(f">ENST{i:03d}.1\n" + ('A' * 50 + 'U' * 200 + 'G' * 150) + "\n")  # different CDS content -- this IS inside the extracted windows, unlike the 5'UTR

summary = build_priesstess_inputs(csv_path, stable_path, unstable_path,
                                   stable_motifs_dir, unstable_motifs_dir, min_fg_seqs=3)
check("writes both directions for all 3 sites (6 total)", len(summary) == 6)
check("all keys are 'direction/site' format",
      all(k.count('/') == 1 for k in summary))
check("stable_motifs/cds_start: fg=5 (stable count), bg=5 (unstable count)",
      summary['stable_motifs/cds_start'] == (5, 5))
check("unstable_motifs/cds_start: fg/bg swapped relative to the other direction",
      summary['unstable_motifs/cds_start'] == (5, 5))

fg_file = os.path.join(stable_motifs_dir, 'cds_start', 'fg.txt')
bg_file = os.path.join(stable_motifs_dir, 'cds_start', 'bg.txt')
swapped_fg_file = os.path.join(unstable_motifs_dir, 'cds_start', 'fg.txt')
check("stable_motifs is its own independent top-level directory (not nested under a shared parent)",
      os.path.exists(fg_file) and not os.path.exists(os.path.join(unstable_motifs_dir, 'stable_motifs')))
check("unstable_motifs/cds_start/fg.txt written (separate top-level directory)", os.path.exists(swapped_fg_file))
with open(fg_file) as f:
    fg_lines = f.read().splitlines()
with open(swapped_fg_file) as f:
    swapped_fg_lines = f.read().splitlines()
check("stable_motifs's fg and unstable_motifs's fg are DIFFERENT sequences (fg/bg genuinely swapped, not duplicated)",
      set(fg_lines) != set(swapped_fg_lines))
with open(bg_file) as f:
    bg_lines = f.read().splitlines()
check("stable_motifs's bg matches unstable_motifs's fg (same underlying unstable sequences)",
      set(bg_lines) == set(swapped_fg_lines))
check("fg.txt lines are 100nt RNA sequences", all(len(l) == 100 and set(l) <= set('ACGUN') for l in fg_lines))

print("build_priesstess_inputs - min_fg_seqs skip logic")
summary2 = build_priesstess_inputs(csv_path, stable_path, unstable_path,
                                    stable_motifs_dir + '_2', unstable_motifs_dir + '_2', min_fg_seqs=6)
check("skips everything when min_fg_seqs > available seqs", len(summary2) == 0)

# ---------------------------------------------------------------------------
print()
print(f"TOTAL: {n_pass} passed, {n_fail} failed")
sys.exit(1 if n_fail else 0)
