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
    to_rna, sliding_windows, build_priesstess_inputs,
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
check("single-line sequence", seqs['ENST002'] == 'TTTT')
os.unlink(fa_path)

# ---------------------------------------------------------------------------
print("load_orf_transcripts")
csv_header = "sequence_id,seq_len,has_orf,utr5_start,utr5_end,utr5_len,cds_start,cds_end,cds_len,utr3_start,utr3_end,utr3_len\n"
csv_rows = [
    "ENST001,300,True,0,50,50,50,250,200,250,300,50\n",   # has ORF, keep
    "ENST002,300,False,0,0,0,0,0,0,0,300,300\n",           # no ORF, drop
    "ENST003,100,True,0,10,10,10,90,80,90,100,10\n",       # short UTRs, but has ORF -> keep (padding handles it)
]
with tempfile.NamedTemporaryFile('w', suffix='.csv', delete=False) as f:
    f.write(csv_header)
    f.writelines(csv_rows)
    csv_path = f.name
orf_t = load_orf_transcripts(csv_path)
check("keeps has_orf==True", 'ENST001' in orf_t)
check("drops has_orf==False", 'ENST002' not in orf_t)
check("keeps short-UTR transcript (no length filter anymore)", 'ENST003' in orf_t)
check("coords parsed correctly", orf_t['ENST001'] == {'cds_start': 50, 'cds_end': 250, 'utr3_end': 300, 'seq_len': 300})
os.unlink(csv_path)

# ---------------------------------------------------------------------------
print("extract_regions - fully in-bounds transcript")
seq = 'A' * 150 + 'C' * 200 + 'G' * 150  # 500nt: 5'UTR(150) CDS(200) 3'UTR(150) -- all >=100nt
coords = {'cds_start': 150, 'cds_end': 350, 'utr3_end': 500, 'seq_len': 500}
regions = extract_regions(seq, coords)
check("returns all 3 sites", set(regions) == {'cds_start', 'cds_end', 'utr3_end'})
check("cds_start window is 200nt", len(regions['cds_start']) == 200)
check("cds_start window has no padding when fully in bounds", 'N' not in regions['cds_start'])
check("cds_start window centered correctly (100 A's then 100 C's)",
      regions['cds_start'] == 'A' * 100 + 'C' * 100)
check("cds_end window centered correctly (100 C's then 100 G's)",
      regions['cds_end'] == 'C' * 100 + 'G' * 100)
check("utr3_end window is last 200nt (50 C's then 150 G's)",
      regions['utr3_end'] == 'C' * 50 + 'G' * 150)

print("extract_regions - short 5'UTR needs upstream padding")
seq2 = 'A' * 20 + 'C' * 200 + 'G' * 50  # only 20nt 5'UTR, need 100 upstream of cds_start
coords2 = {'cds_start': 20, 'cds_end': 220, 'utr3_end': 270, 'seq_len': 270}
regions2 = extract_regions(seq2, coords2)
cds_start_window = regions2['cds_start']
check("pads exactly the missing 80nt with N", cds_start_window[:80] == 'N' * 80)
check("real 5'UTR sequence preserved right after padding", cds_start_window[80:100] == 'A' * 20)
check("downstream side (real CDS) untouched", cds_start_window[100:200] == 'C' * 100)
check("total length still 200", len(cds_start_window) == 200)

print("extract_regions - short 3'UTR needs downstream padding on cds_end window")
seq3 = 'A' * 100 + 'C' * 100 + 'G' * 30  # only 30nt 3'UTR, need 100 downstream of cds_end
coords3 = {'cds_start': 100, 'cds_end': 200, 'utr3_end': 230, 'seq_len': 230}
regions3 = extract_regions(seq3, coords3)
cds_end_window = regions3['cds_end']
check("real CDS tail + real 3'UTR preserved", cds_end_window[:130] == 'C' * 100 + 'G' * 30)
check("pads exactly the missing 70nt with N at the far/downstream end", cds_end_window[130:] == 'N' * 70)

print("extract_regions - transcript shorter than 200nt needs padding on utr3_end window")
seq4 = 'A' * 50 + 'C' * 50 + 'G' * 30  # 130nt total, less than the 200nt utr3_end window needs
coords4 = {'cds_start': 50, 'cds_end': 100, 'utr3_end': 130, 'seq_len': 130}
regions4 = extract_regions(seq4, coords4)
utr3_window = regions4['utr3_end']
check("pads the missing 70nt with N at the start (furthest from the 3' end site)",
      utr3_window[:70] == 'N' * 70)
check("real sequence (all of it) preserved right up to the site",
      utr3_window[70:] == seq4)

print("extract_regions - genuine embedded N in real sequence is dropped")
seq5 = 'A' * 50 + 'N' + 'C' * 199 + 'G' * 50  # N inside the real 5'UTR
coords5 = {'cds_start': 51, 'cds_end': 250, 'utr3_end': 300, 'seq_len': 300}
regions5 = extract_regions(seq5, coords5)
check("dropped because of genuine N in real sequence", regions5 is None)

print("extract_regions - CSV/FASTA length mismatch is dropped")
seq6 = 'A' * 300
coords6 = {'cds_start': 50, 'cds_end': 250, 'utr3_end': 300, 'seq_len': 999}  # seq_len disagrees with len(seq6)
regions6 = extract_regions(seq6, coords6)
check("dropped because seq_len != len(seq)", regions6 is None)

# ---------------------------------------------------------------------------
print("to_rna")
check("T -> U", to_rna("ACGT") == "ACGU")
check("uppercases", to_rna("acgt") == "ACGU")
check("leaves N alone", to_rna("ACGTN") == "ACGUN")

# ---------------------------------------------------------------------------
print("sliding_windows")
region = 'X' * 200
windows = list(sliding_windows(region, window=100, step=10))
check("produces 11 windows for 200nt region, 100nt window, 10nt step", len(windows) == 11)
check("first window starts at 0", windows[0][0] == 0)
check("last window starts at 100", windows[-1][0] == 100)
check("each window is 100nt", all(len(w[1]) == 100 for w in windows))

# ---------------------------------------------------------------------------
print("build_priesstess_inputs - end to end with a tiny synthetic dataset")
tmpdir = tempfile.mkdtemp()
csv_path = os.path.join(tmpdir, 'boundaries.csv')
stable_path = os.path.join(tmpdir, 'stable.fa')
unstable_path = os.path.join(tmpdir, 'unstable.fa')
out_dir = os.path.join(tmpdir, 'out')

with open(csv_path, 'w') as f:
    f.write(csv_header)
    for i in range(10):
        f.write(f"ENST{i:03d},300,True,0,50,50,50,250,200,250,300,50\n")

with open(stable_path, 'w') as f:
    for i in range(5):
        f.write(f">ENST{i:03d}.1\n" + ('A' * 50 + 'C' * 200 + 'G' * 50) + "\n")
with open(unstable_path, 'w') as f:
    for i in range(5, 10):
        f.write(f">ENST{i:03d}.1\n" + ('A' * 50 + 'C' * 200 + 'G' * 50) + "\n")

summary = build_priesstess_inputs(csv_path, stable_path, unstable_path, out_dir, min_fg_seqs=3)
check("writes all 33 site/window combos", len(summary) == 33)
check("fg count matches input (5 stable seqs)", all(v[0] == 5 for v in summary.values()))
check("bg count matches input (5 unstable seqs)", all(v[1] == 5 for v in summary.values()))

fg_file = os.path.join(out_dir, 'cds_start', 'rel_-100', 'fg.txt')
check("fg.txt actually written to disk", os.path.exists(fg_file))
with open(fg_file) as f:
    lines = f.read().splitlines()
check("fg.txt has 5 lines", len(lines) == 5)
check("fg.txt lines are 100nt RNA sequences", all(len(l) == 100 and set(l) <= set('ACGUN') for l in lines))

print("build_priesstess_inputs - min_fg_seqs skip logic")
summary2 = build_priesstess_inputs(csv_path, stable_path, unstable_path, out_dir + '_2', min_fg_seqs=6)
check("skips everything when min_fg_seqs > available seqs", len(summary2) == 0)

# ---------------------------------------------------------------------------
print()
print(f"TOTAL: {n_pass} passed, {n_fail} failed")
sys.exit(1 if n_fail else 0)
