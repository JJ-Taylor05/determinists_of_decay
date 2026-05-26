#!/usr/bin/env python3
"""
visualise_fimo.py
─────────────────
Takes a FIMO output file (fimo.tsv) and the original transcript FASTA,
finds the longest ORF in each sequence to define 5'UTR / CDS / 3'UTR,
annotates every FIMO hit by region, and writes:

  <prefix>_hits.tsv           — annotated hit table
  <prefix>_region_summary.tsv — hit counts per motif per region
  <prefix>_plots.pdf          — visualisations

Usage
─────
  # 1. Run FIMO yourself (example):
  fimo --oc fimo_out streme.txt transcripts.fa

  # 2. Run this script on the output:
  python visualise_fimo.py --fimo fimo_out/fimo.tsv --fasta transcripts.fa

Dependencies: numpy, matplotlib  (pip install numpy matplotlib)
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


# ── FASTA parser ───────────────────────────────────────────────

def parse_fasta(path):
    seqs, order = {}, []
    cur_id, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(buf).upper()
                cur_id = line[1:].split()[0]
                order.append(cur_id)
                buf = []
            elif line:
                buf.append(line)
    if cur_id:
        seqs[cur_id] = "".join(buf).upper()
    return seqs, order


# ── FIMO tsv parser ────────────────────────────────────────────

def parse_fimo(path):
    """
    Parse fimo.tsv. Returns list of hit dicts with 0-based [start, end) coords.
    FIMO uses 1-based inclusive coords, so we convert: start-1, stop unchanged.
    """
    hits = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "motif_id":     # header
                continue
            if len(fields) < 9:
                continue
            try:
                hits.append({
                    "motif_id":    fields[0],
                    "consensus":   fields[1],
                    "seq_id":      fields[2],
                    "start":       int(fields[3]) - 1,   # → 0-based
                    "end":         int(fields[4]),
                    "strand":      fields[5],
                    "score":       float(fields[6]),
                    "pvalue":      float(fields[7]),
                    "qvalue":      fields[8] if fields[8] not in (".", "") else "NA",
                    "matched_seq": fields[9] if len(fields) > 9 else "",
                })
            except (ValueError, IndexError):
                continue
    return hits


# ── Longest-ORF finder ─────────────────────────────────────────

STOPS = {"TAA", "TAG", "TGA"}

def longest_orf(seq):
    best = (None, None, 0)
    for frame in range(3):
        start = None
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon == "ATG" and start is None:
                start = i
            elif codon in STOPS and start is not None:
                end = i + 3
                if end - start > best[2]:
                    best = (start, end, end - start)
                start = None
    return best[0], best[1]


def region_boundaries(seq_len, orf_s, orf_e):
    if orf_s is None:
        return {"utr5": (0, 0), "cds": (0, 0), "utr3": (0, seq_len)}
    return {"utr5": (0, orf_s), "cds": (orf_s, orf_e), "utr3": (orf_e, seq_len)}


def classify(hit_s, hit_e, regions):
    best, best_ov = "unclassified", 0
    for name, (r0, r1) in regions.items():
        ov = max(0, min(hit_e, r1) - max(hit_s, r0))
        if ov > best_ov:
            best, best_ov = name, ov
    return best


# ── Normalised position (each region = equal third of [0,1]) ──

def norm_pos(pos, orf_s, orf_e, seq_len):
    if orf_s is None:
        return pos / seq_len if seq_len else 0.5
    utr5, cds, utr3 = orf_s, orf_e - orf_s, seq_len - orf_e
    if pos < orf_s:
        return (pos / utr5 if utr5 else 0) / 3
    elif pos < orf_e:
        return 1/3 + ((pos - orf_s) / cds if cds else 0) / 3
    else:
        return 2/3 + ((pos - orf_e) / utr3 if utr3 else 0) / 3


# ── Colours / labels ───────────────────────────────────────────

COL = {"utr5": "#4C72B0", "cds": "#55A868", "utr3": "#C44E52",
       "unclassified": "#999999"}
LAB = {"utr5": "5′UTR", "cds": "CDS", "utr3": "3′UTR",
       "unclassified": "Unclassified"}
REGIONS = ["utr5", "cds", "utr3", "unclassified"]


# ── Plots ──────────────────────────────────────────────────────

def shade(ax):
    ax.axvspan(0,   1/3, color=COL["utr5"], alpha=0.10)
    ax.axvspan(1/3, 2/3, color=COL["cds"],  alpha=0.10)
    ax.axvspan(2/3, 1,   color=COL["utr3"], alpha=0.10)
    for x in (1/3, 2/3):
        ax.axvline(x, color="#aaaaaa", lw=0.9, ls="--")


def make_plots(hits, orf_info, seq_order, out_path):
    motif_ids  = sorted({h["motif_id"] for h in hits})
    consensus  = {h["motif_id"]: h["consensus"] for h in hits}
    tab10      = plt.cm.tab10.colors

    # Per-region counts
    rc = defaultdict(lambda: defaultdict(int))
    for h in hits:
        rc[h["motif_id"]][h["region"]] += 1

    hits_by_seq = defaultdict(list)
    for h in hits:
        hits_by_seq[h["seq_id"]].append(h)

    n = len(motif_ids)
    x = np.arange(n)

    with PdfPages(out_path) as pdf:

        # Page 1 — grouped bar: raw counts
        fig, ax = plt.subplots(figsize=(max(7, n * 2 + 2), 5))
        w = 0.18
        for i, region in enumerate(REGIONS):
            vals = [rc[mid].get(region, 0) for mid in motif_ids]
            ax.bar(x + i * w, vals, w, label=LAB[region],
                   color=COL[region], edgecolor="white", lw=0.6)
        ax.set_xticks(x + w * 1.5)
        ax.set_xticklabels([f"{mid}\n({consensus.get(mid,'')})"
                            for mid in motif_ids], rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("FIMO hits")
        ax.set_title("Hit counts per region", fontsize=13)
        ax.legend(frameon=False, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # Page 2 — stacked % bar
        fig, ax = plt.subplots(figsize=(max(7, n * 2 + 2), 5))
        totals  = np.array([sum(rc[mid].values()) or 1 for mid in motif_ids])
        bottoms = np.zeros(n)
        for region in REGIONS:
            vals  = np.array([rc[mid].get(region, 0) for mid in motif_ids])
            fracs = vals / totals * 100
            ax.bar(x, fracs, 0.55, bottom=bottoms, label=LAB[region],
                   color=COL[region], edgecolor="white", lw=0.5)
            for xi, (f, b) in enumerate(zip(fracs, bottoms)):
                if f >= 5:
                    ax.text(xi, b + f / 2, f"{f:.0f}%",
                            ha="center", va="center",
                            fontsize=8, color="white", fontweight="bold")
            bottoms += fracs
        ax.set_xticks(x)
        ax.set_xticklabels([f"{mid}\n({consensus.get(mid,'')})"
                            for mid in motif_ids], rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("% of hits"); ax.set_ylim(0, 100)
        ax.set_title("Regional distribution (%)", fontsize=13)
        ax.legend(frameon=False, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # Pages 3+ — per-motif positional density
        for mid in motif_ids:
            mhits = [h for h in hits if h["motif_id"] == mid]
            if not mhits:
                continue

            npos, regs, neg_logp = [], [], []
            for h in mhits:
                info = orf_info[h["seq_id"]]
                centre = (h["start"] + h["end"]) / 2
                npos.append(norm_pos(centre, info["orf_s"], info["orf_e"], info["seq_len"]))
                regs.append(h["region"])
                neg_logp.append(-np.log10(max(h["pvalue"], 1e-300)))

            fig, (ax_h, ax_s, ax_r) = plt.subplots(
                3, 1, figsize=(11, 7),
                gridspec_kw={"height_ratios": [4, 2, 1]}, sharex=True)

            for ax in (ax_h, ax_s, ax_r):
                shade(ax)

            bins = np.linspace(0, 1, 31)
            for region in REGIONS:
                rp = [p for p, r in zip(npos, regs) if r == region]
                if rp:
                    ax_h.hist(rp, bins=bins, color=COL[region], alpha=0.85,
                              edgecolor="white", lw=0.4, label=LAB[region])
            ax_h.set_ylabel("Hit count", fontsize=10)
            ax_h.set_title(
                f"{mid}  ({consensus.get(mid, '')})   n = {len(mhits):,} hits",
                fontsize=12)
            ax_h.legend(frameon=False, fontsize=8)
            ax_h.spines[["top", "right"]].set_visible(False)

            ax_s.scatter(npos, neg_logp,
                         c=[COL.get(r, "#444") for r in regs],
                         s=12, alpha=0.6, linewidths=0)
            ax_s.set_ylabel("−log₁₀(p)", fontsize=9)
            ax_s.spines[["top", "right"]].set_visible(False)

            for p, r in zip(npos, regs):
                ax_r.axvline(p, color=COL.get(r, "#444"), alpha=0.45, lw=0.8)
            ax_r.set_yticks([])
            ax_r.set_xlim(0, 1)
            ax_r.set_xlabel("Normalised position", fontsize=10)
            ax_r.spines[["top", "right", "left"]].set_visible(False)
            for xc, lbl in ((1/6, "5′UTR"), (1/2, "CDS"), (5/6, "3′UTR")):
                ax_r.text(xc, -0.55, lbl, ha="center", va="top",
                          transform=ax_r.get_xaxis_transform(),
                          fontsize=9, color="#555")

            plt.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # Last page — transcript map (first 60 seqs)
        display = seq_order[:60]
        fig, ax = plt.subplots(figsize=(14, max(5, len(display) * 0.38 + 1.5)))

        for yi, sid in enumerate(display):
            info = orf_info.get(sid)
            if not info:
                continue
            L, orf_s, orf_e = info["seq_len"], info["orf_s"], info["orf_e"]

            if orf_s is not None:
                if orf_s > 0:
                    ax.barh(yi, orf_s, height=0.55, color=COL["utr5"], align="center")
                ax.barh(yi, orf_e - orf_s, left=orf_s, height=0.55,
                        color=COL["cds"], align="center")
                if L - orf_e > 0:
                    ax.barh(yi, L - orf_e, left=orf_e, height=0.55,
                            color=COL["utr3"], align="center")
            else:
                ax.barh(yi, L, height=0.55, color=COL["unclassified"], align="center")

            for h in hits_by_seq.get(sid, []):
                midpt = (h["start"] + h["end"]) / 2
                mi    = motif_ids.index(h["motif_id"])
                ax.plot(midpt, yi, "|", color=tab10[mi % len(tab10)],
                        markersize=8, markeredgewidth=1.5, alpha=0.85)

        ax.set_yticks(range(len(display)))
        ax.set_yticklabels(display, fontsize=6.5)
        ax.set_xlabel("Position (nt)", fontsize=10)
        suffix = f" (first {len(display)})" if len(display) < len(seq_order) else ""
        ax.set_title(f"Transcript map{suffix}", fontsize=12)
        ax.spines[["top", "right"]].set_visible(False)

        reg_patches = [mpatches.Patch(color=COL[r], label=LAB[r])
                       for r in ("utr5", "cds", "utr3", "unclassified")]
        mot_patches = [mpatches.Patch(color=tab10[i % len(tab10)],
                                      label=f"{mid} ({consensus.get(mid,'')})")
                       for i, mid in enumerate(motif_ids)]
        ax.legend(handles=reg_patches + mot_patches,
                  frameon=False, fontsize=7.5, loc="lower right", ncol=2)
        plt.tight_layout(); pdf.savefig(fig); plt.close(fig)


# ── Main ───────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fimo",       required=True, help="fimo.tsv from FIMO output folder")
    ap.add_argument("--fasta",      required=True, help="FASTA file of full transcripts")
    ap.add_argument("--out_prefix", default="fimo_annotated",
                    help="Prefix for output files (default: fimo_annotated)")
    args = ap.parse_args()

    for p, lbl in [(args.fimo, "FIMO tsv"), (args.fasta, "FASTA")]:
        if not Path(p).exists():
            sys.exit(f"ERROR: {lbl} not found: {p}")

    # 1. Load sequences
    print(f"Loading sequences … ", end="")
    seqs, order = parse_fasta(args.fasta)
    print(f"{len(seqs):,} sequences")

    # 2. Load FIMO hits
    print(f"Loading FIMO hits … ", end="")
    hits = parse_fimo(args.fimo)
    print(f"{len(hits):,} hits across "
          f"{len({h['motif_id'] for h in hits})} motif(s)")

    # 3. Find ORFs
    print("Finding longest ORFs … ", end="")
    orf_info = {}
    for sid, seq in seqs.items():
        s, e = longest_orf(seq)
        orf_info[sid] = {
            "orf_s": s, "orf_e": e, "seq_len": len(seq),
            "regions": region_boundaries(len(seq), s, e),
        }
    no_orf = sum(1 for v in orf_info.values() if v["orf_s"] is None)
    print(f"{len(seqs)-no_orf:,}/{len(seqs):,} have ORFs "
          f"({no_orf} labelled 'unclassified')")

    # 4. Annotate hits
    print("Annotating hits … ", end="")
    rc = defaultdict(lambda: defaultdict(int))
    skipped = 0
    for h in hits:
        info = orf_info.get(h["seq_id"])
        if not info:
            skipped += 1
            continue
        h["region"] = classify(h["start"], h["end"], info["regions"])
        rc[h["motif_id"]][h["region"]] += 1
    if skipped:
        print(f"\n  WARNING: {skipped} hit(s) skipped — seq ID not in FASTA")
    print(f"{len(hits)-skipped:,} hits annotated")

    # 5. Write hits TSV
    tsv = f"{args.out_prefix}_hits.tsv"
    with open(tsv, "w", newline="") as f:
        cols = ["seq_id", "motif_id", "consensus", "start", "end", "strand",
                "score", "pvalue", "qvalue", "matched_seq", "region",
                "seq_len", "orf_start", "orf_end"]
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for h in hits:
            if "region" not in h:
                continue
            info = orf_info[h["seq_id"]]
            w.writerow([h["seq_id"], h["motif_id"], h["consensus"],
                        h["start"], h["end"], h["strand"],
                        h["score"], h["pvalue"], h["qvalue"], h["matched_seq"],
                        h["region"], info["seq_len"],
                        info["orf_s"] if info["orf_s"] is not None else "NA",
                        info["orf_e"] if info["orf_e"] is not None else "NA"])
    print(f"Hits table       → {tsv}")

    # 6. Write region summary TSV
    summary = f"{args.out_prefix}_region_summary.tsv"
    consensus_map = {h["motif_id"]: h["consensus"] for h in hits}
    with open(summary, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["motif_id", "consensus", "utr5", "cds", "utr3",
                    "unclassified", "total"])
        for mid in sorted(rc):
            c = rc[mid]
            w.writerow([mid, consensus_map.get(mid, ""),
                        c.get("utr5", 0), c.get("cds", 0),
                        c.get("utr3", 0), c.get("unclassified", 0),
                        sum(c.values())])
    print(f"Region summary   → {summary}")

    # 7. Plots
    plot_path = f"{args.out_prefix}_plots.pdf"
    make_plots(hits, orf_info, order, plot_path)
    print(f"Plots            → {plot_path}")
    print("\n✓ Done.")


if __name__ == "__main__":
    main()
