#!/bin/bash
# run_rnaprot_pipeline.sh
#
# Runs RNAProt on the six FASTA files produced by extract_regions_rnaprot.py:
#   <in-dir>/stable_start_cds.fa    <in-dir>/unstable_start_cds.fa
#   <in-dir>/stable_end_cds.fa      <in-dir>/unstable_end_cds.fa
#   <in-dir>/stable_end_utr3.fa     <in-dir>/unstable_end_utr3.fa
#
# For each of the 3 sites, runs RNAProt in BOTH directions (6 runs total):
#   direction "unstableFg": unstable = positive (fg), stable = negative (bg)
#   direction "stableFg":   stable   = positive (fg), unstable = negative (bg)
# Both directions are needed because RNAProt's model discriminates
# foreground FROM background -- a single direction only reveals what's
# enriched in one class relative to the other, not the reverse.
#
# Usage:
#   ./run_rnaprot_pipeline.sh <in-dir> [out-dir]
#
#   in-dir  : directory containing the six FASTA files from
#             extract_regions_rnaprot.py (required)
#   out-dir : where to write RNAProt's gt/train/eval output
#             (default: rnaprot_runs)
#
# Example:
#   python3 extract_regions_rnaprot.py --csv transcript_boundaries.csv \
#       --stable stable_train.fa --unstable unstable_train.fa \
#       --out-dir rnaprot_inputs
#   ./run_rnaprot_pipeline.sh rnaprot_inputs rnaprot_runs

set -euo pipefail

IN_DIR=${1:?Usage: $0 <in-dir-from-extract_regions_rnaprot.py> [out-dir]}
OUT_DIR=${2:-rnaprot_runs}
mkdir -p "$OUT_DIR"

SITES=("start_cds" "end_cds" "end_utr3")

run_direction () {
    local tag=$1       # e.g. start_cds_unstableFg
    local pos_fa=$2     # positive class fasta (fg)
    local neg_fa=$3     # negative class fasta (bg)
    local outdir="$OUT_DIR/$tag"

    echo "=== $tag ==="
    echo "    fg (positive): $pos_fa"
    echo "    bg (negative): $neg_fa"

    rnaprot gt --in "$pos_fa" --neg-in "$neg_fa" \
        --out "${outdir}_gt_out" --str --report --seed 1

    rnaprot train --in "${outdir}_gt_out" --out "${outdir}_train_out" \
        --verbose-train --seed 1
    # NOTE: default --epochs 200 / --patience 30 (full run).

    rnaprot eval --gt-in "${outdir}_gt_out" --train-in "${outdir}_train_out" \
        --out "${outdir}_eval_out" --report \
        --nr-top-sites 100 200 500 --motif-size 5 7 9
}

for site in "${SITES[@]}"; do
    stable_fa="$IN_DIR/stable_${site}.fa"
    unstable_fa="$IN_DIR/unstable_${site}.fa"

    for f in "$stable_fa" "$unstable_fa"; do
        if [[ ! -s "$f" ]]; then
            echo "!! Missing or empty: $f -- skipping ${site}" >&2
            continue 2
        fi
    done

    # Direction A: unstable = fg (positive), stable = bg (negative)
    run_direction "${site}_unstableFg" "$unstable_fa" "$stable_fa"

    # Direction B: stable = fg (positive), unstable = bg (negative)
    run_direction "${site}_stableFg" "$stable_fa" "$unstable_fa"
done

echo
echo "All runs complete. Reports are at:"
echo "  $OUT_DIR/<tag>_gt_out/report.rnaprot_gt.html    (dataset-level k-mer/structure comparison)"
echo "  $OUT_DIR/<tag>_eval_out/report.rnaprot_eval.html (model + motif logos)"
