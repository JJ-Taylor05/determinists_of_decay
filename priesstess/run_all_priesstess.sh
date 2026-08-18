#!/bin/bash
# run_all_priesstess.sh
#
# Loops over every site/rel_<offset> folder produced by extract_regions.py
# and runs PRIESSTESS on the fg.txt/bg.txt pair inside it.
#
# Usage:
#   ./run_all_priesstess.sh /path/to/priesstess_input
#
# Requires: PRIESSTESS on your PATH (i.e. the PRIESSTESS repo's main
# directory has been added to your bash profile, per the repo's install
# instructions), plus its dependencies (RNAfold, STREME, sklearn, skopt).

set -euo pipefail

INPUT_DIR="${1:?Usage: $0 /path/to/priesstess_input}"

# Loop over every site (cds_start, cds_end, utr3_end)...
for site_dir in "$INPUT_DIR"/*/; do
    site=$(basename "$site_dir")

    # ...and every sliding sub-window within that site
    for win_dir in "$site_dir"rel_*/; do
        offset=$(basename "$win_dir")   # e.g. "rel_-100"
        fg="$win_dir/fg.txt"
        bg="$win_dir/bg.txt"

        run_name="${site}_${offset}"
        echo "=== Running PRIESSTESS: $run_name ==="

        # -o points PRIESSTESS at a directory to create its own
        # PRIESSTESS_output/ folder inside; we give each run its own
        # subfolder here so the 33 runs don't collide.
        mkdir -p "$INPUT_DIR/results/$run_name"

        PRIESSTESS \
            -fg "$fg" \
            -bg "$bg" \
            -o "$INPUT_DIR/results/$run_name"

        echo "=== Done: $run_name ==="
    done
done

echo "All 33 PRIESSTESS runs complete. Results under $INPUT_DIR/results/<site>_rel_<offset>/PRIESSTESS_output/"
