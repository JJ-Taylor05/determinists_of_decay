#!/bin/bash
#SBATCH --job-name=streme_null
#SBATCH --partition=aoraki_short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-100%25
#SBATCH --output=logs/null_%A_%a.out
#SBATCH --error=logs/null_%A_%a.err

PROJECT_DIR="${HOME}/null_dist/streme_null_pipeline"
DATA_DIR="${PROJECT_DIR}/data"
CONDA_ENV="meme_suite"

set -euo pipefail

source ~/miniforge3/bin/activate
conda activate "${CONDA_ENV}"

SEED="${SLURM_ARRAY_TASK_ID}"
OUTDIR="${PROJECT_DIR}/null_runs/iter_${SEED}"
mkdir -p "${OUTDIR}"   # per-iteration output directory - script-made

echo "[$(date)] Iteration ${SEED} starting on $(hostname)"

# --- 1. Dinucleotide-shuffle, same seed for both ---
fasta-shuffle-letters -kmer 2 -rna -seed "${SEED}" \
    "${DATA_DIR}/stable_train.fa"   "${OUTDIR}/shuf_stable.fa"

fasta-shuffle-letters -kmer 2 -rna -seed "${SEED}" \
    "${DATA_DIR}/unstable_train.fa" "${OUTDIR}/shuf_unstable.fa"

# --- 2. Discriminative STREME: shuffled foreground vs unshuffled background ---
# stable-direction null: shuffled stable (p) vs real unstable (n)
streme --rna --evalue \
    -p "${OUTDIR}/shuf_stable.fa" \
    -n "${DATA_DIR}/unstable_train.fa" \
    -o "${OUTDIR}/streme_stable_null"

# unstable-direction null: shuffled unstable (p) vs real stable (n)
streme --rna --evalue \
    -p "${OUTDIR}/shuf_unstable.fa" \
    -n "${DATA_DIR}/stable_train.fa" \
    -o "${OUTDIR}/streme_unstable_null"

# --- 3. Tomtom: your REAL motifs (query) vs THIS iteration's null motifs (target) ---
if [[ -f "${OUTDIR}/streme_stable_null/streme.txt" ]]; then
    tomtom -oc "${OUTDIR}/tomtom_stable" \
        -evalue -thresh 0.5 \
        "${DATA_DIR}/streme_stable.txt" \
        "${OUTDIR}/streme_stable_null/streme.txt" || \
        echo "  (no null stable motifs to compare - streme.txt empty or absent)"
fi

if [[ -f "${OUTDIR}/streme_unstable_null/streme.txt" ]]; then
    tomtom -oc "${OUTDIR}/tomtom_unstable" \
        -evalue -thresh 0.5 \
        "${DATA_DIR}/streme_unstable.txt" \
        "${OUTDIR}/streme_unstable_null/streme.txt" || \
        echo "  (no null unstable motifs to compare - streme.txt empty or absent)"
fi

# --- 4. Drop the shuffled fastas to save disk ---
rm -f "${OUTDIR}/shuf_stable.fa" "${OUTDIR}/shuf_unstable.fa"

echo "[$(date)] Iteration ${SEED} done"
