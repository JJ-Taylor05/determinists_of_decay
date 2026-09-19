#!/bin/bash
#SBATCH --job-name=streme_null_analyze
#SBATCH --partition=aoraki_short
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/analyze_%j.out
#SBATCH --error=logs/analyze_%j.err

PROJECT_DIR="${HOME}/null_dist/streme_null_pipeline"
DATA_DIR="${PROJECT_DIR}/data"
CONDA_ENV="meme_suite"
N_ITER=100

set -euo pipefail

source ~/miniforge3/bin/activate
conda activate "${CONDA_ENV}"

python analyze_null.py \
    --project-dir "${PROJECT_DIR}" \
    --data-dir "${DATA_DIR}" \
    --n-iter "${N_ITER}"
