#!/bin/bash
#SBATCH --job-name=gp_stable_transcriptome
#SBATCH --account=tayem348
#SBATCH --partition=aoraki
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --output=gp_stable_transcriptome.log

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaprot

rnaprot gp \
    --in all_human_halflife_data.fa \
    --train-in unstable_train_train_out/ \
    --out transcriptome_unstable_gp_out \
    --report
