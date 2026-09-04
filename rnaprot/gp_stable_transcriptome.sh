#!/bin/bash
#SBATCH --job-name=gp_stable_transcriptome
#SBATCH --account=tayem348
#SBATCH --partition=aoraki
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --array=1-5
#SBATCH --output=gp_stable_transcriptome_part%a.log

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaprot

rnaprot gp \
    --in human_halflife_data_sorted_part${SLURM_ARRAY_TASK_ID}.fa \
    --train-in stable_train_train_out/ \
    --out transcriptome_stable_gp_out_part${SLURM_ARRAY_TASK_ID}  \
    --report
