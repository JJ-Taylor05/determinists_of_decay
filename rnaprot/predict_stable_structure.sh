#!/bin/bash
#SBATCH --job-name=predict_stable_transcriptome
#SBATCH --account=tayem348
#SBATCH --partition=aoraki_gpu_L40
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --time=48:00:00
#SBATCH --output=predict_stable_transcriptome.log

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaprot

module load CUDA

rnaprot predict \
    --in transcriptome_stable_gp_out/ \
    --train-in stable_train_train_out/ \
    --out transcriptome_stable_predict_out \
    --mode 2 \
    --thr 2
