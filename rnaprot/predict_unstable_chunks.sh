#!/bin/bash
#SBATCH --job-name=predict_unstable_chunk
#SBATCH --account=tayem348
#SBATCH --partition=aoraki_gpu_L40
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150GB
#SBATCH --time=02:00:00
#SBATCH --output=predict_unstable_chunk%a.log

source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaprot

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" manifest.txt)
PART=$(echo $LINE | cut -d' ' -f1)
CHUNKFILE=$(echo $LINE | cut -d' ' -f2)

rnaprot predict \
    --in transcriptome_unstable_gp_out_part${PART}/ \
    --train-in unstable_train_train_out/ \
    --out predict_unstable_chunk_${SLURM_ARRAY_TASK_ID} \
    --mode 2 \
    --thr 2 \
    --site-id $(cat ${CHUNKFILE})
