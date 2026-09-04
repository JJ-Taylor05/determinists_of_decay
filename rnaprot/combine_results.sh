#!/bin/bash
set -e
N=$(wc -l < manifest.txt)

head -1 predict_stable_chunk_1/peak_regions.tsv > transcriptome_stable_predict_combined.tsv
for i in $(seq 1 $N); do
    tail -n +2 predict_stable_chunk_${i}/peak_regions.tsv >> transcriptome_stable_predict_combined.tsv
done

head -1 predict_unstable_chunk_1/peak_regions.tsv > transcriptome_unstable_predict_combined.tsv
for i in $(seq 1 $N); do
    tail -n +2 predict_unstable_chunk_${i}/peak_regions.tsv >> transcriptome_unstable_predict_combined.tsv
done

echo "Done. Combined files:"
wc -l transcriptome_stable_predict_combined.tsv transcriptome_unstable_predict_combined.tsv
