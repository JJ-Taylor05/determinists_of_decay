#!/bin/bash
set -e
mkdir -p chunks
> manifest.txt
for N in 1 2 3 4 5; do
    grep "^>" human_halflife_data_sorted_part${N}.fa | sed 's/^>//' > chunks/part${N}_ids.txt
    split -l 200 -d -a 3 chunks/part${N}_ids.txt chunks/part${N}_chunk_
    for f in chunks/part${N}_chunk_*; do
        echo "${N} ${f}" >> manifest.txt
    done
done
echo "Total chunks created:"
wc -l manifest.txt
