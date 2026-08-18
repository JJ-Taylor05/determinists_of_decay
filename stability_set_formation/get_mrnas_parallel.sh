#!/bin/bash
INPUT="$1"
OUTPUT="${INPUT%.csv}.fa"

process_row() {
    ensembl_ID="$1"
    canonical=$(curl -sf "https://rest.ensembl.org/lookup/id/${ensembl_ID}?expand=1" \
        -H 'Content-type: application/json' | \
        grep -o '"canonical_transcript":"[^"]*"' | \
        grep -o 'ENST[^"]*')
    canonical_trimmed=$(echo "$canonical" | cut -d'.' -f1)

    if [ -z "$canonical_trimmed" ]; then
        echo "Warning: could not find canonical transcript for $ensembl_ID" >&2
    else
        result=$(curl -sf "https://rest.ensembl.org/sequence/id/${canonical_trimmed}?type=cdna" \
            -H 'Content-type:text/x-fasta')
        if [ -z "$result" ]; then
            echo "Warning: sequence not found for $canonical_trimmed" >&2
        else
            echo "$result"
        fi
    fi
}

export -f process_row

# Extract just the ensembl IDs, run 10 workers in parallel, collect output
cut -d',' -f1 "$INPUT" | \
    xargs -P 10 -I {} bash -c 'process_row "$@"' _ {} \
    >> "$OUTPUT"
