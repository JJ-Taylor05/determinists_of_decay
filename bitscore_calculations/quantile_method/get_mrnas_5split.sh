#!/bin/bash
INPUT="$1"
BASENAME="${INPUT%.csv}"

# Count lines and calculate chunk size
TOTAL=$(wc -l < "$INPUT")
CHUNK=$(( (TOTAL + 4) / 5 ))   # ceiling division: (n + 4) / 5

# Define the fetch function
process_row() {
    linenum="$1"
    ensembl_ID="$2"

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
            # Print line number, a tab, then the FASTA (with internal newlines replaced
            # temporarily by a placeholder so the whole entry is one line for sorting)
            echo -e "${linenum}\t${result//$'\n'/___NEWLINE___}"
        fi
    fi
}

export -f process_row

# Loop over 5 chunks

for i in 1 2 3 4 5; do
    START=$(( (i - 1) * CHUNK + 1 ))
    END=$(( i * CHUNK ))
    # Clamp END to TOTAL so the last chunk doesn't overshoot
    [ $END -gt $TOTAL ] && END=$TOTAL

    OUTFILE="${BASENAME}_part${i}.fa"
    echo "Processing chunk $i (lines $START–$END) → $OUTFILE"

    # adds line numbers; filters to only the lines in this chunk's range;
    # pass both the original line number and Ensembl ID to process_row
    nl -ba -nrz "$INPUT" | \
        sed -n "${START},${END}p" | \
        awk '{print $1, $2}' FS='\t|,' OFS='\t' | \
        xargs -P 10 -I {} bash -c '
            args=({});
            process_row "${args[0]}" "${args[1]}"
        ' | \
        sort -n | \
        awk -F'\t' '{
            # Remove the line number column, restore the newline placeholders
            gsub(/___NEWLINE___/, "\n", $2); printf "%s\n", $2
        }' > "$OUTFILE"
done
