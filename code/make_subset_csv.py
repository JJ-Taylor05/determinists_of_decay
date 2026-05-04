import pandas as pd
import mygene
import argparse
import sys

def process_files(fasta_path, csv_path, output_path):
    print(f"1. Extracting Transcript IDs from: {fasta_path}")
    transcript_ids = []
    
    try:
        with open(fasta_path, 'r') as file:
            for line in file:
                # Look for header lines
                if line.startswith('>'):
                    # Extract the ID and remove the version number
                    raw_id = line.strip().split()[0][1:] 
                    base_id = raw_id.split('.')[0]
                    transcript_ids.append(base_id)
    except FileNotFoundError:
        print(f"Error: Could not find the FASTA file at '{fasta_path}'")
        sys.exit(1)
                
    transcript_ids = list(set(transcript_ids))
    print(f"   Found {len(transcript_ids)} unique transcripts.")

    print("\n2. Mapping Transcript IDs (ENST) to Gene IDs (ENSG) via MyGene.info API...")
    mg = mygene.MyGeneInfo()
    mapping_results = mg.querymany(
        transcript_ids, 
        scopes='ensembl.transcript', 
        fields='ensembl.gene', 
        species='human',
        verbose=False
    )

    valid_gene_ids = set()
    for item in mapping_results:
        if 'ensembl' in item:
            ensembl_data = item['ensembl']
            if isinstance(ensembl_data, list):
                for entry in ensembl_data:
                    if 'gene' in entry:
                        valid_gene_ids.add(entry['gene'])
            elif isinstance(ensembl_data, dict):
                if 'gene' in ensembl_data:
                    valid_gene_ids.add(ensembl_data['gene'])

    print(f"   Successfully mapped to {len(valid_gene_ids)} unique Gene IDs.")

    print(f"\n3. Filtering the CSV data from: {csv_path}")
    try:
        df = pd.read_csv(csv_path, header=None, names=['Gene_ID', 'Gene_Symbol', 'Half_Life'])
    except FileNotFoundError:
        print(f"Error: Could not find the CSV file at '{csv_path}'")
        sys.exit(1)
    
    filtered_df = df[df['Gene_ID'].isin(valid_gene_ids)]
    print(f"   Found {len(filtered_df)} matching records in the CSV.")

    print("\n4. Saving the new CSV...")
    filtered_df.to_csv(output_path, index=False, header=False)
    print(f"   Done! Saved to: {output_path}")

if __name__ == "__main__":
    # Set up the command-line argument parser
    parser = argparse.ArgumentParser(description="Filter a CSV file based on Transcript IDs from a FASTA file.")
    
    # Define the arguments
    parser.add_argument("-f", "--fasta", required=True, help="Path to the input FASTA file (.fa)")
    parser.add_argument("-c", "--csv", required=True, help="Path to the input CSV file (.csv)")
    parser.add_argument("-o", "--output", required=True, help="Path for the filtered output CSV file")
    
    # Parse the arguments provided by the user
    args = parser.parse_args()
    
    # Run the main function with the provided arguments
    process_files(args.fasta, args.csv, args.output)