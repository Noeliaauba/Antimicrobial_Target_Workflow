#Standard libraries
import argparse
from pathlib import Path
from Bio import SeqIO
import json

#Local imports form utils.py to reuse declared functions
from utils import build_stats, write_jsonl,ensure_dir,validate_dir


# Extract sequences IDs that pass the filtering criteria (identity, e-value)
def get_ids(tsv_file, id_col=0, identity_col=2, evalue_col=10,min_identity=35.0, max_evalue=1e-3):
    selected_ids = set()
    with open(tsv_file) as f:
        for line in f:
            if not line.strip(): # Skip empty lines
                continue
            # Split the data and extract the relevant columns
            cols = line.rstrip("\n").split("\t")
            protein_id = cols[id_col]
            identity = float(cols[identity_col])
            evalue = float(cols[evalue_col])
            # Filtering threshold
            if identity >= min_identity and evalue <= max_evalue:
                selected_ids.add(protein_id)
    return selected_ids


# Filter FASTA sequences using the selected IDs
def filter_by_ids(input_fasta, output_fasta, selected_ids, mode):
    total = 0
    kept = 0
    removed = 0
    records_out = []
    for record in SeqIO.parse(str(input_fasta), "fasta"):
        total += 1
        record_id = record.id.split()[0] # Extract sequence ID without description
        #Exclude filtering criteria: keep sequences that do not match the selected IDs
        if mode == "exclude":
            if record_id not in selected_ids:
                kept += 1
                records_out.append(record)
            else:
                removed += 1
        #Include filtering criteria: keep sequences that match the selected IDs
        elif mode == "include":
            if record_id in selected_ids:
                kept += 1
                records_out.append(record)
            else:
                removed += 1
        # Raise error if the filtering mode is invalid
        else:
            raise ValueError("mode must be 'include' or 'exclude'")
    # Write filtered FASTA sequences
    SeqIO.write(records_out, str(output_fasta), "fasta")
    kept = len(records_out)
    # Builds the statistics for the filtering process 
    return build_stats(input_fasta.stem, total, kept)


# Process all FASTA files in a directory 
def process_directory(input_dir: Path, tsv_dir: Path, output_dir: Path, blastp_filtered: Path, mode: str):
    ensure_dir(output_dir)
    ensure_dir(blastp_filtered.parent)
    # Opens in writing mode the statistics file in JSONL format
    with open(blastp_filtered, "w") as stats_out:
        for fasta in sorted(input_dir.rglob("*.faa")):
            base = fasta.stem
            tsv_file = tsv_dir / f"{base}.tsv"
            # Computes the statistics for the missing BLAST results 
            if not tsv_file.exists():
                stat = {"genome": base,"total": 0,"kept": 0,"removed": 0,"percentage": 0, "status": "missing_tsv"}
                write_jsonl(stats_out, stat)
                continue
            ids = get_ids(tsv_file)
            output_fasta = output_dir / f"{base}.faa"
            stat = filter_by_ids(fasta, output_fasta, ids, mode)
            # Saves metrics
            write_jsonl(stats_out, stat)


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Filter FASTA sequences based on BLASTp results")
    parser.add_argument("--tsv_dir", type=Path, required=True)
    parser.add_argument("--input_dir", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    parser.add_argument("--stats_file", type=Path, required=True)
    parser.add_argument("--mode", required=True, choices=["include", "exclude"])
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    validate_dir(args.input_dir)
    validate_dir(args.tsv_dir)
    # Run FASTA filtering workflow
    process_directory(args.input_dir,args.tsv_dir,args.output_dir,args.stats_file, args.mode)

if __name__ == "__main__":
    main()
        

