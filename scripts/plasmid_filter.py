#Standard libraries
from pathlib import Path
import argparse
import json

#Local imports form utils.py to reuse declared functions
from utils import build_stats,write_jsonl,ensure_dir,validate_dir

# Filter plasmid-associated sequences '(plasmid)'
def filter_plasmid(input_fasta: Path, output_fasta: Path):
    keep = True
    total = 0
    removed = 0
    kept = 0   
    # Opens the compressed file in read mode and the output file in write mode
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for line in infile: # Read FASTA file line-by-line
            # FASTA header line
            if line.startswith(">"):
                total += 1
                # Detect plasmid-associated sequences when the header contains the term "(plasmid)"
                if "(plasmid)" in line.lower():
                    # Exclude plasmid sequences from the output file 
                    keep = False
                    removed += 1
                else:
                    # Retain non-plasmid sequence in the output file
                    keep = True
                    kept += 1 
            if keep:
                outfile.write(line)
    # Builds the statistics for the filtering process 
    return build_stats(input_fasta.stem, total, kept)


# Filter plasmid sequences from all FASTA files
def filter_directory(input_dir: Path, output_dir: Path, plasmid_stats: Path):
    ensure_dir(output_dir)
    ensure_dir(plasmid_stats.parent)
    # Opens in writing mode the statistics file in JSONL format
    with open(plasmid_stats, "w") as stats_out:
        for faa_path in input_dir.glob("*.faa"):
            output_path = output_dir / faa_path.name
            stat = filter_plasmid(faa_path, output_path)
            # Saves metrics
            write_jsonl(stats_out, stat)


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser( description="Filter plasmid sequences from .faa files" )
    parser.add_argument("--input_dir", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    parser.add_argument("--plasmid_filtered", type=Path, required=True)
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    validate_dir(args.input_dir)
    # Run plasmid filtering pipeline
    filter_directory(args.input_dir, args.output_dir, args.plasmid_filtered)

if __name__ == "__main__":
    main()








