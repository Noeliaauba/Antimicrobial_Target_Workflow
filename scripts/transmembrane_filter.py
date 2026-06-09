#Standard libraries
from pathlib import Path
import argparse
import json

#Local imports form utils.py to reuse declared functions
from utils import build_stats, write_jsonl, ensure_dir, validate_dir

# Extract transmembrane (TM) protein IDs 
def extract_tm_ids(predicted_file: Path):
    tm_ids = set()
    with predicted_file.open() as infile:
        for line in infile:
            if line.startswith(">"): # >protein_id | TM
                # Split the line into identifier and category (TM, BETA, GLOB...)
                left, right = line.strip().split("|", 1)
                protein_id = left.lstrip(">").strip()
                tag = right.strip().upper()  # Normalize the category tag
                # Store transmembrane proteins only
                if tag == "TM":
                    tm_ids.add(protein_id)
    return tm_ids


# Filters FASTA enries whose ID is in tm_ids
def filter_by_id(faa_file: Path, output_file: Path, valid_ids):
    kept = 0
    total = 0
    removed = 0
    # Opens in writing mode the statistics file in JSONL format
    with faa_file.open() as infile, output_file.open("w") as outfile:
        write = False
        # Check line by line if header matches tm_ids
        for line in infile:
            if line.startswith(">"):
                total += 1
                protein_id = line.split()[0].lstrip(">")
                write = protein_id in valid_ids
                if write:
                    kept += 1
                else:
                    removed += 1   
            # Write retained sequences              
            if write:
                outfile.write(line)
    # Builds the statistics for the filtering process 
    return build_stats(faa_file.stem, total, kept)
    
    
# Process all plasmid FASTA files using TM topology prediction
def process_paths(plasmid_dir: Path, topology_dir: Path, output_dir: Path, transmembrane_filtered: Path):
    ensure_dir(output_dir)
    ensure_dir(transmembrane_filtered.parent)
    # Opens in writing mode the statistics file in JSONL format
    with open(transmembrane_filtered, "w") as stats_out:
        plasmid_files = list(plasmid_dir.glob("*.faa"))
        for p_file in plasmid_files:
            sample_name = p_file.stem
            sample_topology_dir = topology_dir / sample_name
            # Warns if no topology directory exists for the sample and skip to the next one
            if not sample_topology_dir.exists():
                print(f"[WARNING] No topology folder for {sample_name}")
                continue
            topology_files = list(sample_topology_dir.glob("*.3line"))
            # Warns about missing topology prediction files
            if not topology_files:
                print(f"[WARNING] No .3line files for {sample_name}")
                continue
            # Extracts all transmembrane IDs
            tm_ids = set()
            for t_file in topology_files:
                tm_ids.update(extract_tm_ids(t_file))
            output_file = output_dir / p_file.name
            stat = filter_by_id(p_file, output_file, tm_ids)
            # Saves metrics
            write_jsonl(stats_out, stat)


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plasmid_dir", type=Path, required=True)
    parser.add_argument("--topology_dir", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    parser.add_argument("--transmembrane_filtered", type=Path, required=True)
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    validate_dir(args.plasmid_dir)
    validate_dir(args.topology_dir)
    # Runs transmembrane filtering pipeline
    process_paths(args.plasmid_dir, args.topology_dir, args.output_dir, args.transmembrane_filtered)

if __name__ == "__main__":
    main()











