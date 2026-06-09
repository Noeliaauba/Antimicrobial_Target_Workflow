#Standard libraries
import argparse
from pathlib import Path
from Bio import SeqIO

#Local imports form utils.py to reuse declared functions
from utils import ensure_dir, validate_dir,find_fastas

# Get all FASTA files in the input directory
def get_fastas(input_dir: Path):
    validate_dir(input_dir)
    fastas = find_fastas(input_dir)
    return fastas


# Generate unique FASTA identifier genomeXXX|proteinXXX
def make_unique_id(file_stem: str, record_id: str):
    return f"{file_stem}|{record_id}"


# All multiple FASTA files are merged into a single file 
def merge_fastas(input_dir: Path, output_dir: Path):
    ensure_dir(output_dir)
    output_fasta = output_dir / "merged.faa"
    total = 0
    # Open merged FASTA file to write the combined sequences
    with open(output_fasta, "w") as out_f:
        for fasta in get_fastas(input_dir):
            prefix = fasta.stem
            print(f"[INFO] Processing {fasta.name}")
            for rec in SeqIO.parse(fasta, "fasta"):
                rec.id = make_unique_id(prefix, rec.id)  # Create unique sequence ID
                rec.description = rec.id # Update description to match the new ID
                SeqIO.write(rec, out_f, "fasta")
                total += 1
    print(f"[INFO] Total sequences merged: {total}")
    print(f"[INFO] Output: {output_fasta}")


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Merge FASTA files from a directory into a single FASTA")
    parser.add_argument("--input_dir", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    merge_fastas(args.input_dir, args.output_dir)
    
if __name__ == "__main__":
    main()