#Standard libraries
import argparse
import subprocess
from pathlib import Path

#Local imports form utils.py to reuse declared functions
from utils import ensure_dir,validate_dir,find_fastas,already_exists,required_paths_exist, run_command

# Searches for FASTA files inside a directory and returns the first candidate 
def resolve_db_fasta(path: Path) -> Path:
    fastas = find_fastas(path)
     # Warn if multiple FASTA candidates exist
    if len(fastas) > 1:
        print(f"[WARN] Multiple candidates found, using first: {fastas[0]}")
    return fastas[0]


# Creates BLAST proteins database from a FASTA file
def run_makeblastdb(db_fasta: Path, db_dir: Path) -> Path:
    ensure_dir(db_dir)
    db_prefix = db_dir / db_fasta.stem
    # Required BLAST database index files to check if db already exists
    required_files= [db_prefix.with_suffix(".pin"), db_prefix.with_suffix(".psq"), db_prefix.with_suffix(".phr"),]
    if required_paths_exist(required_files):
        print(f"[SKIP] BLAST DB already exists: {db_prefix}")
        return db_prefix
    # Builds the command to run makeblastdb
    cmd = ["makeblastdb","-in", str(db_fasta),"-dbtype", "prot","-out", str(db_prefix),"-parse_seqids"]  # Run the command to create db
    run_command(cmd)
    return db_prefix


# Runs BLASTp: bacteria_fasta (query) vs db_fasta (db)
def run_blastp(query_fasta: Path, db_prefix: Path, output_tsv: Path, evalue: float, threads: int):
    ensure_dir(output_tsv.parent)
    if already_exists(output_tsv, "BLAST output"):
        return
    # Builds the command to run blastp
    cmd = ["blastp","-query", str(query_fasta),"-db", str(db_prefix),"-evalue", str(evalue),"-outfmt", "6","-num_threads", str(threads),"-out", str(output_tsv)]
    run_command(cmd) #before query_fasta.parent.name

# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Run BLASTp for multiple query files against a  database")
    parser.add_argument("--query_dir", type=Path, required=True)
    parser.add_argument("--db_fasta", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    parser.add_argument("--db_dir", type=Path, required=True)
    parser.add_argument("--evalue", type=float, default=0.001)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    validate_dir(args.query_dir)
    # Validates the existence of the database FASTA file
    db_fasta = resolve_db_fasta(args.db_fasta)
    if not db_fasta.exists():
        raise FileNotFoundError(f"DB FASTA not found: {db_fasta}")
    # Search query FASTA files from a directory
    fastas = list(args.query_dir.rglob("*.faa"))
    if not fastas:
        raise FileNotFoundError(f"No .faa files found in: {args.query_dir}")
    # Creates BLAST database from the provided FASTA file 
    db_prefix = run_makeblastdb(db_fasta, args.db_dir)
    for fasta in fastas:
        base = fasta.stem # Uses the filename as base for output
        output_tsv = args.output_dir / f"{base}.tsv"
        run_blastp(query_fasta=fasta, db_prefix=db_prefix, output_tsv=output_tsv, evalue=args.evalue, threads=args.threads)
if __name__ == "__main__":
    main()