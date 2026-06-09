#Standard libraries
import argparse
import subprocess
from pathlib import Path

#Local imports form utils.py to reuse declared functions
from utils import ensure_dir,validate_dir, required_paths_exist, run_command

# Verifies that required eggNOG database files exist
def eggnog_db_ready(db_dir: Path):
    required = [
        db_dir / "eggnog.db", # Main annotation database
        db_dir / "eggnog_proteins.dmnd", # DIAMOND protein database 
        db_dir / "eggnog.taxa.db" # Taxonomic database 
    ]
    return required_paths_exist(required)


#Runs the internal script emapper.py to generate annotation tables (GO terms, functionality, Enzymatic activity and metabolic pathways)
def run_eggnog(input_fasta: Path, output_dir: Path, db_dir: Path, db_type: str, emapper_path: Path):
    ensure_dir(output_dir)
    output_prefix = output_dir / "eggnog_annotation"
    # Builds the command to run eggNOG mapper (using DIAMOND search mode)
    cmd = [
        str(emapper_path.resolve()),
        "-i", str(input_fasta),
        "--itype", "proteins",
        "-m", "diamond",
        "--data_dir", str(db_dir),
        "--output", "eggnog_annotation",
        "--output_dir", str(output_dir),
        "--cpu", "8",
        "--override",
        "--target_taxa", db_type,]
    run_command(cmd)


# Parse the input arguments   
def parse_args():
    parser = argparse.ArgumentParser(description="Run eggNOG mapper for functional annotation of protein sequences")
    parser.add_argument("--input_fasta", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    parser.add_argument("--db_dir",type=Path, required=True)
    parser.add_argument("--db_type", type=str, default="2") #2= bacteria (taxaID)
    parser.add_argument("--emapper",type=Path,required=True)                   
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    validate_dir(args.db_dir)
    #Verifies that required eggNOG database files exist before running eggNOG mapper
    if not eggnog_db_ready(args.db_dir):
        raise RuntimeError(f"[ERROR] eggNOG DB incomplete in {args.db_dir}")
    run_eggnog(args.input_fasta, args.output_dir, args.db_dir, args.db_type, args.emapper)

if __name__ == "__main__":
    main()