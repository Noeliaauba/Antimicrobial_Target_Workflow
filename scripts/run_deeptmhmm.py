# Standard libraries
import subprocess
from pathlib import Path
import argparse
import sys
import shutil

# Local imports form utils.py to reuse declared functions
from utils import ensure_dir,validate_dir,already_exists,ensure_dir, run_command
CHUNK_SIZE = 150


# Splits a FASTA file into 150-size chunks 
def split_fasta(fasta_path, chunk_size, out_dir):
    ensure_dir(out_dir)
    # Counters for chunking
    chunks = []
    current_chunk = []
    count = 0
    chunk_idx = 0

    with open(fasta_path) as f:
        for line in f:
            # Detects header lines of the file
            if line.startswith(">"):
                # Save the current chunk when the maximum number of proteins is reached
                if count == chunk_size:
                    chunk_file = out_dir / f"chunk_{chunk_idx}.fasta"
                    # Avoids re-writing existing chunk files
                    if not chunk_file.exists():
                        chunk_file.write_text("".join(current_chunk))
                    chunks.append(chunk_file)
                    current_chunk = []
                    count = 0
                    chunk_idx += 1
                count += 1
            current_chunk.append(line)
        # Saves the last chunk if there is any left sequence
        if current_chunk:
            chunk_file = out_dir / f"chunk_{chunk_idx}.fasta"
            if not chunk_file.exists():
                chunk_file.write_text("".join(current_chunk))    
            chunks.append(chunk_file)
    return chunks

# Run DeepTMHMM prediction on a single FASTA chunked file
def run_chunk(chunk_file, topology_dir, deeptmhmm_path):
    chunk_output_dir = topology_dir / chunk_file.stem
    results_file = chunk_output_dir / "predicted_topologies.3line"
    # Remove incomplete DeepTMHMM outputs (in case of previous failed executions)
    if chunk_output_dir.exists() and not results_file.exists():
        print(f"[CLEAN] Removing incomplete output: {chunk_output_dir}")
        shutil.rmtree(chunk_output_dir)
    # Avoids re-execution
    if results_file.exists():
        print(f"[SKIP] {chunk_file.name} already processed")
        return
    # Build the command to run DeepTMHMM 
    cmd = [sys.executable, deeptmhmm_path.name, "--fasta", str(chunk_file.resolve()),"--output-dir", str(chunk_output_dir.resolve())]
    # Execute DeepTMHMM from its installation directory
    run_command(cmd, cwd=deeptmhmm_path.parent)
    print(f"Finished {chunk_file.name}\n")


# Merge each DeepTMHMM chunk predictions from the same bacteria file
def merge_outputs(output_dir, merged_file):
    with open(merged_file, "w") as outfile:
        for subdir in sorted(output_dir.iterdir()):
            if not subdir.is_dir():
                continue
            # Use each predicted_topologies.3line file to build the merged output    
            results_file = subdir / "predicted_topologies.3line"
            if results_file.exists():
                outfile.write(results_file.read_text())
                outfile.write("\n")

# Coordinate chunking, topology prediction through DeepTMHMM and merging of files
def process_genome(fasta_file,chunk_dir, topology_dir,deeptmhmm_path):
    bacteria_name = fasta_file.stem
    genome_chunk_dir = chunk_dir / bacteria_name
    genome_topology_dir = topology_dir / bacteria_name
    ensure_dir(genome_chunk_dir)
    ensure_dir(genome_topology_dir)
    merged_file = genome_topology_dir / f"{bacteria_name}_merged.3line"
    # Skip genomes that already have the merged prediction 
    if already_exists(merged_file, "Merged topology"):
        return
    # Avoids the re-execution of the prediction topology chunks
    existing_chunks = sorted(genome_chunk_dir.glob("chunk_*.fasta"),key=lambda x: int(x.stem.split("_")[1]))
    if existing_chunks:
        print(f"[SKIP] Using existing chunks for {bacteria_name}")
        chunks = existing_chunks
    else:
        print("Splitting FASTA...")
        chunks = split_fasta(fasta_file, CHUNK_SIZE, genome_chunk_dir)
        print(f"Created {len(chunks)} chunks\n")
    # Runs DeepTMHMM on each chunk and merge outputs
    for i, chunk in enumerate(chunks, start=1):
        print(f"[{bacteria_name}] "f"Chunk {i}/{len(chunks)} -> {chunk.name}")
        try:
            run_chunk(chunk, genome_topology_dir, deeptmhmm_path)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Failed on {chunk.name}")
            print(e)
    # Merges the .3line file of each chunk into a whole genome file
    merge_outputs(genome_topology_dir, merged_file)
    print(f"[DONE] {merged_file}\n")


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser( description="Run DeepTMHMM topology prediction on chunked FASTA genomes.")
    parser.add_argument("--input_dir", required=True,type=Path) #The plasmid_filter
    parser.add_argument("--chunks_dir",required=True,type=Path)
    parser.add_argument("--output_dir",required=True,type=Path) #The execution of DeepTMHMM predicted_topologies
    parser.add_argument("--deeptmhmm",required=True,type=Path) #Path of predict.py
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    validate_dir(args.input_dir)
    fasta_files = sorted(args.input_dir.glob("*.faa"))
    print(f"Found {len(fasta_files)} genomes\n")
    #Process each genome file through chunking, topology prediction and merging
    for fasta_file in fasta_files:
        try:
            process_genome(fasta_file, args.chunks_dir, args.output_dir, args.deeptmhmm)
        except Exception as e:
            print(f"[FATAL ERROR] {fasta_file.name}")
            print(e)


if __name__ == "__main__":
    main()