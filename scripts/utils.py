# Standard library imports
import json
from pathlib import Path
import subprocess
import time
import requests

# File extensions used for protein FASTA files
FASTA_EXTENSIONS = ("*.faa", "*.fa", "*.fasta", "*.aa", "*.dat")

# Search for FASTA files recursively
def find_fastas(path: Path):
    # Single FASTA input
    if path.is_file():
        return [path]
    fastas = []
    for pattern in FASTA_EXTENSIONS:
        fastas.extend(path.rglob(pattern))
    #Get FASTA files for multiple processing
    fastas = sorted(fastas)
    if not fastas:
        raise FileNotFoundError(f"No FASTA files found in {path}")
    return fastas

# Create directory if missing
def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


# Ensure the path exists and if it is a directory.
def validate_dir(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    if not path.is_dir():
        raise NotADirectoryError(path)


# Check if a fileX already exists
def already_exists(path: Path, label="FILE"):
    if path.exists():
        print(f"[SKIP] {label} exists: {path}")
        return True
    return False


# Check if multiple files already exist
def required_paths_exist(paths) -> bool:
    return all(path.exists() for path in paths)


# Download of external files with a retry mechanism
def download_file(url: str, out_path: Path, retries=3):
    for attempt in range(retries):
        # Temporary partial download file
        tmp_path = out_path.with_suffix(".tmp")
        try:
            # Opens streaming connection for the download and writes it in chunks to avoid memory issues with large files
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()  # Raise exception for HTTP errors
            with open(tmp_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            # Rename the temporary file after successful download
            tmp_path.rename(out_path)
            return out_path
        # Catches exceptions and retries the download 
        except Exception as e:
            print(f"[Retry {attempt + 1}] Error: {e}")
            if tmp_path.exists():
                tmp_path.unlink()  # Waits before the next retry
            time.sleep(2)
    raise RuntimeError(f"Download failed: {url}")


# Execute the command
def run_command(cmd, cwd=None):
    print(f"[CMD] {' '.join(map(str, cmd))}")
    subprocess.run(cmd, check=True, cwd=cwd)


# Build filtering statistics file 
def build_stats(genome, total, kept):
    removed = total - kept
    percentage = (round((kept / total) * 100, 2) if total > 0 else 0)
    return {"genome": genome,"total": total,"kept": kept,"removed": removed,"percentage": percentage}


# Write a single JSON object as JSONL line.
def write_jsonl(handle, data):
    handle.write(json.dumps(data) + "\n")
