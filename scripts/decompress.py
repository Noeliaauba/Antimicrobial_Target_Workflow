#Standard libraries  
import gzip
import shutil
import argparse
import zipfile
import tarfile
from pathlib import Path

#Local imports form utils.py to reuse declared functions
from utils import ensure_dir, validate_dir


# Check if the file is a FASTA based on its extension
def is_fasta(file_path: Path):
    return file_path.suffix in (".fa", ".faa", ".fasta")


# Count the number of sequences by counting header lines
def count_sequences(fasta_file: Path):
    count = 0
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


#Decompress a single file with .targ.gz, .tar, .gz,.zip extensions
def decompress(file_path: Path, out_dir: Path, force=False):
    ensure_dir(out_dir)
    outputs = []
    suffixes = file_path.suffixes
    # handles .tar.gz extension
    if len(suffixes) >= 2 and suffixes[-2:] == [".tar", ".gz"]:
        print(f"[TAR.GZ] Extracting {file_path.name}")
        # Open the archive in read mode
        with tarfile.open(file_path, "r:gz") as tar:
            members = tar.getmembers()
            #Extract all contents
            tar.extractall(out_dir)
            for m in members:
                # Store the decompressed file and verify it is FASTA format
                extracted = out_dir / m.name
                outputs.append(extracted)
                if is_fasta(extracted):
                    print(f"{count_sequences(extracted)} sequences")
    
    # handles .tar extension
    elif file_path.suffix == ".tar":
        print(f"[TAR] Extracting {file_path.name}")
        with tarfile.open(file_path, "r") as tar:
            members = tar.getmembers()
            tar.extractall(out_dir)
            for m in members:
                extracted = out_dir / m.name
                outputs.append(extracted)
                if is_fasta(extracted):
                    print(f"{count_sequences(extracted)} sequences")
    
    # handles .gz extension
    elif file_path.suffix == ".gz":
        #Removes .gz from the filename
        out_file = out_dir / file_path.with_suffix("").name
        # Skip decompression if output file already exists
        if out_file.exists() and not force:
            print(f"[SKIP] {out_file.name}")
            return [out_file]
        
        print(f"[GZ] Extracting {file_path.name}")
        # Opens the compressed file in read mode and the output file in write mode
        with gzip.open(file_path, "rb") as f_in, open(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out) # copy on stream

        outputs.append(out_file)
        if is_fasta(out_file):
            print(f"{count_sequences(out_file)} sequences")

    # handles .zip extension
    elif file_path.suffix == ".zip":
        print(f"[ZIP] Extracting {file_path.name}")
        # Opens the compressed file in read mode
        with zipfile.ZipFile(file_path, "r") as z:
            #Extract all contents and verify it is FASTA format
            names = z.namelist()
            z.extractall(out_dir)
            for name in names:
                extracted = out_dir / name
                outputs.append(extracted)
                if is_fasta(extracted):
                    print(f"{count_sequences(extracted)} sequences")
    else:
        print(f"Unsupported format: {file_path}")
    return outputs



# Decompress all files in a directory using the decompress function 
def decompress_directory(input_dir: Path, output_dir: Path) -> list[Path]:
    ensure_dir(output_dir)
    outputs = []
    for file_path in input_dir.iterdir():
        if file_path.is_file():
            outputs.extend(decompress(file_path, output_dir))
    return outputs


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser( description="Decompress all .gz files in a directory")
    parser.add_argument("--input_dir", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    # Validates the existence of the input file
    validate_dir(args.input_dir)
    outputs = decompress_directory(args.input_dir, args.output_dir)
    print(f"Decompressed {len(outputs)} files")

if __name__ == "__main__":
    main()