#Standard libraries  
from pathlib import Path
import requests
import argparse
import time

#Local imports form utils.py to reuse declared functions
from utils import ensure_dir, download_file

# Download DEG database file from a remote URL with 3 retry attempts if download fails
def download_deg(url: str, output_dir: Path, retries=3):
    ensure_dir(output_dir)
    # Extract filename from URL
    filename = url.split("/")[-1] or "deg_download.zip"
    output_path = output_dir / filename
    # Avoid re-downloading if file already exists
    if output_path.exists():
        print(f"Using existing file: {output_path}")
        return output_path
    print(f"[REQUEST] {url}")
    # Download the file into the desired directory
    download_file(url, output_path)
    return output_path


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Download DEG database")
    parser.add_argument("--url", required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    # Calls the function to download DEG file
    download_deg(args.url, args.output_dir)
   
if __name__ == "__main__":
    main()
