#Standard libraries  
import argparse
import time                
import requests            
from Bio import Entrez     # module from Biopython to communicate with NCBI’s system (GenBank, Genome, Nucleotide, Taxonomy, Protein…)
from pathlib import Path   

#Local imports form utils.py to reuse declared functions
from utils import ensure_dir, download_file, already_exists 


# Recives a terminology (Accession number) and searches for the corresponding assembly in NCBI (UID) 
def search_assembly(terminology): 
    # sends a request for search in assembly database (eSearch) and get the first result 
    handle = Entrez.esearch(db="assembly", term=terminology, retmax=1) 
    record = Entrez.read(handle) 
    handle.close() # close the connection for liberate socket, reading buffers, memory...
    # If the list is not empty returns the first (UID). Otherwise, returns None
    if record["IdList"]:
        return record["IdList"][0]
    return None


# Recives an UID (Internal Assembly ID from NCBI) and gets the corresponding FASTA information
def get_fasta(uid):
    handle = Entrez.esummary(db="assembly", id=uid, report="full") # request the complete metadata report 
    summary = Entrez.read(handle)
    handle.close()
    # Object docsum is a dictionary that contains:[Organism name, Accession number, FTP paths, Assembly version, # Metadata fields]
    docsum = summary['DocumentSummarySet']['DocumentSummary'][0] 
    ftp_path = docsum.get('FtpPath_RefSeq') # Gets ONLY the RefSeq (reviewed, standarized and curated version by NCBI)   
    organism = docsum.get('Organism')
    assembly_accession = docsum.get('AssemblyAccession')
    # If the RefSeq path doesn't exist, skip it
    if not ftp_path:
        print(f"Skipping {organism} ({assembly_accession}) — no RefSeq FTP link found.")
        return None, None, None

    # If RefsEQ exists, build URL for the FASTA file by {ftp_path}/{base_name}_protein.faa.gz
    base_name = Path(ftp_path).name                 
    fasta_url = f"{ftp_path}/{base_name}_protein.faa.gz"  
    # Convert FTP to HTTPS for requests
    fasta_url = fasta_url.replace("ftp://", "https://")
    return fasta_url, organism, assembly_accession


# Download the oragnism FASTA file by the Accession Number 
def download_by_accession(accession, bacteria_dir, human_dir):
    print(f"Searching by accession: {accession}")
    # Search UID by the Accession Number
    UID = search_assembly(accession)
    if not UID:
        print(f"NO ASSEMBLY FOUND FOR ACCESSION {accession}\n")
        return False

    # Gets FASTA metadata: URL, organism name, assembly name
    fasta_url, organism, assembly = get_fasta(UID)
    if not fasta_url:
        print(f"NO FASTA LINK FOUND FOR {accession}\n")
        return False
    # If the organism is a human save it in the human_dir. Otherwise, save it in the bacteria_dir
    organism_lower = organism.lower()
    if "homo sapiens" in organism_lower or "human" in organism_lower:
        out_dir = human_dir
    else:
        out_dir = bacteria_dir

    out_path = out_dir / f"{accession}_protein.faa.gz"
    # Avoid re-downloading if file already exists
    if already_exists(out_path, "Downloaded FASTA"):
        return True
    # Download the FASTA file
    try:
        download_file(fasta_url, out_path)
        print(f"Downloaded: {out_path}\n")
        return True
    except Exception as e:
        print(f"Could not download by accession {accession}: {e}\n")
        return False


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Download protein FASTA files from NCBI assemblies")
    parser.add_argument("email")
    parser.add_argument("accessions")
    parser.add_argument("--bacteria_dir", required=True, type=Path)
    parser.add_argument("--human_dir", required=True, type=Path)
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    # Creates the folders
    for d in [args.bacteria_dir, args.human_dir]:
        ensure_dir(d)
    
    # Configure Entrez
    Entrez.email = args.email
    Entrez.tool = "SmartGenomeDownloader"
    # Process each accession number provided in the comma-separated list
    accessions = [acc.strip() for acc in args.accessions.split(",") if acc.strip()]
    for acc in accessions:   
        download_by_accession(acc, args.bacteria_dir, args.human_dir)
        time.sleep(2)  # Respect NCBI requests limits

if __name__ == "__main__":
    main()










