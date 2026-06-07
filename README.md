# Computational workflow for the Identification of Antimicrobial Targets 
This project implements a computational workflow for the identification and prioritization of potential antimicrobial targets over complete bacterial genomes.

The workflow combines:
- Bacterial proteome retrieval
- Transmembrane topology prediction using DeepTMHMM
- Essential and human-homology identification using BLASTp
- Sequences similarity clustering using CD-HIT
- Functional annotation using eggNOG-mapper
- Prioritization of potential proteins targets based on functional and structural criteria 

The workflow was developed under Ubuntu 22.04 (WSL) using Python 3.10 and Bash scripts.
The workflow uses three isolated  Python virtual environments (pipeline_env, deeptmhmm, eggnog_env)

## Workflow overview
1. Bacterial proteome retrieval from NCBI
2. Plasmid filtering
3. Transmembrane protein prediction using DeepTMHMM
4. Human non-homology analysis using BLASTp
5. Essentiality analysis based on DEG and using BLASTp
6. Protein clustering using CD-HIT
7. Functional annotation by eggNOG-mapper
8. Candidate prioritization
9. Visualization of results

# Tested Genomes
The workflow has been validated using the following bacterial proteome and human reference proteome (as it is required in the input) :
NZ_CP015121.1,NC_002516.2,NZ_CP009756.1,NC_002695.2,NC_016845.1,NZ_CP034944.1,NC_010554.1,NZ_CP029736.1,NC_017731.1,NZ_CP063354.1,NZ_CP039729.1,NC_007795.1,NC_017379.1,NC_002163.1,NC_003197.2,NZ_CP012028.1,NZ_CP007593.1,NZ_CP009610.1,NC_004337.2,NZ_CP061527.1,NC_000962.3,GCF_000001405.40


# Installation
## System Requirements
Ubuntu 22.04 LTS 
WSL2 or Linux environment
Python 3.10+ 

## System Dependencies
sudo apt update
sudo apt install unzip

Install the required external tools:
sudo apt install ncbi-blast+
sudo apt install cd-hit
sudo apt install diamond-aligner

Verify installation:
blastp -version
makeblastdb -version
cd-hit -h
diamond version

## Python Virtual Setup
This project uses isolated Python virtual environments (venv), each environment configuration files is provided to ensure reproducibility across different systems.
1. Main Environment 
Create and activate the environment:
python3 -m venv envs/pipeline_env
source envs/pipeline_env/bin/activate

Install dependencies:
pip install -r pipelinebase_requirements.txt

2. DeepTMHMM Environment
Create and activate the environment:
python3 -m venv envs/deeptmhmm
source envs/deeptmhmm/bin/activate

Install dependencies:
pip install -r deeptmhmm_requirements.txt

3. eggNOG Environment
Create and activate the environment:
python3 -m venv envs/eggnog_env
source envs/eggnog_env/bin/activate

Install dependencies:
pip install -r eggnog_requirements.txt
emapper.py --version

## DeepTMHMM Installation
DeepTMHMM is not distributed with this repository and must be obtained from the official website, the retrieve is done with academic license.
https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/ 
1. Complete the form to receive the DeepTMHMM package
2. Place the compressed package in the following directory: project_root/external/
3. Then run:
        cd external
        unzip DeepTMHMM*.zip
        mv DeepTMHMM-Academic-License-v1.0 DeepTMHMM

4. The final structure should be:
   external/
   └── DeepTMHMM/
       └── predict.py
5. DeepTMHMM relies on loading complete serialized objects, recent versions of PyTorch (≥ 2.6) introduced stricter rules for deserialization models, to ensure compatibility the source code must be modified to explicitly disable the weights_only restriction.
FILE TO EDIT: external/DeepTMHMM/utils.py
Replace:
torch.load('esm_model_args.pt')
torch.load('esm_model_alphabet.pt')
torch.load('esm_model_state_dict.pt')
torch.load(path, map_location=device)
with:
torch.load('esm_model_args.pt', weights_only=False) 
torch.load('esm_model_alphabet.pt', weights_only=False) 
torch.load('esm_model_state_dict.pt', weights_only=False)
torch.load(path, map_location=device, weights_only=False)
6. After applying these modifications, test the installation to validate that predictions run successfully:
python external/DeepTMHMM/predict.py \
    --fasta test.fasta \
    --output-dir test_results


## eggNOG Database Installation
The official download_eggnog_data.py utility was not used because the download endpoints referenced by eggNOG-mapper v2.1.13 were unavailable. Direct download links to eggNOG database version 5.0.2 are therefore provided to guarantee reproducibility of the reported results.
cd external/eggnog
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz


## Running the workflow
The complete pipeline can be executed through run_pipeline.sh
The workflow will ask the user to provide:
* Email address (required by NCBI Entrez)
* Genome accession numbers (comma-separated)
* Minimum number of proteins per cluster
* Minimum number of organisms per cluster
* GO term filtering threshold
  
bash run_pipeline.sh
