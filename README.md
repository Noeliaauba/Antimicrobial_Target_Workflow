# Computational workflow for the Identification of Antimicrobial Targets 
This project implements a computational workflow for the identification and prioritization of potential antimicrobial targets over complete bacterial genomes.
The workflow was developed under Ubuntu 22.04 (WSL) using Python 3.10 and Bash scripts.
The workflow uses three isolated  Python virtual environments (pipeline_env, deeptmhmm, eggnog_env)

The workflow combines:
- Bacterial proteome retrieval
- Transmembrane topology prediction using DeepTMHMM
- Essential and human-homology identification using BLASTp
- Sequences similarity clustering using CD-HIT
- Functional annotation using eggNOG-mapper
- Prioritization of potential proteins targets based on functional and structural criteria 


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

The complete pipeline can be executed through run_pipeline.sh
The workflow will ask the user to provide:
* **Email address:** (required by NCBI Entrez)
* **Genome accession numbers:** (comma-separated)
* **Minimum number of proteins per cluster**
* **Minimum number of organisms per cluster**
* **GO term filtering threshold**

## Tested Genomes
The workflow has been validated using the following bacterial proteome and human reference proteome (as it is required in the input) :
**NZ_CP015121.1,NC_002516.2,NZ_CP009756.1,NC_002695.2,NC_016845.1,NZ_CP034944.1,NC_010554.1,NZ_CP029736.1,NC_017731.1,NZ_CP063354.1,NZ_CP039729.1,NC_007795.1,NC_017379.1,NC_002163.1,NC_003197.2,NZ_CP012028.1,NZ_CP007593.1,NZ_CP009610.1,NC_004337.2,NZ_CP061527.1,NC_000962.3,GCF_000001405.40**


# Installation
## System Requirements
* Ubuntu 22.04 LTS (Native or WSL2)
* Python 3.10+

## System Dependencies
Run the following commands to install required system packages and bioinformatic utilities:
```bash
sudo apt update
sudo apt install unzip
```

Install the required external tools:
```bash
sudo apt install ncbi-blast+
sudo apt install cd-hit
sudo apt install diamond-aligner
```

Verify installation:
```bash
blastp -version
makeblastdb -version
cd-hit -h
diamond version
```

## Python Virtual Environment Setup
This project uses isolated Python virtual environments (venv), each environment configuration files is provided to ensure reproducibility across different systems.

### 1. Main Pipeline Environment
Create and activate the environment:
```bash
python3 -m venv envs/pipeline_env
source envs/pipeline_env/bin/activate
```
Install dependencies:
```bash
pip install -r pipelinebase_requirements.txt
deactivate
```

### 2. DeepTMHMM Environment
Create and activate the environment:
```bash
python3 -m venv envs/deeptmhmm
source envs/deeptmhmm/bin/activate
```
Install dependencies:
```bash
pip install -r deeptmhmm_requirements.txt
deactivate
```

### 3. eggNOG Environment
Create and activate the environment:
```bash
python3 -m venv envs/eggnog_env
source envs/eggnog_env/bin/activate
```
Install dependencies:
```bash
pip install -r eggnog_requirements.txt
emapper.py --version
deactivate
```
### Specific External Tools Installation
#### DeepTMHMM Installation
DeepTMHMM is not distributed with this repository and must be obtained from the official website, the retrieve is done with academic license. 
1. Complete the form to receive the DeepTMHMM package from [DTU Health Tech portal](https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/).
2. Place the compressed package in the following directory: project_root/external/
3. Then run:
```bash
   cd external
   unzip DeepTMHMM*.zip
   mv DeepTMHMM-Academic-License-v1.0 DeepTMHMM
   ```

4. The final structure should be:
   ```text
   external/
   └── DeepTMHMM/
       └── predict.py
   ```
       
5. DeepTMHMM relies on loading complete serialized objects, recent versions of PyTorch (≥ 2.6) introduced stricter rules for deserialization models, to ensure compatibility the source code must be modified to explicitly disable the weights_only restriction.
 **FILE TO EDIT** `external/DeepTMHMM/utils.py`
Replace:
```python
  torch.load('esm_model_args.pt')
  torch.load('esm_model_alphabet.pt')
  torch.load('esm_model_state_dict.pt')
  torch.load(path, map_location=lambda storage, loc: storage)

```
New lines:
```python
  torch.load('esm_model_args.pt', weights_only=False) 
  torch.load('esm_model_alphabet.pt', weights_only=False) 
  torch.load('esm_model_state_dict.pt', weights_only=False)
  torch.load(path, map_location=lambda storage, loc: storage, weights_only=False)
 ```

6. After applying these modifications, test the installation to validate that predictions run successfully:
```bash
source envs/deeptmhmm/bin/activate
python external/DeepTMHMM/predict.py --fasta test.fasta --output-dir test_results
deactivate
```

#### eggNOG Database Installation
The official download_eggnog_data.py utility was not used because the download endpoints referenced by eggNOG-mapper v2.1.13 were unavailable. Direct download links to eggNOG database version 5.0.2 are therefore provided to guarantee reproducibility of the reported results.
```bash
mkdir -p external/eggnog
cd external/eggnog
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
```

## Running the workflow
Launch the main automated script:
```bash
bash run_pipeline.sh
```

### Final treemap after the workflow execution
```text
.
├── envs/                               
│   ├── deeptmhmm/
│   ├── eggnog_env/
│   └── pipeline_env/
│
├── external/                           
│   ├── DeepTMHMM/
│   └── eggnog/
│       ├── eggnog.db.gz
│       ├── eggnog.taxa.tar.gz
│       └── eggnog_proteins.dmnd.gz
│
├── data/                               
│   ├── input/                         
│   │   ├── bacterial
│   │   └── human
│   │
│   ├── intermediate/                   
│   │   ├── decompressed/   
│   │   ├── plasmid_filtered/  
│   │   ├── transmembrane_proteins/  
│   │   ├── topology_prediction/  
│   │   ├── chunks/     
│   │   ├── blast_nonhomologous/             
│   │   ├── blast_essentiality/  
│   │   └── cdhit/            
│   │
│   └── results/                        
│       ├── eggnog_mapper/       
│       └── prioritized/                      
│          
├── analysis/                                
│   └── plot.py               
│ 
├── scripts/                               
│   ├── download_ncbi.py
│   ├── decompressed.py
│   ├── plasmid_filter.py
│   ├── run_deeptmhmm.py
│   ├── transmembrane_filter.py
│   ├── transmembrane_phobius.py
│   ├── blastp.py
│   ├── blastp_filter.py
│   ├── download_deg.py
│   ├── merged.py
│   ├── cd_hit.py
│   ├── cluster_filter.py
│   ├── eggNOG_mapper.py
│   ├── prioritization.py
│   └── utils.py
│
├── stats/ 
├── plots/  
│
├── deeptmhmm_requirements.txt
├── eggnog_requirements.txt
├── pipelinebase_requirements.txt
│
├── run_pipeline.sh  
├── config.sh   
└── README.md
```

