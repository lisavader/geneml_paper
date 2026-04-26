This directory contains all data and scripts necessary for reproducing the data presented in the geneML paper.  

## Repository structure
### Data
The 'data' folder contains all data needed to run the pipeline (aside from any publicly available data / databases). 
Because some of the steps take a long time to run, the intermediate data that are produced by these steps are available under 'intermediates'.  
By default, notebooks producing figures and tables take the intermediate data as input, and can be run independently from the other steps.

### Scripts
All bash scripts of the pipeline are found under 'scripts' and are divided into 5 steps:
1. model_training: Training of geneML CNNs
2. benchmarking: Benchmarking of geneML against AUGUSTUS, BRAKER3, Helixer and ANNEVO
3. alt_transcripts: Analysis of alternative transcripts
4. reannotation: Reannotation of training genomes and comparison against original annotations
5. figures: Notebooks for the creation of figures and tables  

A setup script is also included, which includes cloning the geneML repository and downloading of databases.  
All custom python scripts that the pipeline depends on are in 'scripts/python'.  
Docker files are included under 'docker_files' and the geneML model used for benchmarking is under 'models'.

## Dependencies
### External dependencies (install manually)
- Docker (tested with version 29.4.1)
- gffcompare v0.12.10
- gffread v0.12.8
- GNU parallel (tested with version 20210822)
- HMMER 3.3.2
- NCBI datasets 18.18.0

### Python dependencies
Tested with Python 3.12.3.  
All dependencies are in requirements.txt, install with:
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### Databases
- OrthoDB 12v2
- PFAM 38.2
  
Download with:
```bash
bash scripts/00_setup/00_setup.sh
```

### Hardware
The 'model_training' step requires a GPU with CUDA support.  
Models were trained on NVIDIA A100 80GB PCIe, CUDA 13.1.

## Setup
1. Clone this repository
2. Make sure all dependencies are satisfied
3. Run the setup script 'scripts/00_setup/00_01_setup.sh'

## How to run
You can run the whole pipeline with:
```bash
bash run_all.sh
```
It's also possible to run individual steps, for example:
```bash
bash run_all.sh figures
```
