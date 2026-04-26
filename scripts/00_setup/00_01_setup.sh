#!/bin/bash

# Clone geneML v1.1.0
git clone --branch v1.1.0 https://github.com/hexagonbio/geneML.git

# Download OrthoDB 12v2
DATA_PATH=https://data.orthodb.org/v12/download/odb_data_dump

if [ ! -f "databases/orthodb/odb12v2_og_aa_fasta" ]; then
    wget -P databases/orthodb \
    $DATA_PATH/odb12v2_species.tab.gz \
    $DATA_PATH/odb12v2_levels.tab.gz \
    $DATA_PATH/odb12v2_level2species.tab.gz \
    $DATA_PATH/odb12v2_og_aa_fasta.gz
    gunzip databases/orthodb/*.gz
fi

# Download and prepare PFAM v38.2
if [ -f "databases/pfam/Pfam-A.hmm" ]; then
    wget -P databases/pfam https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.hmm.gz
    gunzip databases/pfam/Pfam-A.hmm.gz
    hmmpress databases/pfam/Pfam-A.hmm
fi
