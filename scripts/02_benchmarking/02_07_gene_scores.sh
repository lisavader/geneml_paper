#!/bin/bash

# Use different branch for writing scores
cd geneML
git checkout write_scores
cd ..

python3 -m venv env_geneml_scores
. env_geneml_scores/bin/activate
pip install geneML

mkdir geneml_scores

rm -r geneml_scores_scripts
for a in $(ls benchmarking_dataset); do
    echo "geneml benchmarking_dataset/$a/chromosomes.fna -o geneml_scores/$a.gff3 \
--model models/benchmarking_final.keras --context-length 800 --cores 8 --cpu-only" \
    >> geneml_scores_scripts
done

parallel < geneml_scores_scripts

# Reset to default environment
cd geneML
git checkout main
cd ..

deactivate
source venv/bin/activate
