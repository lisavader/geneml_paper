#!/bin/bash

mkdir -p reannotation/geneml/original

# geneML scripts
rm -f geneml_reannotation_scripts
for g in training_genomes/*.fna; do
    echo geneml $g -o reannotation/geneml/original/$(basename $g .fna).gff --cpu-only --cores 8 \
    >> geneml_reannotation_scripts
done

parallel --joblog geneml_reannotation_log -j4 < geneml_reannotation_scripts
