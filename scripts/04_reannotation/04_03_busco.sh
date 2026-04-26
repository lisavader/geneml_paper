#!/bin/bash

mkdir -p reannotation/reference/busco reannotation/geneml/busco

# Run BUSCO on single isoform data
for f in training_genomes/*.gff; do
    name=$(basename $f .gff)
    for dir in reference geneml; do
        busco --mode protein --lineage fungi_odb12 -i reannotation/$dir/single_isoform/$name.faa \
        -o reannotation/$dir/busco/$name -c8
    done
done
