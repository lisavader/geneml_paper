#!/bin/bash

# Prepare directories
mkdir -p \
reannotation/reference/original \
reannotation/geneml/single_isoform \
reannotation/reference/single_isoform

# Prepare input data
for f in training_genomes/*.gff; do
    name=$(basename $f .gff)
    ln -s $(pwd)/$f reannotation/reference/original/$name.gff
    for dir in reference geneml; do
        gffread reannotation/$dir/original/$name.gff -g training_genomes/$name.fna \
        -y reannotation/$dir/original/$name.faa
        python3 scripts/python/get_longest_isoform_gff.py \
        reannotation/$dir/original/$name.gff \
        reannotation/$dir/single_isoform/$name.gff
        gffread reannotation/$dir/single_isoform/$name.gff -g training_genomes/$name.fna \
        -y reannotation/$dir/single_isoform/$name.faa
        python3 scripts/python/write_sampled_lengths.py reannotation/$dir/single_isoform/$name.faa
        python3 scripts/python/filter_longest_proteins.py reannotation/$dir/original/$name.faa
    done
done
