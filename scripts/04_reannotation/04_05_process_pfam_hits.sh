#!/bin/bash

for f in training_genomes/*.gff; do
    name=$(basename $f .gff)
    python3 scripts/python/get_annotated_transcripts.py \
    --hmmscan reannotation/reference/pfam/${name}_long.domtblout.tsv \
    --gff training_genomes/$name.gff \
    --output reannotation/reference/pfam/${name}_annotations.tsv
    python3 scripts/python/get_annotated_transcripts.py \
    --hmmscan reannotation/geneml/pfam/${name}_long.domtblout.tsv \
    --gff reannotation/geneml/original/$name.gff \
    --output reannotation/geneml/pfam/${name}_annotations.tsv
done

for f in training_genomes/*.gff; do
    name=$(basename $f .gff)
    for dir in reference geneml; do
        python3 scripts/python/extract_gene_functions.py \
        --transcript_annotations reannotation/$dir/pfam/${name}_annotations.tsv \
        --output reannotation/$dir/pfam/${name}_functions.tsv \
        --go-obo data/go.obo --pfam2go data/pfam2go --go-bins data/go_bins.tsv \
        --pfams-mge data/pfams_mge.txt
    done
done
