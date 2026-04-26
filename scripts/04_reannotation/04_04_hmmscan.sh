#!/bin/bash

mkdir -p reannotation/reference/pfam reannotation/geneml/pfam

# Run PFAM hmmscan on proteins of sufficient length
for f in training_genomes/*.gff; do
    name=$(basename $f .gff)
    for dir in reference geneml; do
        echo reannotation/$dir/original/${name}_long.faa >> hmmscan_inputs_$dir
        python3 scripts/python/run_pyhmmer_hmmscan.py \
        databases/pfam/Pfam-A.hmm --cut_tc \
        --batch-tsv hmmscan_inputs_$dir \
        --output-dir reannotation/$dir/pfam \
        --cpu 2 --batch-jobs 12
    done
done

