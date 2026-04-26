#!/bin/bash

mkdir benchmarking_results/reference

# Extract proteins
for a in GCF_021901695.1 GCF_026210795.1; do
    gffread benchmarking_dataset/$a/chromosomes_simplified.gff \
    -g benchmarking_dataset/$a/chromosomes.fna \
    -y benchmarking_results/reference/$a.faa
    for tool in annevo augustus braker3 geneml helixer; do
        gffread benchmarking_results/$tool/$a.gff3 \
        -g benchmarking_dataset/$a/chromosomes.fna \
        -y benchmarking_results/$tool/$a.faa
        python3 scripts/python/filter_longest_proteins.py benchmarking_results/$tool/$a.faa
    done
done

# Run hmmscan
for out in annevo augustus braker3 geneml helixer reference; do
    for a in GCF_021901695.1 GCF_026210795.1; do
        echo benchmarking_results/$out/${a}_long.faa >> hmmscan_inputs_$out
    mkdir benchmarking_results/$out/pfam
    python3 scripts/python/run_pyhmmer_hmmscan.py \
    databases/pfam/Pfam-A.hmm --cut_tc \
    --batch-tsv hmmscan_inputs_$out \
    --output-dir benchmarking_results/$out/pfam \
    --cpu 8 --batch-jobs 2
    done
done

# Process hits
for a in GCF_021901695.1 GCF_026210795.1; do
    python3 scripts/python/get_annotated_transcripts.py \
    --hmmscan benchmarking_results/reference/pfam/${a}.domtblout.tsv \
    --gff benchmarking_dataset/$a/chromosomes_simplified.gff \
    --output benchmarking_results/reference/pfam/${a}_annotations.tsv
    for out in annevo augustus braker3 geneml helixer; do
        python3 scripts/python/get_annotated_transcripts.py \
        --hmmscan benchmarking_results/$out/pfam/${a}.domtblout.tsv \
        --gff benchmarking_results/$out/$a.gff3 \
        --output benchmarking_results/$out/pfam/${a}_annotations.tsv
    done
done
