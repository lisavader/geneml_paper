#!/bin/bash

# Post-processing (BRAKER3)
for a in $(ls benchmarking_results/braker3/); do
    cp benchmarking_results/braker3/$a/braker.gff3 benchmarking_results/braker3/$a.gff3
done

# Post-processing (Helixer)
mkdir benchmarking_results/helixer
for r in $(ls benchmarking_results/helixer_with_UTR/); do
    python3 scripts/python/simplify_gff.py benchmarking_results/helixer_with_UTR/$r \
    benchmarking_results/helixer/${r} --accession $(basename $r .gff)
done
