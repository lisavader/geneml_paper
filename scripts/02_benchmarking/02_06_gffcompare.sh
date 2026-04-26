#!/bin/bash

for acc in $(ls benchmarking_dataset); do
    for out in annevo augustus braker3 geneml helixer; do
        gffcompare -r benchmarking_dataset/$acc/chromosomes_simplified.gff \
        benchmarking_results/$out/$acc.gff3 -T --no-exon-merge --strict-match -e 0
        mkdir -p gffcompare_benchmarking/$acc/
        mv gffcmp.stats gffcompare_benchmarking/$acc/$out.stats
    done
done

#clean up
rm gffcmp.*
